#!/usr/bin/env python3
"""
Flag Severus SV calls that match a Severus panel of normals (PON).

Severus 1.6 applies its PON inside `annotate_mut_type` (breakpoint_finder.py), and when a matched
normal is also given the PON check overwrites the normal-based somatic/germline decision. This
script re-applies the same PON test after the fact, to a VCF Severus already wrote, so matched
tumour/normal runs can keep their normal-based calls and still get the PON as an extra filter.

The matching reproduces Severus 1.6 `extract_pon` + `add_pon` (breakpoint_finder.py):

* PON lines are `chr1,pos1,chr2,pos2,ci1,ci2,svtype,vaf` (comma-separated, no header; lines whose
  chromosomes are not usable are skipped, which also skips the `V1,...` header of the CHM13 PON).
  INS lines are stored as (pos1, ci1, pos1 + pos2, ci2, chr1): their pos2 column is the length.
  Other intra-chromosomal lines are stored as (pos1, ci1, pos2, ci2, chr2) under chr1, kept in
  file order (Severus does not sort them; the published PONs are sorted). Inter-chromosomal lines
  are stored under both (chr1, chr2) and (chr2, chr1) and sorted, as Severus sorts them.
* For an SV with breakpoints bp_1 / bp_2 (bp_1 is POS of the record, or of the `_1` breakend):
  candidates are PON entries whose first position lies in
  [bp_1 - CI_1 - 2000, bp_1 + CI_1 + 2000] (VNTR SVs: [vntr_start - 25, vntr_end + 25]) and whose
  second chromosome equals bp_2's chromosome. SV type and strands are not compared.
* len_diff = |(pon_pos2 - pon_pos1) - SVLEN| for intra-chromosomal SVs, |pon_pos2 - bp_2| for
  translocations, 0 for single breakends; sum_diff = |pon_pos1 - bp_1| + len_diff;
  sum_CI = max(pon_ci1, pon_ci2) + max(CI_1, CI_2).
  Match when (len_diff <= 50 and sum_diff <= 1000) or sum_diff <= sum_CI + max_diff, where
  max_diff is 150 (1000 for VNTR SVs).
* The decision is per breakpoint pair, so both breakends of a BND pair get the same verdict.

Severus does not write the per-breakpoint CI to the VCF. `--ci` sets the value assumed for both
breakpoints (default 0). Every term grows with the CI, so with the default every flagged call is
one Severus would also have flagged; Severus can flag a few more when its CI was larger.
"""

import argparse
import bisect
import gzip
import re
import sys
from array import array
from collections import defaultdict

BUFF = 2000
CLUST_LEN = 150
VNTR_BUFF = 25
VNTR_CLUST_LEN = 1000
MAX_LEN_DIFF = 50
VNTR_BP_TOL = 25  # resolve_vntr.read_vntr_file

FILTER_ID = "PON"
INFO_ID = "PON_MATCH"

ALT_MATE = re.compile(r"[\[\]]([^\[\]:]+):(\d+)[\[\]]")


def open_text(path):
    with open(path, "rb") as fh:
        magic = fh.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt")
    return open(path)


def parse_info(info):
    out = {}
    if info in ("", "."):
        return out
    for item in info.split(";"):
        key, _, value = item.partition("=")
        out[key] = value if _ else True
    return out


class Sv:
    """One Severus breakpoint pair, built from one VCF record or a BND mate pair."""

    __slots__ = ("key", "chr1", "pos1", "chr2", "pos2", "length", "single", "vntr_flag", "vntr")

    def __init__(self, key, chr1, pos1, chr2, pos2, length, single, vntr_flag):
        self.key = key
        self.chr1, self.pos1 = chr1, pos1
        self.chr2, self.pos2 = chr2, pos2
        self.length = length
        self.single = single
        self.vntr_flag = vntr_flag
        self.vntr = None


def record_to_sv(fields):
    chrom, pos, vid, alt = fields[0], int(fields[1]), fields[2], fields[4]
    info = parse_info(fields[7])
    svtype = info.get("SVTYPE", "")
    vntr_flag = "INSIDE_VNTR" in info
    svlen = info.get("SVLEN")
    svlen = abs(int(svlen)) if svlen not in (None, True, "") else None

    if svtype == "sBND" or alt in (".N", "N."):
        # single breakend: bp_2 is bp_1 (DoubleBreak(s_bp, ..., s_bp, ...)), len_diff = 0
        return Sv(vid, chrom, pos, chrom, pos, 0, True, False), None

    mate = ALT_MATE.search(alt)
    if mate:
        chr2, pos2 = mate.group(1), int(mate.group(2))
        base, sep, side = vid.rpartition("_")
        if not sep or side not in ("1", "2"):
            base, side = vid, "1"
        if side == "2":
            # bp_1 is the `_1` breakend; evaluate from it so the window sits where Severus put it
            chrom, pos, chr2, pos2 = chr2, pos2, chrom, pos
        length = svlen if svlen is not None else (abs(pos2 - pos) if chr2 == chrom else 0)
        return Sv(base, chrom, pos, chr2, pos2, length, False, vntr_flag), side

    end = info.get("END")
    end = int(end) if end not in (None, True, "") else pos
    if svtype == "INS":
        # the insertion's second breakpoint sits at the same position (bp_3 in extract_insertions)
        end = pos
    length = svlen if svlen is not None else abs(end - pos)
    return Sv(vid, chrom, pos, chrom, end, length, False, vntr_flag), None


class Pon:
    """PON entries for the keys the calls need, laid out as Severus' extract_pon builds them."""

    def __init__(self):
        self.lists = {}

    @staticmethod
    def _new():
        return [array("q"), array("q"), array("q"), array("q"), [], [], array("d")]

    def load(self, path, needed_keys):
        intra_needed = {k for k in needed_keys if isinstance(k, str)}
        inter_needed = {k for k in needed_keys if isinstance(k, tuple)}
        for k in needed_keys:
            self.lists[k] = self._new()
        unsorted = set()
        last = {}
        n = 0
        with open_text(path) as fh:
            for line in fh:
                parts = line.strip().split(",")
                if len(parts) != 8:
                    parts = line.strip().split("\t")
                if len(parts) != 8:
                    continue
                chr1, pos1, chr2, pos2, ci1, ci2, svtype, vaf = parts
                try:
                    pos1, pos2, ci1, ci2 = int(pos1), int(pos2), int(ci1), int(ci2)
                except ValueError:
                    continue  # header line (e.g. V1,V2,...)
                try:
                    vaf = float(vaf)
                except ValueError:
                    vaf = float("nan")
                n += 1
                if svtype == "INS" or chr1 == chr2:
                    if chr1 not in intra_needed:
                        continue
                    end = pos1 + pos2 if svtype == "INS" else pos2
                    second = chr1 if svtype == "INS" else chr2
                    self._append(chr1, pos1, ci1, end, ci2, second, svtype, vaf)
                    if pos1 < last.get(chr1, pos1):
                        unsorted.add(chr1)
                    last[chr1] = pos1
                else:
                    if (chr1, chr2) in inter_needed:
                        self._append((chr1, chr2), pos1, ci1, pos2, ci2, chr2, svtype, vaf)
                    if (chr2, chr1) in inter_needed:
                        self._append((chr2, chr1), pos2, ci2, pos1, ci1, chr1, svtype, vaf)
        for key in inter_needed:
            ls = self.lists[key]
            if ls[0]:
                rows = sorted(zip(ls[0], ls[1], ls[2], ls[3], ls[4], ls[5], ls[6]),
                              key=lambda r: (r[0], r[1], r[2], r[3], r[4]))
                cols = list(zip(*rows))
                self.lists[key] = [array("q", cols[0]), array("q", cols[1]), array("q", cols[2]),
                                   array("q", cols[3]), list(cols[4]), list(cols[5]), array("d", cols[6])]
        if unsorted:
            sys.stderr.write("WARNING: PON first positions are not sorted on %s; matching follows "
                             "Severus and bisects the file order\n" % ",".join(sorted(unsorted)))
        return n

    def _append(self, key, pos1, ci1, pos2, ci2, chr2, svtype, vaf):
        ls = self.lists[key]
        ls[0].append(pos1)
        ls[1].append(ci1)
        ls[2].append(pos2)
        ls[3].append(ci2)
        ls[4].append(chr2)
        ls[5].append(svtype)
        ls[6].append(vaf)

    def match(self, sv, ci):
        """Severus add_pon for one breakpoint pair; returns the matching PON entry or None."""
        if sv.chr1 == sv.chr2:
            ls = self.lists.get(sv.chr1)
        else:
            ls = self.lists.get((sv.chr1, sv.chr2))
        if not ls or not ls[0]:
            return None
        lo, hi = sv.pos1 - ci - BUFF, sv.pos1 + ci + BUFF
        if sv.vntr:
            lo, hi = sv.vntr[0] - VNTR_BUFF, sv.vntr[1] + VNTR_BUFF
        i0 = bisect.bisect_left(ls[0], lo)
        i1 = bisect.bisect_right(ls[0], hi)
        max_diff = VNTR_CLUST_LEN if sv.vntr else CLUST_LEN
        for i in range(i0, i1):
            if ls[4][i] != sv.chr2:
                continue
            if sv.single:
                len_diff = 0
            elif sv.chr1 == sv.chr2:
                len_diff = abs(ls[2][i] - ls[0][i] - sv.length)
            else:
                len_diff = abs(ls[2][i] - sv.pos2)
            sum_diff = abs(ls[0][i] - sv.pos1) + len_diff
            sum_ci = max(ls[1][i], ls[3][i]) + ci
            if (len_diff <= MAX_LEN_DIFF and sum_diff <= VNTR_CLUST_LEN) or sum_diff <= sum_ci + max_diff:
                return i, ls
        return None


def read_vntr(path):
    """resolve_vntr.read_vntr_file: per chromosome, starts/ends widened by 25 bp, file order."""
    vntr = defaultdict(lambda: (array("q"), array("q")))
    with open_text(path) as fh:
        for line in fh:
            parts = line.strip().split()
            if len(parts) < 3:
                continue
            try:
                start, end = int(parts[1]), int(parts[2])
            except ValueError:
                continue
            starts, ends = vntr[parts[0]]
            starts.append(start - VNTR_BP_TOL)
            ends.append(end + VNTR_BP_TOL)
    return vntr


def check_vntr(sv, vntr):
    """breakpoint_finder.check_vntr: the VNTR interval containing both breakpoints."""
    reg = vntr.get(sv.chr1)
    if not reg or not reg[0]:
        return None
    strt = bisect.bisect_right(reg[0], sv.pos1)
    end = bisect.bisect_left(reg[1], sv.pos2)
    if strt - end == 1:
        return (reg[0][strt - 1], reg[1][strt - 1])
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vcf", required=True, help="Severus VCF (plain or gzipped)")
    ap.add_argument("--pon", required=True, help="Severus PON (chr1,pos1,chr2,pos2,ci1,ci2,svtype,vaf)")
    ap.add_argument("--vntr-bed", help="VNTR BED given to Severus (--vntr-bed); needed for INSIDE_VNTR calls")
    ap.add_argument("--ci", type=int, default=0, help="breakpoint CI assumed for every call [0]")
    ap.add_argument("--flagged", required=True, help="output VCF with every record, PON matches flagged")
    ap.add_argument("--passed", required=True, help="output VCF without the PON matches")
    ap.add_argument("--stats", required=True, help="output stats TSV")
    ap.add_argument("--sample", required=True, help="sample name for the stats TSV")
    args = ap.parse_args()

    header, records = [], []
    with open_text(args.vcf) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("#"):
                header.append(line)
            else:
                records.append(line.split("\t"))
    if not header or not header[-1].startswith("#CHROM"):
        header = [h for h in header if not h.startswith("#CHROM")]
        if not header:
            header = ["##fileformat=VCFv4.2"]
        header.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")

    # one breakpoint pair per record, or per BND pair
    svs = {}
    rec_key = []
    for f in records:
        sv, side = record_to_sv(f)
        if sv.key in svs:
            # keep the `_1` view of a BND pair (it carries bp_1)
            if side == "1":
                svs[sv.key] = sv
        else:
            svs[sv.key] = sv
        rec_key.append(sv.key)

    vntr = read_vntr(args.vntr_bed) if args.vntr_bed else None
    vntr_missed = 0
    for sv in svs.values():
        if sv.vntr_flag:
            sv.vntr = check_vntr(sv, vntr) if vntr is not None else None
            if sv.vntr is None:
                vntr_missed += 1
    if vntr_missed:
        sys.stderr.write("WARNING: %d INSIDE_VNTR call(s) had no VNTR interval%s; matched with the "
                         "non-VNTR window\n" % (vntr_missed, "" if vntr is not None else " (no --vntr-bed)"))

    needed = set()
    for sv in svs.values():
        needed.add(sv.chr1 if sv.chr1 == sv.chr2 else (sv.chr1, sv.chr2))
    pon = Pon()
    n_pon = pon.load(args.pon, needed) if needed else 0

    hits = {}
    for key, sv in svs.items():
        m = pon.match(sv, args.ci)
        if m:
            i, ls = m
            vaf = ls[6][i]
            vaf_s = "NA" if vaf != vaf else "%g" % vaf
            hits[key] = "%s:%d|%s:%d|%s|%s" % (sv.chr1, ls[0][i], ls[4][i], ls[2][i], ls[5][i], vaf_s)

    head = [h for h in header[:-1]
            if not h.startswith("##FILTER=<ID=%s," % FILTER_ID) and not h.startswith("##INFO=<ID=%s," % INFO_ID)]
    head.append('##FILTER=<ID=%s,Description="Breakpoints match the Severus panel of normals '
                '(Severus add_pon criteria, applied after calling)">' % FILTER_ID)
    head.append('##INFO=<ID=%s,Number=1,Type=String,Description="Matching panel-of-normals entry: '
                'chr1:pos1|chr2:pos2|SVTYPE|VAF (INS entries give pos1+length as pos2)">' % INFO_ID)
    head.append('##severus_pon_filter=<PON=%s,CI=%d,VNTR_BED=%s>'
                % (args.pon.split("/")[-1], args.ci, (args.vntr_bed or "none").split("/")[-1]))
    head.append(header[-1])

    flagged = 0
    with open(args.flagged, "x") as out_f, open(args.passed, "x") as out_p:
        out_f.write("\n".join(head) + "\n")
        out_p.write("\n".join(head) + "\n")
        for f, key in zip(records, rec_key):
            if key in hits:
                flagged += 1
                f = list(f)
                f[6] = FILTER_ID if f[6] in ("PASS", ".", "") else f[6] + ";" + FILTER_ID
                tag = "%s=%s" % (INFO_ID, hits[key])
                f[7] = tag if f[7] in (".", "") else f[7] + ";" + tag
                out_f.write("\t".join(f) + "\n")
            else:
                line = "\t".join(f) + "\n"
                out_f.write(line)
                out_p.write(line)

    with open(args.stats, "x") as out_s:
        out_s.write("sample\ttotal\tpon_flagged\n%s\t%d\t%d\n" % (args.sample, len(records), flagged))

    sys.stderr.write("severus_pon_filter: %d PON entries read, %d record(s), %d breakpoint pair(s), "
                     "%d record(s) flagged\n" % (n_pon, len(records), len(svs), flagged))


if __name__ == "__main__":
    main()
