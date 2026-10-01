process CH_VARIANTS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0':
        'biocontainers/bcftools:1.20--h8b25389_0' }"

    input:
    // Each arm is staged into its own directory so that same-named somatic and germline VCFs
    // (e.g. two VEP runs without distinct prefixes) cannot collide in the task directory.
    tuple val(meta), path(som_vcf, stageAs: 'somatic/*'), path(som_tbi, stageAs: 'somatic/*'), path(germ_vcf, stageAs: 'germline/*'), path(germ_tbi, stageAs: 'germline/*')
    path(gene_list)

    output:
    tuple val(meta), path("${prefix}.vcf.gz"), path("${prefix}.vcf.gz.tbi"), emit: vcf
    tuple val(meta), path("${prefix}.tsv")                                , emit: tsv

    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version |& sed '1!d ; s/bcftools //'"), topic: versions, emit: versions_bcftools

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}_CH_variants"
    """
    # Normalise the gene list: drop CRs, surrounding whitespace, blank lines and '#' comments.
    awk '{ sub(/\\r\$/, ""); gsub(/^[ \\t]+|[ \\t]+\$/, "") } \$0 != "" && \$0 !~ /^#/' ${gene_list} \\
        | sort -u > ch_genes.norm.txt
    if [ ! -s ch_genes.norm.txt ]; then
        echo "ERROR: CH gene list '${gene_list}' contains no gene symbols" >&2
        exit 1
    fi

    # Header lines for the two INFO fields stamped on every reported record. The double quotes are
    # written by printf so that no escaped quote has to survive the Nextflow script block.
    printf '##INFO=<ID=CH_GENE,Number=.,Type=String,Description="%s">\\n' 'CH gene(s) from VEP SYMBOL' >  ch.hdr
    printf '##INFO=<ID=CH_ORIGIN,Number=1,Type=String,Description="%s">\\n' 'somatic or germline call set' >> ch.hdr

    # ch_arm <origin> <vcf>
    # Writes <origin>.hits.tsv (one row per CSQ entry in a CH gene) and <origin>.ch.vcf.gz(.tbi)
    # (the input records hitting a CH gene, annotated with CH_GENE/CH_ORIGIN, sample renamed).
    ch_arm() {
        local origin=\$1 vcf=\$2

        # bcftools concat needs identical sample sets, and the two arms come from different
        # callers whose sample columns are not named alike, so each arm is renamed to meta.id.
        local nsamples
        nsamples=\$(bcftools query -l "\$vcf" | wc -l)
        if [ "\$nsamples" -ne 1 ]; then
            echo "ERROR: \$origin VCF '\$vcf' has \$nsamples samples; expected exactly one" >&2
            exit 1
        fi

        # Keep CSQ entries whose SYMBOL is a CH gene; the row layout matches the TSV output. This is
        # streamed rather than staged, since a whole-genome germline VCF expands to millions of CSQ
        # rows. The pipeline runs bash with -C (noclobber), so each file is written exactly once.
        # Header and field list go to files first: grep -q exiting early under pipefail could
        # otherwise SIGPIPE bcftools and turn a match into a false negative.
        bcftools view -h "\$vcf" > \$origin.header.txt
        if grep -q '^##INFO=<ID=CSQ,' \$origin.header.txt; then
            bcftools +split-vep -l "\$vcf" | cut -f2 > \$origin.csq_fields.txt
            local field
            for field in SYMBOL Consequence IMPACT; do
                if ! grep -qx "\$field" \$origin.csq_fields.txt; then
                    echo "ERROR: CSQ header of \$origin VCF '\$vcf' has no '\$field' field" >&2
                    exit 1
                fi
            done
            # -d: one output line per CSQ entry, so every gene annotated at a site is seen.
            # Records without a CSQ value print '.', which is never in the gene list.
            bcftools +split-vep -d \\
                -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%SYMBOL\\t%Consequence\\t%IMPACT\\n' \\
                "\$vcf" \\
                | awk -F'\\t' -v OFS='\\t' -v origin="\$origin" \\
                    'NR == FNR { g[\$1] = 1; next } (\$5 in g) { print \$1, \$2, \$3, \$4, origin, \$5, \$6, \$7 }' \\
                    ch_genes.norm.txt - \\
                | sort -u > \$origin.hits.tsv
        else
            echo "WARNING: \$origin VCF '\$vcf' has no VEP CSQ annotation; no \$origin CH variants reported" >&2
            : > \$origin.hits.tsv
        fi

        # One row per site, CH genes comma-joined (sorted, deduplicated), as a bcftools annotate source.
        awk -F'\\t' -v OFS='\\t' '
            {
                k = \$1 OFS \$2 OFS \$3 OFS \$4
                if (!(k in genes)) { keys[++n] = k; genes[k] = \$6; origin[k] = \$5; seen[k OFS \$6] = 1 }
                else if (!((k OFS \$6) in seen)) { genes[k] = genes[k] "," \$6; seen[k OFS \$6] = 1 }
            }
            END { for (i = 1; i <= n; i++) print keys[i], genes[keys[i]], origin[keys[i]] }
        ' \$origin.hits.tsv | sort -k1,1 -k2,2n > \$origin.sites.tsv

        if [ -s \$origin.sites.tsv ]; then
            bgzip -c \$origin.sites.tsv > \$origin.sites.tsv.gz
            tabix -s1 -b2 -e2 \$origin.sites.tsv.gz
            # REF and ALT are matched too, so a CH annotation never leaks onto another allele.
            bcftools annotate \\
                -a \$origin.sites.tsv.gz \\
                -c CHROM,POS,REF,ALT,INFO/CH_GENE,INFO/CH_ORIGIN \\
                -h ch.hdr \\
                -Ou "\$vcf" \\
                | bcftools view -i 'INFO/CH_ORIGIN!="."' -Oz -o \$origin.ch_tmp.vcf.gz
        else
            # No hits: a header-only VCF that still declares CH_GENE/CH_ORIGIN.
            bcftools annotate -h ch.hdr -Ou "\$vcf" \\
                | bcftools view -h -Oz -o \$origin.ch_tmp.vcf.gz
        fi

        echo "${meta.id}" > \$origin.samples.txt
        bcftools reheader -s \$origin.samples.txt -o \$origin.ch.vcf.gz \$origin.ch_tmp.vcf.gz
        tabix -p vcf \$origin.ch.vcf.gz
    }

    ch_arm somatic  ${som_vcf}
    ch_arm germline ${germ_vcf}

    # -a merges the two coordinate-sorted arms; a site called in both arms is kept twice, once per
    # CH_ORIGIN. concat -a requires indexed inputs, hence the tabix calls above.
    bcftools concat -a -Ou somatic.ch.vcf.gz germline.ch.vcf.gz \\
        | bcftools sort -T ./bcftools_sort_tmp -Oz -o ${prefix}.vcf.gz -W=tbi

    printf 'CHROM\\tPOS\\tREF\\tALT\\tCH_ORIGIN\\tSYMBOL\\tConsequence\\tIMPACT\\n' > ${prefix}.tsv
    # Natural chromosome order (1..22, X, Y, M, then others) via a numeric sort key, since the
    # container's BusyBox sort has no per-key version sort.
    cat somatic.hits.tsv germline.hits.tsv \\
        | awk -F'\\t' -v OFS='\\t' '{
            c = \$1; sub(/^chr/, "", c)
            r = (c ~ /^[0-9]+\$/) ? c + 0 : (c == "X" ? 1000 : (c == "Y" ? 1001 : ((c == "M" || c == "MT") ? 1002 : 2000)))
            print r, \$0
        }' \\
        | sort -k1,1n -k2,2 -k3,3n \\
        | cut -f2- >> ${prefix}.tsv
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_CH_variants"
    """
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    printf 'CHROM\\tPOS\\tREF\\tALT\\tCH_ORIGIN\\tSYMBOL\\tConsequence\\tIMPACT\\n' > ${prefix}.tsv
    """
}
