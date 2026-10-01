process ASAP_PON_FILTER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0':
        'biocontainers/bcftools:1.20--h8b25389_0' }"

    input:
    tuple val(meta) , path(vcf), path(tbi)
    tuple val(meta2), path(asap_vcf), path(asap_tbi)
    tuple val(meta3), path(fasta)
    tuple val(meta4), path(fai)

    output:
    tuple val(meta), path("${prefix}.asap_flagged.vcf.gz"), path("${prefix}.asap_flagged.vcf.gz.tbi"), emit: flagged
    tuple val(meta), path("${prefix}.asap_pass.vcf.gz")   , path("${prefix}.asap_pass.vcf.gz.tbi")   , emit: vcf
    tuple val(meta), path("${prefix}.asap_stats.tsv")                                                , emit: stats

    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version |& sed '1!d ; s/bcftools //'"), topic: versions, emit: versions_bcftools

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix     = task.ext.prefix ?: "${meta.id}_somatic"
    def min_af = task.ext.min_af
    // With a minimum AF the panel's INFO/AF is carried over as INFO/ASAP_AF and a match only counts
    // when the variant is common enough in the panel; without one any panel match is flagged.
    def copy_af   = min_af != null ? "-c INFO/ASAP_AF:=INFO/AF -h asap_af.hdr" : ""
    def expr      = min_af != null ? "INFO/ASAP_PON_MATCH=1 && INFO/ASAP_AF>${min_af}" : "INFO/ASAP_PON_MATCH=1"
    def af_header = min_af != null ? "printf '##INFO=<ID=ASAP_AF,Number=A,Type=Float,Description=\"Allele frequency of the matching allele in the ASAP panel of normals\">\\n' > asap_af.hdr" : ""
    """
    # annotate -a reads the panel through an index; build one when none was supplied.
    if [ ! -e ${asap_vcf}.tbi ] && [ ! -e ${asap_vcf}.csi ]; then
        tabix -p vcf ${asap_vcf}
    fi

    # Split multiallelic sites and left-align so each record carries one ALT allele that can be
    # compared against the panel allele by allele. annotate needs an indexed target to stream the
    # panel alongside it, so the normalised calls are written and indexed rather than piped.
    bcftools norm -m -any -f ${fasta} ${vcf} -Ou \\
        | bcftools sort --temp-dir . -Oz -o normalised.vcf.gz
    tabix -p vcf normalised.vcf.gz

    ${af_header}

    # -m +ASAP_PON_MATCH sets the flag (and declares it in the header) only on records the panel
    # matches. --pair-logic some: the panel keeps multiallelic records (norm +both), so a split call
    # matches when its ALT is one of the panel record's ALTs at the same POS and REF; a different
    # ALT at the same position is not a match. Allele comparison is case-insensitive, which matters
    # because the panel REF/ALT carry soft-masked lowercase bases.
    bcftools annotate \\
        -a ${asap_vcf} \\
        --pair-logic some \\
        ${copy_af} \\
        -m +ASAP_PON_MATCH \\
        normalised.vcf.gz -Ou \\
        | bcftools filter -m + -s ASAP_PON -e '${expr}' -Oz -o ${prefix}.asap_flagged.vcf.gz
    tabix -p vcf ${prefix}.asap_flagged.vcf.gz

    bcftools view -e 'FILTER~"ASAP_PON"' ${prefix}.asap_flagged.vcf.gz -Oz -o ${prefix}.asap_pass.vcf.gz
    tabix -p vcf ${prefix}.asap_pass.vcf.gz

    total=\$(bcftools view -H ${prefix}.asap_flagged.vcf.gz | wc -l)
    flagged=\$(bcftools view -H -i 'FILTER~"ASAP_PON"' ${prefix}.asap_flagged.vcf.gz | wc -l)
    printf 'sample\\ttotal\\tasap_flagged\\n%s\\t%s\\t%s\\n' "${meta.id}" "\$total" "\$flagged" > ${prefix}.asap_stats.tsv
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_somatic"
    """
    echo "" | gzip > ${prefix}.asap_flagged.vcf.gz
    touch ${prefix}.asap_flagged.vcf.gz.tbi
    echo "" | gzip > ${prefix}.asap_pass.vcf.gz
    touch ${prefix}.asap_pass.vcf.gz.tbi
    printf 'sample\\ttotal\\tasap_flagged\\n${meta.id}\\t0\\t0\\n' > ${prefix}.asap_stats.tsv
    """
}
