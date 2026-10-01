process SEVERUS_PON_FILTER {
    tag "$meta.id"
    label 'process_single'

    // The Severus image: python3 for bin/severus_pon_filter.py, plus bgzip and tabix
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/severus:1.6--pyhdfd78af_0':
        'biocontainers/severus:1.6--pyhdfd78af_0' }"

    input:
    // Severus SV calls (SEVERUS.out.somatic_vcf); no index needed, the VCF is read as a stream
    tuple val(meta) , path(vcf)
    // Severus PON (chr1,pos1,chr2,pos2,ci1,ci2,svtype,vaf) and, optionally ([]), the VNTR BED that
    // Severus was run with; without it INSIDE_VNTR calls are matched with the non-VNTR window
    tuple val(meta2), path(pon), path(vntr_bed)

    output:
    tuple val(meta), path("${prefix}.pon_flagged.vcf.gz"), path("${prefix}.pon_flagged.vcf.gz.tbi"), emit: flagged
    tuple val(meta), path("${prefix}.pon_pass.vcf.gz")   , path("${prefix}.pon_pass.vcf.gz.tbi")   , emit: vcf
    tuple val(meta), path("${prefix}.pon_stats.tsv")                                               , emit: stats

    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('tabix') , eval("tabix --version | sed '1!d; s/.* //'")  , topic: versions, emit: versions_tabix

    when:
    task.ext.when == null || task.ext.when

    script:
    // ext.args: e.g. "--ci 25" to assume a breakpoint CI (Severus does not write it to the VCF)
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}_severus_somatic"
    def vntr = vntr_bed ? "--vntr-bed ${vntr_bed}" : ''
    """
    # Re-applies Severus 1.6's own PON test (extract_pon + add_pon) to calls Severus already made,
    # so matched tumour/normal runs keep the normal-based somatic calls and still get the PON.
    severus_pon_filter.py \\
        --vcf ${vcf} \\
        --pon ${pon} \\
        ${vntr} \\
        --sample ${meta.id} \\
        ${args} \\
        --flagged ${prefix}.pon_flagged.vcf \\
        --passed ${prefix}.pon_pass.vcf \\
        --stats ${prefix}.pon_stats.tsv

    bgzip ${prefix}.pon_flagged.vcf
    tabix -p vcf ${prefix}.pon_flagged.vcf.gz
    bgzip ${prefix}.pon_pass.vcf
    tabix -p vcf ${prefix}.pon_pass.vcf.gz
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_severus_somatic"
    """
    echo "" | bgzip -c > ${prefix}.pon_flagged.vcf.gz
    touch ${prefix}.pon_flagged.vcf.gz.tbi
    echo "" | bgzip -c > ${prefix}.pon_pass.vcf.gz
    touch ${prefix}.pon_pass.vcf.gz.tbi
    printf 'sample\\ttotal\\tpon_flagged\\n${meta.id}\\t0\\t0\\n' > ${prefix}.pon_stats.tsv
    """
}
