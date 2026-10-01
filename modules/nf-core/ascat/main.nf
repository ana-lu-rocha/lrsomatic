process ASCAT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/4c/4cf02c7911ee5e974ce7db978810770efbd8d872ff5ab3462d2a11bcf022fab5/data':
        'community.wave.seqera.io/library/ascat_cancerit-allelecount:c3e8749fa4af0e99' }"

    input:
    tuple val(meta), path(input_normal), path(index_normal), path(input_tumor), path(index_tumor)
    val(genomeVersion)
    path(allele_files)
    path(loci_files)
    path(bed_file)  // optional
    path(fasta)     // optional
    path(gc_file)   // optional
    path(rt_file)   // optional

    output:
    tuple val(meta), path("*alleleFrequencies_chr*.txt"),      emit: allelefreqs
    tuple val(meta), path("*BAF.txt"),                         emit: bafs
    tuple val(meta), path("*cnvs.txt"),                        emit: cnvs
    tuple val(meta), path("*LogR.txt"),                        emit: logrs
    tuple val(meta), path("*metrics.txt"),                     emit: metrics
    tuple val(meta), path("*png"),                             emit: png
    tuple val(meta), path("*pdf"),                             emit: pdf, optional: true
    tuple val(meta), path("*purityploidy.txt"),                emit: purityploidy
    tuple val(meta), path("*segments.txt"),                    emit: segments
    tuple val(meta), path("*segments_raw.txt"),                emit: segments_raw, optional: true
    path "versions.yml",                                       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args           = task.ext.args        ?: ''
    def prefix         = task.ext.prefix      ?: "${meta.id}"
    def gender         = args.gender          ?  "$args.gender" :        "NULL"
    def purity         = args.purity          ?  "$args.purity" :        "NULL"
    def ploidy         = args.ploidy          ?  "$args.ploidy" :        "NULL"
    def penalty        = args.penalty         ?  "$args.penalty" :       "NULL"
    def gc_input       = gc_file              ?  "$gc_file" :            "NULL"
    def rt_input       = rt_file              ?  "$rt_file" :            "NULL"
    def pdf_plots      = (args.pdf_plots ?: false) ? "TRUE" : "FALSE"
    def minCounts_arg                    = args.minCounts                     ?  ",minCounts = $args.minCounts" : ""
    def bed_file_arg                     = bed_file                           ?  ",BED_file = '$bed_file'": ""
    def chrom_names_arg                  = args.chrom_names                   ?  ",chrom_names = $args.chrom_names" : ""
    def min_base_qual_arg                = args.min_base_qual                 ?  ",min_base_qual = $args.min_base_qual" : ""
    def min_map_qual_arg                 = args.min_map_qual                  ?  ",min_map_qual = $args.min_map_qual" : ""
    def fasta_arg                        = fasta                              ?  ",ref.fasta = '$fasta'" : ""
    def skip_allele_counting_tumour_arg  = args.skip_allele_counting_tumour   ?  ",skip_allele_counting_tumour = $args.skip_allele_counting_tumour" : ""
    def skip_allele_counting_normal_arg  = args.skip_allele_counting_normal   ?  ",skip_allele_counting_normal = $args.skip_allele_counting_normal" : ""
    
    def normal_exists                    = input_normal                       ? 'TRUE' : 'FALSE'
    def normal_bam                       = input_normal                       ? ",normalseqfile = '$input_normal'" : ""
    def normal_name                      = input_normal                       ? ",normalname = '${prefix}.normal'" : ""
    def longread_bins                    = args.longread_bins                 ? ",loci_binsize = $args.longread_bins" : ""
    def allele_counter_flags             = args.allele_counter_flags          ? ",additional_allelecounter_flags = '$args.allele_counter_flags'" : ""
    """
    #!/usr/bin/env Rscript
    library(RColorBrewer)
    library(ASCAT)
    options(bitmapType='cairo')

    # ASCAT draws every plot with an unqualified png() call. Returns a copy of f (and of the
    # helpers it calls that are named in 'swap') whose png() opens a PDF of the same size instead,
    # so each plot can additionally be written as PDF while the PNGs stay untouched.
    pdf_instead_of_png <- function(f, swap = character()) {
        e <- new.env(parent = environment(f))
        e\$png <- function(filename = "Rplot%03d.png", width = 480, height = 480, units = "px", pointsize = 12, res = NA, ...) {
            if (is.na(res)) res <- 72
            inches <- switch(units, px = 1 / res, `in` = 1, cm = 1 / 2.54, mm = 1 / 25.4)
            grDevices::pdf(file = sub("[.]png\$", ".pdf", filename), width = width * inches, height = height * inches, pointsize = pointsize)
        }
        for (name in swap) {
            g <- get(name, envir = environment(f))
            environment(g) <- e
            assign(name, g, envir = e)
        }
        environment(f) <- e
        f
    }

    # Run an ASCAT plotting function as is (PNGs) and, with pdf_plots, once more to get the PDFs
    # too. Returns the result of the first run.
    with_pdf_plots <- function(f, ..., swap = character()) {
        out <- f(...)
        if ($pdf_plots) invisible(pdf_instead_of_png(f, swap)(...))
        out
    }

    #build prefixes: <abspath_to_files/prefix_chr>
    allele_path = normalizePath("$allele_files")
    allele_prefix = paste0(allele_path, "/", "$allele_files", "_chr")

    loci_path = normalizePath("$loci_files")
    loci_prefix = paste0(loci_path, "/", "$loci_files", "_chr")

    #prepare from BAM files
    ascat.prepareHTS(
        tumourseqfile = "$input_tumor",
        tumourname = paste0("$prefix", ".tumour"),
        allelecounter_exe = "alleleCounter",
        alleles.prefix = allele_prefix,
        loci.prefix = loci_prefix,
        gender = "$gender",
        genomeVersion = "$genomeVersion",
        nthreads = $task.cpus
        $normal_bam
        $normal_name
        $minCounts_arg
        $bed_file_arg
        $chrom_names_arg
        $min_base_qual_arg
        $min_map_qual_arg
        $longread_bins
        $fasta_arg
        $allele_counter_flags
        $skip_allele_counting_tumour_arg
        $skip_allele_counting_normal_arg,
        seed = 42
    )


    #Load the data
    if($normal_exists) {
        print("normal exists")
        ascat.bc = ascat.loadData(
            Tumor_LogR_file = paste0("$prefix", ".tumour_tumourLogR.txt"),
            Tumor_BAF_file = paste0("$prefix", ".tumour_tumourBAF.txt"),
            Germline_LogR_file = paste0("$prefix", ".tumour_normalLogR.txt"),
            Germline_BAF_file = paste0("$prefix", ".tumour_normalBAF.txt"),
            genomeVersion = "$genomeVersion",
            gender = "$gender"
        )
    } else {
        print("normal does not exist")
        ascat.bc = ascat.loadData(
            Tumor_LogR_file = paste0("$prefix", ".tumour_tumourLogR.txt"),
            Tumor_BAF_file = paste0("$prefix", ".tumour_tumourBAF.txt"),
            genomeVersion = "$genomeVersion",
            gender = "$gender")
        gg = ascat.predictGermlineGenotypes(ascat.bc, platform = "WGS_hg38_50X")
        
    }
    print("printing ascat.bc")
    print(ascat.bc)

    #Plot the raw data
    with_pdf_plots(ascat.plotRawData, ascat.bc, img.prefix = paste0("$prefix", ".before_correction."))

    # optional LogRCorrection
    if("$gc_input" != "NULL") {
        gc_input = paste0(normalizePath("$gc_input"))

        if("$rt_input" != "NULL"){
            rt_input = paste0(normalizePath("$rt_input"))
            ascat.bc = ascat.correctLogR(ascat.bc, GCcontentfile = gc_input, replictimingfile = rt_input)
            #Plot raw data after correction
            with_pdf_plots(ascat.plotRawData, ascat.bc, img.prefix = paste0("$prefix", ".after_correction_gc_rt."))
        }
        else {
            ascat.bc = ascat.correctLogR(ascat.bc, GCcontentfile = gc_input)
            #Plot raw data after correction
            with_pdf_plots(ascat.plotRawData, ascat.bc, img.prefix = paste0("$prefix", ".after_correction_gc."))
        }
    }

    #Segment the data
    if($normal_exists) {
        ascat.bc = ascat.aspcf(ascat.bc, seed=42, penalty = $penalty)
    } else {
        ascat.bc = ascat.aspcf(ascat.bc, seed=42, penalty = $penalty, ascat.gg = gg)
    }

    #Plot the segmented data
    with_pdf_plots(ascat.plotSegmentedData, ascat.bc)

    #Run ASCAT to fit every tumor to a model, inferring ploidy, normal cell contamination, and discrete copy numbers
    #rho (purity) and psi (ploidy) are only passed when manually set
    #The sunrise and profile plots are drawn by the internal runASCAT, so its png() is swapped too.
    #runAscat is deterministic, so the PDF re-run draws the same fit; pdfPlot stays FALSE so the PNGs are kept.
    runAscat_args <- list(ascat.bc, gamma = 1)
    if (!is.null($purity)) runAscat_args\$rho_manual <- $purity
    if (!is.null($ploidy)) runAscat_args\$psi_manual <- $ploidy
    ascat.output <- do.call(with_pdf_plots, c(list(ascat.runAscat), runAscat_args, list(swap = "runASCAT")))

    #Extract metrics from ASCAT profiles
    QC = ascat.metrics(ascat.bc,ascat.output)

    #Write out segmented regions (including regions with one copy of each allele)
    write.table(ascat.output[["segments"]], file=paste0("$prefix", ".segments.txt"), sep="\t", quote=F, row.names=F)
    
    #Write out raw segmented regions (including regions with one copy of each allele)
    tryCatch({ # In case segments_raw is not selected
      write.table(
        ascat.output[["segments_raw"]],
        file = paste0("$prefix", ".segments_raw.txt"),
        sep = "\t", quote = FALSE, row.names = FALSE
      )
    }, error = function(e) {
      message("Error in writing segments_raw: ", conditionMessage(e))
    })

    #Write out CNVs in bed format
    cnvs=ascat.output[["segments"]][2:6]
    write.table(cnvs, file=paste0("$prefix",".cnvs.txt"), sep="\t", quote=F, row.names=F, col.names=T)

    #Write out purity and ploidy info
    summary <- tryCatch({
            matrix(c(ascat.output[["aberrantcellfraction"]], ascat.output[["ploidy"]]), ncol=2, byrow=TRUE)}, error = function(err) {
                # error handler picks up where error was generated
                print(paste("Could not find optimal solution:  ",err))
                return(matrix(c(0,0),nrow=1,ncol=2,byrow = TRUE))
        }
    )
    colnames(summary) <- c("AberrantCellFraction","Ploidy")
    write.table(summary, file=paste0("$prefix",".purityploidy.txt"), sep="\t", quote=F, row.names=F, col.names=T)

    write.table(QC, file=paste0("$prefix", ".metrics.txt"), sep="\t", quote=F, row.names=F)

    # version export
    f <- file("versions.yml","w")
    alleleCounter_version = system(paste("alleleCounter --version"), intern = T)
    ascat_version = sessionInfo()\$otherPkgs\$ASCAT\$Version
    writeLines(paste0('"', "$task.process", '"', ":"), f)
    writeLines(paste("    alleleCounter:", alleleCounter_version), f)
    writeLines(paste("    ascat:", ascat_version), f)
    close(f)

    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def pdf_plots = task.ext.args instanceof Map && task.ext.args.pdf_plots
    """
    echo stub > ${prefix}.after_correction.gc_rt.test.tumour.germline.png
    echo stub > ${prefix}.after_correction.gc_rt.test.tumour.tumour.png
    echo stub > ${prefix}.before_correction.test.tumour.germline.png
    echo stub > ${prefix}.before_correction.test.tumour.tumour.png
    echo stub > ${prefix}.cnvs.txt
    echo stub > ${prefix}.metrics.txt
    echo stub > ${prefix}.normal_alleleFrequencies_chr21.txt
    echo stub > ${prefix}.normal_alleleFrequencies_chr22.txt
    echo stub > ${prefix}.purityploidy.txt
    echo stub > ${prefix}.segments.txt
    echo stub > ${prefix}.segments_raw.txt
    echo stub > ${prefix}.tumour.ASPCF.png
    echo stub > ${prefix}.tumour.sunrise.png
    echo stub > ${prefix}.tumour_alleleFrequencies_chr21.txt
    echo stub > ${prefix}.tumour_alleleFrequencies_chr22.txt
    echo stub > ${prefix}.tumour_normalBAF.txt
    echo stub > ${prefix}.tumour_normalLogR.txt
    echo stub > ${prefix}.tumour_tumourBAF.txt
    echo stub > ${prefix}.tumour_tumourLogR.txt
    if [ "${pdf_plots}" = "true" ]; then
        for png in *.png; do echo stub > "\${png%.png}.pdf"; done
    fi

    echo "${task.process}:" > versions.yml
    echo ' alleleCounter: 4.3.0' >> versions.yml
    echo ' ascat: 3.2.0' >> versions.yml

    """

}
