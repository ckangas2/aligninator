process STARSOLO {
    tag "${sample}"
    label 'process_star_align'
    publishDir "${params.outdir}/samples/${sample}", mode: params.publish_dir_mode,
        saveAs: { f -> f.endsWith('.bam') && !params.save_bam ? null : f }

    input:
    tuple val(sample), path(r1), path(r2)
    path index
    path whitelist

    output:
    tuple val(sample), path('Solo.out'),                 emit: solo
    tuple val(sample), path("${sample}.star_summary.csv"), emit: star_summary
    path 'Log.final.out',                                emit: log
    path '*.bam', optional: true,                        emit: bam

    script:
    def no_wl = whitelist.name == 'NO_WHITELIST'
    // --soloCBmatchWLtype and the CellRanger-flavoured UMI options are only
    // valid with a whitelist. With --soloCBwhitelist None they are an error.
    def wl_args = no_wl
        ? "--soloCBwhitelist None --soloUMIdedup 1MM_All --soloUMIfiltering MultiGeneUMI"
        : "--soloCBwhitelist ${whitelist} --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIdedup 1MM_CR --soloUMIfiltering MultiGeneUMI_CR"
    def cb  = params.cb_len as int
    def umi = params.umi_len as int
    def bam_args = params.save_bam
        ? "--outSAMtype BAM SortedByCoordinate --outSAMattributes NH HI nM AS CR UR CB UB GX GN"
        : "--outSAMtype None"
    """
    set -euo pipefail
    EXTRA=""
    if [ "${params.save_bam}" = "true" ]; then
      # Sorting opens runThreadN x outBAMsortingBinsN files at once; 32 threads
      # x 50 bins exceeds a 1,024 soft limit and STAR dies partway through.
      need=\$(( ${task.cpus} * 50 + 256 ))
      cur=\$(ulimit -n); hard=\$(ulimit -Hn)
      if [ "\$cur" != "unlimited" ] && [ "\$cur" -lt "\$need" ]; then
        if [ "\$hard" = "unlimited" ] || [ "\$hard" -ge "\$need" ]; then
          ulimit -n "\$need"
        else
          bins=\$(( (hard - 256) / ${task.cpus} )); [ "\$bins" -lt 1 ] && bins=1
          EXTRA="--outBAMsortingBinsN \$bins"
        fi
      fi
    fi

    # readFilesIn takes cDNA FIRST, then the barcode read.
    STAR --runThreadN ${task.cpus} \\
         --genomeDir ${index} \\
         --readFilesIn ${r2} ${r1} \\
         --readFilesCommand zcat \\
         --soloType CB_UMI_Simple \\
         --soloCBstart 1 --soloCBlen ${cb} \\
         --soloUMIstart ${cb + 1} --soloUMIlen ${umi} \\
         --soloBarcodeReadLength 0 \\
         --soloStrand ${params.solo_strand} \\
         ${wl_args} \\
         --soloFeatures ${params.solo_features} \\
         --soloCellFilter ${params.solo_cell_filter} \\
         --soloMultiMappers Unique EM \\
         ${bam_args} \$EXTRA \\
         ${task.ext.args ?: ''}

    if [ -f Aligned.sortedByCoord.out.bam ]; then mv Aligned.sortedByCoord.out.bam ${sample}.bam; fi
    cp Solo.out/${params.solo_features.tokenize(' ')[0]}/Summary.csv ${sample}.star_summary.csv
    rm -rf _STARtmp
    """
}
