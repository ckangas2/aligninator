process SOLO_TO_H5AD {
    tag "${sample}:${feature}"
    label 'process_low'
    publishDir path: { "${params.outdir}/samples/${sample}" }, mode: params.publish_dir_mode

    input:
    tuple val(sample), val(feature), path(filtered, stageAs: 'filtered')
    path viral_genes

    output:
    path "${sample}.${feature}.{h5ad,csv}",      emit: matrix
    path "${sample}.${feature}.summary.tsv",     emit: summary
    path "${sample}.${feature}.viral.log",       emit: log

    script:
    // Viral features come from viral_genes.txt written alongside the index,
    // never from a name prefix. --strict fails if the two disagree.
    """
    set -o pipefail
    solo_to_h5ad.py \\
        --solo-dir filtered \\
        --viral-genes-file ${viral_genes} \\
        --sample ${sample} \\
        --feature-set ${feature} \\
        --exclude-from-viral '${params.exclude_from_viral}' \\
        --summary-out ${sample}.${feature}.summary.tsv \\
        --strict \\
        --out ${sample}.${feature}.h5ad \\
        | tee ${sample}.${feature}.viral.log
    """
}
