process BARCODE_STRUCTURE {
    tag "${sample}"
    label 'process_single'
    publishDir "${params.outdir}/barcode_structure", mode: params.publish_dir_mode

    input:
    tuple val(sample), path(r1), path(r2)

    output:
    path "${sample}.barcode_structure.log", emit: log

    script:
    """
    barcode_structure.py --r1 ${r1} --top-k ${params.barcode_top_k} | tee ${sample}.barcode_structure.log
    """
}
