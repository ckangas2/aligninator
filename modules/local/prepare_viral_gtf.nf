process PREPARE_VIRAL_GTF {
    tag "${gtf.name}"
    label 'process_single'
    publishDir "${params.outdir}/reference", mode: params.publish_dir_mode

    input:
    path gtf

    output:
    path 'viral.normalised.gtf', emit: gtf
    path 'viral_genes.txt',      emit: genes
    path 'prepare_viral_gtf.log', emit: log

    script:
    def rename = (params.viral_contig_rename && params.viral_contig_rename != true) ? "--rename-contig '${params.viral_contig_rename}'" : ''
    """
    prepare_viral_gtf.py \\
        --in-gtf ${gtf} \\
        --out-gtf viral.normalised.gtf \\
        --gene-list-out viral_genes.txt \\
        ${rename} \\
        --prefix '${params.viral_prefix}' \\
        --add-gene-rows --copy-id-to-name \\
        | tee prepare_viral_gtf.log
    """
}
