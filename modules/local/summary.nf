process SUMMARY {
    label 'process_single'
    publishDir "${params.outdir}", mode: params.publish_dir_mode

    input:
    path viral_rows
    path star_summaries

    output:
    path 'summary.tsv', emit: tsv

    script:
    """
    collect_summary.py --viral ${viral_rows} --star ${star_summaries} --out summary.tsv
    """
}
