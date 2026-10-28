// Fails the run before genomeGenerate spends an hour on a broken reference:
// contig mismatch, missing exon rows, host/virus symbol collisions.
process VALIDATE_REFERENCE {
    tag "${viral_fasta.name}"
    label 'process_single'
    publishDir "${params.outdir}/reference", mode: params.publish_dir_mode

    input:
    path viral_fasta, stageAs: 'in/viral/*'
    path viral_gtf,   stageAs: 'in/viral_gtf/*'
    path host_gtf,    stageAs: 'in/host_gtf/*'

    output:
    val true,                      emit: ok
    path 'validate_reference.log', emit: log

    script:
    """
    set -o pipefail
    validate_reference.py \\
        --viral-fasta ${viral_fasta} \\
        --viral-gtf   ${viral_gtf} \\
        --host-gtf    ${host_gtf} \\
        --viral-prefix '${params.viral_prefix}' \\
        | tee validate_reference.log
    """
}
