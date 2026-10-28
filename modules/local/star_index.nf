process STAR_INDEX {
    tag "host+virus"
    label 'process_star_index'
    publishDir "${params.outdir}/reference", mode: params.publish_dir_mode,
        saveAs: { f -> (f == 'star_index' && !params.save_reference) ? null : f }

    // Inputs are staged into their own folders. They arrive as symlinks to
    // the user's real files, so an output that shares an input's name (the
    // 10x reference ships as genome.fa / genes.gtf) would be written THROUGH
    // the link and destroy the original. Separate folders make that
    // impossible, and noclobber below refuses to overwrite anything.
    input:
    path host_fasta,  stageAs: 'in/host/*'
    path host_gtf,    stageAs: 'in/host_gtf/*'
    path viral_fasta, stageAs: 'in/viral/*'
    path viral_gtf,   stageAs: 'in/viral_gtf/*'
    val  validated      // ordering only: never index an unvalidated reference

    output:
    path 'star_index',            emit: index
    path 'star_index_check.log',  emit: log

    script:
    // genomeSAindexNbases must shrink for small genomes or STAR segfaults /
    // hangs; min(14, log2(L)/2 - 1) is the formula from the STAR manual.
    """
    set -euo pipefail
    set -o noclobber    # `>` fails rather than overwrite an existing file

    # Host first so contig order matches the host-only reference people compare to.
    zcat -f ${host_fasta} ${viral_fasta} > merged.genome.fa
    zcat -f ${host_gtf} > merged.genes.gtf
    cat ${viral_gtf} >> merged.genes.gtf

    # A host FASTA with no contigs (empty, truncated, or the wrong file)
    # would otherwise index the virus alone and look like success.
    n_host=\$(zcat -f ${host_fasta} | grep -c '^>' || true)
    if [ "\$n_host" -eq 0 ]; then
      echo "FATAL: host FASTA ${host_fasta} contains no contigs." >&2; exit 1
    fi

    n_contigs=\$(grep -c '^>' merged.genome.fa)
    glen=\$(grep -v '^>' merged.genome.fa | tr -d '\\n' | wc -c)
    sa=\$(python3 -c "import math; print(min(14, int(math.log2(\$glen)/2 - 1)))")

    mkdir star_index
    STAR --runMode genomeGenerate \\
         --runThreadN ${task.cpus} \\
         --genomeDir star_index \\
         --genomeFastaFiles merged.genome.fa \\
         --sjdbGTFfile merged.genes.gtf \\
         --sjdbOverhang ${params.sjdb_overhang} \\
         --genomeSAindexNbases \$sa \\
         --limitGenomeGenerateRAM ${task.memory.toBytes()} \\
         --outFileNamePrefix star_index/

    rm -f merged.genome.fa

    # If this is 0 the reference is broken regardless of what STAR printed.
    n_viral=\$(grep -c "^${params.viral_prefix}" star_index/geneInfo.tab || true)
    {
      echo "host contigs            : \$n_host"
      echo "contigs in merged FASTA : \$n_contigs"
      echo "genome length           : \$glen"
      echo "genomeSAindexNbases     : \$sa"
      echo "sjdbOverhang            : ${params.sjdb_overhang}"
      echo "viral genes in index    : \$n_viral"
    } | tee star_index_check.log
    if [ "\$n_viral" -eq 0 ]; then
      echo "FATAL: no viral genes in the index. Do not align against this." >&2
      exit 1
    fi
    """
}
