process STAR_INDEX {
    tag "host+virus"
    label 'process_star_index'
    publishDir "${params.outdir}/reference", mode: params.publish_dir_mode,
        saveAs: { f -> (f == 'star_index' && !params.save_reference) ? null : f }

    input:
    path host_fasta
    path host_gtf
    path viral_fasta
    path viral_gtf
    val  validated      // ordering only: never index an unvalidated reference

    output:
    path 'star_index',            emit: index
    path 'star_index_check.log',  emit: log

    script:
    // genomeSAindexNbases must shrink for small genomes or STAR segfaults /
    // hangs; min(14, log2(L)/2 - 1) is the formula from the STAR manual.
    """
    set -euo pipefail

    # Host first so contig order matches the host-only reference people compare to.
    zcat -f ${host_fasta} ${viral_fasta} > genome.fa
    zcat -f ${host_gtf} > genes.gtf
    cat ${viral_gtf} >> genes.gtf

    n_contigs=\$(grep -c '^>' genome.fa)
    glen=\$(grep -v '^>' genome.fa | tr -d '\\n' | wc -c)
    sa=\$(python3 -c "import math; print(min(14, int(math.log2(\$glen)/2 - 1)))")

    mkdir star_index
    STAR --runMode genomeGenerate \\
         --runThreadN ${task.cpus} \\
         --genomeDir star_index \\
         --genomeFastaFiles genome.fa \\
         --sjdbGTFfile genes.gtf \\
         --sjdbOverhang ${params.sjdb_overhang} \\
         --genomeSAindexNbases \$sa \\
         --limitGenomeGenerateRAM ${task.memory.toBytes()} \\
         --outFileNamePrefix star_index/

    rm -f genome.fa

    # If this is 0 the reference is broken regardless of what STAR printed.
    n_viral=\$(grep -c "^${params.viral_prefix}" star_index/geneInfo.tab || true)
    {
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
