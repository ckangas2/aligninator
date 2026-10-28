// Cheap structural check before alignment. A swapped R1/R2 or a barcode read
// shorter than CB+UMI does not make STAR fail; it makes ~0 valid barcodes.
process CHECK_READS {
    tag "${sample}"
    label 'process_single'
    publishDir path: { "${params.outdir}/samples/${sample}" }, mode: params.publish_dir_mode,
        pattern: '*.log'

    input:
    tuple val(sample), path(r1), path(r2)

    output:
    tuple val(sample), path(r1), path(r2), emit: reads
    path "${sample}.read_check.log",       emit: log

    script:
    def need = (params.cb_len as int) + (params.umi_len as int)
    """
    # No pipefail here: head closing the pipe early would SIGPIPE zcat.
    set -u
    mode_len() { zcat -f "\$1" | head -40000 | awk 'NR%4==2 {print length(\$0)}' \\
                 | sort -n | uniq -c | sort -rn | head -1 | awk '{print \$2}'; }
    l1=\$(mode_len ${r1}); l2=\$(mode_len ${r2})
    {
      echo "sample              : ${sample}"
      echo "fastq_1 (barcode)   : ${r1.name}  mode length \$l1 bp"
      echo "fastq_2 (cDNA)      : ${r2.name}  mode length \$l2 bp"
      echo "CB + UMI configured : ${params.cb_len} + ${params.umi_len} = ${need} bp"
      echo "sjdbOverhang in use : ${params.sjdb_overhang}   (ideal for this cDNA read: \$((l2 - 1)))"
    } | tee ${sample}.read_check.log

    if [ -z "\$l1" ] || [ -z "\$l2" ]; then
      echo "FATAL: could not read sequences from one of the FASTQs." >&2; exit 1
    fi
    if [ "\$l1" -gt "\$l2" ]; then
      echo "FATAL: fastq_1 is longer than fastq_2. fastq_1 must be the barcode read" >&2
      echo "and fastq_2 the cDNA read; these look swapped in the samplesheet." >&2
      exit 1
    fi
    if [ "\$l1" -lt ${need} ]; then
      echo "FATAL: barcode read is \$l1 bp but CB+UMI needs ${need} bp." >&2
      echo "Wrong chemistry, or the barcode read was trimmed." >&2
      exit 1
    fi
    """
}
