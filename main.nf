#!/usr/bin/env nextflow
/*
 * aligninator — host + virus single-cell alignment with STARsolo.
 *
 *   reference:  normalise viral GTF -> validate -> concatenate -> STAR index
 *               -> confirm viral genes survived into the index
 *   per sample: check read structure -> STARsolo -> matrix + per-cell viral_frac
 *   run:        one summary.tsv across every sample and feature set
 *
 * Every check that the shell pipeline makes is kept. Each exists because the
 * failure it guards against produces a successful run and a wrong matrix.
 */
nextflow.enable.dsl = 2

include { PREPARE_VIRAL_GTF  } from './modules/local/prepare_viral_gtf'
include { VALIDATE_REFERENCE } from './modules/local/validate_reference'
include { STAR_INDEX         } from './modules/local/star_index'
include { CHECK_READS        } from './modules/local/check_reads'
include { BARCODE_STRUCTURE  } from './modules/local/barcode_structure'
include { STARSOLO           } from './modules/local/starsolo'
include { SOLO_TO_H5AD       } from './modules/local/solo_to_h5ad'
include { SUMMARY            } from './modules/local/summary'

def helpMessage() {
    log.info """
    aligninator v${workflow.manifest.version}

    Usage:
      nextflow run ckangas2/aligninator -profile test,docker
      nextflow run ckangas2/aligninator -profile docker --input samplesheet.csv \\
          --host_fasta genome.fa --host_gtf genes.gtf \\
          --viral_fasta virus.fa --viral_gtf virus.gtf \\
          --cb_len 12 --umi_len 8 --whitelist None

    Input
      --input              CSV with columns sample,fastq_1,fastq_2
                           fastq_1 = barcode read (CB+UMI), fastq_2 = cDNA read
                           relative paths resolve against the samplesheet's folder

    Reference (build from scratch)
      --host_fasta, --host_gtf
      --viral_fasta, --viral_gtf
      --viral_contig_rename  OLD:NEW if the viral GTF contig name differs from the FASTA
      --viral_prefix         prefix for viral gene names (default: ${params.viral_prefix})
      --sjdb_overhang        cDNA read length - 1 (default: ${params.sjdb_overhang})
    Reference (reuse)
      --star_index           prebuilt index directory
      --viral_genes          viral_genes.txt written when that index was built

    Chemistry (explicit; nothing is guessed)
      --cb_len, --umi_len    leave --cb_len unset to run barcode detection and stop
      --barcode_top_k        detection only: roughly the expected cell count (default: ${params.barcode_top_k})
      --whitelist            path to a barcode whitelist, or None (Drop-seq, Seq-Well)
      --solo_strand          Forward | Reverse | Unstranded (default: ${params.solo_strand})
      --solo_cell_filter     STARsolo --soloCellFilter (default: ${params.solo_cell_filter})

    Output
      --outdir               (default: ${params.outdir})
      --save_bam             write coordinate-sorted BAMs (default: ${params.save_bam})
      --save_reference       publish the STAR index (default: ${params.save_reference})
    """.stripIndent()
}

// ---------------------------------------------------------------- params
// On the command line `--cb_len null` arrives as the string "null"; treat it,
// false and '' as unset rather than letting them reach STAR.
def isSet(v) {
    return v != null && !(v instanceof Boolean) && !(v.toString().trim() in ['', 'null', 'false', 'true'])
}

// Collects every parameter problem, reports them all, then stops.
def validateParams() {
    def problems = []
    if (!params.input) problems << "--input is required (CSV: sample,fastq_1,fastq_2)"
    if (params.star_index) {
        if (!params.viral_genes)
            problems << "--star_index needs --viral_genes: the viral feature list must come from the same build as the index"
    } else {
        ['host_fasta', 'host_gtf', 'viral_fasta', 'viral_gtf'].each { p ->
            if (!params[p]) problems << "--${p} is required when --star_index is not given"
        }
    }
    if (isSet(params.cb_len)) {
        if (!(params.cb_len.toString() ==~ /\d+/))  problems << "--cb_len must be a whole number of bases, got '${params.cb_len}'"
        if (!(params.umi_len.toString() ==~ /\d+/)) problems << "--umi_len is required alongside --cb_len, as a whole number of bases"
    }
    if (params.whitelist == null) problems << "--whitelist is required: a file path, or None for whitelist-free chemistries"
    if (isSet(params.viral_contig_rename) && !(params.viral_contig_rename.toString() ==~ /[^:]+:[^:]+/))
        problems << "--viral_contig_rename must be OLD:NEW"
    if (problems) {
        problems.each { msg -> log.error msg }
        error "Fix the parameters above (see --help)."
    }
}

def resolvePath(base, String p) {
    return (p.contains('://') || p.startsWith('/')) ? p : "${base}/${p}"
}

// Relative FASTQ paths resolve against the samplesheet's own directory, so a
// samplesheet and its FASTQs can move together.
def readSamplesheet(path) {
    def seen = [] as Set
    def base = file(path, checkIfExists: true).parent
    channel.fromPath(path, checkIfExists: true)
        .splitCsv(header: true, strip: true)
        .map { row ->
            ['sample', 'fastq_1', 'fastq_2'].each { k ->
                if (!row[k]) error "samplesheet row missing '${k}': ${row}"
            }
            if (!(row.sample ==~ /[A-Za-z0-9._-]+/))
                error "sample name '${row.sample}' may only contain letters, digits, '.', '_' and '-'"
            if (!seen.add(row.sample))
                error "duplicate sample name in samplesheet: ${row.sample}"
            tuple(row.sample,
                  file(resolvePath(base, row.fastq_1), checkIfExists: true),
                  file(resolvePath(base, row.fastq_2), checkIfExists: true))
        }
}

// ---------------------------------------------------------------- workflow
workflow {

    if (params.help) {
        helpMessage()
    } else {
        validateParams()
        def chemistry_set = isSet(params.cb_len)
        def reads = readSamplesheet(params.input)

        if (!chemistry_set) {
            // Chemistry is never inferred and then used silently. With no
            // --cb_len, run the whitelist-free detector on every sample,
            // report, and stop.
            BARCODE_STRUCTURE(reads)
            BARCODE_STRUCTURE.out.log
                .collect()
                .subscribe {
                    log.warn """
                    --cb_len was not set. Barcode structure was detected and written to
                    ${params.outdir}/barcode_structure/. Read the verdict, then re-run with
                    --cb_len, --umi_len and --whitelist set explicitly. The pipeline does not
                    align on an auto-detected value: chemistry belongs in the command line
                    where it is visible and reviewable.
                    """.stripIndent()
                }
        } else {
            align(reads)
        }
    }
}

workflow align {
    take:
    reads

    main:
    // ------------------------------------------------ reference
    def index       = null
    def viral_genes = null
    if (params.star_index) {
        index       = channel.value(file(params.star_index, checkIfExists: true))
        viral_genes = channel.value(file(params.viral_genes, checkIfExists: true))
    } else {
        def host_fasta  = file(params.host_fasta,  checkIfExists: true)
        def host_gtf    = file(params.host_gtf,    checkIfExists: true)
        def viral_fasta = file(params.viral_fasta, checkIfExists: true)
        def viral_gtf   = file(params.viral_gtf,   checkIfExists: true)

        PREPARE_VIRAL_GTF(viral_gtf)
        VALIDATE_REFERENCE(viral_fasta, PREPARE_VIRAL_GTF.out.gtf, host_gtf)
        STAR_INDEX(host_fasta, host_gtf, viral_fasta,
                   PREPARE_VIRAL_GTF.out.gtf, VALIDATE_REFERENCE.out.ok)
        index       = STAR_INDEX.out.index
        viral_genes = PREPARE_VIRAL_GTF.out.genes
    }

    def whitelist = (params.whitelist.toString() == 'None')
        ? file("${projectDir}/assets/NO_WHITELIST")
        : file(params.whitelist, checkIfExists: true)

    // ------------------------------------------------ per sample
    CHECK_READS(reads)
    STARSOLO(CHECK_READS.out.reads, index, whitelist)

    def matrices = STARSOLO.out.solo
        .flatMap { sample, solo ->
            params.solo_features.tokenize(' ').collect { feat ->
                tuple(sample, feat, file("${solo}/${feat}/filtered"))
            }
        }
    SOLO_TO_H5AD(matrices, viral_genes)

    SUMMARY(SOLO_TO_H5AD.out.summary.collect(),
            STARSOLO.out.star_summary.map { _sample, f -> f }.collect())
    SUMMARY.out.tsv.subscribe {
        log.info "aligninator: per-sample viral load written to ${params.outdir}/summary.tsv"
    }
}
