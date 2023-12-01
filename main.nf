// Use newest nextflow dsl
nextflow.enable.dsl = 2

// Print pipeline configuration
log.info """\
    ============================================
            Simple DNASeq Pipeline Configuration
    ============================================
    reads           : ${params.reads}
    samplesheet     : ${params.samplesheet}
    genome          : ${params.genome_file}
    genome index    : ${params.genome_index_files}
    index directory : ${params.indexDir}
    output directory: ${params.outdir}
    results directory: ${params.resDir}
    fractions       : ${params.fractions}
""".stripIndent()

// Include modules in sequence
include { indexGenome } from './modules/indexGenome'
include { FASTQC } from './modules/FASTQC'
include { alignReads } from './modules/alignReads'
include { downsampleBam } from './modules/downsampleBam'
include { sortBam } from './modules/sortBam'
include { markDuplicates } from './modules/markDuplicates'
include { indexBam } from './modules/indexBam'
include { haplotypeCaller } from './modules/haplotypeCaller'
include { filterVCF } from './modules/filterVCF'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    // Set up channel for the genome fasta file including indexes if already created
    genome_ch = file(params.genome_file)

    // Set up a channel for the genome index files including indexes if already created. If empty, send a placeholder
    //filepath to trigger the next process (using the ifEmpty syntax)
    genome_index_ch = Channel.fromPath(params.genome_index_files).collect().ifEmpty {
        // Provide a filepath to a placeholder file when the channel is empty
        file("$params.indexDir/placeholder.fasta.phd")
    }

    // Invoke indexGenome process with the genome_file parameter and an alternative empty channel
    indexed_genome_ch = indexGenome(genome_ch, genome_index_ch)

    // Set channel to gather read_pairs
    // read_pairs_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)
    // read_pairs_ch.view()


    // Set channel to gather read_pairs
    read_pairs_ch = Channel
        .fromPath(params.samplesheet)
        .splitCsv(sep: '\t')
        .map { row -> tuple(row[0], [row[1], row[2]]) }

    read_pairs_ch.view()

    // Run FASTQC on read pairs
    FASTQC(read_pairs_ch)

    // Align reads to the indexed genome
    align_ch = alignReads(read_pairs_ch)

    // Downsample BAM files
    downsample_ch = downsampleBam(align_ch)

    // Sort BAM files
    sort_ch = sortBam(downsample_ch.collect())

    // Mark duplicates in BAM files
    mark_ch = markDuplicates(sort_ch.collect())

    // Index the BAM files and collect the output channel
    indexed_bam_ch = indexBam(mark_ch.collect())

    // Call Variants
    vcf_call_ch = haplotypeCaller(mark_ch.collect(), genome_ch, Channel.fromPath(params.genome_index_files).collect(), indexed_bam_ch.collect())

    // Filter Variants
    filterVCF(vcf_call_ch.collect(), genome_ch, Channel.fromPath(params.genome_index_files).collect())
}

workflow.onComplete {
    log.info ( workflow.success ? "\nworkflow is done!\n" : "Oops .. something went wrong" )
}