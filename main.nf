// Disable buffering for stdout and stderr
System.setProperty("org.gradle.native", "true")
System.setProperty("org.gradle.console", "plain")

// Define project parameters
params.reads = "$projectDir/data/samples/SRR1518253_{1,2}.fastq"
params.genome_file = "$projectDir/data/genome/Homo_sapiens_assembly38.fasta"
params.genome_index_files = "$projectDir/data/genome/*.fasta.*"
params.indexDir = "$projectDir/data/genome"
params.outdir = "$projectDir/files"
params.resDir = "$projectDir/results"
params.fractions = [0.001, 0.005, 0.01, 0.05, 0.1, 0.5]

// Print pipeline configuration
log.info """\
    Simple DNASeq Pipeline Without BQSR and VQSR
    ============================================
    genome          : ${params.genome_file}
    reads           : ${params.reads}
    downsample      : ${params.fractions}
    output_directory: ${params.outdir}
    results         : ${params.resDir}
""".stripIndent()

// Include modules in sequence
include { indexGenome } from './modules/indexGenome.nf'
include { FASTQC } from './modules/FASTQC.nf'
include { alignReads } from './modules/alignReads.nf'
include { downsampleBam } from './modules/downsampleBam.nf'
include { sortBam } from './modules/sortBam.nf'
include { markDuplicates } from './modules/markDuplicates.nf'
include { indexBam } from './modules/indexBam.nf'
include { haplotypeCaller } from './modules/haplotypeCaller.nf'
include { filterVCF } from './modules/filterVCF.nf'

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
    read_pairs_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)

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
