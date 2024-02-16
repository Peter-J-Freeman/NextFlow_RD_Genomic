// Use newest nextflow dsl
nextflow.enable.dsl = 2

// Print pipeline configuration
//removed indexDir from params
log.info """\
    ============================================
            Simple DNASeq Pipeline Configuration
    ============================================
    samplesheet     : ${params.samplesheet}
    genome          : ${params.genome_file}
    genome index    : ${params.genome_index_files}
    index_genome    : ${params.index_genome}
    output directory: ${params.outdir}
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
    // genome_ch = file(params.genome_file)
    // can just use params.genome_file directly here as it will auto create a channel for it if not declared

    //SCADD - create a fractions channel
    fractions_ch = Channel.fromList(params.fractions)

    //SCADD - user decides to index genome or not
    if (params.index_genome){
        //flatten as is of format [fasta, [rest of files..]]
        indexed_genome_ch = indexGenome(params.genome_file).flatten()
        indexed_genome_ch.view()
    }
    else {
        indexed_genome_ch = Channel.fromPath(params.genome_index_files)
        indexed_genome_ch.view()
    }


    // Set channel to gather read_pairs
    read_pairs_ch = Channel
        .fromPath(params.samplesheet)
        .splitCsv(sep: '\t')
        .map { row -> tuple(row[0], [row[1], row[2]]) }

    read_pairs_ch.view()

    // Run FASTQC on read pairs
    FASTQC(read_pairs_ch)

    // Align reads to the indexed genome
    // SCADD- added index_genpme_ch channel to make it wait for index if one not already created
    // read_pairs_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)
    align_ch = alignReads(read_pairs_ch, indexed_genome_ch.collect())

    // Downsample BAM files
    // SCADD- first combine fractions with bam channel https://www.nextflow.io/docs/latest/operator.html#combine
    // this means every downsample process will get its own instance so downsampling can run in parallel
    bam_fracs_ch = fractions_ch.combine(align_ch)
    downsample_ch = downsampleBam(bam_fracs_ch)

    // Sort BAM files
    // SCADD - removed collect to run sorting in parallel
    sort_ch = sortBam(downsample_ch)

    // Mark duplicates in BAM files
    // SCADD - removed collect to run sorting in parallel
    mark_ch = markDuplicates(sort_ch)

    // Index the BAM files and collect the output channel
    // SCADD - removed collect to run sorting in parallel
    // changed indexbam to output sample_id, bam, bai file in tuple so that correct bai file stays with bamfile
    indexed_bam_ch = indexBam(mark_ch)

    // Call Variants
    // SCADD- removed collect from vcf bam ch
    // altered haplotype caller to take in the index bam ch
    // replaced params.genome_index_files with indexed_genome_ch
    // vcf_call_ch = haplotypeCaller(indexed_bam_ch, params.genome_file, indexed_genome_ch)
    vcf_call_ch = haplotypeCaller(indexed_bam_ch, indexed_genome_ch.collect())
    
    // Filter Variants
    // SCADD - removed collect from vcf channel
    // replaced params.genome_index_files with indexed_genome_ch
    // filterVCF(vcf_call_ch, params.genome_file, indexed_genome_ch)
    filterVCF(vcf_call_ch, indexed_genome_ch.collect())
}

workflow.onComplete {
    log.info ( workflow.success ? "\nworkflow is done!\n" : "Oops .. something went wrong" )
}