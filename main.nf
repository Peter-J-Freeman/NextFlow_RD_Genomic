// Disable buffering for stdout and stderr
System.setProperty("org.gradle.native", "true")
System.setProperty("org.gradle.console", "plain")

// Define project parameters
params.reads = "$projectDir/data/samples/NA12878_WES_R{1,2}.fastq"
params.genome_file = "$projectDir/data/genome/hg38_v0_Homo_sapiens_assembly38.fasta"
params.genome_index_files = "$projectDir/data/genome/*.fasta.*"
params.indexDir = "$projectDir/data/genome"
params.multiqc = "$projectDir/multiqc"
params.outdir = "$projectDir/files"
params.resDir = "$projectDir/results"
params.logDir = "$projectDir/logs"
params.fractions = [0.001, 0.005, 0.01, 0.05, 0.1, 0.5]

// Set up log directory
logDir params.logDir

// Print pipeline configuration
log.info """\
    Simple DNASeq Pipeline Without BQSR and VQSR
    ============================================
    genome          : ${params.genome_file}
    reads           : ${params.reads}
    downsample      : ${params.fractions}
    output_directory: ${params.outdir}
    multiqc         : ${params.multiqc}
    results         : ${params.resDir}
""".stripIndent()

/*
 * Define the indexGenome process that creates a BWA index
 * given the genome fasta file
 */
process indexGenome {

    // Publish indexed files to the specified directory
    publishDir(params.indexDir, mode: "copy")

    input:
    path genomeFasta
    path indexFiles

    output:
    path "${genomeFasta}.amb"
    path "${genomeFasta}.ann"
    path "${genomeFasta}.bwt"
    path "${genomeFasta}.pac"
    path "${genomeFasta}.sa"
    path "${genomeFasta}.fai"
    path "${genomeFasta}.dict"

    script:
    """
    echo "Running Index Genome"
    set -euo pipefail

    # Use the 'genomeFasta' variable to refer to the genome file path
    genomeFasta="${genomeFasta}"

    # Use the 'indexFiles' variable to refer to the list of existing index files
    indexFiles=(${indexFiles})

    requiredIndexFiles=("\${genomeFasta}.amb" "\${genomeFasta}.ann" "\${genomeFasta}.bwt" "\${genomeFasta}.pac"
                        "\${genomeFasta}.sa" "\${genomeFasta}.fai" "\${genomeFasta}.dict")

    echo "Required index files:"
    echo "\${requiredIndexFiles[@]}"

    echo "Existing index files:"
    ls -l "\${indexFiles[@]}"

    echo "Checking for the presence of required index files in the indexFiles directory..."
    missingFiles=()
    for file in "\${requiredIndexFiles[@]}"; do
        found=false
        for idxFile in "\${indexFiles[@]}"; do
            if [ "\$idxFile" == "\$file" ]; then
                found=true
                break
            fi
        done

        if [ "\$found" == false ]; then
            missingFiles+=("\$file")
        fi
    done

    if [ \${#missingFiles[@]} -gt 0 ]; then
        echo "Some or all of the required index files are missing. Generating BWA index and samtools faidx..."

        # Generate BWA index
        bwa index "\${genomeFasta}"

        # Generate samtools faidx
        samtools faidx "\${genomeFasta}"

        # Generate Fasta dict
        picard CreateSequenceDictionary R="\${genomeFasta}" O="\${genomeFasta}.dict"
    fi

    echo "Genome Indexing complete."
    """
}

/*
 * Run fastq on the read fastq files
 */
process FASTQC {
    // Add a tag to identify the process
    tag "FASTQC on $sample_id"

    // Specify the output directory for the FASTQC results
    publishDir(params.resDir, mode: "copy", mkdirs: true)

    input:
    tuple val(sample_id), path(reads)

    output:
    path "fastqc_${sample_id}_logs/*"

    script:
    """
    echo "Running FASTQC"
    mkdir -p fastqc_${sample_id}_logs
    fastqc ${reads[0]} ${reads[1]} -o fastqc_${sample_id}_logs
    echo "FASTQC Complete"
    """
}

/*
 * Align reads to the indexed genome
 */
process alignReads {

    // Input channels: sample_id
    input:
    tuple val(sample_id), path(reads)

    output:
    file("${sample_id}.bam")

    script:
    """
    echo "Running Align Reads"
    # Use the 'genomeFasta' variable to refer to the genome file path
    genomeFasta="${params.genome_file}"

    requiredIndexFiles=("\${genomeFasta}.amb" "\${genomeFasta}.ann" "\${genomeFasta}.bwt" "\${genomeFasta}.pac"
                        "\${genomeFasta}.sa" "\${genomeFasta}.fai" "\${genomeFasta}.dict")

    echo "Required index files:"
    echo "\${requiredIndexFiles[@]}"

    # Wait until all the required index files are present
    while true; do
        found=true
        for file in "\${requiredIndexFiles[@]}"; do
            if [ ! -f "\$file" ]; then
                found=false
                break
            fi
        done

        if [ "\$found" == true ]; then
            break
        fi

        echo "Some of the required index files are missing. Waiting for 30 sec..."
        sleep 30
    done

    echo "All required index files are present. Proceed with the alignment."

    # Align reads using BWA and generate BAM file
    bwa mem -t ${task.cpus} "\${genomeFasta}" ${reads[0]} ${reads[1]} |
    samtools view -b - |
    samtools addreplacerg -r "@RG\\tID:${sample_id}_1.fastq\\tSM:${sample_id}_2.fastq\\tPL:illumina" - |
    samtools addreplacerg -r "@RG\\tID:${sample_id}_2.fastq\\tSM:${sample_id}_1.fastq\\tPL:illumina" - > ${sample_id}.bam
    echo "Alignment complete"
    """
}

/*
 * Downsample the BAM files
 */
process downsampleBam {

    input:
    file bamFile

    output:
    tuple val(bamFile.baseName), file("*_downsampled_*.bam")

    script:
    """
    echo "Running Downsampling"
    fractions=(${params.fractions.join(' ')})
    for fraction in "\${fractions[@]}"; do
        outputBam="${bamFile.baseName}_downsampled_\${fraction//./}.bam"
        picard DownsampleSam I=${bamFile} O="\${outputBam}" P="\${fraction}"
        mv "\${outputBam}" "\${outputBam}"
    done
    mv ${bamFile} ${bamFile.baseName}_downsampled_100.bam
    echo "Downsampling complete"
    """
}

/*
 * Sort the BAM files
 */
process sortBam {

    input:
    tuple val(sample_id), file(bamFiles)

    output:
    tuple val(sample_id), file("${sample_id}_*_sorted.bam")

    script:
    """
    echo "Running Sort Bam"
    # Define an array to store the input BAM files
    inputBams=(${bamFiles})

    # Sort each input BAM file and generate sorted BAM files
    for inputBam in "\${inputBams[@]}"; do
        outputBam="\$(basename \${inputBam} .bam)_sorted.bam"

        # Use samtools to sort the input BAM file and save the output to the sorted BAM file
        samtools sort -o \${outputBam} \${inputBam}
        echo "\${outputBam}"
    done
    echo "BAM Sorting complete"
    """
}

/*
 * Mark duplicates in BAM files
 */
process markDuplicates {

    // Publish deduplicated BAM files to the specified directory
    publishDir(params.outdir, mode: "copy")

    input:
    tuple val(sample_id), file(bamFiles)

    output:
    tuple val(sample_id), file("${sample_id}_downsampled_*_dedup.bam")

    script:
    """
    echo "Running Mark Duplicates"
    # Define an array to store the input BAM files
    inputBams=(${bamFiles})

    # Loop through each input BAM file and mark duplicates
    for inputBam in "\${inputBams[@]}"; do
        outputBam="\$(basename \${inputBam} _sorted.bam)_dedup.bam"
        metricsFile="\$(basename \${inputBam} _sorted.bam)_dedup_metrics.txt"

        # Use Picard tools to mark duplicates in the input BAM file
        picard MarkDuplicates I=\${inputBam} O="\${outputBam}" M="\${metricsFile}"

        echo "\${outputBam}"
    done
    echo "Mark Duplicates Complete"
    """
}

/*
 * Index the BAM files
 */
process indexBam {

    // Publish indexed BAM files to the specified directory
    publishDir(params.outdir, mode: "copy")

    input:
    tuple val(sample_id), file(bamFiles)

    output:
    tuple val(sample_id), file("${sample_id}_downsampled_*_dedup.bam.bai")

    script:
    """
    echo "Running Index Bam"
    # Define an array to store the input BAM files
    inputBams=(${bamFiles})

    # Loop through each input BAM file and create the index
    for inputBam in "\${inputBams[@]}"; do
        outputBai="\$(basename \${inputBam} _dedup.bam)_dedup.bam.bai"

        # Use samtools to create the index for the input BAM file
        samtools index \${inputBam} \${outputBai}

    done
    echo "Bam Indexing Complete"
    """
}

/*
 * Call variants using HaplotypeCaller
 */
process haplotypeCaller {

    input:
    tuple val(sample_id), file(bamFiles)
    path genomeFasta
    path indexFiles
    tuple val(sample_id), file(bamIndexFiles)

    output:
    tuple val(sample_id), file("${sample_id}_downsampled_*.vcf")

    script:
    """
    echo "Running HaplotypeCaller for Sample: ${sample_id}"

    # Use the 'genomeFasta' variable passed as an input
    genomeFasta="${genomeFasta}"

    # Define an array to store the input BAM files
    inputBams=(${bamFiles})

    # Rename the dictionary file to the expected name
    mv "\${genomeFasta}.dict" "\${genomeFasta%.*}.dict"

    # Loop through each input BAM file and call variants
    echo "Calling VCFs for Sample: ${sample_id}..."
    for inputBam in "\${inputBams[@]}"; do
        outputVcf="\$(basename \${inputBam} _dedup.bam).vcf"

        # Use GATK HaplotypeCaller to call variants
        gatk HaplotypeCaller -R \${genomeFasta} -I \${inputBam} -O "\${outputVcf}"

        echo "Sample: ${sample_id} VCF: \${outputVcf}"
    done

    # Move the dictionary file back to its original name
    mv "\${genomeFasta%.*}.dict" "\${genomeFasta}.dict"

    echo "Variant Calling for Sample: ${sample_id} Complete"
    """
}

/*
 * Filter the raw VCF files using GATK VariantFiltration
 */
process filterVCF {

    // Publish VCF files to the specified directory
    publishDir(params.resDir, mode: "copy")

    input:
    tuple val(sample_id), file(vcfFiles)
    path genomeFasta
    path indexFiles

    // Output channel for sample_id and filtered VCF files
    output:
    tuple val(sample_id), file("${sample_id}_downsampled_*.filtered.vcf")

    // Script section to run the process
    script:
    """
    # Print a message indicating the start of the process for the current sample
    echo "Running Variant Filtration for Sample: ${sample_id}"

    # Use the 'genomeFasta' variable passed as an input
    genomeFasta="${genomeFasta}"

    # Define an array to store the input VCF files
    inputVcfs=(${vcfFiles})

    # Rename the dictionary file to the expected name
    mv "\${genomeFasta}.dict" "\${genomeFasta%.*}.dict"

    # Loop through each input VCF file and perform variant filtration
    for inputVcf in "\${inputVcfs[@]}"; do
        outputVcf="\$(basename \${inputVcf} .vcf).filtered.vcf"

        # Use GATK VariantFiltration to filter the input VCF file
        gatk VariantFiltration -R "\${genomeFasta}" -V "\${inputVcf}" -O "\${outputVcf}"

        # Print a message indicating the completion of variant filtration for the current sample
        echo "Sample: ${sample_id} VCF: \${outputVcf}"
    done

    # Move the dictionary file back to its original name
    mv "\${genomeFasta%.*}.dict" "\${genomeFasta}.dict"

    # Print a message indicating the completion of variant filtration for the current sample
    echo "Variant Filtering for Sample: ${sample_id} Complete"
    """
}

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
