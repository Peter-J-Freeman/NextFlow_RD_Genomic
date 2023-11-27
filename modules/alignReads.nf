/*
 * Align reads to the indexed genome
 */
process alignReads {

    container = 'variantvalidator/indexgenome:1.1.0'

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
