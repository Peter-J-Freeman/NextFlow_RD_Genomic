/*
 * Mark duplicates in BAM files
 */
process markDuplicates {

    container = 'variantvalidator/indexgenome:1.1.0'

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
