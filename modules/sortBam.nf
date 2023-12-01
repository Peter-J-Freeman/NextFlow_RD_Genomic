/*
 * Sort the BAM files
 */
process sortBam {

    label 'process_single'
    container 'variantvalidator/indexgenome:1.1.0'

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
