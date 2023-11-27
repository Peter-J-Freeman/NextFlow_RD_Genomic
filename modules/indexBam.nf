/*
 * Index the BAM files
 */
process indexBam {

    container = 'variantvalidator/indexgenome:1.1.0'

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
