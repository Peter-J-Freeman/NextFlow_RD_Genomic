/*
 * Call variants using HaplotypeCaller
 */
process haplotypeCaller {

    label 'process_medium'
    container 'variantvalidator/gatk4:4.3.0.0'

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
