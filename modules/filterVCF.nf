/*
 * Filter the raw VCF files using GATK VariantFiltration
 */
process filterVCF {

    label 'process_medium'
    container 'variantvalidator/gatk4:4.3.0.0'

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
