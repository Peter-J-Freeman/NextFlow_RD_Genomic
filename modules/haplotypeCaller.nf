/*
 * Call variants using HaplotypeCaller
 */
process haplotypeCaller {

    label 'process_medium'
    container 'variantvalidator/gatk4:4.3.0.0'

    tag "$bamFile"

    input:
    tuple val(sample_id), file(bamFile), file(bamIndex)
    path indexFiles

    output:
    tuple val(sample_id), file("*.vcf")

    script:
    """
    echo "Running HaplotypeCaller for Sample: ${bamFile}"

    if [[ -n params.genome_file ]]; then
        genomeFasta=\$(basename ${params.genome_file})
    else
        genomeFasta=\$(find -L . -name '*.fasta')
    fi

    echo "Genome File: \${genomeFasta}"

    # Rename the dictionary file to the expected name if it exists
    if [[ -e "\${genomeFasta}.dict" ]]; then
        mv "\${genomeFasta}.dict" "\${genomeFasta%.*}.dict"
    fi

    outputVcf="\$(basename ${bamFile} _dedup.bam).vcf"

    # Use GATK HaplotypeCaller to call variants
    gatk HaplotypeCaller -R "\${genomeFasta}" -I ${bamFile} -O "\${outputVcf}"

    echo "Sample: ${sample_id} VCF: \${outputVcf}"

    echo "Variant Calling for Sample: ${sample_id} Complete"
    """
}