/*
 * Downsample the BAM files
 */
process downsampleBam {

    label 'process_single'
    container 'variantvalidator/indexgenome:1.1.0'

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
        picard DownsampleSam I=${bamFile} O="\${outputBam}" P="\${fraction}" QUIET=true VALIDATION_STRINGENCY=SILENT > "\${outputBam}.log" 2>&1
    done
    if [ -e "${bamFile}" ]; then
        mv "${bamFile}" "${bamFile.baseName}_downsampled_100.bam"
    fi
    echo "Downsampling complete"
    """
}
