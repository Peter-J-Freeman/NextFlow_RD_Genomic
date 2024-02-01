/*
 * Define the indexGenome process that creates a BWA index
 * given the genome fasta file
 */
process indexGenome{

    label 'process_medium'
    container 'variantvalidator/indexgenome:1.1.0'

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
