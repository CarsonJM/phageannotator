process NUCMER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mummer:3.23--pl5262h1b792b2_12' :
        'biocontainers/mummer:3.23--pl5262h1b792b2_12' }"

    input:
    tuple val(meta), path(query)

    output:
    tuple val(meta), path("*.delta")        , emit: delta
    tuple val(meta), path("*.coords.tsv")   , emit: coords
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def is_compressed_query = query.getName().endsWith(".gz") ? true : false
    def fasta_name_query    = query.getName().replace(".gz", "")
    """
    if [ "${is_compressed_query}" == "true" ]; then
        gzip -c -d ${query} > ${fasta_name_query}
    fi

    nucmer \\
        -p ${prefix} \\
        ${args} \\
        ${fasta_name_query} \\
        ${fasta_name_query}
    
    show-coords -c -q -T -H ${prefix}.delta > ${prefix}.coords.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nucmer: \$( nucmer --version 2>&1  | grep "version" | sed -e "s/NUCmer (NUCleotide MUMmer) version //g; s/nucmer//g;" )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.delta
    touch ${prefix}.coords

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nucmer: \$( nucmer --version 2>&1  | grep "version" | sed -e "s/NUCmer (NUCleotide MUMmer) version //g; s/nucmer//g;" )
    END_VERSIONS
    """
}
