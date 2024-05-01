process PROPAGATE_EXTRACTPROVIRUSES {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-80c23cbcd32e2891421c54d1899665046feb07ef:77a31e289d22068839533bf21f8c4248ad274b60-0' :
        'biocontainers/mulled-v2-80c23cbcd32e2891421c54d1899665046feb07ef:77a31e289d22068839533bf21f8c4248ad274b60-0' }"

    input:
    tuple val(meta) , path(fasta)
    tuple val(meta2), path(virus_summary)
    tuple val(meta3), path(contamination_summary)

    output:
    tuple val(meta), path("${prefix}_provirus_scaffolds.fasta.gz")  , emit: provirus_scaffolds
    tuple val(meta), path("${prefix}_provirus_coords.tsv")          , emit: provirus_coords
    path "versions.yml"                                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    extract_proviruses.py \\
        --fasta ${fasta} \\
        --genomad ${virus_summary} \\
        --checkv ${contamination_summary} \\
        --output ${prefix}_provirus_scaffolds.fasta \\
        --tsv ${prefix}_provirus_coords.tsv \\
        ${args}

    gzip ${prefix}_provirus_scaffolds.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed 's/Python //' )
        biopython: \$(echo \$(biopython_version.py 2>&1))
        numpy: \$(echo \$(numpy_version.py 2>&1))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}_provirus_scaffolds.fasta.gz
    touch ${prefix}_provirus_coords.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed 's/Python //' )
        biopython: \$(echo \$(biopython_version.py 2>&1))
        numpy: \$(echo \$(numpy_version.py 2>&1))
    END_VERSIONS
    """
}
