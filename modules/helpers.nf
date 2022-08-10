
process VERSIONS {
    label "process_low"
    
    output:
    file("versions.txt")

    script:
    """
    fastp --version 2>&1 > versions.txt
    seqfu version 2>&1 >> versions.txt
    humann --version 2>&1 >> versions.txt
    """
}

process JSON_STATS {
    input:
    path "*"

    output:
    path "stats.tsv"

    script:
    """
    getCountsFromFastp.py --fields 2 --output stats.tsv *.fastp.json
    """
}