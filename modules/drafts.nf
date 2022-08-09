process JOIN_STATS {
    // UNUSED
    input:
    path "*"

    output:
    path "stats.tsv"

    script:
    """
    cat *.tsv | head -n 1 > stats.tsv
    grep -v 'File' *.tsv | cut -f 2 -d ":" | sort >> stats.tsv
    """
}