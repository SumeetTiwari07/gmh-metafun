process HUMANN {
    tag "$sample_id"
    publishDir "$params.outdir/humann/", mode:'copy'
    label "humann"
    
    input:
    tuple val(sample_id), path(reads) 
    path(chocophlan)
    path(uniref)
    path(metaphlandb)
    
    output:
    tuple val(sample_id), path("${sample_id}_genefamilies.tsv"), emit: genefamilies
    tuple val(sample_id), path("${sample_id}_pathabundance.tsv"), emit: pathabundance
    tuple val(sample_id), path("${sample_id}_pathcoverage.tsv"), emit: pathcoverage
    tuple val(sample_id), path("${sample_id}_metaphlan_bugs_list.tsv"), emit: metaphlan
    
  
    script:
    """
    INDEX=\$(cat mpa/mpa_latest)
    humann -i ${reads} -o $sample_id --threads ${task.cpus} --output-basename $sample_id \
        --nucleotide-database ${chocophlan} \
        --protein-database ${uniref} \
        --metaphlan-options="-x \$INDEX --bowtie2db ${metaphlandb} --unknown_estimation"
    mv ${sample_id}/*tsv .
    mv ${sample_id}/${sample_id}_humann_temp/${sample_id}_metaphlan_bugs_list.tsv ${sample_id}_metaphlan_bugs_list.tsv
    """  
    stub:
    """
    getoutput.py $sample_id
    """
}  
 

process FASTP {
    /* 
       fastp process to remove adapters and low quality sequences
    */
    tag "filt $sample_id"
    label "process_low"

    input:
    tuple val(sample_id), path(reads) 
    path("versions.lock")
    
    output:
    tuple val(sample_id), path("${sample_id}.fq.gz"), emit: reads
    tuple val(sample_id), path("${sample_id}.fastp.json"), emit: json
    tuple val(sample_id), path("${sample_id}.fastp.html"), emit: html
    
  
    script:
    """
    fastp -i ${reads[0]} -I ${reads[1]} --stdout \
      --detect_adapter_for_pe --length_required 75 --thread ${task.cpus} \
      -h ${sample_id}.fastp.html -j ${sample_id}.fastp.json \
      | seqfu cat --strip-name --strip-comments  --prefix "${sample_id}." \
      | gzip -c > ${sample_id}.fq.gz
    """     
    stub:
    """
    cat ${reads[0]} | head -n 4000 | gzip -c > ${sample_id}.fq.gz
    echo '{ "summary": { "after_filtering": { "total_reads":1090016, "total_bases":159862544,"read1_mean_length":146}}}' > ${sample_id}.fastp.json
    touch ${sample_id}.fastp.html
    """
}


process SEQFU_STATS {
    /* 
       fastp process to remove adapters and low quality sequences
    */
    tag "$sample_id"
    label "process_low"

    input:
    tuple val(sample_id), path(reads) 
    
    
    output:
    tuple val(sample_id), path("${sample_id}.tsv")
  
    script:
    """
    seqfu stats ${reads[0]} > ${sample_id}.tsv
    """     
    stub:
    """
    echo 'File,#Seq,Total bp,Avg,N50,N75,N90,auN,Min,Max' > tmp
    echo -n '${sample_id}' >> tmp
    echo ',16434,2465100,150.0,150,150,150,0.009,150,150' >> tmp
    sed 's/,/\\t/g' tmp > ${sample_id}.tsv
    """
}


process JOIN_HUMANN {
    publishDir "$params.outdir", mode:'copy'
    
    input:
    path '*'

    output:
    path 'tables/*.tsv'

    script:
    """
    mkdir -p tables
    humann_join_tables -i ./ --file_name _genefamilies.tsv -s -o tables/GeneFamilies.tsv
    humann_join_tables -i ./ --file_name _pathabundance.tsv -s -o tables/PathAbundance.tsv
    humann_join_tables -i ./ --file_name _pathcoverage.tsv -s -o tables/PathCoverage.tsv
    """
}
