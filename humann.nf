

nextflow.enable.dsl = 2
params.reads = "$baseDir/reads/*_R{1,2}.fastq.gz"
params.outdir = "nf-humann"

params.uniref = "/qib/platforms/Informatics/transfer/outgoing/databases/humann_db/uniref/"
params.chocophlan = "/qib/platforms/Informatics/transfer/outgoing/databases/humann_db/chocophlan/"
params.metaphlandb = "/qib/platforms/Informatics/transfer/outgoing/databases/humann_db/mpa/"

reads = Channel
        .fromFilePairs(params.reads, checkIfExists: true)

log.info """
         GMH Humann (version 0.2)
         ===================================
         input reads  : ${params.reads}
         outdir       : ${params.outdir}
         
         metaphlan    : ${params.metaphlandb}
         chocophlan   : ${params.chocophlan}
         uniref       : ${params.uniref}

         Running on   : ${params.max_memory}, ${params.max_cpus} cores
         """
         .stripIndent()
         
def UNIREF = file(params.uniref, checkIfExists: true)
def METAPHLANDB = file(params.metaphlandb, checkIfExists: true)
def CHOCOPHLAN = file(params.chocophlan, checkIfExists: true)

include  { VERSIONS; JSON_STATS } from "./modules/helpers.nf"
include  { HUMANN; JOIN_HUMANN; FASTP }   from "./modules/modules.nf"




process INTERLEAVE {
    tag "ilv $sample_id"
    label "process_low"
    
    input:
    tuple val(sample_id), path(reads) 
    
    output:
    tuple val(sample_id), path("${sample_id}.fastq")
  
    script:
    """
    seqfu interleave -1 ${reads[0]} -2 ${reads[1]} | seqfu cat --strip-name --strip-comments  --prefix "${sample_id}." > ${sample_id}.fastq
    """  
    stub:
    """
    gzip -dc ${reads[0]} > ${sample_id}.fastq
    """
}  

process multiqc {
    publishDir params.outdir, mode:'copy'
       
    input:
    path '*'  
    
    output:
    path 'multiqc_*'
     
    script:
    """
    multiqc . 
    """
} 



process TOP_TAXA {
    label 'process_filtering'
    publishDir "$params.outdir/top_taxa/", mode:'copy'

    input:
    path "merged-taxonomy.tsv"

    output:
    path '*.csv'

    script:
    """
    #Estimate the top N genefamilies and pathways in the study [update to: extract_topNtaxaabundance]
    extract_taxaabundance.py -i merged-taxonomy.tsv -o Species -r S -t 20
    extract_taxaabundance.py -i merged-taxonomy.tsv -o Genus -r G -t 20
    """
}
process JOIN_TAXONOMY {
    publishDir "$params.outdir", mode:'copy'
    input:
    path "*"

    output:
    path "tables/Taxa_merged.tsv"

    script:
    """
    mkdir -p tables
    merge_metaphlan_tables.py *_metaphlan_bugs_list.tsv > tables/Taxa_merged.tsv
    """
}

process TOP_PATHWAYS {
    publishDir "$params.outdir/", mode:'copy'
    input:
    path tables
    path 'stats.tsv'

    output:
    path 'summary'
    
    script:
    """
    #Estimate the top N genefamilies and pathways in the study
    mkdir -p summary
    extract_genefamilies.py -i GeneFamilies.tsv -st stats.tsv -qc summary/qc -o summary/top20 -t 20
    extract_pathabundance.py -i PathAbundance.tsv -o summary/top20 -t 20
    """
}

process BUBBLE_PLOTS {
    publishDir "$params.outdir/", mode:'copy'

    input:
    path summary

    output:
    path 'plots-pathways'

    script:
    """
    #Plot the top N genefamilies and pathways in the study
    bubble_plot.R summary/qc_report.csv summary/top20_GFabundance-rel.csv top20gf False ./
    bubble_plot.R summary/qc_report.csv summary/top20_PWabundance-rel.csv top20pa False ./
    mv plots plots-pathways
    """
}
process BUBBLE_TAXA {
    publishDir "$params.outdir/", mode:'copy'
    input:
    path top_taxa

    output:
    path 'plots-taxa'

    script:
    """
    #Plot the top N genefamilies and pathways in the study
    mkdir -p plots
    heatmap.R Species-rel.csv ./plots-taxa/
    heatmap.R Genus-rel.csv ./plots-taxa
    """
}
workflow {
   // Collect versions (and ensures that the programs are available)
   VERSIONS()
   // Clean reads, interleave (versions is used as lock to run VERSION first)
   FASTP(reads, VERSIONS.out)

   // Count reads, merge counts
   //STATS(FASTP.out.reads)
   //JOIN_STATS( STATS.out.map{it -> it[1]}.collect() )
   JSON_STATS(FASTP.out.json.map{it -> it[1]}.collect())

   // RUN HUMANN and use its tools to merge tables
   HUMANN(FASTP.out.reads, CHOCOPHLAN, UNIREF, METAPHLANDB)
   JOIN_TAXONOMY((HUMANN.out.metaphlan).map{it -> it[1]}.collect())
   JOIN_HUMANN( HUMANN.out.genefamilies.mix( HUMANN.out.pathabundance, HUMANN.out.pathcoverage).map{it -> it[1]}.collect())

   // Sumeet Tiwari scripts to summarise results
   TOP_PATHWAYS(JOIN_HUMANN.out, JSON_STATS.out)
   BUBBLE_PLOTS(TOP_PATHWAYS.out)
   TOP_TAXA(JOIN_TAXONOMY.out)
   BUBBLE_TAXA(TOP_TAXA.out)
}
