#!/usr/bin/env python3

import pandas as pd
import argparse
import sys

"""
This script will do the following things:
    1. Extract the genefamilies abundance from the output of the Humann2/3.
    2. Genefamilies abundance at species level are discarded.
    3. Genefamilies abundance at community level are kept.
    Note: For detailied description of the input file format please refer to the https://github.com/biobakery/humann#1-gene-families-file
    3. TopN Gene families in whole experiment was estimated.
    4. Abundance of TopN Gene families across the samples is extracted.
    5. Proportion of TopN Gene families contributing each sample is estimated.

Final outputs of the scripts:
Case 1: For genefamilies summary file.
  Inputs are:
    --input = Sample_Genefamilies.tsv
    --stats = stats_10.tsv.
    --topN = 5
    --output = test
    
  Outputs will be:
    test-rel.csv
    test-abs.csv
    qc_report.csv

Case 2: When pathabudance summary file is provided no qc file will be generated
  Inputs are:
    --input = Sample_Genefamilies.tsv
    --topN = 5
    --output = test

Outputs will be:
    test-rel.csv
    test-abs.csv
"""

def extract_GenefamiliesAbundance(gf, gstring):
    """
    Extract the community-level genefamilies abundance per sample.
    Store the unmapped reads per sample.
    """
    gf.columns = gf.columns.str.replace(gstring,"") # Renaming columnames by replacing based on string default: "_Abundance-RPKs"
    Ureads = gf.loc[gf[gf.columns[0]] == 'UNMAPPED'] # Storing Unmapped reads
    # Discarding unmapped reads.
    gf = gf.drop(Ureads.index)
    return gf, Ureads

def extract_PathwaysAbundance(input_file, pstring):
    """
    Extract the community-level pathways abundance per samples.
    Discrad the UNINTEGRATED and UNMAPPED Reads per sample.
    """
    pa = input_file[~input_file[input_file.columns[0]].str.contains('\|')] # Extracting community-level pathways abundance.
    pa.columns = pa.columns.str.replace(pstring,"") # Renaming columnames by replacing based on string default: "_Abundance"

    # Discard UNINTEGRATED
    pa = pa.loc[pa[pa.columns[0]] != 'UNINTEGRATED']

    # Discard UNMAPPED reads.
    pa = pa.loc[pa[pa.columns[0]] != 'UNMAPPED']

    # Extracting the pathway names
    pa[pa.columns[0]] = pa[pa.columns[0]].str.split(":", expand=True)[0]

    return pa

def topN_gfORpa_across_samples(df, N):
    """
    Estimate topN genefamilies or pathways 
    """
    df.set_index(df.columns[0], inplace = True) #Setting index at column name: ID.
    
    # Calculate total abundance of each gene family/pathways
    total_abundance = df.sum(axis=1)
    
    # Get top N genefamily/pathways based on abundance
    topN = df.loc[total_abundance.nlargest(N).index]
    topN = topN.T # tranpose
    # Estimate the abundance of others.
    others = df.loc[total_abundance.nsmallest(len(df)-N).index].sum(axis=0)
    
    # Add others to topN
    topN = pd.concat([topN,others],axis=1)
    topN.set_axis([*topN.columns[:-1], 'Others'], axis=1, inplace=True)
    
    # Estimate relative abundance
    topN_rel = topN.div(topN.sum(axis = 1), axis = 0)

    return topN, topN_rel

def qc_data(stats, unmap):
    """
    Generate qc file with: Samplenames, # Total_reads, # Unmapped, # Mapped
    """
    stats.set_index(stats.columns[0], inplace = True) # Set index on first column
    unmap.set_index(unmap.columns[0], inplace = True) # Set index on first column
    unmap = unmap.T

    # Check if the sample-names in stats matched with the sample-names in unmap.
    samplenames_in_stats = set(stats.index.to_list())
    samplenames_in_unmap = set(unmap.index.to_list())

    if samplenames_in_stats == samplenames_in_unmap:
        qc_data = pd.concat([stats,unmap], axis = 1, join = "inner") # Concat to dataframes
        qc_data['Mapped'] = qc_data[qc_data.columns[0]]-qc_data[qc_data.columns[1]] # Estimate mapped reads
        qc_data = qc_data.rename(columns={qc_data.columns[0]:"Total_reads",qc_data.columns[1]:"Unmapped"}) # Rename column names.
        qc_data = qc_data.astype(int)

    else:
        qc_data = pd.DataFrame()
        sys.exit("Sample name miss match or missing\n")

    return qc_data

def generateOutput(abs, rel, metadata):
    if metadata.empty:
        # Write absoulute and relative abundance.
        abs.to_csv(f'{args.output}-abs.csv')
        rel.to_csv(f'{args.output}-rel.csv')
    else:
        # Adding metadata of the samples if provided.
        abs_meta=pd.concat([metadata,abs],axis=1,join="inner")
        rel_meta=pd.concat([metadata,rel],axis=1,join="inner")
        # Writing  output to the file
        abs_meta.to_csv(f'{args.output}-abs.csv')
        rel_meta.to_csv(f'{args.output}-rel.csv')


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__ , formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-i','--input', type=str, required=True, help="Gene Families/Path(way) Abundance summary file across multiple samples")
    parser.add_argument('-o','--output', type=str, required=True, help="Prefix for output file")
    parser.add_argument('-st','--stats', type=str, help="A tsv file containing sample names and total number of reads per sample.")
    parser.add_argument('-gs','--gstring', type=str, default="_Abundance-RPKs", help="String suffix to trim from sample names in genefamilies abundance summary file (default: %(default)s)")
    parser.add_argument('-ps','--pstring', type=str, default="_Abundance", help="String suffix to trim from sample names in PathAbundance summary file (default: %(default)s)")
    parser.add_argument('-qc','--qc',type=str, default="qc", help="Save the qc report (default: %(default)s)")
    parser.add_argument('-m','--metadata', type=str, help="Metadata about the samples in csv format with colum name: 'SampleId' containing the sample names")
    parser.add_argument('-t','--topN', type=int, default=20, help="top N genefamilies (default: %(default)d)")
    parser.add_argument('--version', action='version', version='%(prog)s 1.1')
    args = parser.parse_args()

    # Parsing input file and checking if its a genefamily abundance or pathways abundance summary file.
    input_file=pd.DataFrame()
    for chunk in pd.read_table(args.input, header = 0, chunksize=100000):
        if input_file.empty:
            input_file = chunk[~chunk[chunk.columns[0]].str.contains('\|')] # Extracting abundance at community-level.
        else:
            chunk = chunk[~chunk[chunk.columns[0]].str.contains('\|')]
            input_file = pd.concat([input_file,chunk])
    # Count the occurrence of UNINTEGRATED. This is only suppose to find in pathways abundance file.
    count = input_file[input_file.columns[0]].str.count("UNINTEGRATED").sum()
    
    if count > 0:

        pathways_abundance = extract_PathwaysAbundance(input_file, args.pstring)

        # Estimating topN Gene families in the whole experiment and their abundance across the samples.        
        topNpa_abs, topNpa_rel = topN_gfORpa_across_samples(pathways_abundance, args.topN)
        
        # Metadata file is provided.
        if args.metadata:
            metadata=pd.read_csv(args.metadata,header=0,index_col="SampleId")
        else:
            metadata = pd.DataFrame() # empty dataframe
        
        generateOutput(topNpa_abs, topNpa_rel, metadata)
                
    else:

        genefamilies_abundance, unmapped_reads = extract_GenefamiliesAbundance(input_file, args.gstring)

        # Estimating topN Gene families in the whole experiment and their abundance across the samples.        
        topNgf_abs, topNgf_rel = topN_gfORpa_across_samples(genefamilies_abundance, args.topN)
                
        # If metadata file is provided.
        if args.metadata:
            metadata=pd.read_csv(args.metadata,header=0,index_col="SampleId")
        else:
            metadata = pd.DataFrame() # empty dataframe
        
        generateOutput(topNgf_abs, topNgf_rel, metadata)
        
        # Check if the qc statistic file is required. This will only be used along with Genefamilies file
        if args.stats is not None :

            stats = pd.read_table(args.stats, header = None) # Reads stats file.
            qc_stats = qc_data(stats, unmapped_reads) # Generate QC file.
            
            if not qc_stats.empty:
                # Write a final qc report of each samples: Samplenames,Total_reads,Unmapped,Mapped 
                qc_stats.to_csv(f'{args.qc}-report.csv')


  