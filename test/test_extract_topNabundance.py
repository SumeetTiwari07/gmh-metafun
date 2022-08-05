#!/usr/bin/env python3
import extract_topNabundance as et
import pandas as pd
import compare
import sys

"""
This script will perform small tests and compare results with an expected outcome.
Aim is to check if functions in the scripts are producing expected results.
To achive this different functions are imported from extract_topNabundance.py

"""


if __name__ == '__main__':

    # Inputs Genefamilies abundance:
    print("\033[1m"+"a) Testing extract_GenefamiliesAbundance ...!!!"+"\033[0;0m")
    sys.stdout.flush()
    gf_summary_file = "./input/Sample_Genefamilies.tsv"
    StringToTrim  = "_Abundance-RPKs"
    gf_expected_abundance = "./output/gf_abundance.csv"
    
    # Get genefamilies abundnace both extract and expected
    gf_file = pd.read_table(gf_summary_file, header=0) # Test genefamilies abundance file 
    expected_abundance = pd.read_csv(gf_expected_abundance, header = 0, index_col = 0) # Expected output
    extract_gf, extract_UR = et.extract_GenefamiliesAbundance(gf_file, StringToTrim) # Extracted ouput
    extract_gf.set_index(extract_gf.columns[0], inplace=True)
    compare.compare_dataframe(extract_gf, expected_abundance) # Check if the expected and extracted outputs are same.
    
    # Input Pathways abundance
    print("\033[1m"+"b) Testing extract_PathwaysAbundance function...!!!"+"\033[0;0m")
    sys.stdout.flush()
    pa_summary_file = "./input/Sample_Pathabundance.tsv"
    StringToTrim  = "_Abundance"
    pa_expected_abundance = "./output/pa_abundance.csv"

    # Get pathways abundance both extracted and expected
    expected_abundance = pd.read_csv(pa_expected_abundance, header=0, index_col = 0) # Expected output file.
    pa_file = pd.read_table(pa_summary_file, header=0) # Test pathways abundance file
    extract_pa = et.extract_PathwaysAbundance(pa_file, StringToTrim)
    extract_pa.set_index(extract_pa.columns[0], inplace=True)
    compare.compare_dataframe(extract_pa, expected_abundance)

    # Inputs top N genefamilies/pathways abundance
    N = 5 # Top 5
    print("\033[1m"+"c) Testing topN_gfORpa_across_samples function...!!!"+"\033[0;0m")
    sys.stdout.flush()
    abundance_file = "./input/abundance.csv"
    expected_output = "./output/Top5_gf_abundance.csv"

    # Get top N genefamilies estimated and expected.
    estimated_df = pd.read_csv(abundance_file, header = 0) # Absolute Genefamilies abundance file accross the samples.
    estimated_abs, estimated_rel = et.topN_gfORpa_across_samples(estimated_df, N) # Get top 5 genfamilies and their absolute abundance across the samples estimated by using the function.
    expected_abs = pd.read_csv(expected_output,header = 0, index_col=0) # Expected top 5 genefamilies and their absolute abundance across the sample.
    compare.compare_dataframe(estimated_abs, expected_abs) # Compare the two dataframes if the results are same.
    
    # Compare the stats file
    print("\033[1m"+"d) Testing qc_data function...!!!"+"\033[0;0m")
    total_reads = pd.read_table("./input/stats_10.tsv", header = 0) # Input stats file per sample (samplename "\t" Total_reads)
    unmapped_reads = pd.read_table("./input/unmapped_count.tsv", header = 0) # Unmapped reads per sample
    estimated_qcstats = et.qc_data(total_reads, unmapped_reads) # Estimated QC table
    expected_qcstats = pd.read_csv("./output/qc-report.csv", header = 0) # Expected QC table
    expected_qcstats.set_index(expected_qcstats.columns[0], inplace = True)
    compare.compare_dataframe(estimated_qcstats, expected_qcstats) # Compare the two dataframes.
    # Inputs for generating the stats file (expected vs estimated)
    print("\033[1m"+"All functions were successfully tested"+"\033[0;0m")
    

    