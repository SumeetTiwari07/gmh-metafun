#!/usr/bin/env python3
import extract_topNtaxaabundance as eta
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
    print("\033[1m"+"a) Testing topN_species_across_samples ...!!!"+"\033[0;0m")
    sys.stdout.flush()
    sample_taxaabundance = "./input/Sample_taxaabundance.tsv"
    suffix  = "_metaphlan_bugs_list"
    topN=5 # Top 5 genus/species

    # Estimating top 5 genus
    taxa_abundance,unknown=eta.extract_taxa(sample_taxaabundance,suffix)
    ranks_lvl=eta.create_ranks()
    r=ranks_lvl["G"] # Genus
    estimated_top5genusabundance=eta.topN_species_across_samples(taxa_abundance,unknown,r,topN) # Estimated top 5 genus
    
    # Loading expected top 5 genus
    expected_top5genusabundance = pd.read_csv("./output/top5genus-rel.csv",header=0, index_col="Unnamed: 0") # Expected top 5 genus
    
    # Check if both results are same.
    compare.compare_dataframe(estimated_top5genusabundance, expected_top5genusabundance) # Check if the expected and extracted outputs are same.