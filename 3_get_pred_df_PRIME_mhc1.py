#! /python3

# 2023-11-14 (modified 23-11-23)
# author: Barbara Walkowiak bw450
# script to obtain the dataframe with predictions for different CH variants-allele combinations from PRIME2.0

# INPUT: 
# out/csv file (from PRIME_specify_hla.sh)

# OUTPUT: 
# output: a csv file which has binding affinity scores for each (specified) HLA for a list of desired variants

# HOW TO RUN THIS SCRIPT:
# cd ~/Desktop/msc_thesis/PRIME_out/csv
# python3 ../../task1_predict_binding_to_HLA/scripts02/3_get_pred_df_PRIME_mhc1.py /Users/barbarawalkowiak/Desktop/msc_thesis/PRIME_out/csv/csv_20240220

# %%

# IMPORTS 
import sys
import pandas as pd
import os
import time

# %%
timestr = time.strftime("%Y%m%d") # get current date 

# path to the directory with csv files generated from previous step  
directory = sys.argv[1] 

# %%

# create an empty dataframe to append to at the end of the loop
df_all = pd.DataFrame() 

# for each HLA allele, PRIME outputs 3 parameters: %rank, score, %rank_binding 
# %Rank indicates the fraction of random peptides from the human proteome (length 8 to 14) that would have a score higher or equal to the peptide given in input

for file in os.listdir(directory):

    filename = os.fsdecode(file)

    # for PRIME, I obtained csv files first so will want to analyse files in those format to get to the dataframe 
    if filename.endswith(".csv"): 

        # read in the file 
        csv = pd.read_csv(file)
        
        # identify the gene and variant, and whether predictions for wt or mutant kmers 
        gene = filename.split('_')[2]
        variant = filename.split('_')[3] 
        genotype = filename.split('_')[4] 

        # Obtain predictions of the %Rank 
        rank = '%Rank_'
        ranks = csv.filter(regex=rank, axis=1)
        stacked_ranks = ranks.stack().reset_index()
        stacked_ranks.columns = ['index', 'allele', '%Rank']
        stacked_ranks['allele'] = stacked_ranks['allele'].str.split('_').str[1]

        # remove the "bestAllele" column
        mask = stacked_ranks['allele'] != 'bestAllele'
        stacked_ranks = stacked_ranks[mask]
    
        # additional parameters that may be useful (nr of peptides with %EL_rank < 0.5 and nr of peptides with %EL_rank < 2)
        min_values_rank = stacked_ranks.groupby('allele')['%Rank'].min().reset_index()
        sum_05_rank = stacked_ranks.groupby('allele')['%Rank'].apply(lambda x: (x < 0.5).sum()).reset_index(name='sum_peptides_below_05')
        sum_2_rank = stacked_ranks.groupby('allele')['%Rank'].apply(lambda x: (x < 2).sum()).reset_index(name='sum_peptides_below_2')

        # merge dfs by shared column (allele)
        merged_ranks = pd.merge(sum_05_rank, sum_2_rank, on='allele', how="inner")
        merged_ranks = pd.merge(min_values_rank, merged_ranks, on='allele', how="inner")

        merged_ranks['gene'] = gene
        merged_ranks['variant'] = variant 
        merged_ranks['genotype'] = genotype
        
        # make allele names consistent across different methdos 
        string_to_add = 'HLA-' # add "HLA-" to make consistent with other lists  
        merged_ranks['allele'] = string_to_add + merged_ranks['allele']
        character_to_add = ':'
        merged_ranks['allele'] = merged_ranks['allele'].apply(lambda x: x[:-2] + character_to_add + x[-2:])

        # specify column names 
        merged_ranks.columns = ['allele', 'min_rank', 'sum_peptides_below_05', 'sum_peptides_below_2', 'gene', 'variant', 'genotype']
        
        df_all = pd.concat([df_all, merged_ranks])

    else:
         continue

# # write the complete file to csv
df_all.to_csv("../../scores/" + timestr + "_percent_ranks_for_each_variant_by_HLA.csv", index = False)
