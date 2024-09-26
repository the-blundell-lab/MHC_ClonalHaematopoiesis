#! /python3

# 20234-02-10
# cd Desktop/msc_thesis/task1_predict_binding_to_HLA
# script to obtain the dataframe with predictions for different CH variants-allele combinations from NetMHCIIpan-4.3 
# author: Barbara Walkowiak bw450

# INPUT: 
# out/txt file (from netMHCIIpan_specify_hla_kmers_add_affinity.sh)

# OUTPUT: 
# csv file which has binding affinity scores for each (specified) HLA for a list of desired variants

# HOW TO RUN THIS SCRIPT:
# cd ~/Desktop/msc_thesis/netMHCII_out/txt
# python3 ../../task1_predict_binding_to_HLA/scripts02/3_get_pred_df_netMHCIIpan_mhc2.py /Users/barbarawalkowiak/Desktop/msc_thesis/netMHCII_out/txt

# %%

# IMPORTS 
import sys
import pandas as pd
import os
import time

# %%
timestr = time.strftime("%Y%m%d") # get current date 

# path to the directory with txt files generated from previous step
directory = sys.argv[1] 

# %%
# create an empty dataframe in which to store scores
dfs_all = pd.DataFrame()

# for each txt file in the directory provided
for file in os.listdir(directory):

    filename = os.fsdecode(file)
    if filename.endswith(".txt"): 

        input_file_path = os.path.join(directory, filename)
        
        out_filename = filename.split('.')[0]
        out_filename = out_filename + "_scores.csv"
        output_csv_path = os.path.join(directory, out_filename)

        # Initialize variables to store header and data
        header = None
        data = []

        # extract gene and variant data from file name
        print(filename)

        gene = filename.split('_')[2]
        variant = filename.split('_')[3]
        genotype = filename.split('_')[4]

        # read in the input file 
        with open(input_file_path, 'r') as file:
            lines = file.readlines()

        # Process lines to extract header and data
        for line in lines:
            if line.startswith('#'):
                # ignore comment lines 
                continue
            elif '---' in line:
                # if there is a ---, the next line is going to be the header
                header = lines[lines.index(line) + 1].strip().split()
            elif line.strip():
                # if there is sth, treat as data 
                data.append(line.strip().split())

        # check these columns for a specific value 
        column_to_check1 = 3 

        # specify the value to exclude rows with in the specified column
        value_to_exclude1 = "binders:"  # we don't want to write rows which specify the number of strong / weak binders

        # remove rows with the specified value in the specified column
        data = [row for row in data if row[column_to_check1] != value_to_exclude1]

        # find and keep only one instance of each duplicate row as the header
        seen_rows = set()
        unique_rows = []
        for row in data:
            row_tuple = tuple(row)
            if row_tuple not in seen_rows:
                seen_rows.add(row_tuple)
                unique_rows.append(row)

        df = pd.DataFrame(unique_rows)
        df.columns = df.iloc[0]
        df = df[1:].reset_index(drop=True)

        # some dfs have 15 columns and some have 16
        # the last column does not contain useful data and can be dropped 
        if df.shape[1] == 16:
            df = df.iloc[:, :-1]
        
        df.columns = ('Pos', 'MHC', 'Peptide', 'Of', 'Core', 'Core_Rel', 'Inverted', 'Identity', 
        'Score_EL', '%Rank_EL', 'Exp_Bind', 'Score_BA', '%Rank_BA', 'Affinity (nM)', 'BindLevel')

        df["gene"] = gene
        df["variant"] = variant 
        df["genotype"] = genotype
        df["gene_var"] = df["gene"] + "_" + df["variant"] # add gene_var column
        
        # note: NOT ALL COLUMNS HAVE BA SCORE / %RANK_BA / AFFINITY DUE TO restrictions on this software performance (it can't predict affinities for all combinations)
        dfs_all = pd.concat([dfs_all, df])

    else:
        continue

# clean up column names 
dfs_all = dfs_all[dfs_all["Pos"] != "Pos"] # clean up headers
dfs_all.rename(columns = {'MHC':'HLA'}, inplace = True) # rename column to match desired name going forwad

# save the dataframe to scores (assumes there is a scores folder created: NetMHC_out/scores)
dfs_all.to_csv("../scores/" + timestr + '_NetMHC_HLA_UKBB_with_affinities.csv', index=False)
