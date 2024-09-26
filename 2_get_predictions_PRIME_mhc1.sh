#!/bin/bash

# 2023-11-14 
# This is a script to automate analysis of peptide affinity using PRIME (see https://github.com/GfellerLab/PRIME)
# author: Barbara Walkowiak bw450

# INPUT 
# list of HLA alleles to generate predictions for 
# folder with files (.txt) that contain sequences of peptides of interest
# each contains sequences of peptides to examine - kmers generated with get_kmers_from_csv_all_updated_stop.py script 

# OUTPUT 
# txt files (in PRIME_out/txt folder) and xls files (in PRIME_out/xls folder)
# which will contain predictions (EL)for all peptide kmers in files contained in the kmer folders in the input 
# each file will contain predictions for a specific set of peptides that all contain the same variant (one txt file of kmers > one txt file of predictions)

# HOW TO RUN THIS SCRIPT  
# cd ~/Desktop/msc_thesis/task1_predict_binding_to_HLA
# run from command line: bash scripts02/2_get_predictions_PRIME_mhc1.sh data/hla_15_most_common.csv kmers/kmers_20231116

# CREATE DIRECTORIES TO SAVE OUTPUT TO 
mkdir -p ~/Desktop/msc_thesis/PRIME_out/csv
mkdir -p ~/Desktop/msc_thesis/PRIME_out/txt
mkdir -p ~/Desktop/msc_thesis/logs

# IDENTIFY THE FOLDER with sequences of peptides (8- to 11-mers) which can be analysed for presentation 
directory=$2

# RUN PREDICTIONS
for fasta in "$directory"/*

do
    base=`basename $fasta .txt`
    string='STOP'
    base="${base/"*"/$string}"
    hla=$1
    base_hla=`basename $hla .txt`
    hla_list=$(<$hla)

    # PRIME does not accept * in output 
    # 'stars in output paths are not supported'
    # find the star and modify base name (change * to STOP)
    log=~/Desktop/msc_thesis/logs/${base}_${base_hla}_PRIME_log.txt
    out_txt=~/Desktop/msc_thesis/PRIME_out/txt/${base}_${base_hla}_PRIME_out.txt
    out_txt2=~/Desktop/msc_thesis/PRIME_out/txt/${base}_${base_hla}_PRIME_out_data.txt
    out_csv=~/Desktop/msc_thesis/PRIME_out/csv/${base}_${base_hla}_PRIME_out.csv
    # exec >$log 2>&1

    echo "Sample name is $base, HLA alleles examined are $hla_list"
    
    PRIME -i $fasta -o $out_txt -a $hla_list 

    # convert .txt files to .csv
    grep -v '^#' $out_txt > $out_txt2
    tr -s " " < $out_txt2 | sed 's/\t/,/g' > $out_csv

done
