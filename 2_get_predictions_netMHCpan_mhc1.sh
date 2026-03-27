#!/bin/bash

# 2023-11-24
# This is a script to automate analysis of peptide affinity using netMHCpan-4.1

# INPUT: 
# list of HLA alleles to generate predictions for 
# folder with files (.txt) that contain sequences of peptides of interest
# each contains sequences of peptides to examine - kmers generated with get_kmers_from_csv_all_updated_stop.py script 

# OUTPUT: 
# txt files (in netMHC_out/txt folder) and xls files (in NetMHC_out/xls folder)
# which will contain predictions (binding affinity + EL) for all peptide kmers in files contained in the kmer folders in the input 
# each file will contain predictions for a specific set of peptides that all contain the same variant (one txt file of kmers > one txt file of predictions)

# HOW TO RUN THIS SCRIPTS 
# cd ~/Desktop/msc_thesis (move the the correct directory)
# example usage: bash task1_predict_binding_to_hla/scripts02/2_get_predictions_netMHCpan_mhc1.sh data/20231108_list_of_HLA_alleles.txt kmers/kmers_231116
# provide path to the alleles for which predictions are needed and folder with kmer sequences 

# CREATE DIRECTORIES TO SAVE OUTPUT TO 
mkdir -p ~/Desktop/msc_thesis/netMHC_out/xls
mkdir -p ~/Desktop/msc_thesis/netMHC_out/txt
mkdir -p ~/Desktop/msc_thesis/logs

# IDENTIFY THE FOLDER with sequences of peptides (8- to 11-mers) which can be analysed for presentation
directory=$2

# RUN PREDICTIONS
for fasta in "$directory"/*
do

    base=`basename $fasta .txt`
    hla=$1 # LIST OF HLA ALLELES OF INTEREST 
    base_hla=`basename $hla .txt`
    hla_list=$(<$hla)

    log=~/Desktop/msc_thesis/logs/${base}_${base_hla}_netMHC4.1_log.txt
    out=~/Desktop/msc_thesis/netMHC_out/txt/${base}_${base_hla}_netMHC_out_affinities.txt
    out_xls=~/Desktop/msc_thesis/netMHC_out/xls/${base}_${base_hla}_netMHC_out_affinities.xls
    out_xlsx=~/Desktop/msc_thesis/netMHC_out/xls/${base}_${base_hla}_netMHC_out.xlsx
    
    echo "Sample name is $base, HLA alleles examined are $hla_list"
    
    ../netMHCpan-4.1/netMHCpan -p -a $hla_list $fasta -BA > $out # run prediction 
    ../netMHCpan-4.1/netMHCpan -p $fasta -BA -xls -a $hla_list -xlsfile $out_xls # save as .xls file 
    ssconvert $out_xls $out_xlsx # convert to xlsx file 

done


