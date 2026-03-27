#!/bin/bash

# 2023-11-24
# This is a script to automate analysis of peptide affinity using netMHCIIpan-4.3 for DR class MHC II alleles 
# author: Barbara Walkowiak bw450

# INPUT 
# list of HLA alleles to generate predictions for 
# folder with files (.txt) that contain sequences of peptides of interest
# each contains sequences of peptides to examine - kmers generated with get_kmers_from_csv_all_updated_stop.py script 

# OUTPUT 
# txt files (in netMHCII_out/txt folder) and xls files (in netMHCII_out/xls folder)
# which will contain predictions (EL)for all peptide kmers in files contained in the kmer folders in the input 
# each file will contain predictions for a specific set of peptides that all contain the same variant (one txt file of kmers > one txt file of predictions)

# HOW TO RUN THIS SCRIPT 
# cd ~/Desktop/msc_thesis
# example usage: bash task1_predict_binding_to_HLA/scripts02/2_get_predictions_netMHCIIpan_dr.sh data/20231108_list_of_HLA_alleles.txt kmers/kmers_231116

mkdir -p ~/Desktop/msc_thesis/netMHCII_out/xls
mkdir -p ~/Desktop/msc_thesis/netMHCII_out/txt_dr
mkdir -p ~/Desktop/msc_thesis/logs

# IDENTIFY THE FOLDER with sequences of peptides (15-mers) which can be analysed for presentation 
directory=$2

# RUN PREDICTIONS
for fasta in "$directory"/*
do

    echo $fasta
    base=`basename $fasta .txt`
    string='STOP'
    base="${base/"*"/$string}" # replace * with STOP (in case there is a * in the name of the file)
    hla=$1
    base_hla=`basename $hla .txt`
    hla_list=$(<$hla)

    log=~/Desktop/msc_thesis/logs/${base}_${base_hla}_netMHCII_log.txt
    out=~/Desktop/msc_thesis/netMHCII_out/20240613_txt_dr/${base}_${base_hla}_netMHCII_out_affinities.txt
    echo "Sample name is $base, HLA alleles examined are $hla_list"
    
    for hla in "${hla_list[@]}"; do
    
    netMHCIIpan-4.3/netMHCIIpan -a $hla -BA -f $fasta >> $out # run prediction "
    
    done

done


