#!/bin/bash

# 2023-11-02 
# This is a test script to automate analysis of peptide affinity using netMHCpan-4.3 for DP-DQ class MHC II allele combinations 
# author: Barbara Walkowiak bw450

# INPUT 
# list of HLA alleles to generate predictions for (DP-DQ combinations)
# folder with files (.txt) that contain sequences of peptides of interest
# each contains sequences of peptides to examine - kmers generated with get_kmers_from_csv_all_updated_stop.py script 

# OUTPUT 
# txt files (in netMHCII_out/txt folder) and xls files (in netMHCII_out/xls folder)
# which will contain predictions (EL)for all peptide kmers in files contained in the kmer folders in the input 
# each file will contain predictions for a specific set of peptides that all contain the same variant (one txt file of kmers > one txt file of predictions)

# HOW TO RUN THIS SCRIPT 
# cd ~/Desktop/msc_thesis
# example usage: bash task1_predict_binding_to_HLA/scripts02/2_get_predictions_netMHCIIpan_dpq.sh data/20231108_list_of_HLA_alleles.txt kmers/kmers_231116

mkdir -p ~/Desktop/msc_thesis/netMHCII_out/20240613_txt_dpq
mkdir -p ~/Desktop/msc_thesis/logs

# IDENTIFY THE FOLDER with sequences of peptides (15-mers) which can be analysed for presentation 
directory=$2

# RUN PREDICTIONS
for fasta in "$directory"/* 

do
    
    base=`basename $fasta .txt`
    string='STOP'
    base="${base/"*"/$string}" # change to stop if there is a stop codon (easier to run, * throw this thing off)
    hla=$1
    base_hla=`basename $hla .txt`
    hla_list=$(<$hla)
    output_txt=~/Desktop/msc_thesis/netMHCII_out/20240613_txt_dpq/${base}_${base_hla}_netMHCII_out.txt
    # note that I am completely ignoring excel files and will only write txt

    echo "Sample name is $base, HLA alleles examined are $hla_list"

    # define function to split HLA-DPAxxx-DPByyy into DPA and DPB in correct formats 
    split_arguments() {

        # loop over each argument 
        
        for arg in "$@"; do

            # split (because netMHC wants alpha and beta separately)
            arg1="${arg%%-*}"
            arg2="${arg#*-}"
            netMHCIIpan-4.3/netMHCIIpan -BA -choose ON -cha $arg1 -chb $arg2 -f $fasta >> $output_txt # write everything to a single txt file 
            
        done
    }

    # split HLA list into separate alpha / beta chain list 
    IFS=',' read -ra arguments <<< "$hla_list"

    # call the function 
    split_arguments "${arguments[@]}"

done
