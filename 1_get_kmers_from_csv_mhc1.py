#!/usr/bin/python3

# 20230-11-02 / 2023-11-03 updated in February 2024
# script to automate generation of 8-11 k-mers (MHC I presentation) from fasta sequences 

# INPUT: 
# this script requires a dataframe with 3 columns
# the first column = gene name, the second column = variant (e.g., R882H), the third column = URL to UniProt fasta seq
# example URL to UniProt fasta seq: https://rest.uniprot.org/uniprotkb/P22681.fasta where P22681 is a UniProt ID for the proteoform of interest

# OUTPUT: this script outputs a folder (named kmers_timestr where timestr is the current date) which contains kmers generated from each variant
# for a given variant, files with kmer sequences for wt (reference) and mutant peptides are saved separately

# HOW TO RUN THIS SCRIPT:
# cd ~/Desktop/msc_thesis/task1_predict_binding_to_HLA (move to the path to the folder which contains the folder with scripts)
# python3 scripts02/1_get_kmers_from_csv_mhc1.py data/ch_var_test.csv (run the script, providing the path to the dataframe as the argument)
# test_ch_var.csv is a dataframe that contains three columns: gene, variant, URL link to UniProt sequence (protein)

# %%
# IMPORTS 
from urllib.request import urlopen
import sys
import pandas as pd
import time
import re 
import os

# %%

# read the file with variants and sequences 
file = sys.argv[1]
data = pd.read_csv(file, sep=",", chunksize=1) 

# Get current date (used to name folder with kmers to track when kmers were generated)
timestr = time.strftime("%Y%m%d")  

# %% 
# FUNCTIONS

# Define function to build kmers 
def build_kmers(fasta, ksize):
    kmers = [] # initialise list 
    for i in range(ksize): 
        kmer = fasta[11-ksize+i : 11+i] 
        kmers.append(kmer)

    return kmers

# %%
# GENERATE KMERS 

# for each CH variant included in the dataframe: 
for line in data:

# link = str(sys.argv[1]) # use if you want to provide this as command line 
    url = line.iloc[:,2].values[0] # get the URL on the spreadsheet 
    print(url)
    url = urlopen(url)
    fasta = url.read().decode("utf-8", "ignore") # import the url sequence 

    # identify which gene is being looked at 
    desc = fasta.split("\n",1)[0] # get the description of the sequence 
    gene = re.sub(".+GN=", "", desc) # find the gene
    gene = re.sub(r'\s* PE.*', '', gene) 
    print("The name of the gene being investigated is", gene)

    # get the protein sequence 
    fasta = fasta.split("\n",1)[1] # get only the protein sequence 
    fasta = "".join(fasta.split()) # remove spaces and lines 
    print("The sequence of the", gene, "protein is", fasta)
    print("The sequence is", len(fasta), "amino acids long.") 

    # identify the specific residue that is mutated 
    mut = line.iloc[:,1].values[0] # you need to provide the ID of the mutation (in the format R882H)
    nr_sub = mut[1:-1] # determine the number of the residue the mutation occurs in (remove the first - original, last - mutated)
    nr_sub = int(nr_sub) - 1 # python starts indexing at 0, so your n turns into n-1 (this will be used to identify the position in the fasta sequence)
    original = mut[0] # the reference sequence at this position
    ch_variant = mut[-1:] # what is the sequence at that site in the CH variant? (last entry e.g., if R882H, prints H)
    print("At position:", nr_sub+1, "the reference residue is", original, ".")
    print("In CH, this residue is mutated to", ch_variant)

    # check: is the aa in the hotspot position what it should be?
    if (fasta[nr_sub] == original) == True:
        print("All good, can proceed!")
    else: 
        print("Error! There is an issue in how the fasta sequence was processed")
        exit() # stop executing the script if the residue at the position does not match what is expected from the variant  

    # find aa +/- 10 from that residue (our max length of peptide is 11 and need to include the mutant residue)
    fasta_sub = fasta[nr_sub-10 : nr_sub+11]
    print("The sequence containing the original variant and the residues 10 aa away from it is:", fasta_sub)

    # set path to current directory 
    path = "/Users/barbarawalkowiak/Desktop/msc_thesis"
    os.chdir(path)
    
    # create new folder to save k-mers into if it does not exist yet
    if not os.path.exists("/Users/barbarawalkowiak/Desktop/msc_thesis/kmers/kmers_" + timestr): 
        os.mkdir("/Users/barbarawalkowiak/Desktop/msc_thesis/kmers/kmers_" + timestr)
    
    # generate kmers 
    kmers = []
    for i in range(8,12):
        ki = build_kmers(fasta_sub, i)
        kmers.append(ki)
    
    unlisted_kmers = [item for sublist in kmers for item in sublist]
    print(unlisted_kmers)
    with open('kmers/kmers_' + timestr + '/' + timestr + '_kmers_' + gene + '_' + mut + '_refseq.txt', 'w') as outfile: # refseq indicates I am using the "wild-type" sequence of the protein 
        outfile.write('\n'.join(unlisted_kmers))
    
    print("Completed writing reference k-mer sequences from", gene, "to a file")

    # modify the fasta_sub seqeunce so that it carries the MUTANT CH variant 
    fasta_ch = list(fasta_sub)

    # check 
    if (fasta_sub[10] == original) == True:
        print("All good, can proceed!")
    else: 
        print("Error! There is an issue in subsetting the fasta sequence")
        exit() # stop executing the script if the residue does not match what should be in the variant 

    # STOP CODON VARIANTS
        
    # I want to modify the kmers written depending on if there is a stop codon or not 
    if ch_variant == '*':

        # create all kmers as usual 
        fasta_ch[10] = str(ch_variant) # replace with variant identified in CH
        fasta_ch = ''.join(fasta_ch)
        print("The sequence containing the CH variant and the residues 10 aa away from it is:", fasta_ch)
        kmers = []
        
        for i in range(8,12):
            ki = build_kmers(fasta_ch, i)
            kmers.append(ki)
        
        unlisted_kmers = [item for sublist in kmers for item in sublist]
        
        kmers_stop = [seq for seq in unlisted_kmers if seq.endswith('*')]
        kmers_stop_aa = [seq.replace('*', '') for seq in kmers_stop if seq.endswith('*')]
        kmers_stop_aa = [seq for seq in kmers_stop_aa if len(seq) >= 8] # I don't want stuff shorter than 8 aa

        with open('kmers/kmers_' + timestr + '/' + timestr + '_kmers_' + gene + '_' + mut +'_ch_variant.txt', 'w') as outfile: # refseq indicates I am using the "wild-type" sequence of the protein 
            outfile.write('\n'.join(kmers_stop_aa))
        
        print("Completed writing ch-variant", mut, "k-mer sequences from", gene, "to a file")

        # note: this only outputs 3 sequences of length 8, 9, 10 amino acids ('back-tracked' from the end)
        # note: truncating variants are not expected to be immunogenic. However, we can use them as a sort of negative control where we do not expect a difference. 

    # SUBSTITUTIONS  
    else:
    
        fasta_ch[10] = str(ch_variant) # replace with variant identified in CH
        fasta_ch = ''.join(fasta_ch)
        print("The sequence containing the ch variant and the residues 10 aa away from it is:", fasta_ch)
        kmers = []
        
        for i in range(8,12):
            ki = build_kmers(fasta_ch, i)
            kmers.append(ki)
        
        unlisted_kmers = [item for sublist in kmers for item in sublist]
        print(unlisted_kmers)

        with open('kmers/kmers_' + timestr + '/' + timestr + '_kmers_' + gene + '_' + mut +'_ch_variant.txt', 'w') as outfile: # refseq indicates I am using the "wild-type" sequence of the protein 
            outfile.write('\n'.join(unlisted_kmers))
        
        print("Completed writing ch-variant", mut, "k-mer sequences from", gene, "to a file")
