# MHC_ClonalHaematopoiesis
Code associated with "No evidence of immunosurveillance in mutation-hotspot driven clonal haematopoiesis"

PACKAGE VERSIONS / SOFTWARE

- NetMHCpan-4.1
- NEtMHCIIpan-4.3
- PRIME2.0

- Jupyter Notebook 6.5.4
- Python 3.11.5
- Pandas 2.0.3
- Numpy 1.24.3
- Seaborn 0.13.2
- Matplotlib 3.7.2

SCRIPTS 
1. Generate kmers of proteins with CH variants of appropriate length.
2. Obtain predictions for kmer-MHC allele binding
3. Obtain datasframes with prediction scores
4. Obtain summary dataframes (choose best variant-MHC combination)
5. Obtain dataframe with MHC genotype and CH status for UKB participants
6. Add MHC binding scores to each participant examined
7. (A) Generate main figures (1-5), supplementary figures (1, 2, 4, 5, 7, 9, 10-13 (run with PRIME predictions), 18, 19), (B) Generate supplementary figures 14-17
8. Generate supplementary figure 6
9. Run and analyse power simulations (supplementary figures 3 and 8) 

SCRIPTS added in response to review comments 
- response_1_6_15_17.py - analysis relevant for comments 1 (age association), 6 (de novo variant calling QC), 15 (depth distribution), 17 (confounding by associations with sex / ethnicity)
- response_7.py - analysis relevant for comment 7 (germline contamination)
- response_9_1.ipynb and response_9_2.py - analysis relevant for comment 9 (analysis by quartiles). Dataframe with quartiles created in the script response_9_1.ipynb, analysis performed in the script response_9_2.py
- response_11.py - analysis relevant for comment 11 (top MHC allele from binding predictions)
  
