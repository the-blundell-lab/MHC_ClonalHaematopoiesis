# %%

# This script is intended to reproduce the main figures for analysis (1-3 + SF 1-3) for MHC I
# Input = dataframes with binding predictions and labels (top / bottom half of binding) 
# Output = figures (1-3) + supplementary figures (1-3)
# This script is based on jupyter notebook used to analyse predictions from NetMHC I (20240620_script3_mhci_netmhc.ipynb, version from 09/07/2024)

# EXAMPLE: To run python3 scripts02/7_main_figures_mhc2.py '/Users/barbarawalkowiak/Desktop/msc_thesis/results/20240624_netmhc2_scores_for_all_var.csv' '/Users/barbarawalkowiak/Desktop/msc_thesis/results/20240708_netmhc2_scores_for_all_var_with_labels.csv' 'netmhc_2reads'
# sys.arg[0] = name of the script to run
# sys.arg[1] = name of the file with scores 
# sys.arg[2] = name of the file with labels 
# sys.arg[3] = name of the file (what you want on figure titles etc so it should be descriptive enough to indicate software for predictions + nr of reads) 

# %%
# IMPORTS
import random
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, FuncFormatter, MaxNLocator
from matplotlib.font_manager import FontProperties
from matplotlib.lines import Line2D
from scipy.stats import ks_2samp 
import glob, os
import re
import pandas as pd
from decimal import *
import time 
import csv
import seaborn as sns 
import scipy.stats as stats
from scipy.stats import mannwhitneyu
from scipy.stats import chi2_contingency
import sys 

# %% 
# SET UP

# specify font for plotting 
plt.rcParams.update({'font.sans-serif':'Helvetica'})
# this is to make things editable in Illustrator
plt.rcParams['pdf.fonttype'] = 42 
# stop printing warnings 
import warnings
warnings.filterwarnings("ignore")
# get current date 
timestr = time.strftime("%Y%m%d") 

# COLORS 
# colors for relative binding 
col0r = '#0AAE37' # this is used to plot the top group (highest score = most immunogenic / strongest binding) (orange is binding?)
col1r = '#8FE8A7' # this is used to plot the middle group 
col2r = '#B018ED' # this is used to plot the bottom group (turquoise = does not bind)

# colors for absolute binding (I think we decided we want to do the same colors for the moment)
col0a = '#0AAE37' # this is used to plot the top group (strong binding based on threshold) 
col1a = '#8FE8A7' # this is used to plot the middle group (weak binding based on threshold)
col2a = '#B018ED' # this is used to plot the bottom group (no binding based on threshold)

# for CH specifically for absolute binding I am just using slightly darker colors but overall not changing it much
col0a1 = '#024F18' # this is used to plot the top group (strong binding based on threshold) for CH-neg in SF2 
col1a1 = '#078A2B' # this is used to plot the middle group (weak binding based on threshold) for CH-neg in SF2 
col2a1 = '#580479' # this is used to plot the bottom group (no binding based on threshold) for CH-neg in SF2 

# colors for CH-positive / CH-negative 
col_pos = '#BB0733' # color for CH-positive individuals (~reddish)
col_pos2 = '#fb9893'
col_neg = '#6892ED' # color for CH-negative individuals (blue)
col_neg2 = '#8AAAEF' # lighter color for background (dots)

# light grey for background (alternating variants)
col_background = '#EEEEEE'

# FONTS
title_font = 14
xaxis_font = 12
yaxis_font = 12
xticks_font = 9
yticks_font = 9
legend_title = 12
legend_font = 11
text_font = 9

# %% 
# IMPORT REQUIRED DATAFRAMES 

# # file with scores 
# netmhc2_df = pd.read_csv(sys.argv[1])
# # melted file with scores + labels 
# netmhc2_df_labels = pd.read_csv(sys.argv[2])
# # specify df name (for saving plots etc)
# df = sys.argv[3]

# # # file with scores 
netmhc2_df = pd.read_csv('/Users/barbarawalkowiak/Desktop/msc_thesis/results/dataframes/20240907_netmhc2_scores_for_all_var.csv')
# melted file with scores + labels 
netmhc2_df_labels = pd.read_csv('/Users/barbarawalkowiak/Desktop/msc_thesis/results/dataframes/20240907_netmhc2_scores_for_all_var_with_labels.csv')
# specify df name (for saving plots etc)
df = 'netmhc2_2reads'

# %%
# FORMATTING AND SETTING UP THE DATAFRAME

# only use variants for which there are enough carriers 
variants_identified = netmhc2_df_labels.gene_var.unique().tolist()
# make sure that you are only looking at variants for which you have more than 10 individuals
variants_identified = [var for var in variants_identified if netmhc2_df_labels[netmhc2_df_labels['gene_var'] == var].shape[0] >= 10]
variants_with_scores = netmhc2_df_labels.CH_variant.unique().tolist()
vars_with_no_carriers = [var for var in variants_with_scores if var not in variants_identified]
vars_with_carriers = [var for var in variants_with_scores if var in variants_identified]

# subset the dataframe 
netmhc2_df_labels_sub = netmhc2_df_labels[~netmhc2_df_labels['CH_variant'].isin(vars_with_no_carriers)]

# specify order of variants to plot (by median score)
order = netmhc2_df_labels_sub.sort_values(by = 'median_score', ascending = False).CH_variant.unique()
order2 = [var.replace('_', '\n') for var in order] # make it prettier to plot 

# specify levels for CH status (carrier, non-carrier)
netmhc2_df_labels_sub['CH_status'] = netmhc2_df_labels_sub['CH_status'].astype('category') # make sure you convert this to category first of all

# format the names of genetic variants so they are a bit nicer for plotting 
netmhc2_df_labels_sub['gene_var2'] = netmhc2_df_labels_sub['gene_var'].str.replace('_', '\n')
netmhc2_df_labels_sub['CH_variant2'] = netmhc2_df_labels_sub['CH_variant'].str.replace('_', '\n')

# rename CH_status column descriptions (we want CH-positive / negative to avoid associations with genetic carriers eg germline carriers)
netmhc2_df_labels_sub['CH_status2'] = np.where(netmhc2_df_labels_sub['CH_status']==1, 'CH-positive', 'CH-negative')
netmhc2_df_labels_sub['CH_status2'] = pd.Categorical(netmhc2_df_labels_sub['CH_status2'], categories = ['CH-positive', 'CH-negative']) # make sure you convert this to category first of all

# %%
# IDENTIFY VARIANTS WHICH CAN BE BOUND WELL AND POORLY (NON-UNIFORM ACROSS THE POPULATION)

# identify the variants values of which span different binding thresholds (ie some people bind it well, some people don't bind it etc)
# maybe choose variants where there are at least 5 carriers > -log10(1) and at leat 5 carriers is < -log10(5) (thresholds for MHC II)

variants_thresh = [var for var in vars_with_carriers if sorted(netmhc2_df_labels_sub[(netmhc2_df_labels_sub['variable']==f'score_{var}') & (netmhc2_df_labels_sub['gene_var']==var)].log_score.tolist())[-5] > -np.log10(1)]
variants_thresh = [var for var in variants_thresh if sorted(netmhc2_df_labels_sub[(netmhc2_df_labels_sub['variable']==f'score_{var}') & (netmhc2_df_labels_sub['gene_var']==var)].log_score.tolist())[4] < -np.log10(5)]
variants_thresh2 = [var.replace('_', '\n') for var in variants_thresh]

# %%
# CREATE NEW DATAFRAME TO STORE COUNT OF CARRIERS IN DIFFERENT GROUPS
df_counts_carriers = pd.DataFrame()

# loop for all variants
for var in order:

    # first, determine how many people are classified as top / bottom binding in general 
    var_df = netmhc2_df_labels_sub[netmhc2_df_labels_sub['CH_variant']==var] # df with scores for the variant for everyone
    var_df_carriers = var_df[var_df['gene_var']==var] # df with scores for carriers only
    
    n_all_top = var_df[var_df['group']=='top half'].shape[0]
    n_all_bottom = var_df[var_df['group']=='bottom half'].shape[0]
    n_all_total = var_df.shape[0]

    # count carriers in top binding group
    n_carriers_top = var_df_carriers[var_df_carriers['group']=='top half'].shape[0]
    n_carriers_bottom = var_df_carriers[var_df_carriers['group']=='bottom half'].shape[0]
    n_carriers_total = var_df_carriers.shape[0]

    # add to a new dataframe 
    df_counts_carriers = pd.concat([df_counts_carriers, 
        pd.DataFrame([var, n_all_top, n_all_bottom, n_all_total, n_carriers_top, n_carriers_bottom, n_carriers_total]).transpose()], axis = 0)

df_counts_carriers.columns = ['gene_var', 'n_all_top', 'n_all_bottom', 'n_all_total', 'n_carriers_top', 'n_carriers_bottom', 'n_carriers_total']
df_counts_carriers['p_carriers_top'] = df_counts_carriers['n_carriers_top'] / df_counts_carriers['n_all_top'] * 100
df_counts_carriers['p_carriers_bottom'] = df_counts_carriers['n_carriers_bottom'] / df_counts_carriers['n_all_bottom'] * 100

# MELT THE DF
# melt the dataframe to make this easier to plot 
df_counts_carriers_melted = pd.melt(df_counts_carriers, id_vars = 'gene_var')
df_counts_carriers_melted[['param', 'status', 'group']] = df_counts_carriers_melted.variable.str.split('_', expand = True)
df_counts_carriers_melted['group'] = df_counts_carriers_melted['group'] + ' half' # add more explicit name for plotting 

# specify order (by total nr of cases)
df_counts_carriers['gene_var2'] = df_counts_carriers['gene_var'].str.replace('_', '\n')
order_by_total = df_counts_carriers.sort_values(by = 'n_carriers_total', ascending = False)['gene_var'].tolist()
order_by_total2 = df_counts_carriers.sort_values(by = 'n_carriers_total', ascending = False)['gene_var2'].tolist()

# ERROR BARS 
# given that carriers are rare and independent of each other, it may make sense to use the Poisson distribution 
# then, standard error can be estimated as sqrt of the count (so I presume st error of carriers would be n carriers ± sqrt(n carriers))
# calculating std error as sqrt(count)
# that does return values which look more reasonable but not sure how I know this is correctly calculated at the end of the day
df_counts_carriers['std_error_count_top'] = np.sqrt(pd.to_numeric(df_counts_carriers['n_carriers_top']))
df_counts_carriers['std_error_count_bottom'] = np.sqrt(pd.to_numeric(df_counts_carriers['n_carriers_bottom']))

# FORMATTING
# format gene var name to something nicer
df_counts_carriers_melted['gene_var2'] = df_counts_carriers_melted['gene_var'].str.replace('_', '\n')
df_counts_carriers_melted['gene_var3'] = np.where(df_counts_carriers_melted['gene_var2'].isin(variants_thresh), '*' + df_counts_carriers_melted['gene_var2'], df_counts_carriers_melted['gene_var2'])
# subset dataframe to only have counts or percentages 
df_counts_carriers_melted_counts = df_counts_carriers_melted[df_counts_carriers_melted['param']=='n'] # only select counts 
df_counts_carriers_melted_percents = df_counts_carriers_melted[df_counts_carriers_melted['param']=='p'] # only select counts 

df_counts_carriers_melted_counts_tb = df_counts_carriers_melted_counts[df_counts_carriers_melted_counts['group'].isin(['top half', 'bottom half'])] # only compare top / bottom
df_counts_carriers_melted_percents_tb = df_counts_carriers_melted_percents[df_counts_carriers_melted_percents['group'].isin(['top half', 'bottom half'])] # only compare top / bottom

# %%
#############################################################################
#############################################################################
#############################################################################
# SUPPLEMENTARY FIGURE 1
# show the number of carriers of each variant
# this is just for me to help me keep track of who I have and how many people to expect
netmhc2_df_labels_sub_carriers = netmhc2_df_labels_sub[(~netmhc2_df_labels_sub['gene_var'].isna()) & (netmhc2_df_labels_sub['CH_status'] == 1)] # select CH-positive individuals
netmhc2_df_labels_carriers_counts = pd.DataFrame(netmhc2_df_labels_sub_carriers.gene_var2.value_counts()).reset_index()

plt.figure(figsize = (16,4))
sns.barplot(data = netmhc2_df_labels_carriers_counts, x = 'gene_var2', y = 'count', order = order_by_total2, color = col_pos, width = 0.6)
plt.xlim(-0.5, len(order_by_total)-0.5)
plt.ylim(0, 1700)
plt.xlabel(f'CH variant', fontsize = xaxis_font)
plt.ylabel('Number of CH-positive individuals', fontsize = yaxis_font)
plt.xticks(fontsize = xticks_font, rotation = 90)
plt.yticks(fontsize = yticks_font)

for i, var in enumerate(order_by_total2):
    dfx = netmhc2_df_labels_carriers_counts[netmhc2_df_labels_carriers_counts['gene_var2']==var]
    nr = dfx['count'].iloc[0]
    plt.text(i, nr+20, nr, ha = 'center', va = 'bottom', fontsize = text_font)
# this is not dependent 
plt.savefig(f'/Users/barbarawalkowiak/Desktop/msc_thesis/results/supplementary_files/{df}/sf1/sf1_nr_of_carriers_{df}.pdf', bbox_inches='tight')

# %% 
# may be worth combining this with what I have for NetMHC I
# # # file with scores 
netmhc1_df = pd.read_csv('/Users/barbarawalkowiak/Desktop/msc_thesis/results/dataframes/20240907_netmhc1_scores_for_all_var.csv')
# melted file with scores + labels 
netmhc1_df_labels = pd.read_csv('/Users/barbarawalkowiak/Desktop/msc_thesis/results/dataframes/20240907_netmhc1_scores_for_all_var_with_labels.csv')
# specify df name (for saving plots etc)
df1 = 'netmhc_2reads'

# subset the dataframe 
netmhc1_df_labels_sub = netmhc1_df_labels[~netmhc1_df_labels['CH_variant'].isin(vars_with_no_carriers)]

# specify levels for CH status (carrier, non-carrier)
netmhc1_df_labels_sub['CH_status'] = netmhc1_df_labels_sub['CH_status'].astype('category') # make sure you convert this to category first of all

# format the names of genetic variants so they are a bit nicer for plotting 
netmhc1_df_labels_sub['gene_var2'] = netmhc1_df_labels_sub['gene_var'].str.replace('_', '\n')
netmhc1_df_labels_sub['CH_variant2'] = netmhc1_df_labels_sub['CH_variant'].str.replace('_', '\n')

# rename CH_status column descriptions (we want CH-positive / negative to avoid associations with genetic carriers eg germline carriers)
netmhc1_df_labels_sub['CH_status2'] = np.where(netmhc1_df_labels_sub['CH_status']==1, 'CH-positive', 'CH-negative')
netmhc1_df_labels_sub['CH_status2'] = pd.Categorical(netmhc1_df_labels_sub['CH_status2'], categories = ['CH-positive', 'CH-negative']) # make sure you convert this to category first of all

netmhc1_df_labels_sub_carriers = netmhc1_df_labels_sub[(~netmhc1_df_labels_sub['gene_var'].isna()) & (netmhc1_df_labels_sub['CH_status'] == 1)] # select CH-positive individuals
netmhc1_df_labels_carriers_counts = pd.DataFrame(netmhc1_df_labels_sub_carriers.gene_var2.value_counts()).reset_index()

netmhc1_df_labels_carriers_counts['MHC'] = 'MHC I'
netmhc2_df_labels_carriers_counts['MHC'] = 'MHC II'

netmhc_df_labels_carriers_counts = pd.concat([netmhc1_df_labels_carriers_counts, netmhc2_df_labels_carriers_counts])

# %%
##### SUPPLEMENTARY FIG1 EDITION 2
colors = [col_pos, col_pos2]
plt.figure(figsize = (16,4))
sns.barplot(data = netmhc_df_labels_carriers_counts, x = 'gene_var2', y = 'count', hue = 'MHC', order = order_by_total2, palette = colors, width = 0.6)
plt.xlim(-0.5, len(order_by_total)-0.5)
plt.ylim(0, 1700)
plt.xlabel(f'CH variant', fontsize = xaxis_font)
plt.ylabel('Number of CH-positive individuals', fontsize = yaxis_font)
plt.xticks(fontsize = xticks_font, rotation = 90)
plt.yticks(fontsize = yticks_font)
plt.legend(frameon = False)
plt.savefig(f'/Users/barbarawalkowiak/Desktop/msc_thesis/results/supplementary_files/{df}/sf1/sf1_nr_of_carriers_{df}_mhci_ii.pdf', bbox_inches='tight')

# %%
#############################################################################
#############################################################################
#############################################################################
# FIGURE 1b
# NR OF CARRIERS VS PREDICTED RELATIVE BINDING TO MHC I (SPLIT INTO TOP / BOTTOM HALF)
# REQUIRES USE OF THE DATAFRAEME WITH LABELS

colors = [col0r, col2r]

plt.figure(figsize = (16, 4))
ax = sns.stripplot(data = df_counts_carriers_melted_counts_tb[df_counts_carriers_melted_counts_tb['status']=='carriers'], 
                   dodge = True, jitter = False, x = 'gene_var2', y = 'value', hue = 'group', palette = colors, size = 4.7, edgecolor = 'black', order = order_by_total2)

# GREY BACKGROUND (every other variant)
for i in range(1, len(order_by_total2), 2):
    ax.axvspan(i-0.5, i+0.5, color=col_background, alpha = 0.8)

# ADJUST HOW DODGED THE HUE IS (ie how separate are ppl who are better vs worse in MHC peptide binding)
# (Doing this to make sure error bars are placed correctly)
for i, artist in enumerate(ax.collections):
    # Get the current positions of the points
    offsets = artist.get_offsets()
    dodge_extent = 0
    offsets[:, 0] += (i % 2) * dodge_extent - dodge_extent / 2
    # Update the positions
    artist.set_offsets(offsets)

# ERROR BARS
for i, var in enumerate(order_by_total):

    small_df = df_counts_carriers[df_counts_carriers['gene_var'] == var]
    # extract values
    std_error_top = small_df['std_error_count_top'].iloc[0]
    std_error_bottom = small_df['std_error_count_bottom'].iloc[0]
    n_carriers_top = small_df['n_carriers_top'].iloc[0]
    n_carriers_bottom = small_df['n_carriers_bottom'].iloc[0]
    
    # add error bars 
    plt.errorbar(x=i-0.2, y=n_carriers_top, yerr=[[std_error_top], [std_error_top]], fmt='none', capsize = 0.5, color='black', capthick=0)
    plt.errorbar(x=i+0.21, y=n_carriers_bottom, yerr=[[std_error_bottom], [std_error_bottom]], fmt='none', capsize = 0.5, color='black', capthick=0)

# CHI SQUARED TEST 
p_values = {} # you want to save the p values and write them into a table
p_values_list = []

for i, var in enumerate(order_by_total):

    # find max value (nr carriers + standard error)
    small_df = df_counts_carriers[df_counts_carriers['gene_var'] == var]
    std_error_top = small_df['std_error_count_top'].iloc[0]
    std_error_top = small_df['std_error_count_top'].iloc[0]
    std_error_bottom = small_df['std_error_count_bottom'].iloc[0]
    n_carriers_top = small_df['n_carriers_top'].iloc[0]
    n_carriers_bottom = small_df['n_carriers_bottom'].iloc[0]
    max_value = max(n_carriers_top + std_error_top, n_carriers_bottom + std_error_bottom) # I need to have this value to determine where to place the indication of significance

    # determine significance (contingency table + chi-squared test)

    sel_df = df_counts_carriers_melted[df_counts_carriers_melted['gene_var'] == var]
    observed = [[sel_df[(sel_df['group'] == 'top half') & (sel_df['status'] == 'carriers')].value.iloc[0], netmhc2_df.shape[0]/2  - sel_df[(sel_df['group'] == 'top half') & (sel_df['status'] == 'carriers')].value.iloc[0]], 
                [sel_df[(sel_df['group'] == 'bottom half') & (sel_df['status'] == 'carriers')].value.iloc[0], netmhc2_df.shape[0]/2 - sel_df[(sel_df['group'] == 'bottom half') & (sel_df['status'] == 'carriers')].value.iloc[0]]]
    
    # Perform chi-square test
    chi2, p_value, dof, expected = chi2_contingency(observed)
    p_values_list.append(p_value) # save to dictionary
    p_values[var] = p_value

    significance = ''
    if p_value > 0.05 / len(order_by_total2): # correction for multiple testing
        significance = 'ns'
    elif p_value < 0.01 / len(order_by_total2):
        significance = '**'
    else:
        significance = '*'
    plt.text(i, max_value + 20, significance, ha='center', va='bottom', fontsize=text_font)

# axes labels 
plt.xlabel(f'CH variant', fontsize = xaxis_font)
plt.ylabel('Number of CH-positive individuals', fontsize = yaxis_font)

plt.xlim(-0.5, len(order_by_total2)-0.5)
plt.ylim(0, df_counts_carriers_melted[(df_counts_carriers_melted['gene_var'] == 'DNMT3A_R882H') & (df_counts_carriers_melted['status'] == 'carriers') & (df_counts_carriers_melted['group'].isin(['top half', 'bottom half']))].value.max() + 100) # I know this will be the max value for all variants bc R882H is most common

# BOLD X TICKS for which you have binding predictions that span both no binding and good binding
labels = []
for label in order_by_total2:
    if label in variants_thresh2:
        part1 = label.split('\n')[0]
        part2 = label.split('\n')[1]
        label = '$\mathbf{' + part1 + r'}$' + '\n' + '$\mathbf{' + part2 + r'}$'
        labels.append(label)
    else:
        labels.append(label)

ax.set_xticklabels(labels, fontsize = xticks_font, rotation = 90)

plt.yticks(fontsize = yticks_font)
plt.legend(['stronger binding\n(top half)', 'weaker binding\n(bottom half)'], 
           markerscale = 1.5, loc = 'upper right', fontsize = legend_font, frameon = False, handletextpad = 0.2)

# save the main figure 
plt.savefig(f'/Users/barbarawalkowiak/Desktop/msc_thesis/results/supplementary_files/{df}/fig1/1b_carriers_better_worse_binding_mhcii_{df}.pdf', bbox_inches='tight')

# %%
#############################################################################
#############################################################################
#############################################################################
# SUPPLEMENTARY FIGURE 1: P-VALUES (Fisher's exact test)
       
# p values need to be saved into supp table + plot distribution
p_values_df = pd.DataFrame([order_by_total, p_values_list]).T
p_values_df.columns = ['variant', 'p_value (not normalized)']
# sort values by how low the p value is (lowest p value at the top)
p_values_df = p_values_df.sort_values(by = 'p_value (not normalized)', ascending = True)
p_values_df.to_csv(f'/Users/barbarawalkowiak/Desktop/msc_thesis/results/supplementary_files/{df}/sf1/sf1b_carriers_better_worse_binding_p_value_{df}.csv')

# plot the distribution of p values and save to supplementary files 
plt.hist(p_values_list, bins = 5, edgecolor = 'black')
plt.xlabel('p-value', fontsize = xaxis_font)
plt.ylabel('Frequency', fontsize = yaxis_font)
plt.ylim(0, 15)
plt.title("Distribution of p-values from Fisher's exact test (Fig. 1B)")
plt.savefig(f'/Users/barbarawalkowiak/Desktop/msc_thesis/results/supplementary_files/{df}/sf1/sf1b_carriers_better_worse_binding_pval_distribution_{df}.pdf', bbox_inches='tight')

# %%
#############################################################################
#############################################################################
#############################################################################
# FIGURE 1c
# COMPARE PREDICTED BINDING B/N CH-POSITIVE AND CH-NEGATIVE INDIVIDUALS
# REQUIRES THE USE OF DATAFRAME WITH LABELS 

# identify carriers 
netmhc2_df_labels_sub_carriers = netmhc2_df_labels_sub[netmhc2_df_labels_sub['CH_status']==1]

# SAMPLE 2,000 individuals for CH-negative section
np.random.seed(1) # set seed to ensure reproducibility
sample_size = 2000 # select 2000-non-carriers individuals 

# for each variant, I need to sample 2,000 rows which contain scores for this variant
netmhc2_df_labels_sub_fig1c = pd.DataFrame() # empty dataframe to write rows to 

for var in variants_identified:

    netmhc2_df_labels_sub_noncarriers = netmhc2_df_labels_sub[netmhc2_df_labels_sub['gene_var']!=var] # anyone who does not carry this variant (they may carry a different one)
    df_scores = netmhc2_df_labels_sub_noncarriers[netmhc2_df_labels_sub_noncarriers['variable']==f'score_{var}'] # select dataframe with scores for the variant you are looking at (in people who are negative for this variant)
    
    # sample 2,000 rows (each for different individual) at random
    # replace = False so that each row is chosen only once 
    # random state is specified for reproducibilty
    netmhc2_df_labels_sub_sample = df_scores.sample(n = sample_size, replace = False, random_state = 10) 
    netmhc2_df_labels_sub_fig1c = pd.concat([netmhc2_df_labels_sub_fig1c, netmhc2_df_labels_sub_sample], axis = 0)

netmhc2_df_labels_sub_fig1c = pd.concat([netmhc2_df_labels_sub_fig1c, netmhc2_df_labels_sub_carriers], axis = 0)

# %%
#####################################################################
#####################################################################
#####################################################################
# PLOT 1c 
# plot for all variants together, with indicated median and binding thresholds
colors = [col_pos, col_neg2]

plt.figure(figsize=(16,4)) # set figure size
ax = sns.stripplot(y='log_score', x='CH_variant2', hue = 'CH_status2', data=netmhc2_df_labels_sub_fig1c, dodge = True, jitter = False, palette = colors, size = 1.5, legend = True, order = order_by_total2, alpha = 0.4)

# add grey background
for i in range(1, len(order_by_total2), 2):
    ax.axvspan(i - 0.5, i + 0.5, color=col_background, alpha=0.8)

# offset so that the middle of the median bar aligns with the datapoints
for i, artist in enumerate(ax.collections):
    # Get the current positions of the points
    offsets = artist.get_offsets()
    dodge_extent = 0.06
    offsets[:, 0] += (i % 2) * dodge_extent - dodge_extent / 2
    # Update the positions
    artist.set_offsets(offsets)

# find median score and add onto the plot 
for i, category in enumerate(order_by_total):
            
    median_carrier = netmhc2_df_labels_sub_fig1c[(netmhc2_df_labels_sub_fig1c['CH_variant'] == f'{category}') & (netmhc2_df_labels_sub_fig1c['CH_status'] == 1)].log_score.median()
    median_noncarrier = netmhc2_df_labels_sub_fig1c[(netmhc2_df_labels_sub_fig1c['CH_variant'] == f'{category}') & (netmhc2_df_labels_sub_fig1c['CH_status'] == 0)].log_score.median()

    # Plot text for each hue group
    plt.text(i, median_carrier, '-', ha='right', va='center', fontsize=30, fontweight='bold', color = col_pos)
    plt.text(i, median_noncarrier, '-', ha='left', va='center', fontsize=30, fontweight='bold', color = col_neg)

# add Mann-Whitney U test between groups (non-parametric t-test alternative)
p_values = []    

for i, category in enumerate(order_by_total):
    
    category_data = netmhc2_df_labels_sub_fig1c[netmhc2_df_labels_sub_fig1c['CH_variant'] == f'{category}']
    max_value = category_data['log_score'].max()
    
    scores_carrier = netmhc2_df_labels_sub_fig1c[(netmhc2_df_labels_sub_fig1c['CH_variant'] == f'{category}') & (netmhc2_df_labels_sub_fig1c['CH_status'] == 1)].log_score.tolist()
    scores_noncarrier = netmhc2_df_labels_sub_fig1c[(netmhc2_df_labels_sub_fig1c['CH_variant'] == f'{category}') & (netmhc2_df_labels_sub_fig1c['CH_status'] == 0)].log_score.tolist()
    
    statistic, p_value = mannwhitneyu(scores_carrier, scores_noncarrier)
    p_values.append(p_value)

    # adjust to the number of tests performed 
    significance = ''
    if p_value > (0.05 / len(order_by_total2)):
        significance = 'ns'
    elif p_value < (0.01 / len(order_by_total2)):
        significance = '**'
    else:
        significance = '*'
    plt.text(i, 0.2+max_value, significance, ha='center', va='bottom', fontsize=text_font)

# Add labels to the right of the markers

plt.text(0.6, -1*np.log10(1)+2, 'strong\nbinding', color='black', verticalalignment='center', fontsize = 11, ha='center')
plt.text(0.6, -1*np.log10(5)-2, 'no\nbinding', color='black', verticalalignment='center', fontsize = 11, ha='center')

# add lines to indicate thresholds for strong and weak binding
plt.axhline(-1*np.log10(1), color = '#3D3D3D', linewidth = 0.75, linestyle='--')
plt.axhline(-1*np.log10(5), color = '#3D3D3D', linewidth = 0.75, linestyle='--')

# axes labelling
plt.xlabel(f'CH variant', fontsize = xaxis_font)
plt.ylabel('MHC-variant binding', fontsize = yaxis_font)

labels = []
for label in order_by_total2:
    if label in variants_thresh2:
        part1 = label.split('\n')[0]
        part2 = label.split('\n')[1]
        label = '$\mathbf{' + part1 + r'}$' + '\n' + '$\mathbf{' + part2 + r'}$'
        labels.append(label)
    else:
        labels.append(label)
ax.set_xticklabels(labels, fontsize = xticks_font, rotation = 90)

plt.yticks(fontsize = yticks_font)

plt.xlim(-0.5, len(order_by_total2)-0.5) # we don't have scores which are this high / low but this gives a bit more space 
plt.ylim(-3.2, 2.5)

# specify legend 
legend = plt.legend(markerscale = 4.7, loc = 'upper right', fontsize = legend_font, frameon = False, handletextpad = 0.2)
for legend_handle in legend.legendHandles:
    legend_handle.set_alpha(1)
plt.savefig(f'/Users/barbarawalkowiak/Desktop/msc_thesis/results/supplementary_files/{df}/fig1/1c_scores_carriers_noncarriers_{df}.pdf', bbox_inches='tight')

# %%
#############################################################################
#############################################################################
#############################################################################
# SUPPLEMENTARY FIGURE 1: P-VALUES (Mann-Whitney U test)
       
# p values need to be saved into supp table + plot distribution
p_values_df = pd.DataFrame([order_by_total, p_values]).T
p_values_df.columns = ['variant', 'p_value (not normalized)']
# sort values by how low the p value is (lowest p value at the top)
p_values_df = p_values_df.sort_values(by = 'p_value (not normalized)', ascending = True)
p_values_df.to_csv(f'/Users/barbarawalkowiak/Desktop/msc_thesis/results/supplementary_files/{df}/sf1/sf1c_scores_carriers_noncarriers_p_value_{df}.csv')

# plot the distribution of p values and save to supplementary files 
plt.hist(p_values, bins = 5, edgecolor = 'black')
plt.xlabel('p-value', fontsize = xaxis_font)
plt.ylabel('Frequency', fontsize = yaxis_font)
plt.ylim(0, 15)
plt.title("Distribution of p-values from Mann-Whintey U test (Fig. 1C)")
plt.savefig(f'/Users/barbarawalkowiak/Desktop/msc_thesis/results/supplementary_files/{df}/sf1/sf1c_scores_carriers_noncarriers_pval_distribution_{df}.pdf', bbox_inches='tight')

# %%
#############################################################################
#############################################################################
#############################################################################
# FIGURE 2
# NR OF CARRIERS IN GROUPS PREDICTED ABSOLULTE BINDING TO MHC I (SPLIT INTO TOP / BOTTOM HALF)
# REQUIRES USE OF THE DATAFRAEME WITH LABELS

# SELECT DF + ASSIGN BINDING LABEL BASED ON ABSOLUTE SCORE
# first, select only variants which you know will span different binding thresholds
netmhc2_df_labels_sub_thresh = netmhc2_df_labels_sub[netmhc2_df_labels_sub['CH_variant2'].isin(variants_thresh2)]

# split into binding / no binding (based on threshold for weak binding)
conditions1 = [
    (netmhc2_df_labels_sub_thresh['log_score'] < -1 * np.log10(5)), # does not bind
    (netmhc2_df_labels_sub_thresh['log_score'] >= -1 * np.log10(5)) # higher = binds (at all)
]
values1 = ['no binding', 'binding']
netmhc2_df_labels_sub_thresh['any_binding'] = np.select(conditions1, values1)
netmhc2_df_labels_sub_thresh['any_binding']= pd.Categorical(netmhc2_df_labels_sub_thresh['any_binding'], categories=['no binding', 'binding'], ordered=True)

# split into binding / no binding (based on threshold for strong binding)
conditions2 = [
    (netmhc2_df_labels_sub_thresh['log_score'] < -1 * np.log10(1)), # does not bind
    (netmhc2_df_labels_sub_thresh['log_score'] >= -1 * np.log10(1)) # higher = binds strongly
]
values2 = ['no binding', 'binding']
netmhc2_df_labels_sub_thresh['strong_binding'] = np.select(conditions2, values2)
netmhc2_df_labels_sub_thresh['strong_binding']= pd.Categorical(netmhc2_df_labels_sub_thresh['strong_binding'], categories=['no binding', 'binding'], ordered=True)

# split into strong / weak / no binding
conditions3 = [
    (netmhc2_df_labels_sub_thresh['log_score'] < -1 * np.log10(5)), # does not bind
    ((netmhc2_df_labels_sub_thresh['log_score'] >= -1 * np.log10(5)) & (netmhc2_df_labels_sub_thresh['log_score'] < -1 * np.log10(1))), # binds weakly
    (netmhc2_df_labels_sub_thresh['log_score'] >= -1 * np.log10(1)) # binds strongly
]
values3 = ['no binding', 'weak binding', 'strong binding']
netmhc2_df_labels_sub_thresh['binding'] = np.select(conditions3, values3)
netmhc2_df_labels_sub_thresh['binding']= pd.Categorical(netmhc2_df_labels_sub_thresh['binding'], categories=['no binding', 'weak binding', 'strong binding'], ordered=True)

# %%
#############################################################################
#############################################################################
#############################################################################
### 1 ABSOLUTE BINDING (BINDING VS NON-BINDING, BASED ON THRESHOLD %EL RANK = 2)

# CREATE A NEW DF TO STORE COUNTS FOR ABSOLUTE BINDING (ANY)
df_counts_carriers_binding = pd.DataFrame()

# select variants which we established have enough individuals in binding / non-binding categories
thresh_order = [var for var in order if var in variants_thresh]

# loop for all variants 
for var in thresh_order:

    # first, determine how many people are classified as top / bottom binding in general 
    var_df = netmhc2_df_labels_sub_thresh[netmhc2_df_labels_sub_thresh['CH_variant']==var] # df with scores for the variant for everyone
    var_df_carriers = var_df[var_df['gene_var']==var] # df with scores for carriers only
    
    n_all_top = var_df[var_df['any_binding']=='binding'].shape[0]
    n_all_bottom = var_df[var_df['any_binding']=='no binding'].shape[0]
    n_all_total = var_df.shape[0]

    # count carriers in top binding group
    n_carriers_top = var_df_carriers[var_df_carriers['any_binding']=='binding'].shape[0]
    n_carriers_bottom = var_df_carriers[var_df_carriers['any_binding']=='no binding'].shape[0]
    n_carriers_total = var_df_carriers.shape[0]

    # add to a new dataframe 
    df_counts_carriers_binding = pd.concat([df_counts_carriers_binding, 
    pd.DataFrame([var, n_all_top, n_all_bottom, n_all_total, n_carriers_top, n_carriers_bottom, n_carriers_total]).transpose()], axis = 0)

df_counts_carriers_binding.columns = ['gene_var', 'n_all_binding', 'n_all_nonbinding', 'n_all_total', 'n_carriers_binding', 'n_carriers_nonbinding', 'n_carriers_total']

# calculate what proportion of individuals in each group is CH+
df_counts_carriers_binding['prop_carriers_binding'] = df_counts_carriers_binding['n_carriers_binding'] / df_counts_carriers_binding['n_all_binding']
df_counts_carriers_binding['prop_carriers_nonbinding'] = df_counts_carriers_binding['n_carriers_nonbinding'] / df_counts_carriers_binding['n_all_nonbinding']
# convert to percentage value for display 
df_counts_carriers_binding['p_carriers_binding'] = df_counts_carriers_binding['prop_carriers_binding'] * 100
df_counts_carriers_binding['p_carriers_nonbinding'] = df_counts_carriers_binding['prop_carriers_nonbinding'] * 100

# make sure everything is in num format 
df_counts_carriers_binding['prop_carriers_binding'] = pd.to_numeric(df_counts_carriers_binding['prop_carriers_binding'])
df_counts_carriers_binding['prop_carriers_nonbinding'] = pd.to_numeric(df_counts_carriers_binding['prop_carriers_nonbinding'])
df_counts_carriers_binding['n_all_binding'] = pd.to_numeric(df_counts_carriers_binding['n_all_binding'])
df_counts_carriers_binding['n_all_nonbinding'] = pd.to_numeric(df_counts_carriers_binding['n_all_nonbinding'])

# ERROR BARS
# the idea is that we have a binomial process so the variance is n * p * (1-p) where n = nr people, p = probability of success
# to get error bars we divide by the total nr of people which has no error so we do sqrt(np(1-p)) / n = sqrt(p(1-p)) / sqrt(n)
# this is what I found people do: https://www2.sjsu.edu/faculty/gerstman/StatPrimer/conf-prop.htm#:~:text=The%20standard%20error%20of%20a,in%20the%20population%20proportion%2C%20p.
df_counts_carriers_binding['stderror_carriers_binding'] = np.sqrt(df_counts_carriers_binding['prop_carriers_binding'] * (1 - df_counts_carriers_binding['prop_carriers_binding'])) / np.sqrt(df_counts_carriers_binding['n_all_binding']) 
df_counts_carriers_binding['stderror_carriers_nonbinding'] = np.sqrt(df_counts_carriers_binding['prop_carriers_nonbinding'] * (1 - df_counts_carriers_binding['prop_carriers_nonbinding'])) / np.sqrt(df_counts_carriers_binding['n_all_nonbinding']) 

# MELT THE DF
# melt the dataframe to make this easier to plot 
df_counts_carriers_binding_melted = pd.melt(df_counts_carriers_binding, id_vars = 'gene_var')
df_counts_carriers_binding_melted[['param', 'status', 'group']] = df_counts_carriers_binding_melted.variable.str.split('_', expand = True)

# format gene var name to something nicer
df_counts_carriers_binding_melted['gene_var2'] = df_counts_carriers_binding_melted['gene_var'].str.replace('_', '\n')

# subset dataframe to only have counts or percentages 
df_counts_carriers_binding_melted_counts = df_counts_carriers_binding_melted[df_counts_carriers_binding_melted['param']=='n'] # only select counts 
df_counts_carriers_binding_melted_props = df_counts_carriers_binding_melted[df_counts_carriers_binding_melted['param']=='prop'] # only select proportion
df_counts_carriers_binding_melted_percent = df_counts_carriers_binding_melted[df_counts_carriers_binding_melted['param']=='p'] # only select proportion
df_counts_carriers_binding_melted_stderr = df_counts_carriers_binding_melted[df_counts_carriers_binding_melted['param']=='stderror'] # only select standard error 

# specify order
order_by_total = df_counts_carriers_binding.sort_values(by = 'n_carriers_total', ascending = False)['gene_var'].tolist()
order_by_total2 = [o.replace('_', '\n') for o in order_by_total]

# %%
# CONTINGENCY TABLES (just in case, save this to somewhere as a csv file)
# create a dictionary to store contingency tables and p values 
cont_tables = {}
p_values = {}

# create contingency tables for each variant 
for i, var in enumerate(order_by_total):

    small_df = df_counts_carriers_binding_melted_counts[df_counts_carriers_binding_melted_counts['gene_var'] == var] # subset df for a particular variant
    observed = [[small_df[(small_df['group'] == 'binding') & (small_df['status'] == 'carriers')].value.iloc[0], 
                small_df[(small_df['group'] == 'binding') & (small_df['status'] == 'all')].value.iloc[0] - small_df[(small_df['group'] == 'binding') & (small_df['status'] == 'carriers')].value.iloc[0]], 
                [small_df[(small_df['group'] == 'nonbinding') & (small_df['status'] == 'carriers')].value.iloc[0], 
                small_df[(small_df['group'] == 'nonbinding') & (small_df['status'] == 'all')].value.iloc[0] - small_df[(small_df['group'] == 'nonbinding') & (small_df['status'] == 'carriers')].value.iloc[0]]]
    
    if np.all(observed): # np.all(observed) = True if all values are different to 0, otherwise does not make sense to do this
        chi2, p_value, dof, expected = chi2_contingency(observed)
        cont_tables[var] = expected
        p_values[var] = p_value

# save dictionaries to csv
with open(f'/Users/barbarawalkowiak/Desktop/msc_thesis/results/supplementary_files//{df}/sf2/2a_carriers_mhci_contingency_tables_anybinding2_{df}.csv', "w", newline="") as f:
    w = csv.DictWriter(f, cont_tables.keys())
    w.writeheader()
    w.writerow(cont_tables)

with open(f'/Users/barbarawalkowiak/Desktop/msc_thesis/results/supplementary_files/{df}/sf2/2a_carriers_mhci_p_values_anybinding2_{df}.csv', "w", newline="") as f:
    w = csv.DictWriter(f, p_values.keys())
    w.writeheader()
    w.writerow(p_values)

# %%
# PLOT results (for all variants, indicate significant)
# respecify genes (to now only include the ones for which it makes sense to do this)
df_counts_carriers_binding_melted_props.gene_var2 = df_counts_carriers_binding_melted_props.gene_var2.astype('category')

order_by_total3 = [o for o in order_by_total2 if o in df_counts_carriers_binding_melted_props.gene_var2.tolist()]
order_by_total31 = [o.replace('\n', '_') for o in order_by_total3]

# plot by order of total nr cases 
colors = [col0a, col2a]

# plot percentages 
plt.figure(figsize = (6, 4))
ax = sns.stripplot(data = df_counts_carriers_binding_melted_percent, x = 'gene_var2', y = 'value', hue = 'group', 
                   dodge = True, jitter = False, palette = colors, size = 6, edgecolor = 'black', order = order_by_total3)

# specify offsets
for i, artist in enumerate(ax.collections):
    # Get the current positions of the points
    offsets = artist.get_offsets()
    dodge_extent = 0.05
    offsets[:, 0] += (i % 2) * dodge_extent - dodge_extent / 2
    # Update the positions
    artist.set_offsets(offsets)

# adjust background (grey out every other variant)
for i in range(1, len(order_by_total31), 2):
    ax.axvspan(i-0.5, i+0.5, color=col_background, alpha = 0.8)

plt.xlabel(f'CH hotspot variant', fontsize = xaxis_font)
plt.ylabel('Individuals identified as CH-positive (%)', fontsize = yaxis_font)
plt.xticks(fontsize = xticks_font, rotation = 90)
plt.yticks(fontsize = yticks_font)

# Add standard errors
for i, var in enumerate(order_by_total31):

    # extract values of standard error 
    err_df = df_counts_carriers_binding_melted_stderr[df_counts_carriers_binding_melted_stderr['gene_var'] == var]
    std_error_bind = err_df[err_df['group'] == 'binding'].value.iloc[0] * 100
    std_error_nonbind = err_df[err_df['group'] == 'nonbinding'].value.iloc[0] * 100
    
    # get percentage values 
    p_df = df_counts_carriers_binding_melted_percent[df_counts_carriers_binding_melted_percent['gene_var']==var]
    p_bind = p_df[p_df['group'] == 'binding'].value.iloc[0]
    p_nonbind = p_df[p_df['group'] == 'nonbinding'].value.iloc[0]
    
    # add error bars 
    plt.errorbar(x=i-0.225, y=p_bind, yerr=[[std_error_bind], [std_error_bind]], fmt='none', capsize=1, capthick=0, color='black')
    plt.errorbar(x=i+0.225, y=p_nonbind, yerr=[[std_error_nonbind], [std_error_nonbind]], fmt='none', capsize=1, capthick=0, color='black')

    significance = p_values[var] 
    formatted_significance = f"{significance:.3f}"
    max_value = max(df_counts_carriers_binding_melted_percent[df_counts_carriers_binding_melted_percent['gene_var']==var].value) + max(std_error_bind, std_error_nonbind) + 0.01
    
    if significance < 0.05 / (len(order_by_total31) * 3):
        plt.text(i, max_value, f'$\mathbf{{{formatted_significance}}}$', 
                ha='center', va='bottom', fontsize=text_font)
    else:
        plt.text(i, max_value, f'{formatted_significance}', 
                ha='center', va='bottom', fontsize=text_font)

plt.legend(['binding\n(%EL rank < 5)', 'no binding\n(%EL rank > 5)'], markerscale = 1.2, loc = 'upper right', fontsize = legend_font, frameon = False, handletextpad = 0.2)
plt.xlim(-0.5, len(order_by_total31)-0.5)
plt.ylim(0, 1.45 * max(df_counts_carriers_binding_melted_percent.value))
if df == 'netmhc_2reads':
    plt.savefig(f'/Users/barbarawalkowiak/Desktop/msc_thesis/results/figures/fig2/2a_carriers_absolute_any_binding2_mhci_{df}_percentage.pdf', bbox_inches='tight')
else:
    plt.savefig(f'/Users/barbarawalkowiak/Desktop/msc_thesis/results/supplementary_files/{df}/fig2/2a_carriers_absolute_any_binding2_mhci_{df}_percentage.pdf', bbox_inches='tight')

# %%
#############################################################################
#############################################################################
#############################################################################
### 2 ABSOLUTE BINDING (BINDING VS NON-BINDING, BASED ON THRESHOLD %EL RANK = 1)

# CREATE A NEW DF TO STORE COUNTS FOR ABSOLUTE BINDING 
df_counts_carriers_binding2 = pd.DataFrame()

# select variants which we established have enough individuals in binding / non-binding categories
thresh_order = [var for var in order if var in variants_thresh]

# loop for all variants and fill in the df with values 
for var in thresh_order:

    # first, determine how many people are classified as top / bottom binding in general 
    var_df = netmhc2_df_labels_sub_thresh[netmhc2_df_labels_sub_thresh['CH_variant']==var] # df with scores for the variant for everyone
    var_df_carriers = var_df[var_df['gene_var']==var] # df with scores for carriers only
    
    n_all_top = var_df[var_df['strong_binding']=='binding'].shape[0]
    n_all_bottom = var_df[var_df['strong_binding']=='no binding'].shape[0]
    n_all_total = var_df.shape[0]

    # count carriers in top binding group
    n_carriers_top = var_df_carriers[var_df_carriers['strong_binding']=='binding'].shape[0]
    n_carriers_bottom = var_df_carriers[var_df_carriers['strong_binding']=='no binding'].shape[0]
    n_carriers_total = var_df_carriers.shape[0]

    # add to a new dataframe 
    df_counts_carriers_binding2 = pd.concat([df_counts_carriers_binding2, 
    pd.DataFrame([var, n_all_top, n_all_bottom, n_all_total, n_carriers_top, n_carriers_bottom, n_carriers_total]).transpose()], axis = 0)

# DF FORMATTING
# specify column names 
df_counts_carriers_binding2.columns = ['gene_var', 'n_all_binding', 'n_all_nonbinding', 'n_all_total', 'n_carriers_binding', 'n_carriers_nonbinding', 'n_carriers_total']
    
# calculate what proportion of people in each group are carriers 
df_counts_carriers_binding2['prop_carriers_binding'] = df_counts_carriers_binding2['n_carriers_binding'] / df_counts_carriers_binding2['n_all_binding']
df_counts_carriers_binding2['prop_carriers_nonbinding'] = df_counts_carriers_binding2['n_carriers_nonbinding'] / df_counts_carriers_binding2['n_all_nonbinding']

# get the percentage value 
df_counts_carriers_binding2['p_carriers_binding'] = df_counts_carriers_binding2['prop_carriers_binding'] * 100
df_counts_carriers_binding2['p_carriers_nonbinding'] = df_counts_carriers_binding2['prop_carriers_nonbinding'] * 100

# make sure values are numeric type
df_counts_carriers_binding2['prop_carriers_binding'] = pd.to_numeric(df_counts_carriers_binding2['prop_carriers_binding'])
df_counts_carriers_binding2['prop_carriers_nonbinding'] = pd.to_numeric(df_counts_carriers_binding2['prop_carriers_nonbinding'])
df_counts_carriers_binding2['n_all_binding'] = pd.to_numeric(df_counts_carriers_binding2['n_all_binding'])
df_counts_carriers_binding2['n_all_nonbinding'] = pd.to_numeric(df_counts_carriers_binding2['n_all_nonbinding'])

# ERROR BARS
# the idea is that we have a binomial process so the variance is n * p * (1-p) where n = nr people, p = probability of success
# to get error bars we divide by the total nr of people which has no error so we do sqrt(np(1-p)) / n = sqrt(p(1-p)) / sqrt(n)
# this is what I found people do: https://www2.sjsu.edu/faculty/gerstman/StatPrimer/conf-prop.htm#:~:text=The%20standard%20error%20of%20a,in%20the%20population%20proportion%2C%20p.
df_counts_carriers_binding2['stderror_carriers_binding'] = np.sqrt(df_counts_carriers_binding2['prop_carriers_binding'] * (1 - df_counts_carriers_binding2['prop_carriers_binding'])) / np.sqrt(df_counts_carriers_binding2['n_all_binding']) 
df_counts_carriers_binding2['stderror_carriers_nonbinding'] = np.sqrt(df_counts_carriers_binding2['prop_carriers_nonbinding'] * (1 - df_counts_carriers_binding2['prop_carriers_nonbinding'])) / np.sqrt(df_counts_carriers_binding2['n_all_nonbinding']) 

# MELT THE DATAFRAME 
# melt the dataframe to make this easier to plot 
df_counts_carriers_binding2_melted = pd.melt(df_counts_carriers_binding2, id_vars = 'gene_var')
df_counts_carriers_binding2_melted[['param', 'status', 'group']] = df_counts_carriers_binding2_melted.variable.str.split('_', expand = True)

# format gene var name to something nicer
df_counts_carriers_binding2_melted['gene_var2'] = df_counts_carriers_binding2_melted['gene_var'].str.replace('_', '\n')

# subset dataframe to only have counts or percentages 
df_counts_carriers_binding2_melted_counts = df_counts_carriers_binding2_melted[df_counts_carriers_binding2_melted['param']=='n'] # only select counts 
df_counts_carriers_binding2_melted_props = df_counts_carriers_binding2_melted[df_counts_carriers_binding2_melted['param']=='prop'] # only select proportion
df_counts_carriers_binding2_melted_percent = df_counts_carriers_binding2_melted[df_counts_carriers_binding2_melted['param']=='p'] # only select proportion
df_counts_carriers_binding2_melted_stderr = df_counts_carriers_binding2_melted[df_counts_carriers_binding2_melted['param']=='stderror'] # only select standard error 

# specify order
order_by_total = df_counts_carriers_binding2.sort_values(by = 'n_carriers_total', ascending = False)['gene_var'].tolist()
order_by_total2 = [o.replace('_', '\n') for o in order_by_total]

# %%
# CONTINGENCY TABLES (just in case, save this to somewhere as a csv file)
# create a dictionary to store contingency tables and p values 
cont_tables = {}
p_values = {}

# create contingency tables for each variant 
for i, var in enumerate(order_by_total):

    small_df = df_counts_carriers_binding2_melted_counts[df_counts_carriers_binding2_melted_counts['gene_var'] == var] # subset df for a particular variant
    observed = [[small_df[(small_df['group'] == 'binding') & (small_df['status'] == 'carriers')].value.iloc[0], 
                small_df[(small_df['group'] == 'binding') & (small_df['status'] == 'all')].value.iloc[0] - small_df[(small_df['group'] == 'binding') & (small_df['status'] == 'carriers')].value.iloc[0]], 
                [small_df[(small_df['group'] == 'nonbinding') & (small_df['status'] == 'carriers')].value.iloc[0], 
                small_df[(small_df['group'] == 'nonbinding') & (small_df['status'] == 'all')].value.iloc[0] - small_df[(small_df['group'] == 'nonbinding') & (small_df['status'] == 'carriers')].value.iloc[0]]]
    
    if np.all(observed): # np.all(observed) = True if all values are different to 0, otherwise does not make sense to do this
        chi2, p_value, dof, expected = chi2_contingency(observed)
        cont_tables[var] = expected
        p_values[var] = p_value

# save dictionaries to csv
with open(f'/Users/barbarawalkowiak/Desktop/msc_thesis/results/supplementary_files/{df}/sf2/2a_carriers_mhci_contingency_tables_anybinding05_{df}.csv', "w", newline="") as f:
    w = csv.DictWriter(f, cont_tables.keys())
    w.writeheader()
    w.writerow(cont_tables)

with open(f'/Users/barbarawalkowiak/Desktop/msc_thesis/results/supplementary_files/{df}/sf2/2a_carriers_mhci_p_values_anybinding05_{df}.csv', "w", newline="") as f:
    w = csv.DictWriter(f, p_values.keys())
    w.writeheader()
    w.writerow(p_values)

# %%
# PLOT results (for all variants, indicate significance from chi-squared test)
# respecify genes (to now only include the ones for which it makes sense to do this)
df_counts_carriers_binding2_melted_props.gene_var2 = df_counts_carriers_binding2_melted_props.gene_var2.astype('category')

# plot (percentage)
plt.figure(figsize = (6, 4))
ax = sns.stripplot(data = df_counts_carriers_binding2_melted_percent, x = 'gene_var2', y = 'value', hue = 'group', 
                   dodge = True, jitter = False, palette = colors, size = 6, edgecolor = 'black', order = order_by_total3)

# specify offsets
for i, artist in enumerate(ax.collections):
    # Get the current positions of the points
    offsets = artist.get_offsets()
    dodge_extent = 0.05
    offsets[:, 0] += (i % 2) * dodge_extent - dodge_extent / 2
    # Update the positions
    artist.set_offsets(offsets)

# adjust background (grey out every other variant)
for i in range(1, len(order_by_total31), 2):
    ax.axvspan(i-0.5, i+0.5, color=col_background, alpha = 0.8)

plt.xlabel(f'CH hotspot variant', fontsize = xaxis_font)
plt.ylabel('Individuals identified as CH-positive (%)', fontsize = yaxis_font)
plt.xticks(fontsize = xticks_font, rotation = 90)
plt.yticks(fontsize = yticks_font)

# Add standard errors
for i, var in enumerate(order_by_total31):

    # extract values of standard error 
    err_df = df_counts_carriers_binding2_melted_stderr[df_counts_carriers_binding2_melted_stderr['gene_var'] == var]
    std_error_bind = err_df[err_df['group'] == 'binding'].value.iloc[0] * 100
    std_error_nonbind = err_df[err_df['group'] == 'nonbinding'].value.iloc[0] * 100
    
    # get percentage values 
    p_df = df_counts_carriers_binding2_melted_percent[df_counts_carriers_binding2_melted_percent['gene_var']==var]
    p_bind = p_df[p_df['group'] == 'binding'].value.iloc[0]
    p_nonbind = p_df[p_df['group'] == 'nonbinding'].value.iloc[0]
    
    # add error bars 
    plt.errorbar(x=i-0.225, y=p_bind, yerr=[[std_error_bind], [std_error_bind]], fmt='none', capsize=1, capthick=0, color='black')
    plt.errorbar(x=i+0.225, y=p_nonbind, yerr=[[std_error_nonbind], [std_error_nonbind]], fmt='none', capsize=1, capthick=0, color='black')

    significance = p_values[var] 
    formatted_significance = f"{significance:.3f}"
    max_value = max(df_counts_carriers_binding2_melted_percent[df_counts_carriers_binding2_melted_percent['gene_var']==var].value) + max(std_error_bind, std_error_nonbind) + 0.01
    
    if significance < 0.05 / (len(order_by_total31) * 3):
        plt.text(i, max_value, f'$\mathbf{{{formatted_significance}}}$', 
                ha='center', va='bottom', fontsize=text_font)
    else:
        plt.text(i, max_value, f'{formatted_significance}', 
                ha='center', va='bottom', fontsize=text_font)

plt.legend(['binding\n(%EL rank < 1)', 'no binding\n(%EL rank > 1)'], markerscale = 1.2, loc = 'upper right', fontsize = legend_font, frameon = False, handletextpad = 0.2)
plt.xlim(-0.5, len(order_by_total31)-0.5)
plt.ylim(0, 1.45 * max(df_counts_carriers_binding_melted_percent.value))
plt.savefig(f'/Users/barbarawalkowiak/Desktop/msc_thesis/results/supplementary_files/{df}/fig2/2a_carriers_absolute_any_binding_mhci_{df}_percentage.pdf', bbox_inches='tight')

# %%
#############################################################################
#############################################################################
#############################################################################
### ABSOLUTE BINDING (STRONG / WEAK / NON-BINDING)

# CREATE EMPTY DF
# create a new dataframe where you will store counts of carriers in top and bottom binding groups
df_counts_carriers_binding3 = pd.DataFrame()

# loop for all variants
for var in thresh_order:

    # first, determine how many people are classified as top / bottom binding in general 
    var_df = netmhc2_df_labels_sub_thresh[netmhc2_df_labels_sub_thresh['CH_variant']==var] # df with scores for the variant for everyone
    var_df_carriers = var_df[var_df['gene_var']==var] # df with scores for carriers only
    
    n_all_top = var_df[var_df['binding']=='strong binding'].shape[0]
    n_all_mid = var_df[var_df['binding']=='weak binding'].shape[0]
    n_all_bottom = var_df[var_df['binding']=='no binding'].shape[0]
    n_all_total = var_df.shape[0]

    # count carriers in top binding group
    n_carriers_top = var_df_carriers[var_df_carriers['binding']=='strong binding'].shape[0]
    n_carriers_mid = var_df_carriers[var_df_carriers['binding']=='weak binding'].shape[0]
    n_carriers_bottom = var_df_carriers[var_df_carriers['binding']=='no binding'].shape[0]
    n_carriers_total = var_df_carriers.shape[0]

    # add to a new dataframe 
    df_counts_carriers_binding3 = pd.concat([df_counts_carriers_binding3, 
    pd.DataFrame([var, n_all_top, n_all_mid, n_all_bottom, n_all_total, n_carriers_top, n_carriers_mid, n_carriers_bottom, n_carriers_total]).transpose()], axis = 0)

df_counts_carriers_binding3.columns = ['gene_var', 'n_all_strong', 'n_all_weak', 'n_all_nonbinding', 'n_all_total', 'n_carriers_strong', 'n_carriers_weak', 'n_carriers_nonbinding', 'n_carriers_total']

# calculate the proportion of individuals in each group who are CH-positive
df_counts_carriers_binding3['prop_carriers_strong'] =  df_counts_carriers_binding3['n_carriers_strong'] / df_counts_carriers_binding3['n_all_strong']
df_counts_carriers_binding3['prop_carriers_weak'] = df_counts_carriers_binding3['n_carriers_weak'] / df_counts_carriers_binding3['n_all_weak']
df_counts_carriers_binding3['prop_carriers_nonbinding'] = df_counts_carriers_binding3['n_carriers_nonbinding'] / df_counts_carriers_binding3['n_all_nonbinding']

# get the percentage value 
df_counts_carriers_binding3['p_carriers_strong'] = df_counts_carriers_binding3['prop_carriers_strong'] * 100
df_counts_carriers_binding3['p_carriers_weak'] = df_counts_carriers_binding3['prop_carriers_weak'] * 100
df_counts_carriers_binding3['p_carriers_nonbinding'] = df_counts_carriers_binding3['prop_carriers_nonbinding'] * 100
 
# make sure columns are in numeric format 
df_counts_carriers_binding3['prop_carriers_strong'] = pd.to_numeric(df_counts_carriers_binding3['prop_carriers_strong'])
df_counts_carriers_binding3['prop_carriers_weak'] = pd.to_numeric(df_counts_carriers_binding3['prop_carriers_weak'])
df_counts_carriers_binding3['prop_carriers_nonbinding'] = pd.to_numeric(df_counts_carriers_binding3['prop_carriers_nonbinding'])
df_counts_carriers_binding3['n_all_strong'] = pd.to_numeric(df_counts_carriers_binding3['n_all_strong'])
df_counts_carriers_binding3['n_all_weak'] = pd.to_numeric(df_counts_carriers_binding3['n_all_weak'])
df_counts_carriers_binding3['n_all_nonbinding'] = pd.to_numeric(df_counts_carriers_binding3['n_all_nonbinding'])

df_counts_carriers_binding3['stderror_carriers_strong'] = np.sqrt(df_counts_carriers_binding3['prop_carriers_strong'] * (1 - df_counts_carriers_binding3['prop_carriers_strong']) / df_counts_carriers_binding3['n_all_strong']) 
df_counts_carriers_binding3['stderror_carriers_weak'] = np.sqrt(df_counts_carriers_binding3['prop_carriers_weak'] * (1 - df_counts_carriers_binding3['prop_carriers_weak']) / df_counts_carriers_binding3['n_all_weak']) 
df_counts_carriers_binding3['stderror_carriers_nonbinding'] = np.sqrt(df_counts_carriers_binding3['prop_carriers_nonbinding'] * (1 - df_counts_carriers_binding3['prop_carriers_nonbinding']) / df_counts_carriers_binding3['n_all_nonbinding']) 

# MELT DF
df_counts_carriers_binding3_melted = pd.melt(df_counts_carriers_binding3, id_vars = 'gene_var')
df_counts_carriers_binding3_melted[['param', 'status', 'group']] = df_counts_carriers_binding3_melted.variable.str.split('_', expand = True)

# FORMATTING
# format gene var name to something nicer
df_counts_carriers_binding3_melted['gene_var2'] = df_counts_carriers_binding2_melted['gene_var'].str.replace('_', '\n')

# subset dataframe to only have counts or percentages 
df_counts_carriers_binding3_melted_counts = df_counts_carriers_binding3_melted[df_counts_carriers_binding3_melted['param']=='n'] # only select counts 
df_counts_carriers_binding3_melted_percent = df_counts_carriers_binding3_melted[df_counts_carriers_binding3_melted['param']=='p'] # only select counts 
df_counts_carriers_binding3_melted_stderr = df_counts_carriers_binding3_melted[df_counts_carriers_binding3_melted['param']=='stderror'] # only select counts 

# specify order 
order_by_total = df_counts_carriers_binding2.sort_values(by = 'n_carriers_total', ascending = False)['gene_var'].tolist()
order_by_total2 = [o.replace('_', '\n') for o in order_by_total]

# %%
# CONTINGENCY TABLES
# run contingency table to see if there are any significant results  
# create a dictionary to store contingency tables and p values 
cont_tables = {}
p_values = {}

# create contingency tables for each variant 
for i, var in enumerate(order_by_total):

    small_df = df_counts_carriers_binding3_melted_counts[df_counts_carriers_binding3_melted_counts['gene_var'] == var] # subset df for a particular variant
    
    observed = [[small_df[(small_df['group'] == 'strong') & (small_df['status'] == 'carriers')].value.iloc[0], 
                small_df[(small_df['group'] == 'strong') & (small_df['status'] == 'all')].value.iloc[0] - small_df[(small_df['group'] == 'strong') & (small_df['status'] == 'carriers')].value.iloc[0]], 
                [small_df[(small_df['group'] == 'weak') & (small_df['status'] == 'carriers')].value.iloc[0], 
                small_df[(small_df['group'] == 'weak') & (small_df['status'] == 'all')].value.iloc[0] - small_df[(small_df['group'] == 'weak') & (small_df['status'] == 'carriers')].value.iloc[0]],
                [small_df[(small_df['group'] == 'nonbinding') & (small_df['status'] == 'carriers')].value.iloc[0], 
                small_df[(small_df['group'] == 'nonbinding') & (small_df['status'] == 'all')].value.iloc[0] - small_df[(small_df['group'] == 'nonbinding') & (small_df['status'] == 'carriers')].value.iloc[0]]]
    
    if np.all(observed): # np.all(observed) = True if all values are different to 0
        chi2, p_value, dof, expected = chi2_contingency(observed)
        cont_tables[var] = expected
        p_values[var] = p_value

with open(f'/Users/barbarawalkowiak/Desktop/msc_thesis/results/supplementary_files/{df}/sf2/2b_carriers_mhci_contingency_tables_swbinding_{df}.csv', "w", newline="") as f:
    w = csv.DictWriter(f, cont_tables.keys())
    w.writeheader()
    w.writerow(cont_tables)

with open(f'/Users/barbarawalkowiak/Desktop/msc_thesis/results/supplementary_files/{df}/sf2/2b_carriers_mhci_p_values_swbinding_{df}.csv', "w", newline="") as f:
    w = csv.DictWriter(f, p_values.keys())
    w.writeheader()
    w.writerow(p_values)

# %%
# PLOT (for all variants for which it makes sense to do this)
# re-specify categories in this dataframe
df_counts_carriers_binding3_melted_percent.gene_var2 = df_counts_carriers_binding3_melted_percent.gene_var.str.replace('_', '\n')
df_counts_carriers_binding3_melted_percent.gene_var2 = df_counts_carriers_binding3_melted_percent.gene_var2.astype('category')

# plot by order of total nr cases 
colors = [col0a, col1a, col2a]

plt.figure(figsize = (6,4))
ax = sns.stripplot(data = df_counts_carriers_binding3_melted_percent, x = 'gene_var2', y = 'value', hue = 'group', 
                  dodge = True, jitter = False, size = 6, palette = colors, edgecolor = 'black', order = order_by_total3)

# grey background for every other variant 
for i in range(1, len(order_by_total31), 2):
    ax.axvspan(i-0.5, i+0.5, color=col_background, alpha = 0.8)

plt.xlabel(f'CH hotspot variant', fontsize = xaxis_font)
plt.ylabel('Individuals identified as CH-positive (%)', fontsize = yaxis_font)
plt.xticks(fontsize = xticks_font, rotation = 90)
plt.yticks(fontsize = yticks_font)

# Add error bars
for i, var in enumerate(order_by_total31):

    # extract values of standard error 
    err_df = df_counts_carriers_binding3_melted_stderr[df_counts_carriers_binding3_melted_stderr['gene_var'] == var]
    std_error_strong = err_df[err_df['group'] == 'strong'].value.iloc[0] * 100
    std_error_weak = err_df[err_df['group'] == 'weak'].value.iloc[0] * 100
    std_error_nonbind = err_df[err_df['group'] == 'nonbinding'].value.iloc[0] * 100
    
    # get percentage values 
    p_df = df_counts_carriers_binding3_melted_percent[df_counts_carriers_binding3_melted_percent['gene_var']==var]
    p_strong = p_df[p_df['group'] == 'strong'].value.iloc[0]
    p_weak = p_df[p_df['group'] == 'weak'].value.iloc[0]
    p_nonbind = p_df[p_df['group'] == 'nonbinding'].value.iloc[0]
    
    # add error bars 
    plt.errorbar(x=i-0.266, y=p_strong, yerr=[[std_error_strong], [std_error_strong]], fmt='none', capsize=2, capthick = 0, color='black')
    plt.errorbar(x=i, y=p_weak, yerr=[[std_error_weak], [std_error_weak]], fmt='none', capsize=2, capthick = 0, color='black')
    plt.errorbar(x=i+0.266, y=p_nonbind, yerr=[[std_error_nonbind], [std_error_nonbind]], fmt='none', capsize=2, capthick = 0, color='black')

    significance = p_values[var] 
    formatted_significance = f"{significance:.3f}"
    max_value = max(df_counts_carriers_binding3_melted_percent[df_counts_carriers_binding3_melted_percent['gene_var']==var].value) + max(std_error_weak, std_error_strong, std_error_nonbind) + 0.01
    if significance < 0.05 / (len(order_by_total31) * 3):
        plt.text(i, max_value + 0.005, f'$\mathbf{{{formatted_significance}}}$', 
                ha='center', va='bottom', fontsize=text_font)
    else:
        plt.text(i, max_value + 0.005, f'{formatted_significance}', 
                ha='center', va='bottom', fontsize=text_font)

plt.legend(['strong binding\n(%EL rank < 1)', 'weak binding\n(1 < %EL rank < 5)', 'no binding\n(%EL rank > 5)'], 
           markerscale = 1.2, loc = 'upper right', fontsize = legend_font, frameon = False, handletextpad = 0.2)
plt.xlim(-0.5, len(order_by_total31)-0.5)
plt.ylim(0, 1.6 * max(df_counts_carriers_binding2_melted_percent.value))
if df == 'netmhc_2reads':
    plt.savefig(f'/Users/barbarawalkowiak/Desktop/msc_thesis/results/figures/fig2/2b_carriers_absolute_strong_weak_nonbinding_mhci_{df}.pdf', bbox_inches='tight')
else:
    plt.savefig(f'/Users/barbarawalkowiak/Desktop/msc_thesis/results/supplementary_files/{df}/fig2/2b_carriers_absolute_strong_weak_nonbinding_mhci_{df}.pdf', bbox_inches='tight')

# %% SUPPLEMENTARY FIGURE 2
######################################################
######################################################
######################################################

# FOR EACH VARIANT, SHOW % OF STRONG / WEAK / NO BINDING IN CARRIERS AND NON-CARRIERS
# REQUIRES THE USE OF DATAFRAME WITH LABELS 

# CREATE DF WITH DATA TO PLOT 
percentages_binding_thresholds = pd.DataFrame()

# create a list of all examined variants 
variants_examined = netmhc2_df_labels_sub.CH_variant.unique().tolist()

# calculate percentage values for each variant 
for var in variants_examined:

    # subset the main dataframe so only includes data for this variant 
    df_var = netmhc2_df_labels_sub[netmhc2_df_labels_sub['CH_variant']==var]
    
    # first get the numbers for everyone (carriers + non-carriers)
    N = df_var.shape[0] # number of people we have scores for (should be MHC-complete ppl we screened, and should be the same for each variant)
    
    # count nr of people who bind the variant strong / weak / not
    # NB value so I am not doing log-transformation
    n_strong_binding_all = df_var[df_var['value'] < 1].shape[0]
    n_weak_binding_all = df_var[(df_var['value'] < 5) & (df_var['value'] > 1)].shape[0]
    n_any_binding_all = df_var[df_var['value'] < 5].shape[0]
    n_no_binding_all = df_var[df_var['value'] >= 5].shape[0]
    
    p_strong_binding_all = n_strong_binding_all / N * 100
    p_weak_binding_all = n_weak_binding_all / N * 100
    p_any_binding_all = n_any_binding_all / N * 100
    p_no_binding_all = n_no_binding_all / N * 100

    # now just for carriers
    df_var_carriers = df_var[df_var['CH_status']==1]
    N = df_var_carriers.shape[0] # number of people we have scores for (should be MHC-complete ppl we screened, and should be the same for each variant)
    # count nr of people who bind the variant strong / weak / not
    n_strong_binding_carriers = df_var_carriers[df_var_carriers['value'] < 1].shape[0]
    n_weak_binding_carriers = df_var_carriers[(df_var_carriers['value'] < 5) & (df_var_carriers['value'] > 1)].shape[0]
    n_any_binding_carriers = df_var_carriers[df_var_carriers['value'] < 5].shape[0]
    n_no_binding_carriers = df_var_carriers[df_var_carriers['value'] >= 5].shape[0]
    
    p_strong_binding_carriers = n_strong_binding_carriers / N * 100
    p_weak_binding_carriers = n_weak_binding_carriers / N * 100
    p_any_binding_carriers = n_any_binding_carriers / N * 100
    p_no_binding_carriers = n_no_binding_carriers / N * 100

    # now just for non-carriers 
    df_var_ncarriers = df_var[df_var['CH_status']==0]
    N = df_var_ncarriers.shape[0] # number of people we have scores for (should be MHC-complete ppl we screened, and should be the same for each variant)
    # count nr of people who bind the variant strong / weak / not
    n_strong_binding_noncarriers = df_var_ncarriers[df_var_ncarriers['value'] < 1].shape[0]
    n_weak_binding_noncarriers = df_var_ncarriers[(df_var_ncarriers['value'] < 5) & (df_var_ncarriers['value'] > 1)].shape[0]
    n_any_binding_noncarriers = df_var_ncarriers[df_var_ncarriers['value'] < 5].shape[0]
    n_no_binding_noncarriers = df_var_ncarriers[df_var_ncarriers['value'] >= 5].shape[0]
    
    p_strong_binding_noncarriers = n_strong_binding_noncarriers / N * 100
    p_weak_binding_noncarriers = n_weak_binding_noncarriers / N * 100
    p_any_binding_noncarriers = n_any_binding_noncarriers / N * 100
    p_no_binding_noncarriers = n_no_binding_noncarriers / N * 100

    # write results to the empty dataframe you have created 
    percentages_binding_thresholds = pd.concat([percentages_binding_thresholds, 
        pd.DataFrame([var, n_strong_binding_all, n_weak_binding_all, n_any_binding_all, n_no_binding_all,
        p_strong_binding_all, p_weak_binding_all, p_any_binding_all, p_no_binding_all,
        n_strong_binding_carriers, n_weak_binding_carriers, n_any_binding_carriers, n_no_binding_carriers,
        p_strong_binding_carriers, p_weak_binding_carriers, p_any_binding_carriers, p_no_binding_carriers,
        n_strong_binding_noncarriers, n_weak_binding_noncarriers, n_any_binding_noncarriers, n_no_binding_noncarriers,
        p_strong_binding_noncarriers, p_weak_binding_noncarriers, p_any_binding_noncarriers, p_no_binding_noncarriers]).transpose()], axis = 0)

# FORMAT THE DATAFRAME TO MAKE IT SUITABLE FOR PLOTTING
# specify column names 
percentages_binding_thresholds.columns = ['CH_variant', 
        'n_strong_binding_all', 'n_weak_binding_all', 'n_any_binding_all', 'n_no_binding_all',
        'p_strong_binding_all', 'p_weak_binding_all', 'p_any_binding_all', 'p_no_binding_all',
        'n_strong_binding_carriers', 'n_weak_binding_carriers', 'n_any_binding_carriers', 'n_no_binding_carriers',
        'p_strong_binding_carriers', 'p_weak_binding_carriers', 'p_any_binding_carriers', 'p_no_binding_carriers',
        'n_strong_binding_noncarriers', 'n_weak_binding_noncarriers', 'n_any_binding_noncarriers', 'n_no_binding_noncarriers',
        'p_strong_binding_noncarriers', 'p_weak_binding_noncarriers', 'p_any_binding_noncarriers', 'p_no_binding_noncarriers']

# okay now what we want to do is to create stacked barplot where this is plotted side by side
percentages_binding_thresholds_melted = pd.melt(percentages_binding_thresholds, id_vars = 'CH_variant')

# split variable into 3 columns
percentages_binding_thresholds_melted[['np', 'level', 'level2', 'group']] = percentages_binding_thresholds_melted['variable'].str.split('_', expand=True)
# don't need the second (level2) column but this is just to be quicker with spliting this 

# select only percentage values (not absolute count values)
percentages_binding_thresholds_melted_percents = percentages_binding_thresholds_melted[percentages_binding_thresholds_melted['np'] == 'p']

# %%
# PLOT

colors1 = [col2a, col1a, col0a] # colors for CH-positive
colors2 = [col2a1, col1a1, col0a1] # colors for CH-positive
bar_width = 0.4

plt.figure(figsize=(16, 4))

for i, var in enumerate(order):

    percentage_var = percentages_binding_thresholds_melted_percents[percentages_binding_thresholds_melted_percents['CH_variant']==var]
    percentage_var_carriers = percentage_var[percentage_var['group']=='carriers']
    percentage_var_noncarriers = percentage_var[percentage_var['group']=='noncarriers']

    plt.bar(2*i, percentage_var_carriers[percentage_var_carriers['level']=='no']['value'], color=colors1[0], width=bar_width, label='no binding, CH-positive' if i == 0 else None)
    plt.bar(2*i+0.5, percentage_var_noncarriers[percentage_var_noncarriers['level']=='no']['value'], color=colors2[0], width=bar_width, label='no binding, CH-negative' if i == 0 else None)
    
    no_carriers = percentage_var_carriers[percentage_var_carriers['level']=='no']['value']
    no_noncarriers = percentage_var_noncarriers[percentage_var_noncarriers['level']=='no']['value']
    plt.bar(2*i, percentage_var_carriers[percentage_var_carriers['level']=='weak']['value'], color=colors1[1], width=bar_width, label='weak binding, CH-positive' if i == 0 else None, bottom=no_carriers)
    plt.bar(2*i+0.5, percentage_var_noncarriers[percentage_var_noncarriers['level']=='weak']['value'], color=colors2[1], width=bar_width, label='weak binding, CH-negative' if i == 0 else None, bottom=no_noncarriers)
    
    no_weak_carriers = percentage_var_carriers[percentage_var_carriers['level'].isin(['no', 'weak'])]['value'].sum()
    no_weak_noncarriers = percentage_var_noncarriers[percentage_var_noncarriers['level'].isin(['no', 'weak'])]['value'].sum()
    plt.bar(2*i, percentage_var_carriers[percentage_var_carriers['level']=='strong']['value'], color=colors1[2], width=bar_width, label='strong binding, CH-positive' if i == 0 else None, bottom=no_weak_carriers)
    plt.bar(2*i+0.5, percentage_var_noncarriers[percentage_var_noncarriers['level']=='strong']['value'], color=colors2[2], width=bar_width, label='strong binding, CH-negative' if i == 0 else None, bottom=no_weak_noncarriers)
    
plt.xlim(-0.5, 2*i+2)

plt.xlabel(f'CH variant', fontsize = xaxis_font)
plt.ylabel('Percentage (%)', fontsize = yaxis_font)
plt.legend(['no binding\n(%EL rank > 5)\nCH-positive', 'no binding\n(%EL rank > 5)\nCH-negative',
            'weak binding\n(1 < %EL rank < 5)\nCH-positive', 'weak binding\n(1 < %EL rank < 5)\nCH-negative',
            'strong binding\n(%EL rank < 1)\nCH-positive','strong binding\n(%EL rank < 1)\nCH-negative'], 
            markerscale = 1, loc = 'center left', bbox_to_anchor=(1, 0.5), fontsize = legend_font, frameon = False)

tick_positions = [i*2+0.25 for i, tick in enumerate(order2)]

# bold variants which you are looking at for absolute binding 
labels = []
for label in order2:
    if label in variants_thresh2:
        part1 = label.split('\n')[0]
        part2 = label.split('\n')[1]
        label = '$\mathbf{' + part1 + r'}$' + '\n' + '$\mathbf{' + part2 + r'}$'
        labels.append(label)
    else:
        labels.append(label)

plt.xticks(tick_positions, labels, rotation=90, fontsize = xticks_font)
plt.yticks(fontsize = yticks_font)
plt.tick_params(left = True, bottom = False) 

plt.savefig(f'/Users/barbarawalkowiak/Desktop/msc_thesis/results/supplementary_files/{df}/sf2/2_scores_carriers_noncarriers_percentages_mhci_{df}.pdf', bbox_inches='tight')

# %% 
#################################################################################
#################################################################################
#################################################################################
# FIGURE 3
# VAF (IN CARRIERS) VS RELATIVE BINDING CAPACITY TO MHC I
# REQUIRES USE OF THE DATAFRAME WITH LABELS
    
# %%
# DEFINE FUNCTIONS 
# FUNCTION TO SPLIT INTO GROUPS
def split_into_equal_groups(df, num_groups):

    # first, assign a rank to everyone based on score
    # the highest score = the lowest rank (ie highest score means you get rank 1)
    # if two people have the same score, assign consecutive ranks 
    df['rank'] = df['log_score'].rank(method='first', ascending=False)
    
    # assign groups based on rank 

    # first, determine the number of samples
    total_samples = len(df)

    # now, determine the number of samples in each group
    samples_per_group = total_samples // num_groups
    remainder = total_samples % num_groups
    group_sizes = [samples_per_group + 1 if i < remainder else samples_per_group for i in range(num_groups)] # add one person if there is a reminder 
    
    # now, assign the group based on the rank 
    df_sort = df.sort_values(by = 'rank')
    group_assignments = []
    group_number = 1
    start = 0
    for size in group_sizes:
        end = start + size
        group_assignments.extend([group_number] * (end - start))
        start = end
        group_number += 1

    # Add a new column 'group' to the DataFrame indicating the group assignment for each row
    df_sort[f'group_{num_groups}'] = group_assignments
    print(df_sort)
    
    return df_sort

# %%
# PLOT CUMULATIVE DISTRIBUTION
# define function to fix formatting of the y axis ticks 
def custom_formatter(y, pos):
    if y == 0:
        return '0'
    elif y >= 1:
        return f'{y:.1f}'.rstrip('0').rstrip('.')
    else:
        return f'{y:.2g}'

# PLOT CUMULATIVE DISTRIBUTION (log-log plot)
def plot_cumulative(data, x, y, hue, colors, var, name, df, nr_groups, style, ks = None, legend = None):

    var_name = var.split('_')[0:2]
    var_name = ' '.join(var_name)
    part1 = var_name.split(' ')[0]
    part2 = var_name.split(' ')[1]
    var_name = '$\mathbf{' + part1 + r'}$' + '\n' + '$\mathbf{' + part2 + r'}$'
    
    if var == 'ALL_VARIANTS':
        nr_cases = data.gene_var.value_counts().sum()
    else:
        nr_cases = data.shape[0]

    if var == 'ALL_VARIANTS':
        plot_width = 3.3
        plot_height = 3.3
        left_margin = 0.8
        right_margin = 0.8
        top_margin = 0.8
        bottom_margin = 0.8
    else:
        plot_width = 1.35
        plot_height = 1.35
        left_margin = 0.35
        right_margin = 0.35
        top_margin = 0.35
        bottom_margin = 0.35

    # Calculate the figure size
    fig_width = plot_width + left_margin + right_margin
    fig_height = plot_height + top_margin + bottom_margin
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    plt.subplots_adjust(left=left_margin/fig_width,
                    right=1-right_margin/fig_width,
                    top=1-top_margin/fig_height,
                    bottom=bottom_margin/fig_height)

    ax = sns.lineplot(x=x, y=y, data=data, hue=hue, palette=colors, legend = False)
    scatter = sns.scatterplot(x=x, y=y, data=data, hue=hue, palette=colors, legend = True, size = 11, linewidth=0)

    # specify title and axes labels 
    if style == 'log':
        if var == 'ALL_VARIANTS':
            plt.text(0.8, 0.9, f'{var_name}\n(n={nr_cases})', transform = plt.gca().transAxes,
                fontsize = 10, verticalalignment = 'center', horizontalalignment = 'center')
        else:
            plt.text(0.8, 0.85, f'{var_name}\n(n={nr_cases})', transform = plt.gca().transAxes,
                fontsize = 6, verticalalignment = 'center', horizontalalignment = 'center')
    else:
        if var == 'ALL_VARIANTS':
            plt.text(0.8, 0.85, f'{var_name}\n(n={nr_cases})', transform = plt.gca().transAxes,
                fontsize = 10, verticalalignment = 'center', horizontalalignment = 'center')
        else:
            plt.text(0.8, 0.9, f'{var_name}\n(n={nr_cases})', transform = plt.gca().transAxes,
                fontsize = 6, verticalalignment = 'center', horizontalalignment = 'center')
    
    plt.xlabel('Variant allele frequency (%)', fontsize = xaxis_font-5)
    plt.ylabel('Reverse cumulative probability', fontsize = yaxis_font-5)

    # convert both axes to log scale (if log-log)
    if style == 'log':
        plt.gca().set_xscale('log')
        plt.gca().set_yscale('log')

        x_ticks = [2, 3, 10, 30, 100]
        plt.xticks(x_ticks, x_ticks, fontsize = xticks_font-3)

    # add results of K-test (for log-log)
    if style == 'log':
        if ks == None:
            plt.text(0.8, 0.55, ' ') # don't add any text if KS test not run
        else:
            if var == 'ALL_VARIANTS':
                plt.text(0.8, 0.8, f'p={ks}', transform = plt.gca().transAxes,
                fontsize = 9, verticalalignment = 'center', horizontalalignment = 'center')
            else:
                plt.text(0.2, 0.15, f'p={ks}', transform = plt.gca().transAxes,
                fontsize = 6, verticalalignment = 'center', horizontalalignment = 'center')
    else:
        if ks == None:
            plt.text(0.8, 0.55, ' ') # don't add any text if KS test not run
        else:
            if var == 'ALL_VARIANTS':
                plt.text(0.8, 0.55, f'p={ks}', transform = plt.gca().transAxes,
                fontsize = 8.5, verticalalignment = 'center', horizontalalignment = 'center')
            else:
                plt.text(0.8, 0.65, f'p={ks}', transform = plt.gca().transAxes,
                fontsize = 6, verticalalignment = 'center', horizontalalignment = 'center')

    # specify y axis ticks
    scatter.yaxis.set_major_formatter(FuncFormatter(custom_formatter))
    scatter.tick_params(axis = 'y', labelsize = yticks_font-3) 
    
    # specify limits on both axes 
    if style == 'log':
        if min(data[y].tolist()) > 0.1:
            plt.xlim(2, 100)
            plt.ylim(0.05, 1.2) # show where 0.1 is for each plot
        else:
            plt.xlim(2, 100)
            plt.ylim(None, 1.2)

    elif style == 'linear':
        plt.xlim(2, None)
        plt.ylim(-0.05, 1.05)

# add legend (not on all plots though)
    if legend == True:
        if style == 'log':
            if nr_groups == 3:
                handles, labels = scatter.get_legend_handles_labels()
                plt.legend(labels = ['stronger binding\n(top 1/3)', 'mid binding\n(mid 1/3)', 'weaker binding\n(bottom 1/3)'],
                    handles = handles[:-1], loc='lower left', frameon = False,
                    fontsize = text_font, markerscale = 1.2, handletextpad = 0.2)
            elif nr_groups == 2:
                handles, labels = scatter.get_legend_handles_labels()
                plt.legend(labels = ['stronger binding\n(top half)', 'weaker binding\n(bottom half)'],
                    handles = handles[:-1], loc='lower left', frameon = False,
                    fontsize = text_font, markerscale = 1.2, handletextpad = 0.2)
            else: 
                print('wrong number of groups')
        elif style == 'linear':
            if nr_groups == 3:
                handles, labels = scatter.get_legend_handles_labels()
                plt.legend(labels = ['stronger binding\n(top 1/3)', 'mid binding\n(mid 1/3)', 'weaker binding\n(bottom 1/3)'],
                    handles = handles[:-1], loc='lower right', frameon = False,
                    fontsize = text_font, markerscale = 1.2, handletextpad = 0.2)
            elif nr_groups == 2:
                handles, labels = scatter.get_legend_handles_labels()
                plt.legend(labels = ['stronger binding\n(top half)', 'weaker binding\n(bottom half)'],
                    handles = handles[:-1], loc='lower right', frameon = False,
                    fontsize = text_font, markerscale = 1.2, handletextpad = 0.2)
            else: 
                print('wrong number of groups')
    else:
        plt.legend('',frameon=False) # explicitly do NOT include the legend
    
    if df == 'netmhc_2reads':
        plt.savefig(f'/Users/barbarawalkowiak/Desktop/msc_thesis/results/figures/fig3/3a_cumulative_plot_{var}_mhci_{name}_{df}_{nr_groups}_groups_{style}.pdf', bbox_inches='tight')
    else:
        plt.savefig(f'/Users/barbarawalkowiak/Desktop/msc_thesis/results/supplementary_files/{df}/fig3/3a_cumulative_plot_{var}_mhci_{name}_{df}_{nr_groups}_groups_{style}.pdf', bbox_inches='tight')
    
# %%
#####################################################################
# first subset for carriers data only and only keep scores for variant that is carried by the specific person so filter this out
netmhc2_df_labels_sub_carriers = netmhc2_df_labels_sub[netmhc2_df_labels_sub['CH_variant'] == netmhc2_df_labels_sub['gene_var']].iloc[:, 1:]
# NOTE: there are 8,982 CH-positive individuals in this dataset
# however, two of them have two variants (hence, there are 8,984 rows)
# one individual (Person ID = 4315153) has DNMT3A_R882H and TP53_R175H
# the other individual (Person ID = 2822684) has TP53_R273H and DNMT3A_R882S

#####################################################################
# SPLIT CARRIERS INTO TWO GROUPS (better vs worse binding)
netmhc2_df_labels_sub_groups2_eqsize = netmhc2_df_labels_sub_carriers.groupby('gene_var').apply(lambda x: split_into_equal_groups(x, num_groups=2))
netmhc2_df_labels_sub_groups2_eqsize['group_name'] =  netmhc2_df_labels_sub_groups2_eqsize['group_2'].map({1: 'top', 2: 'bottom'})
netmhc2_df_labels_sub_groups2_eqsize = netmhc2_df_labels_sub_groups2_eqsize.droplevel(0)
netmhc2_df_labels_sub_groups2_eqsize = netmhc2_df_labels_sub_groups2_eqsize.reset_index()
netmhc2_df_labels_sub_groups2_eqsize_sorted = netmhc2_df_labels_sub_groups2_eqsize.sort_values(by = 'rank', ascending = False)

# SPLIT CARRIERS INTO THREE GROUPS (better / medium / worse binding)
netmhc2_df_labels_sub_groups3_eqsize = netmhc2_df_labels_sub_carriers.groupby('gene_var').apply(lambda x: split_into_equal_groups(x, num_groups=3))
netmhc2_df_labels_sub_groups3_eqsize['group_name'] =  netmhc2_df_labels_sub_groups3_eqsize['group_3'].map({1: 'top', 2: 'mid', 3: 'bottom'})
netmhc2_df_labels_sub_groups3_eqsize = netmhc2_df_labels_sub_groups3_eqsize.droplevel(0)
netmhc2_df_labels_sub_groups3_eqsize = netmhc2_df_labels_sub_groups3_eqsize.reset_index()
netmhc2_df_labels_sub_groups3_eqsize_sorted = netmhc2_df_labels_sub_groups3_eqsize.sort_values(by = 'rank', ascending = False)

# this looks fine to me: no changes in the nr of CH-positive individuals identified 
# save dataframes to check that you split carriers correctly, and have equal (±1 nr of people in both groups)
dfx_check2 = netmhc2_df_labels_sub_groups2_eqsize_sorted.groupby('gene_var').group_2.value_counts()
dfx_check3 = netmhc2_df_labels_sub_groups3_eqsize_sorted.groupby('gene_var').group_3.value_counts()
dfx_check2.to_csv(f'/Users/barbarawalkowiak/Desktop/msc_thesis/results/supplementary_files/{df}/sf3/check_splitting_carriers_2groups.csv')
dfx_check3.to_csv(f'/Users/barbarawalkowiak/Desktop/msc_thesis/results/supplementary_files/{df}/sf3/check_splitting_carriers_3groups.csv')

# %%
#########################################################################
# PLOT RESULTS OF SPLITTING (supplementary figure)

# SPLITTING INTO 2 GROUPS
colors2 = [col0r, col2r]
order_by_total = df_counts_carriers.sort_values(by = 'n_carriers_total', ascending = False)['gene_var'].tolist()
order_by_total2 = df_counts_carriers.sort_values(by = 'n_carriers_total', ascending = False)['gene_var2'].tolist()

# PLOT
plt.figure(figsize = (16, 4))
ax = sns.stripplot(data=netmhc2_df_labels_sub_groups2_eqsize_sorted, x='gene_var2', y='log_score', jitter = False, hue='group_2', palette=colors2, order=order_by_total2, size = 4, alpha = 0.2)
# specify offsets
for i, artist in enumerate(ax.collections):
    offsets = artist.get_offsets()
    dodge_extent = 0
    offsets[:, 0] += (i % 3) * dodge_extent - dodge_extent / 3
    artist.set_offsets(offsets)
# grey background for every other variant 
for i in range(1, len(order), 2):
    ax.axvspan(i-0.5, i+0.5, color=col_background, alpha = 0.8)
# labels
plt.xlabel(f'CH variant', fontsize = xaxis_font)
plt.ylabel('MHC-variant binding', fontsize = yaxis_font)
labels = []
for label in order_by_total2:
    if label in variants_thresh2:
        part1 = label.split('\n')[0]
        part2 = label.split('\n')[1]
        label = '$\mathbf{' + part1 + r'}$' + '\n' + '$\mathbf{' + part2 + r'}$'
        labels.append(label)
    else:
        labels.append(label)
ax.set_xticklabels(labels, fontsize = xticks_font, rotation = 90)
plt.yticks(fontsize = yticks_font)
plt.ylim(-3.5, 2.5)
plt.xlim(-0.5, len(order)-0.5)
# legend 
plt.legend(['stronger binding\n(top half)', 'weaker binding\n(bottom half)'], handletextpad = 0.2,
            loc = 'upper right', fontsize = legend_font, markerscale = 2, frameon = False)
legend = plt.gca().get_legend()
for handle, color in zip(legend.legendHandles, colors2):
    handle.set_color(color)
for lh in legend.legendHandles:
    lh.set_alpha(1)
# save plot 
plt.savefig(f'/Users/barbarawalkowiak/Desktop/msc_thesis/results/supplementary_files/{df}/sf3/sf3a_carriers_split_by_score_mhci_2groups_{df}.pdf', bbox_inches='tight')

# SPLITTING INTO 3 GROUPS
colors3 = [col0r, col1r, col2r]

plt.figure(figsize = (16, 4))
ax = sns.stripplot(data=netmhc2_df_labels_sub_groups3_eqsize_sorted, x='gene_var2', y='log_score', jitter = False, hue='group_3', palette=colors3, order=order_by_total2, size = 4, alpha = 0.2)
# specify offsets
for i, artist in enumerate(ax.collections):
    # Get the current positions of the points
    offsets = artist.get_offsets()
    dodge_extent = 0
    offsets[:, 0] += (i % 3) * dodge_extent - dodge_extent / 3
    # Update the positions
    artist.set_offsets(offsets)
# grey background for every other variant 
for i in range(1, len(order), 2):
    ax.axvspan(i-0.5, i+0.5, color=col_background, alpha = 0.8)
# labels
plt.xlabel(f'CH hotspot variant', fontsize = xaxis_font)
plt.ylabel('MHC-variant binding', fontsize = yaxis_font)
labels = []
for label in order_by_total2:
    if label in variants_thresh2:
        part1 = label.split('\n')[0]
        part2 = label.split('\n')[1]
        label = '$\mathbf{' + part1 + r'}$' + '\n' + '$\mathbf{' + part2 + r'}$'
        labels.append(label)
    else:
        labels.append(label)
ax.set_xticklabels(labels, fontsize = xticks_font, rotation = 90)
plt.yticks(fontsize = yticks_font)
plt.ylim(-3.5, 2.5)
plt.xlim(-0.5, len(order)-0.5)
# legend formatting 
plt.legend(['stronger binding', 'mid binding', 'weaker binding'], handletextpad = 0.2,
            loc = 'upper right', fontsize = legend_font, markerscale = 2, frameon = False)
legend = plt.gca().get_legend()
for handle, color in zip(legend.legendHandles, colors3):
    handle.set_color(color)
for lh in legend.legendHandles:
    lh.set_alpha(1)
# save plot 
plt.savefig(f'/Users/barbarawalkowiak/Desktop/msc_thesis/results/supplementary_files/{df}/sf3/sf3a_carriers_split_by_score_mhci_3groups_{df}.pdf', bbox_inches='tight')

# %%
#########################################################################
#########################################################################
#########################################################################
#########################################################################

# CUMULATIVE DISTRIBUTION (2 groups)
# add VAF data 
netmhc2_df_labels_sub_groups2_eqsize_sorted['VAF_percent'] = netmhc2_df_labels_sub_groups2_eqsize_sorted['VAF'] * 100
netmhc2_df_labels_sub_groups2_eqsize_sorted['group_name'] = netmhc2_df_labels_sub_groups2_eqsize_sorted['group_name'].astype('category') # reorder group names
netmhc2_df_labels_sub_groups2_eqsize_sorted['group_name'].cat.reorder_categories(['top', 'bottom'])

# CUMULATIVE DISTRIBUTION (3 groups)
netmhc2_df_labels_sub_groups3_eqsize_sorted['VAF_percent'] = netmhc2_df_labels_sub_groups3_eqsize_sorted['VAF'] * 100
netmhc2_df_labels_sub_groups3_eqsize_sorted['group_name'] = netmhc2_df_labels_sub_groups3_eqsize_sorted['group_name'].astype('category') # reorder group names
netmhc2_df_labels_sub_groups3_eqsize_sorted['group_name'].cat.reorder_categories(['top', 'mid', 'bottom'])

# %%
# CUMULATIVE DISTRIBUTION: PLOTS

# NOTE ON THE STATISTICS: KS test is a non-parametric test to check whether two cumulative distributions are different
# how this works: for each value (VAF value) in the data, you count how many cases you have
# you then calculate the probability of being at each VAF value (nr cases @ this VAF / total nr of cases)
# you do this for both groups
# then you calculate the cumulative probability (basically add the probabilities for lower values)
# the KS statistic is the maximum value of the absolute difference between the two distributions

# PLOTS: SPLITTING INTO 2 GROUPS

# Log-log plot 
colors = {'top': col0r, 'bottom': col2r}
nr_groups = 2

# plot: separate for each variant 
for var in variants_identified: 

    data = netmhc2_df_labels_sub_groups2_eqsize_sorted[netmhc2_df_labels_sub_groups2_eqsize_sorted['gene_var'] == var]
    data = data.sort_values(by='VAF')
    N = data.shape[0] # number of carriers of the variant 
    
    # add number of people in each category 
    data['group_size'] = data.groupby('group_name')['group_name'].transform('count')
    data['index_vaf_group'] = data.groupby('group_name')['VAF'].rank(ascending=False) # index by VAF
    data['fraction_index_presenter_group'] = data['index_vaf_group'] / data['group_size']
    
    # make sure group name is the correct order of levels 
    data['group_name'] = data['group_name'].astype('category')
    data['group_name'].cat.reorder_categories(['top', 'bottom'])
    data['group_name'] = pd.Categorical(data['group_name'], categories=['top','bottom'], ordered=True)
    
    # need to add Kolmogorov-Smirnov test to check if distirbutions are significantly different from one another
    # I have two independent samples so will do the two-sample test 
    vaf1 = data[data['group_name']=='top'].VAF.tolist() # compare VAFs 
    vaf2 = data[data['group_name']=='bottom'].VAF.tolist() 
    stat, ks = ks_2samp(vaf1, vaf2, alternative="greater") # ks will indicate the p-value from the test 
    # your H1 is that vaf1 < vaf2 hence CDF vaf 1 > CDF vaf 2 so you need to use 'greater' as alternative
    ks = ks.round(3) 
    plot_cumulative(data, 'VAF_percent', 'fraction_index_presenter_group', 'group_name', colors, var, f'{nr_groups}_groups_eqsize_allon1', df, 2, 'log', ks = ks)

# %%
# plot: aggregate for all variants
netmhc2_df_labels_sub_groups2_eqsize_vaf = netmhc2_df_labels_sub_groups2_eqsize_sorted.sort_values(by='VAF')
netmhc2_df_labels_sub_groups2_eqsize_vaf['group_size'] = netmhc2_df_labels_sub_groups2_eqsize_vaf.groupby('group_name')['group_name'].transform('count')
netmhc2_df_labels_sub_groups2_eqsize_vaf['index_vaf_group'] = netmhc2_df_labels_sub_groups2_eqsize_vaf.groupby('group_name')['VAF'].rank(ascending=False) # index by VAF
netmhc2_df_labels_sub_groups2_eqsize_vaf['fraction_index_presenter_group'] = netmhc2_df_labels_sub_groups2_eqsize_vaf['index_vaf_group'] / netmhc2_df_labels_sub_groups2_eqsize_vaf['group_size']
netmhc2_df_labels_sub_groups2_eqsize_vaf['group_name'] = pd.Categorical(netmhc2_df_labels_sub_groups2_eqsize_vaf['group_name'], categories=['top', 'bottom'], ordered=True)
vaf1 = netmhc2_df_labels_sub_groups2_eqsize_vaf[netmhc2_df_labels_sub_groups2_eqsize_vaf['group_name']=='top'].fraction_index_presenter_group.tolist() # cumulative probabilities in sample 1 (top group)
vaf2 = netmhc2_df_labels_sub_groups2_eqsize_vaf[netmhc2_df_labels_sub_groups2_eqsize_vaf['group_name']=='bottom'].fraction_index_presenter_group.tolist() # cumulative probabilities in sample 2 (bottom group)
stat, ks = ks_2samp(vaf1, vaf2, alternative = 'greater') # ks will indicate the p-value of the test 
ks = ks.round(3)
plot_cumulative(netmhc2_df_labels_sub_groups2_eqsize_vaf, 'VAF_percent', 'fraction_index_presenter_group', 'group_name', colors, 'ALL_VARIANTS', f'{nr_groups}_groups_eqsize_allon1', df, 2, 'log', ks = ks, legend = True) 

# %%
# PLOT LINEAR VERSION OF THE LOG-LOG PLOT

# plot: separate for each variant
for var in variants_identified: 

    data = netmhc2_df_labels_sub_groups2_eqsize_sorted[netmhc2_df_labels_sub_groups2_eqsize_sorted['gene_var'] == var]
    data = data.sort_values(by='VAF')
    N = data.shape[0]

    # add number of people in each category 
    data['group_size'] = data.groupby('group_name')['group_name'].transform('count')
    data['index_vaf_group'] = data.groupby('group_name')['VAF'].rank(ascending=False) # index by VAF
    data['fraction_index_presenter_group'] = data['index_vaf_group'] / data['group_size']

    # make sure group name is the correct order of levels 
    data['group_name'] = data['group_name'].astype('category')
    data['group_name'].cat.reorder_categories(['top', 'bottom'])
    data['group_name'] = pd.Categorical(data['group_name'], categories=['top','bottom'], ordered=True)

    # KS test
    vaf1 = data[data['group_name']=='top'].VAF.tolist() # compare VAFs 
    vaf2 = data[data['group_name']=='bottom'].VAF.tolist() 
    stat, ks = ks_2samp(vaf1, vaf2, alternative="greater") # ks will indicate the p-value of the test 
    ks = ks.round(3)
    plot_cumulative(data, 'VAF_percent', 'fraction_index_presenter_group', 'group_name', colors, var, f'{nr_groups}_groups_eqsize_allon1', df, 2, 'linear', ks = ks)

# plot: aggregate for all variants
netmhc2_df_labels_sub_groups2_eqsize_vaf = netmhc2_df_labels_sub_groups2_eqsize_vaf.sort_values(by='VAF')
netmhc2_df_labels_sub_groups2_eqsize_vaf['group_size'] = netmhc2_df_labels_sub_groups2_eqsize_vaf.groupby('group_name')['group_name'].transform('count')
netmhc2_df_labels_sub_groups2_eqsize_vaf['index_vaf_group'] = netmhc2_df_labels_sub_groups2_eqsize_vaf.groupby('group_name')['VAF'].rank(ascending=False) # index by VAF
netmhc2_df_labels_sub_groups2_eqsize_vaf['fraction_index_presenter_group'] = netmhc2_df_labels_sub_groups2_eqsize_vaf['index_vaf_group'] / netmhc2_df_labels_sub_groups2_eqsize_vaf['group_size']
netmhc2_df_labels_sub_groups2_eqsize_vaf['group_name'] = pd.Categorical(netmhc2_df_labels_sub_groups2_eqsize_vaf['group_name'], categories=['top', 'bottom'], ordered=True)
# KS test
vaf1 = netmhc2_df_labels_sub_groups2_eqsize_vaf[netmhc2_df_labels_sub_groups2_eqsize_vaf['group_name']=='top'].VAF.tolist() # compare VAFs 
vaf2 = netmhc2_df_labels_sub_groups2_eqsize_vaf[netmhc2_df_labels_sub_groups2_eqsize_vaf['group_name']=='bottom'].VAF.tolist() 
stat, ks = ks_2samp(vaf1, vaf2, alternative = 'greater') # ks will indicate the p-value of the test 
ks = ks.round(3)
plot_cumulative(netmhc2_df_labels_sub_groups2_eqsize_vaf, 'VAF_percent', 'fraction_index_presenter_group', 'group_name', colors, 'ALL_VARIANTS', f'{nr_groups}_groups_eqsize_allon1', df, 2, 'linear', ks = ks, legend = True) 

# %%
# PLOT: SPLITTING INTO 3 GROUPS 
# PLOT CUMULATIVE DISTRIBUTION (LOG-LOG PLOTS)
# specify the colors for plotting 
colors = {'top': col0r, 'mid': col1r, 'bottom': col2r}
nr_groups = 3

# plot: separate for each variant 
for var in variants_identified: 

    data = netmhc2_df_labels_sub_groups3_eqsize_sorted[netmhc2_df_labels_sub_groups3_eqsize_sorted['gene_var'] == var]
    data = data.sort_values(by='VAF')
    N = data.shape[0]

    # add number of people in each category 
    data['group_size'] = data.groupby('group_name')['group_name'].transform('count')
    data['index_vaf_group'] = data.groupby('group_name')['VAF'].rank(ascending=False) # index by VAF
    data['fraction_index_presenter_group'] = data['index_vaf_group'] / data['group_size']

    # make sure group name is the correct order of levels 
    data['group_name'] = data['group_name'].astype('category')
    data['group_name'].cat.reorder_categories(['top', 'mid', 'bottom'])
    data['group_name'] = pd.Categorical(data['group_name'], categories=['top', 'mid', 'bottom'], ordered=True)

    # plot
    plot_cumulative(data, 'VAF_percent', 'fraction_index_presenter_group', 'group_name', colors, var, f'{nr_groups}_groups_eqsize_allon1', df, 3, 'log')

# plot: aggregate for all variants
netmhc2_df_labels_sub_groups3_eqsize_vaf = netmhc2_df_labels_sub_groups3_eqsize_sorted.sort_values(by='VAF')
netmhc2_df_labels_sub_groups3_eqsize_vaf['group_size'] = netmhc2_df_labels_sub_groups3_eqsize_vaf.groupby('group_name')['group_name'].transform('count')
netmhc2_df_labels_sub_groups3_eqsize_vaf['index_vaf_group'] = netmhc2_df_labels_sub_groups3_eqsize_vaf.groupby('group_name')['VAF'].rank(ascending=False) # index by VAF
netmhc2_df_labels_sub_groups3_eqsize_vaf['fraction_index_presenter_group'] = netmhc2_df_labels_sub_groups3_eqsize_vaf['index_vaf_group'] / netmhc2_df_labels_sub_groups3_eqsize_vaf['group_size']
netmhc2_df_labels_sub_groups3_eqsize_vaf['group_name'] = pd.Categorical(netmhc2_df_labels_sub_groups3_eqsize_vaf['group_name'], categories=['top', 'mid', 'bottom'], ordered=True)
plot_cumulative(netmhc2_df_labels_sub_groups3_eqsize_vaf, 'VAF_percent', 'fraction_index_presenter_group', 'group_name', colors, 'ALL_VARIANTS', f'{nr_groups}_groups_eqsize_allon1', df, 3, 'log', legend = True) 

# %%
# PLOT ON A LINEAR SCALE
# specify the colors for plotting 
colors = {'top': col0r, 'mid': col1r, 'bottom': col2r}
nr_groups = 3

# plot: separate for each variant 
for var in variants_identified: 

    data = netmhc2_df_labels_sub_groups3_eqsize_vaf[netmhc2_df_labels_sub_groups3_eqsize_vaf['gene_var'] == var]
    data = data.sort_values(by='VAF')
    N = data.shape[0]

    # add number of people in each category 
    data['group_size'] = data.groupby('group_name')['group_name'].transform('count')
    data['index_vaf_group'] = data.groupby('group_name')['VAF'].rank(ascending=False) # index by VAF
    data['fraction_index_presenter_group'] = data['index_vaf_group'] / data['group_size']

    # make sure group name is the correct order of levels 
    data['group_name'] = data['group_name'].astype('category')
    data['group_name'].cat.reorder_categories(['top', 'mid', 'bottom'])
    data['group_name'] = pd.Categorical(data['group_name'], categories=['top', 'mid', 'bottom'], ordered=True)

    # plot
    plot_cumulative(data, 'VAF_percent', 'fraction_index_presenter_group', 'group_name', colors, var, f'{nr_groups}_groups_eqsize_allon1', df, 3, 'linear')

# plot: aggregate for all variants
netmhc2_df_labels_sub_groups3_eqsize_vaf = netmhc2_df_labels_sub_groups3_eqsize_vaf.sort_values(by='VAF')
netmhc2_df_labels_sub_groups3_eqsize_vaf['group_size'] = netmhc2_df_labels_sub_groups3_eqsize_vaf.groupby('group_name')['group_name'].transform('count')
netmhc2_df_labels_sub_groups3_eqsize_vaf['index_vaf_group'] = netmhc2_df_labels_sub_groups3_eqsize_vaf.groupby('group_name')['VAF'].rank(ascending=False) # index by VAF
netmhc2_df_labels_sub_groups3_eqsize_vaf['fraction_index_presenter_group'] = netmhc2_df_labels_sub_groups3_eqsize_vaf['index_vaf_group'] / netmhc2_df_labels_sub_groups3_eqsize_vaf['group_size']
netmhc2_df_labels_sub_groups3_eqsize_vaf['group_name'] = pd.Categorical(netmhc2_df_labels_sub_groups3_eqsize_vaf['group_name'], categories=['top', 'mid', 'bottom'], ordered=True)
plot_cumulative(netmhc2_df_labels_sub_groups3_eqsize_vaf, 'VAF_percent', 'fraction_index_presenter_group', 'group_name', colors, 'ALL_VARIANTS', f'{nr_groups}_groups_eqsize_allon1', df, 3, 'linear', legend = True) 

# %% 
#################################################################################
#################################################################################
#################################################################################
# FIGURE 4
# VAF (IN CARRIERS) VS ABSOLUTE BINDING CAPACITY TO MHC I
# REQUIRES USE OF THE DATAFRAME WITH LABELS

# %%  
# DEFINE FUNCTION TO PLOT

# PLOT CUMULATIVE DISTRIBUTION (function to plot either on log-log or linear scale)
def plot_cumulative_absolute(data, x, y, hue, colors, var, df, split, style, n1, n2, n3 = None, ks = None, legend = None):

    var_name = var.split('_')[0:2]
    var_name = ' '.join(var_name)
    part1 = var_name.split(' ')[0]
    part2 = var_name.split(' ')[1]
    var_name = '$\mathbf{' + part1 + r'}$' + '\n' + '$\mathbf{' + part2 + r'}$'
    
    if var == 'ALL_VARIANTS':
        nr_cases = data.gene_var.value_counts().sum()
    else:
        nr_cases = data.shape[0]

    plot_width = 1.4
    plot_height = 1.4

    # Margins around the plot area (in inches)
    left_margin = 0.4
    right_margin = 0.4
    top_margin = 0.4
    bottom_margin = 0.4

    # Calculate the figure size
    fig_width = plot_width + left_margin + right_margin
    fig_height = plot_height + top_margin + bottom_margin
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    plt.subplots_adjust(left=left_margin/fig_width,
                    right=1-right_margin/fig_width,
                    top=1-top_margin/fig_height,
                    bottom=bottom_margin/fig_height)

    ax = sns.lineplot(x=x, y=y, data=data, hue=hue, palette=colors, legend = False)
    scatter = sns.scatterplot(x=x, y=y, data=data, hue=hue, palette=colors, legend = True, size = 11, alpha = 0.8, linewidth=0)

    # add text on the plot (variant name and nr of carriers) 
    plt.text(0.8, 0.82, f'{var_name}\n(n={nr_cases})', transform = plt.gca().transAxes,
         fontsize = 7, verticalalignment = 'center', horizontalalignment = 'center')
    
    # specify title and axes labels 
    plt.xlabel('Variant allele frequency (%)', fontsize = xaxis_font-4)
    plt.ylabel('Reverse cumulative probabilty', fontsize = yaxis_font-4)

    # convert both axes to log scale (if log-log)
    if style == 'log':
        plt.gca().set_xscale('log')
        plt.gca().set_yscale('log')

        x_ticks = [2, 3, 10, 30, 100]
        plt.xticks(x_ticks, x_ticks, fontsize = xticks_font-3)
    else:
        plt.xticks(fontsize = xticks_font-3)
   
    # add results of the KS test
    if ks == None:
        plt.text(0.8, 0.55, ' ') # don't add any text if KS test not run
    else:
        plt.text(0.8, 0.65, f'p={ks}', transform = plt.gca().transAxes,
        fontsize = 5.5, verticalalignment = 'center', horizontalalignment = 'center')
        
    # specify y axis ticks
    scatter.yaxis.set_major_formatter(FuncFormatter(custom_formatter))
    scatter.tick_params(axis = 'y', labelsize = yticks_font-3) 
    
    # specify limits on both axes 
    if style == 'log':
        if min(data[y].tolist()) > 0.1:
            plt.xlim(2, 100)
            plt.ylim(0.05, 1.2) # show where 0.1 is for each plot
        else:
            plt.xlim(2, 100)
            plt.ylim(None, 1.2)

    elif style == 'linear':
        plt.xlim(2, None)
        plt.ylim(-0.05, 1.05)

    # add legend (will only use in the plot for all variants)
    if style == 'log':
        if legend == True: # add legend for the color code you are using 
            if split == '2_strong':
                handles, labels = scatter.get_legend_handles_labels()
                plt.legend(labels = [f'binding\n(%EL rank < 1)\nn = {n1}', f'no binding\n(%EL rank > 1)\nn = {n2}'],
                        handles = handles[:-1], loc='lower left', frameon = False,
                        fontsize = text_font-4, markerscale = 1, handletextpad = 0.2)
            elif split == '2_weak':
                handles, labels = scatter.get_legend_handles_labels()
                plt.legend(labels = [f'binding\n(%EL rank < 5)\nn = {n1}', f'no binding\n(%EL rank > 5)\nn = {n2}'],
                        handles = handles[:-1], loc='lower left', frameon = False,
                        fontsize = text_font-4, markerscale = 1, handletextpad = 0.2)
            else:
                handles, labels = scatter.get_legend_handles_labels()
                plt.legend(labels = [f'strong binding\n(%EL rank < 1)\nn = {n1}', f'weak binding\n(1 < %EL rank < 5)\nn = {n2}', f'no binding\n(%EL rank > 5)\nn = {n3}'],
                        handles = handles[:-1], loc='lower left', frameon = False,
                        fontsize = text_font-4, markerscale = 1, handletextpad = 0.2)
        else: # only include the numbers of individuals in different groups
            if split == '2_strong':
                    handles, labels = scatter.get_legend_handles_labels()
                    plt.legend(labels = [f'n = {n1}', f'n = {n2}'],
                            handles = handles[:-1], loc='lower left', frameon = False,
                            fontsize = text_font-4, markerscale = 1, handletextpad = 0.2)
            elif split == '2_weak':
                handles, labels = scatter.get_legend_handles_labels()
                plt.legend(labels = [f'n = {n1}', f'n = {n2}'],
                        handles = handles[:-1], loc='lower left', frameon = False,
                        fontsize = text_font-4, markerscale = 1, handletextpad = 0.2)
            else:
                handles, labels = scatter.get_legend_handles_labels()
                plt.legend(labels = [f'n = {n1}', f'n = {n2}', f'n = {n3}'],
                        handles = handles[:-1], loc='lower left', frameon = False,
                        fontsize = text_font-4, markerscale = 1, handletextpad = 0.2) # explicitly do NOT include the legend
    elif style == 'linear':
        if legend == True: # add legend for the color code you are using 
            if split == '2_strong':
                handles, labels = scatter.get_legend_handles_labels()
                plt.legend(labels = [f'binding\n(%EL rank < 1)\nn = {n1}', f'no binding\n(%EL rank > 1)\nn = {n2}'],
                        handles = handles[:-1], loc='upper left', frameon = False,
                        fontsize = text_font-4, markerscale = 1, handletextpad = 0.2)
            elif split == '2_weak':
                handles, labels = scatter.get_legend_handles_labels()
                plt.legend(labels = [f'binding\n(%EL rank < 5)\nn = {n1}', f'no binding\n(%EL rank > 5)\nn = {n2}'],
                        handles = handles[:-1], loc='upper left', frameon = False,
                        fontsize = text_font-4, markerscale = 1, handletextpad = 0.2)
            else:
                handles, labels = scatter.get_legend_handles_labels()
                plt.legend(labels = [f'strong binding\n(%EL rank < 1)\nn = {n1}', f'weak binding\n(1 < %EL rank < 5)\nn = {n2}', f'no binding\n(%EL rank > 5)\nn = {n3}'],
                        handles = handles[:-1], loc='upper left', frameon = False,
                        fontsize = text_font-4, markerscale = 1, handletextpad = 0.2)
        else: # only include the numbers of individuals in different groups
            if split == '2_strong':
                    handles, labels = scatter.get_legend_handles_labels()
                    plt.legend(labels = [f'n = {n1}', f'n = {n2}'],
                            handles = handles[:-1], loc='upper left', frameon = False,
                            fontsize = text_font-4, markerscale = 1, handletextpad = 0.2)
            elif split == '2_weak':
                handles, labels = scatter.get_legend_handles_labels()
                plt.legend(labels = [f'n = {n1}', f'n = {n2}'],
                        handles = handles[:-1], loc='upper left', frameon = False,
                        fontsize = text_font-4, markerscale = 1, handletextpad = 0.2)
            else:
                handles, labels = scatter.get_legend_handles_labels()
                plt.legend(labels = [f'n = {n1}', f'n = {n2}', f'n = {n3}'],
                        handles = handles[:-1], loc='upper left', frameon = False,
                        fontsize = text_font-4, markerscale = 1, handletextpad = 0.2) # explicitly do NOT include the legend
    
    if df == 'netmhc_2reads':
        plt.savefig(f'/Users/barbarawalkowiak/Desktop/msc_thesis/results/figures/fig4/VAF_threshold_{var}_{split}_{style}_{df}.pdf', bbox_inches='tight')
    else:
        plt.savefig(f'/Users/barbarawalkowiak/Desktop/msc_thesis/results/supplementary_files/{df}/fig4/VAF_threshold_{var}_{split}_{style}_{df}.pdf', bbox_inches='tight')

# %%
# select the dataframe with variants for which it makes sense to do this 
netmhc2_df_labels_sub_carriers_sub = netmhc2_df_labels_sub_carriers[netmhc2_df_labels_sub_carriers['CH_variant'].isin(variants_thresh)]

# indicate level of binding (binding / no binding, threshold = 2)
conditions1 = [
    (netmhc2_df_labels_sub_carriers_sub['log_score'] < -1 * np.log10(5)),
    (netmhc2_df_labels_sub_carriers_sub['log_score'] >= -1 * np.log10(5))
]

values1 = ['no binding', 'binding']
netmhc2_df_labels_sub_carriers_sub['any_binding'] = np.select(conditions1, values1)
netmhc2_df_labels_sub_carriers_sub['any_binding']= pd.Categorical(netmhc2_df_labels_sub_carriers_sub['any_binding'], categories=['binding', 'no binding'], ordered=True)

# indicate level of binding (binding / no binding)
conditions2 = [
    (netmhc2_df_labels_sub_carriers_sub['log_score'] < -1 * np.log10(1)),
    (netmhc2_df_labels_sub_carriers_sub['log_score'] >= -1 * np.log10(1))
]

values2 = ['no binding', 'binding']
netmhc2_df_labels_sub_carriers_sub['strong_binding'] = np.select(conditions2, values2)
netmhc2_df_labels_sub_carriers_sub['strong_binding']= pd.Categorical(netmhc2_df_labels_sub_carriers_sub['strong_binding'], categories=['binding', 'no binding'], ordered=True)

# indicate level of binding (no / weak / strong binding)
conditions3 = [
    (netmhc2_df_labels_sub_carriers_sub['log_score'] < -1 * np.log10(5)),
    ((netmhc2_df_labels_sub_carriers_sub['log_score'] >= -1 * np.log10(5)) & (netmhc2_df_labels_sub_carriers_sub['log_score'] < -1 * np.log10(1))),
    (netmhc2_df_labels_sub_carriers_sub['log_score'] >= -1 * np.log10(1))
]

values3 = ['no binding', 'weak binding', 'strong binding']
netmhc2_df_labels_sub_carriers_sub['binding'] = np.select(conditions3, values3)
netmhc2_df_labels_sub_carriers_sub['binding']= pd.Categorical(netmhc2_df_labels_sub_carriers_sub['binding'], categories=['strong binding', 'weak binding', 'no binding'], ordered=True)
  
# %%
# add VAF data 
netmhc2_df_labels_sub_carriers_subscores_vaf = pd.merge(netmhc2_df_labels_sub_carriers_sub, netmhc2_df[['Person_ID', 'VAF']])
netmhc2_df_labels_sub_carriers_subscores_vaf['VAF_percent'] = netmhc2_df_labels_sub_carriers_subscores_vaf['VAF'] * 100

# PLOT: ANY BINDING
colors1 = {'binding': col0a, 'no binding': col2a}
colors2 = {'strong binding': col0a, 'weak binding': col1a, 'no binding': col2a}

# %%
# for each variant / binding v non-binding at threshold = 5
for i, var in enumerate(order_by_total31): 

    data = netmhc2_df_labels_sub_carriers_subscores_vaf[netmhc2_df_labels_sub_carriers_subscores_vaf['gene_var'] == var]        
    data = data.sort_values(by='VAF')
    N = data.shape[0]
    n1 = data[data['any_binding'] == 'binding'].shape[0]
    n2 = data[data['any_binding'] == 'no binding'].shape[0]

    # add number of people in each category 
    data['group_size'] = data.groupby('any_binding')['any_binding'].transform('count')
    data['index_vaf_group'] = data.groupby('any_binding')['VAF'].rank(ascending=False) # index by VAF
    data['fraction_index_presenter_group'] = data['index_vaf_group'] / data['group_size']

    # KS test
    vaf1 = data[data['any_binding']=='binding'].VAF.tolist() # compare VAFs 
    vaf2 = data[data['any_binding']=='no binding'].VAF.tolist() 
    stat, ks = ks_2samp(vaf1, vaf2, alternative = 'greater') # ks will indicate the p-value of the test 
    ks = ks.round(3)
    if i == 0: # include legend only on the first plot
        plot_cumulative_absolute(data, 'VAF_percent', 'fraction_index_presenter_group', 'any_binding', colors1, var, df, '2_weak', 'log', n1, n2, ks = ks, legend = True)
    else:
        plot_cumulative_absolute(data, 'VAF_percent', 'fraction_index_presenter_group', 'any_binding', colors1, var, df, '2_weak', 'log', n1, n2, ks = ks)

# %%
# for each variant / binding v non-binding at threshold = 1
for i, var in enumerate(order_by_total31): 

    data = netmhc2_df_labels_sub_carriers_subscores_vaf[netmhc2_df_labels_sub_carriers_subscores_vaf['gene_var'] == var]        
    data = data.sort_values(by='VAF')
    N = data.shape[0]
    n1 = data[data['strong_binding'] == 'binding'].shape[0]
    n2 = data[data['strong_binding'] == 'no binding'].shape[0]

    # add number of people in each category 
    data['group_size'] = data.groupby('strong_binding')['strong_binding'].transform('count')
    data['index_vaf_group'] = data.groupby('strong_binding')['VAF'].rank(ascending=False) # index by VAF
    data['fraction_index_presenter_group'] = data['index_vaf_group'] / data['group_size']

    # KS test
    vaf1 = data[data['strong_binding']=='binding'].VAF.tolist() # compare VAFs 
    vaf2 = data[data['strong_binding']=='no binding'].VAF.tolist() 
    stat, ks = ks_2samp(vaf1, vaf2, alternative = 'greater') # ks will indicate the p-value of the test 
    ks = ks.round(3)
    if i == 0: # include legend only on the first plot
        plot_cumulative_absolute(data, 'VAF_percent', 'fraction_index_presenter_group', 'strong_binding', colors1, var, df, '2_strong', 'log', n1, n2, ks = ks, legend = True)
    else:
        plot_cumulative_absolute(data, 'VAF_percent', 'fraction_index_presenter_group', 'strong_binding', colors1, var, df, '2_strong', 'log', n1, n2, ks = ks)

# %%
# PLOT: STRONG / WEAK / NO BINDING
# for each variant / strong vs weak vs non-binding 
for i, var in enumerate(order_by_total31): 

    data = netmhc2_df_labels_sub_carriers_subscores_vaf[netmhc2_df_labels_sub_carriers_subscores_vaf['gene_var'] == var]        
    data = data.sort_values(by='VAF')
    N = data.shape[0]
    n1 = data[data['binding'] == 'strong binding'].shape[0]
    n2 = data[data['binding'] == 'weak binding'].shape[0]
    n3 = data[data['binding'] == 'no binding'].shape[0]

    # add number of people in each category 
    data['group_size'] = data.groupby('binding')['binding'].transform('count')
    data['index_vaf_group'] = data.groupby('binding')['VAF'].rank(ascending=False) # index by VAF
    data['fraction_index_presenter_group'] = data['index_vaf_group'] / data['group_size']

    if i == 0: # include legend only on the first plot
        plot_cumulative_absolute(data, 'VAF_percent', 'fraction_index_presenter_group', 'binding', colors2, var, df, '3_sw', 'log', n1, n2, n3 = n3, legend = True)
    else:
        plot_cumulative_absolute(data, 'VAF_percent', 'fraction_index_presenter_group', 'binding', colors2, var, df, '3_sw', 'log', n1, n2, n3 = n3)

# %%
# NOW PLOT ALL THIS BUT LINEAR PLOT
        
# split into two groups (threshold = 2)
for i, var in enumerate(order_by_total31): 

    data = netmhc2_df_labels_sub_carriers_subscores_vaf[netmhc2_df_labels_sub_carriers_subscores_vaf['gene_var'] == var]        
    data = data.sort_values(by='VAF')
    N = data.shape[0]
    n1 = data[data['any_binding'] == 'binding'].shape[0]
    n2 = data[data['any_binding'] == 'no binding'].shape[0]

    # add number of people in each category 
    data['group_size'] = data.groupby('any_binding')['any_binding'].transform('count')
    data['index_vaf_group'] = data.groupby('any_binding')['VAF'].rank(ascending=False) # index by VAF
    data['fraction_index_presenter_group'] = data['index_vaf_group'] / data['group_size']

    # KS test
    vaf1 = data[data['any_binding']=='binding'].VAF.tolist() # compare VAFs 
    vaf2 = data[data['any_binding']=='no binding'].VAF.tolist() 
    stat, ks = ks_2samp(vaf1, vaf2, alternative = 'greater') # ks will indicate the p-value of the test 
    ks = ks.round(3)
    if i == 0: # include legend only on the first plot
        plot_cumulative_absolute(data, 'VAF_percent', 'fraction_index_presenter_group', 'any_binding', colors1, var, df, '2_weak', 'linear', n1, n2, ks = ks, legend = True)
    else:
        plot_cumulative_absolute(data, 'VAF_percent', 'fraction_index_presenter_group', 'any_binding', colors1, var, df, '2_weak', 'linear', n1, n2, ks = ks)
  
# split into two groups (threshold = 0.5)
for i, var in enumerate(order_by_total31): 

    data = netmhc2_df_labels_sub_carriers_subscores_vaf[netmhc2_df_labels_sub_carriers_subscores_vaf['gene_var'] == var]        
    data = data.sort_values(by='VAF')
    N = data.shape[0]
    n1 = data[data['strong_binding'] == 'binding'].shape[0]
    n2 = data[data['strong_binding'] == 'no binding'].shape[0]

    # add number of people in each category 
    data['group_size'] = data.groupby('strong_binding')['strong_binding'].transform('count')
    data['index_vaf_group'] = data.groupby('strong_binding')['VAF'].rank(ascending=False) # index by VAF
    data['fraction_index_presenter_group'] = data['index_vaf_group'] / data['group_size']

    # KS test
    vaf1 = data[data['strong_binding']=='binding'].VAF.tolist() # compare VAFs 
    vaf2 = data[data['strong_binding']=='no binding'].VAF.tolist() 
    stat, ks = ks_2samp(vaf1, vaf2, alternative = 'greater') # ks will indicate the p-value of the test 
    ks = ks.round(3)
    if i == 0: # include legend only on the first plot
        plot_cumulative_absolute(data, 'VAF_percent', 'fraction_index_presenter_group', 'strong_binding', colors1, var, df, '2_strong', 'linear', n1, n2, ks = ks, legend = True)
    else:
        plot_cumulative_absolute(data, 'VAF_percent', 'fraction_index_presenter_group', 'strong_binding', colors1, var, df, '2_strong', 'linear', n1, n2, ks = ks)

# split into three groups (strong / weak / no binding)
for i, var in enumerate(order_by_total31): 

    data = netmhc2_df_labels_sub_carriers_subscores_vaf[netmhc2_df_labels_sub_carriers_subscores_vaf['gene_var'] == var]        
    data = data.sort_values(by='VAF')
    N = data.shape[0]
    n1 = data[data['binding'] == 'strong binding'].shape[0]
    n2 = data[data['binding'] == 'weak binding'].shape[0]
    n3 = data[data['binding'] == 'no binding'].shape[0]

    # add number of people in each category 
    data['group_size'] = data.groupby('binding')['binding'].transform('count')
    data['index_vaf_group'] = data.groupby('binding')['VAF'].rank(ascending=False) # index by VAF
    data['fraction_index_presenter_group'] = data['index_vaf_group'] / data['group_size']

    if i == 0: # include legend only on the first plot
        plot_cumulative_absolute(data, 'VAF_percent', 'fraction_index_presenter_group', 'binding', colors2, var, df, '3_sw', 'linear', n1, n2, n3 = n3, legend = True)
    else:
        plot_cumulative_absolute(data, 'VAF_percent', 'fraction_index_presenter_group', 'binding', colors2, var, df, '3_sw', 'linear', n1, n2, n3 = n3)

# %%
#################################################################################
#################################################################################
#################################################################################
# SF (discussion): scores for present vs absent variants for carriers
   
# format the dataframe 
# look only at variants you are including in the rest of the analysis
variants_identified = netmhc2_df_labels.gene_var.unique().tolist()
variants_identified = [var for var in variants_identified if netmhc2_df_labels[netmhc2_df_labels['gene_var'] == var].shape[0] >= 10]
variants_with_scores = netmhc2_df_labels.CH_variant.unique().tolist()
vars_with_no_carriers = [var for var in variants_with_scores if var not in variants_identified]
vars_with_carriers = [var for var in variants_with_scores if var in variants_identified]

# only consider people with a single variant 
netmhc2_df_dp =  netmhc2_df[~netmhc2_df.Person_ID.duplicated()]
 
# subset the dataframe 
netmhc2_df_sub = netmhc2_df_dp[~netmhc2_df_dp['gene_var'].isin(vars_with_no_carriers)]

# select only carriers
netmhc2_df_carriers = netmhc2_df_sub[netmhc2_df_sub['ch_status']==1]

# select only columns with genetic variant they carry and scores for other variants 
scores_all = [col for col in netmhc2_df_carriers.columns if 'score_' in col]
scores_selected = [col for col in scores_all if col.replace('score_', '') in vars_with_carriers]
scores_cols_df = netmhc2_df_carriers[scores_selected]
netmhc2_df_carriers_scores = pd.concat([netmhc2_df_carriers[['Person_ID', 'gene_var']], scores_cols_df], axis = 1)

# melt the dataframe
netmhc2_df_carriers_scores_melted = pd.melt(netmhc2_df_carriers_scores, id_vars = ['Person_ID', 'gene_var'])

# add column to indicate if score if for variant present or absent 
netmhc2_df_carriers_scores_melted['is_variant_present'] = np.where(netmhc2_df_carriers_scores_melted['gene_var'] == netmhc2_df_carriers_scores_melted['variable'].str.replace('score_', ''), 'present', 'absent')

# formatting names for plotting 
netmhc2_df_carriers_scores_melted['gene_var2'] = netmhc2_df_carriers_scores_melted['gene_var'].str.replace('_', '\n')
netmhc2_df_carriers_scores_melted['log_score'] = -1 * np.log10(netmhc2_df_carriers_scores_melted['value'])
                                                               
# %%
# plot
colors = [col_pos, col_neg]

plt.figure(figsize = (16, 4))
ax = sns.stripplot(data = netmhc2_df_carriers_scores_melted, 
                   dodge = True, jitter = False, x = 'gene_var2', y = 'log_score', hue = 'is_variant_present', palette = colors, size = 1.5, edgecolor = 'black', order = order_by_total2, alpha = 0.4)

# GREY BACKGROUND (every other variant)
for i in range(1, len(order_by_total2), 2):
    ax.axvspan(i-0.5, i+0.5, color=col_background, alpha = 0.8)

for i, artist in enumerate(ax.collections):
    # Get the current positions of the points
    offsets = artist.get_offsets()
    dodge_extent = 0.01
    offsets[:, 0] += (i % 2) * dodge_extent - dodge_extent / 2
    # Update the positions
    artist.set_offsets(offsets)

# find median score and add onto the plot 
for i, category in enumerate(order_by_total):
            
    median_present = netmhc2_df_carriers_scores_melted[(netmhc2_df_carriers_scores_melted['gene_var'] == f'{category}') & (netmhc2_df_carriers_scores_melted['is_variant_present'] == 'present')].log_score.median()
    median_absent = netmhc2_df_carriers_scores_melted[(netmhc2_df_carriers_scores_melted['gene_var'] == f'{category}') & (netmhc2_df_carriers_scores_melted['is_variant_present'] == 'absent')].log_score.median()

    # Plot text for each hue group
    plt.text(i, median_present, '-', ha='right', va='center', fontsize=30, fontweight='bold', color = col_pos)
    plt.text(i, median_absent, '-', ha='left', va='center', fontsize=30, fontweight='bold', color = col_neg)

# Add labels to the right of the markers
plt.text(0.6, -1*np.log10(0.5)+1.9, 'strong\nbinding', color='black', verticalalignment='center', fontsize = 11, ha = 'center')
plt.text(0.6, -1*np.log10(2)-2.2, 'no\nbinding', color='black', verticalalignment='center', fontsize = 11, ha = 'center')

# add lines to indicate thresholds for strong and weak binding
plt.axhline(-1*np.log10(0.5), color = '#3D3D3D', linewidth = 0.75, linestyle='--')
plt.axhline(-1*np.log10(2.0), color = '#3D3D3D', linewidth = 0.75, linestyle='--')

# axes labelling
plt.xlabel(f'CH variant', fontsize = xaxis_font)
plt.ylabel('MHC-variant binding', fontsize = yaxis_font)

plt.xlim(-0.5, len(order_by_total2)-0.5)
plt.ylim(-3.2, 3) # I know this will be the max value for all variants bc R882H is most common

# BOLD X TICKS for which you have binding predictions that span both no binding and good binding
labels = []
for label in order_by_total2:
    if label in variants_thresh2:
        part1 = label.split('\n')[0]
        part2 = label.split('\n')[1]
        label = '$\mathbf{' + part1 + r'}$' + '\n' + '$\mathbf{' + part2 + r'}$'
        labels.append(label)
    else:
        labels.append(label)

ax.set_xticklabels(labels, fontsize = xticks_font, rotation = 90)

plt.yticks(fontsize = yticks_font)
legend = plt.legend(['variant present', 'variants absent'], markerscale = 4.7, loc = 'upper right', fontsize = legend_font, frameon = False, handletextpad = 0.2)
for legend_handle in legend.legendHandles:
    legend_handle.set_alpha(1)

# save the main figure 
if df == 'netmhc_2reads':    
    plt.savefig(f'/Users/barbarawalkowiak/Desktop/msc_thesis/results/figures/fig5/5a_carriers_present_vs_absent_variants_{df}.pdf', bbox_inches='tight')
else: # for other dataframes / ways of analyzing, save to the folder in supplmentary figures 
    plt.savefig(f'/Users/barbarawalkowiak/Desktop/msc_thesis/results/supplementary_files/{df}/fig5/5a_carriers_present_vs_absent_variants_{df}.pdf', bbox_inches='tight')

# %%
# plot variants when aggregated
colors = [col_pos, col_neg]

plt.figure(figsize = (2, 4))
ax = sns.stripplot(data = netmhc2_df_carriers_scores_melted, 
                   jitter = False, x = 'is_variant_present', y = 'log_score', palette = colors, size = 1.5, edgecolor = 'black', alpha = 0.4)

# find median score and add onto the plot    
median_present = netmhc2_df_carriers_scores_melted[netmhc2_df_carriers_scores_melted['is_variant_present'] == 'present'].log_score.median()
median_absent = netmhc2_df_carriers_scores_melted[netmhc2_df_carriers_scores_melted['is_variant_present'] == 'absent'].log_score.median()

# add median to the plot
plt.text(0, median_present, '--', ha='center', va='center', fontsize=30, color = col_pos)
plt.text(1, median_absent, '--', ha='center', va='center', fontsize=30, color = col_neg)

# add MW test for aggregated scores 
scores_present = netmhc2_df_carriers_scores_melted[netmhc2_df_carriers_scores_melted['is_variant_present'] == 'present'].log_score.tolist()
scores_absent = netmhc2_df_carriers_scores_melted[netmhc2_df_carriers_scores_melted['is_variant_present'] == 'absent'].log_score.tolist()
# run the MW test
statistic, p_value = mannwhitneyu(scores_present, scores_absent)
p_value = f"{p_value:.4g}"
plt.text(0.5, 2.3, f'p-value:\n{p_value}', ha='center', va='bottom', fontsize=text_font)

# add lines to indicate thresholds for strong and weak binding
plt.axhline(-1*np.log10(0.5), color = '#3D3D3D', linewidth = 0.75, linestyle='--')
plt.axhline(-1*np.log10(2.0), color = '#3D3D3D', linewidth = 0.75, linestyle='--')

# axes labelling
plt.xlabel(f'Is variant present?', fontsize = xaxis_font)
plt.ylabel('MHC-variant binding', fontsize = yaxis_font)

plt.ylim(-3.2, 3) # I know this will be the max value for all variants bc R882H is most common

# x ticks 
plt.xticks(fontsize = xticks_font)
plt.yticks(fontsize = yticks_font)

# save the main figure 
if df == 'netmhc_2reads':    
    plt.savefig(f'/Users/barbarawalkowiak/Desktop/msc_thesis/results/figures/fig5/5b_carriers_present_vs_absent_variants_aggregared_{df}.pdf', bbox_inches='tight')
else: # for other dataframes / ways of analyzing, save to the folder in supplmentary figures 
    plt.savefig(f'/Users/barbarawalkowiak/Desktop/msc_thesis/results/supplementary_files/{df}/fig5/5b_carriers_present_vs_absent_variants_aggregated_{df}.pdf', bbox_inches='tight')

# %%
# plot variants when aggregated (removing most common ones)
plt.figure(figsize = (2, 4))
dnmt3a_to_remove = ['DNMT3A_R882H', 'DNMT3A_R882C']
netmhc2_df_carriers_scores_melted_dnmt3r = netmhc2_df_carriers_scores_melted[~netmhc2_df_carriers_scores_melted['gene_var'].isin(dnmt3a_to_remove)]
scores_to_remove = ['score_DNMT3A_R882H', 'score_DNMT3A_R882C']
netmhc2_df_carriers_scores_melted_dnmt3r = netmhc2_df_carriers_scores_melted_dnmt3r[~netmhc2_df_carriers_scores_melted_dnmt3r['variable'].isin(scores_to_remove)]

ax = sns.stripplot(data = netmhc2_df_carriers_scores_melted_dnmt3r, 
                   jitter = False, x = 'is_variant_present', y = 'log_score', palette = colors, size = 1.5, edgecolor = 'black', alpha = 0.4)

# find median score and add onto the plot    
median_present = netmhc2_df_carriers_scores_melted_dnmt3r[netmhc2_df_carriers_scores_melted_dnmt3r['is_variant_present'] == 'present'].log_score.median()
median_absent = netmhc2_df_carriers_scores_melted_dnmt3r[netmhc2_df_carriers_scores_melted_dnmt3r['is_variant_present'] == 'absent'].log_score.median()

# add median to the plot
plt.text(0, median_present, '--', ha='center', va='center', fontsize=30,  color = col_pos)
plt.text(1, median_absent, '--', ha='center', va='center', fontsize=30, color = col_neg)

# add MW test for aggregated scores 
scores_present = netmhc2_df_carriers_scores_melted_dnmt3r[netmhc2_df_carriers_scores_melted_dnmt3r['is_variant_present'] == 'present'].log_score.tolist()
scores_absent = netmhc2_df_carriers_scores_melted_dnmt3r[netmhc2_df_carriers_scores_melted_dnmt3r['is_variant_present'] == 'absent'].log_score.tolist()
# run the MW test
statistic, p_value = mannwhitneyu(scores_present, scores_absent)
p_value = f"{p_value:.4g}"
plt.text(0.5, 2.3, f'p-value:\n{p_value}', ha='center', va='bottom', fontsize=text_font)

# add lines to indicate thresholds for strong and weak binding
plt.axhline(-1*np.log10(0.5), color = '#3D3D3D', linewidth = 0.75, linestyle='--')
plt.axhline(-1*np.log10(2.0), color = '#3D3D3D', linewidth = 0.75, linestyle='--')

# axes labelling
plt.xlabel(f'Is variant present?', fontsize = xaxis_font)
plt.ylabel('MHC-variant binding', fontsize = yaxis_font)

plt.ylim(-3.2, 3) # I know this will be the max value for all variants bc R882H is most common

# x ticks 
plt.xticks(fontsize = xticks_font)
plt.yticks(fontsize = yticks_font)

# save the main figure 
if df == 'netmhc_2reads':    
    plt.savefig(f'/Users/barbarawalkowiak/Desktop/msc_thesis/results/figures/fig5/5b_carriers_present_vs_absent_variants_aggregared_rm_mostcommon{df}.pdf', bbox_inches='tight')
else: # for other dataframes / ways of analyzing, save to the folder in supplmentary figures 
    plt.savefig(f'/Users/barbarawalkowiak/Desktop/msc_thesis/results/supplementary_files/{df}/fig5/5b_carriers_present_vs_absent_variants_aggregated_rm_mostcommon{df}.pdf', bbox_inches='tight')
    

# %%
