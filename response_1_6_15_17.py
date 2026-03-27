# %%
# REVIEWER RESPONSE - RESPONSES 1, 6, 15, 17

# December 2024
# Barbara Walkowiak bw450 

# Script to reproduce the figures included in response to reviewers related to quality control
# responses: 1 (age association), 6 (quality control of de novo calls), 15 (confounding by sex), 17 (distribution of depths)

# input: requires dataframes with scores (NetMHC I)

# %%
# IMPORTS
import warnings
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
from decimal import *
import time 
import csv
import seaborn as sns 
import scipy.stats as stats
from scipy.stats import mannwhitneyu
import sys 
import statsmodels.formula.api as smf
from statsmodels.stats.anova import anova_lm

# %% 
# SET UP

# specify font for plotting 
plt.rcParams.update({'font.sans-serif':'Helvetica'})
# this is to make things editable in Illustrator
plt.rcParams['pdf.fonttype'] = 42 
# stop printing warnings 
warnings.filterwarnings("ignore")
# get current date 
timestr = time.strftime("%Y%m%d") 

# COLORS 
# colors for relative binding 
col0r = '#0AAE37' # this is used to plot the top group (highest score = most immunogenic / strongest binding)
col1r = '#8FE8A7' # this is used to plot the middle group 
col2r = '#B018ED' # this is used to plot the bottom group (does not bind)

# colors for absolute binding (I think we decided we want to do the same colors for the moment)
col0a = '#0AAE37' # this is used to plot the top group (strong binding based on threshold) 
col1a = '#8FE8A7' # this is used to plot the middle group (weak binding based on threshold)
col2a = '#B018ED' # this is used to plot the bottom group (no binding based on threshold)

# for CH specifically for absolute binding I am just using slightly darker colors but overall not changing it much
col0a1 = '#001707' # this is used to plot the top group (strong binding based on threshold) for CH-neg in SF2 
col1a1 = '#037724' # this is used to plot the middle group (weak binding based on threshold) for CH-neg in SF2 
col2a1 = '#580479' # this is used to plot the bottom group (no binding based on threshold) for CH-neg in SF2 

# colors for CH-positive / CH-negative 
col_pos = '#BB0733' # color for CH-positive individuals (~reddish)
col_neg = '#6892ED' # color for CH-negative individuals (blue)
col_neg2 = '#8AAAEF' # lighter color for background (dots)

# MHC I / MHC II
col_mhc1 = '#a0032a'
col_mhc2 = '#e78ea4'

# light grey for background (alternating variants)
col_background = '#EEEEEE'

# FONTS
title_font = 15
xaxis_font = 14
yaxis_font = 14
xticks_font = 12
yticks_font = 12
legend_title = 14
legend_font = 13
text_font = 9

# %% 
# IMPORT REQUIRED DATAFRAMES 

# # file with scores 
netmhc1_df = pd.read_csv('/Users/barbarawalkowiak/Desktop/msc_thesis/results/dataframes/20240907_netmhc1_scores_for_all_var.csv')
# melted file with scores + labels 
netmhc1_df_labels = pd.read_csv('/Users/barbarawalkowiak/Desktop/msc_thesis/results/dataframes/20240907_netmhc1_scores_for_all_var_with_labels.csv')
# specify df name (for saving plots etc)
df = 'netmhc_2reads'

# %%
# FORMATTING AND SETTING UP THE DATAFRAME

# only use variants for which there are enough carriers (10 with MHC)
variants_identified = netmhc1_df_labels.gene_var.unique().tolist() 
variants_identified = [var for var in variants_identified if netmhc1_df_labels[netmhc1_df_labels['gene_var'] == var].shape[0] >= 10]
variants_with_scores = netmhc1_df_labels.CH_variant.unique().tolist()
vars_with_no_carriers = [var for var in variants_with_scores if var not in variants_identified]
vars_with_carriers = [var for var in variants_with_scores if var in variants_identified]
vars_with_carriers = [var for var in vars_with_carriers if len(netmhc1_df[netmhc1_df['gene_var']==var].Person_ID.tolist()) >= 10]

# subset the dataframe 
netmhc1_df_labels_sub = netmhc1_df_labels[netmhc1_df_labels['CH_variant'].isin(vars_with_carriers)]

# specify order of variants to plot (by median score)
order = netmhc1_df_labels_sub.sort_values(by = 'median_score', ascending = False).CH_variant.unique() # order by median score 
order2 = [var.replace('_', '\n') for var in order] # formatting  

# specify levels for CH status (carrier, non-carrier)
netmhc1_df_labels_sub['CH_status'] = netmhc1_df_labels_sub['CH_status'].astype('category') # category 

# formatting 
netmhc1_df_labels_sub['gene_var2'] = netmhc1_df_labels_sub['gene_var'].str.replace('_', '\n')
netmhc1_df_labels_sub['CH_variant2'] = netmhc1_df_labels_sub['CH_variant'].str.replace('_', '\n')

# rename CH_status column descriptions (we want CH-positive / negative to avoid associations with genetic carriers eg germline carriers)
netmhc1_df_labels_sub['CH_status2'] = np.where(netmhc1_df_labels_sub['CH_status']==1, 'CH-positive', 'CH-negative')
netmhc1_df_labels_sub['CH_status2'] = pd.Categorical(netmhc1_df_labels_sub['CH_status2'], categories = ['CH-positive', 'CH-negative']) # category

# %%
# IDENTIFY VARIANTS WHICH CAN BE BOUND WELL AND POORLY (NON-UNIFORM ACROSS THE POPULATION)

# identify the variants values of which span different binding thresholds 
# choose variants where there are at least 5 carriers > -log10(0.5) (strong binding) and at leat 5 carriers is < -log10(2) (no binding)
variants_thresh = [var for var in vars_with_carriers if sorted(netmhc1_df_labels_sub[(netmhc1_df_labels_sub['variable']==f'score_{var}') & (netmhc1_df_labels_sub['gene_var']==var)].log_score.tolist())[-5] > -np.log10(0.5)]
variants_thresh = [var for var in variants_thresh if sorted(netmhc1_df_labels_sub[(netmhc1_df_labels_sub['variable']==f'score_{var}') & (netmhc1_df_labels_sub['gene_var']==var)].log_score.tolist())[4] < -np.log10(2)]
variants_thresh2 = [var.replace('_', '\n') for var in variants_thresh]

# %%
# COUNTS OF CARRIERS + ORDER BY TOTAL NR OF CARRIERS

# create a new dataframe where you will store counts of variants in top and bottom binding groups
df_counts_carriers = pd.DataFrame()

# loop for all variants
for var in order:

    # first, determine how many variants are classified as top / bottom binding in general 
    var_df = netmhc1_df_labels_sub[netmhc1_df_labels_sub['CH_variant']==var] # df with scores for the variant for everyone
    var_df_carriers = var_df[var_df['gene_var']==var] # df with scores for cases with variants only
    
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

# %%
# ANALYSIS ADDRESSING COMMENT 1: AGE ASSOCIATION

# FIGURE 1 OF REPLY 1

# compare the distribution of age b/n carrier predicted to bind a variant well vs poorly 
# this is taking CARRIERS from each group (top / bottom half of the POPULATION by binding) and comparing their scores

# add age to dataframe with labels 
netmhc1_df_labels_age = pd.merge(netmhc1_df_labels, netmhc1_df[['Person_ID', 'age']])
netmhc1_df_labels_age_carriers =netmhc1_df_labels_age[netmhc1_df_labels_age['CH_status']==1].drop_duplicates()
# 9624 cases total 

# %%
# PLOT 
# Compare age distribution in carriers (better vs worse binding)

colors = [col0r, col2r]

plt.figure(figsize = (16, 4))
ax = sns.stripplot(data = netmhc1_df_labels_age_carriers, 
                   dodge = True, jitter = 0.2, x = 'CH_variant', y = 'age', hue = 'group', 
                   palette = colors, size = 2.5, edgecolor = 'black', order = order_by_total, alpha = 0.6)

# GREY BACKGROUND (every other variant)
for i in range(1, len(order_by_total2), 2):
    ax.axvspan(i-0.5, i+0.5, color=col_background, alpha = 0.8)

# ADJUST HOW DODGED THE HUE IS (ie how separate are ppl who are better vs worse in MHC peptide binding)
for i, artist in enumerate(ax.collections):
    # Get the current positions of the points
    offsets = artist.get_offsets()
    dodge_extent = 0
    offsets[:, 0] += (i % 2) * dodge_extent - dodge_extent / 2
    # Update the positions
    artist.set_offsets(offsets)

# find median score and add onto the plot 
for i, category in enumerate(order_by_total):
            
    median_top = netmhc1_df_labels_age_carriers[(netmhc1_df_labels_age_carriers['gene_var'] == f'{category}') & (netmhc1_df_labels_age_carriers['group'] == 'top half')].age.median()
    median_bottom = netmhc1_df_labels_age_carriers[(netmhc1_df_labels_age_carriers['gene_var'] == f'{category}') & (netmhc1_df_labels_age_carriers['group'] == 'bottom half')].age.median()

    # Plot text for each hue group
    plt.text(i, median_top, '-', ha='right', va='center', fontsize=30, fontweight='bold', color = col0r)
    plt.text(i, median_bottom, '-', ha='left', va='center', fontsize=30, fontweight='bold', color = col2r)

# add Mann-Whitney U test between groups (non-parametric t-test alternative)
p_values = []    

for i, category in enumerate(order_by_total):
    
    category_data = netmhc1_df_labels_age_carriers[netmhc1_df_labels_age_carriers['gene_var'] == f'{category}']
    max_value = category_data['age'].max()
    
    age_top = netmhc1_df_labels_age_carriers[(netmhc1_df_labels_age_carriers['gene_var'] == f'{category}') & (netmhc1_df_labels_age_carriers['group'] == 'top half')].age.tolist()
    age_bottom = netmhc1_df_labels_age_carriers[(netmhc1_df_labels_age_carriers['gene_var'] == f'{category}') & (netmhc1_df_labels_age_carriers['group'] == 'bottom half')].age.tolist()
    
    statistic, p_value = mannwhitneyu(age_top, age_bottom)
    p_values.append(p_value)

    # adjust to the number of tests performed 
    significance = ''
    if p_value > (0.05 / len(order_by_total2)):
        significance = 'ns'
    elif p_value < (0.01 / len(order_by_total2)):
        significance = '**'
    else:
        significance = '*'
    plt.text(i, 2+max_value, significance, ha='center', va='center', fontsize=text_font)

# axes axes labels 
plt.xlabel(f'CH hotspot variant', fontsize = xaxis_font)
plt.ylabel('Age', fontsize = yaxis_font)

plt.xlim(-0.5, len(order_by_total2)-0.5)
plt.ylim(20, 85) # I know this will be the max value for all variants bc R882H is most common

# BOLD X TICKS for which you have binding predictions that span both no binding and good binding
new_labels = []
for lab in order_by_total2:
    if '\n' in lab:
        gene, variant = lab.split('\n', 1)
        gene = gene.replace('_', r'\_')
        variant = variant.replace('_', r'\_')

        if lab in variants_thresh2:
            formatted = rf'$\it{{{gene}}}$' + '\n' + variant
        else:
            formatted = rf'$\it{{{gene}}}$' + '\n' + variant
        new_labels.append(formatted)
    else:
        new_labels.append(lab)
ax.set_xticklabels(new_labels, fontsize=xticks_font, rotation=90)

ax = plt.gca()  # Get current axes
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.yticks(fontsize = yticks_font)
plt.xticks(fontsize = xticks_font-2)

# specify legend 
legend = plt.legend(['stronger binding\n(top half)', 'weaker binding\n(bottom half)'],
                    markerscale = 3, loc = 'lower right', fontsize = legend_font, frameon = False, handletextpad = 0.2)
for legend_handle in legend.legendHandles:
    legend_handle.set_alpha(1)

# save the main figure (NB I don't think it helped)
plt.savefig(f'/Users/barbarawalkowiak/Desktop/msc_thesis/results/reply/age/1_age_distribution_topvsbottom_{df}_jitter.pdf', bbox_inches='tight')

# %%
# PLOT THE FRACTION OF CH-POSITIVE INDIVIDUALS 
# FIGURE 2A, RESPONSE 1
    
# plot for all variants together (aggregated)
netmhc1_df_labels_age = netmhc1_df_labels_age.drop_duplicates()

bins = pd.qcut(netmhc1_df_labels_age.age, q=4, labels=False)

# add column based on age bin
netmhc1_df_labels_age['age_bin'] =  bins 

# count CH-positive and CH-negative in different age_bins
values_to_count = [0, 1, 2, 3] # 4 age bins, labelled 0, 1, 2, 3

counts_dt = pd.DataFrame(netmhc1_df_labels_age.groupby(['group', 'CH_status', 'age_bin']).size()).reset_index()
counts_dt.columns.values[3] = 'counts'

# pivot dt 
counts_dt2 = counts_dt.pivot(index=['group', 'age_bin'], columns=['CH_status'], values='counts').reset_index()
counts_dt2.columns.values[2] = 'ch_neg'
counts_dt2.columns.values[3] = 'ch_pos'
counts_dt2['group'] = counts_dt2['group'].astype('category')
counts_dt2['group'] = counts_dt2['group'].cat.reorder_categories(['top half', 'bottom half'])

counts_dt2['fraction_CH'] = counts_dt2['ch_pos'] / (counts_dt2['ch_neg'] + counts_dt2['ch_pos'])
counts_dt2['std_err'] = np.sqrt(counts_dt2['fraction_CH'] * (1 - counts_dt2['fraction_CH'])) / np.sqrt(counts_dt2['ch_neg'] + counts_dt2['ch_pos']) 

counts_dt2['fraction_CH'].fillna(0)
counts_dt2['std_err'].fillna(0)

# create the plot 
ax = sns.stripplot(data = counts_dt2, x = 'age_bin', hue = 'group', y = 'fraction_CH', palette = [col0r, col2r],
                size = 10, alpha = 1, dodge = False, jitter = False)

plt.xlabel(f'Age bin', fontsize = xaxis_font)
plt.ylabel('Fraction CH-positive individuals', fontsize = yaxis_font)
plt.title(f'All variants', fontsize = title_font)

plt.xticks(fontsize = xticks_font)
plt.yticks(fontsize = yticks_font)
plt.ylim(0, counts_dt2['fraction_CH'].max() * 1.2)

# replace ticks with what the actual values for the age bins
x = [0, 1, 2, 3]
min0 = netmhc1_df_labels_age[netmhc1_df_labels_age['age_bin']==0].age.min()
max0 = netmhc1_df_labels_age[netmhc1_df_labels_age['age_bin']==0].age.max()
min1 = netmhc1_df_labels_age[netmhc1_df_labels_age['age_bin']==1].age.min()
max1 = netmhc1_df_labels_age[netmhc1_df_labels_age['age_bin']==1].age.max()
min2 = netmhc1_df_labels_age[netmhc1_df_labels_age['age_bin']==2].age.min()
max2 = netmhc1_df_labels_age[netmhc1_df_labels_age['age_bin']==2].age.max()
min3 = netmhc1_df_labels_age[netmhc1_df_labels_age['age_bin']==3].age.min()
max3 = netmhc1_df_labels_age[netmhc1_df_labels_age['age_bin']==3].age.max()

ax.set_xticks(x)
ax.set_xticklabels([f'{min0}-{max0}', f'{min1}-{max1}', f'{min2}-{max2}', f'{min3}-{max3}'])

# add confidence intervals for each age bin 
for i, ab in enumerate(values_to_count):

    n_top = counts_dt2[(counts_dt2['age_bin']==ab) & (counts_dt2['group']=='top half')]['fraction_CH'].iloc[0]
    n_bottom = counts_dt2[(counts_dt['age_bin']==ab) & (counts_dt2['group']=='bottom half')]['fraction_CH'].iloc[0]
    std_error_top = counts_dt2[(counts_dt2['age_bin']==ab) & (counts_dt2['group']=='top half')]['std_err'].iloc[0]
    std_error_bottom = counts_dt2[(counts_dt2['age_bin']==ab) & (counts_dt2['group']=='bottom half')]['std_err'].iloc[0]

    # add error bars 
    plt.errorbar(x=i, y=n_top, yerr=[[std_error_top], [std_error_top]], fmt='none', capsize = 0.25, color='black', capthick=0)
    plt.errorbar(x=i, y=n_bottom, yerr=[[std_error_bottom], [std_error_bottom]], fmt='none', capsize = 0.25, color='black', capthick=0)

legend = plt.legend(loc = 'upper left', fontsize = legend_font, frameon = False, handletextpad = 0.2)

ax = plt.gca()  # Get current axes
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.savefig(f'/Users/barbarawalkowiak/Desktop/msc_thesis/results/reply/age/3_age_bins_fraction_ch_positive_allVariants.pdf', bbox_inches='tight')

# %%
# PLOT THE FRACTION OF CH-POSITIVE INDIVIDUALS 
# FIGURE 2B, RESPONSE 1
rows = 8
cols = 5

with PdfPages('/Users/barbarawalkowiak/Desktop/msc_thesis/results/reply/age/3_age_bins_fraction_ch_positive.pdf') as pdf:

    fig, axes = plt.subplots(rows, cols, figsize=(28, 40)) 
    axes = axes.flatten()

    for a, ax in enumerate(axes[:len(order_by_total)]):
        
        var = order_by_total[a]
        parts = var.split('_')[0:2]
        part1, part2 = parts[0], parts[1]
        part1 = part1.replace('_', r'\_')
        part2 = part2.replace('_', r'\_')
        if var == "ALL_VARIANTS":
            var_name = part1 + '\n' + part2
        else:
            var_name = rf'$\mathit{{{part1}}}$' + ' ' + part2
        data = netmhc1_df_labels_age[ (netmhc1_df_labels_age['variable']==f'score_{var}')]
        bins = pd.qcut(data.age, q=4, labels=False)

        # add column based on age bin
        data['age_bin'] =  bins 

        # count CH-positive and CH-negative in different age_bins
        values_to_count = [0, 1, 2, 3] # 4 age bins, labelled 0, 1, 2, 3

        counts_dt = pd.DataFrame(data.groupby(['group', 'CH_status', 'age_bin']).size()).reset_index()
        counts_dt.columns.values[3] = 'count'
        counts_dt = counts_dt.pivot(index=['group', 'age_bin'], columns=['CH_status'], values='count').reset_index()
        counts_dt.columns.values[2] = 'ch_neg'
        counts_dt.columns.values[3] = 'ch_pos'
        counts_dt['group'] = counts_dt['group'].astype('category')
        counts_dt['group'] = counts_dt['group'].cat.reorder_categories(['top half', 'bottom half'])

        counts_dt['fraction_CH'] = counts_dt['ch_pos'] / (counts_dt['ch_neg'] + counts_dt['ch_pos'])
        counts_dt['std_err'] = np.sqrt(counts_dt['fraction_CH'] * (1 - counts_dt['fraction_CH'])) / np.sqrt(counts_dt['ch_neg'] + counts_dt['ch_pos']) 

        counts_dt['fraction_CH'].fillna(0)
        counts_dt['std_err'].fillna(0)

        ax = sns.stripplot(data = counts_dt, x = 'age_bin', hue = 'group', y = 'fraction_CH', palette = [col0r, col2r],
                        size = 11, alpha = 1, dodge = False, jitter = False, ax = ax)

        ax.set_xlabel(f'Age bin', fontsize = xaxis_font+3)
        ax.set_ylabel('Fraction CH-positive individuals', fontsize = yaxis_font+3)
        ax.set_title(f'{var_name}', fontsize = title_font+4)

        ax.set_ylim(0, counts_dt['fraction_CH'].max() * 1.2)

        # replace ticks with what the actual values for the age bins
        x = [0, 1, 2, 3]
        min0 = data[data['age_bin']==0].age.min()
        max0 = data[data['age_bin']==0].age.max()
        min1 = data[data['age_bin']==1].age.min()
        max1 = data[data['age_bin']==1].age.max()
        min2 = data[data['age_bin']==2].age.min()
        max2 = data[data['age_bin']==2].age.max()
        min3 = data[data['age_bin']==3].age.min()
        max3 = data[data['age_bin']==3].age.max()

        ax.set_xticks(x)
        ax.set_xticklabels([f'{min0}-{max0}', f'{min1}-{max1}', f'{min2}-{max2}', f'{min3}-{max3}'], fontsize = xticks_font+1)

        # add confidence intervals for each age bin 
        for i, ab in enumerate(values_to_count):

            n_top = counts_dt[(counts_dt['age_bin']==ab) & (counts_dt['group']=='top half')]['fraction_CH'].iloc[0]
            n_bottom = counts_dt[(counts_dt['age_bin']==ab) & (counts_dt['group']=='bottom half')]['fraction_CH'].iloc[0]
            std_error_top = counts_dt[(counts_dt['age_bin']==ab) & (counts_dt['group']=='top half')]['std_err'].iloc[0]
            std_error_bottom = counts_dt[(counts_dt['age_bin']==ab) & (counts_dt['group']=='bottom half')]['std_err'].iloc[0]

            # add error bars 
            ax.errorbar(x=i, y=n_top, yerr=[[std_error_top], [std_error_top]], fmt='none', capsize = 0.25, color='black', capthick=0)
            ax.errorbar(x=i, y=n_bottom, yerr=[[std_error_bottom], [std_error_bottom]], fmt='none', capsize = 0.25, color='black', capthick=0)

        ax.legend(loc = 'upper left', fontsize = legend_font+4, frameon = False, handletextpad = 0.2)
        
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    for j in range(40, len(axes)):
        fig.delaxes(axes[j])

    plt.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)

# %%
# ANALYSIS ADDRESSING RESPONSE 6: QUALITY CONTROL OF DE NOVO CALLS 
# FIGURES 1, 2 AND 3 OF RESPONSE 6

# FIGURE 1 OF RESPONSE 6: COMPARISON OF NR OF READS TO GU ET AL
# comparison to number of calls from Gu et al

variants_compare_gu = ['DNMT3A R882H','DNMT3A R882C', 'DNMT3A Y735C', 'DNMT3A P904L', 'DNMT3A 736H', 
'DNMT3A R771*', 'DNMT3A R736C', 'DNMT3A R598*', 'DNMT3A R326C', 'DNMT3A R729W', 
'DNMT3A R320*', 'GNB1 K57E', 'IDH2 R140Q', 'JAK2 V617F', 'MPL W515L', 'NRAS G12D', 
'SF3B1 K666N', 'SF3B1 K700E', 'SRSF2 P95H', 'SRSF2 P95L', 'SRSF2 P95R']

variant_counts = pd.read_csv('/Users/barbarawalkowiak/Desktop/msc_thesis/results/reply/qc/variant_counts.csv')

# load the data from Gu et al 2023 (subsetted supplementary table 4 from the paper)
variant_counts_gu = pd.read_csv('/Users/barbarawalkowiak/Desktop/msc_thesis/results/reply/qc/varIDcounts_gu2023.csv')
variant_counts_gu['study'] = 'Gu et al'
variant_counts['study'] = 'This study'

# concat both dataframes 
variant_counts_compare = pd.concat([variant_counts, variant_counts_gu], axis = 0)
variant_counts_compare

order_vars = variant_counts_compare[variant_counts_compare['study']=='This study'].sort_values(by = 'total', ascending = False)['varID'].tolist()

# plot 
sns.barplot(data = variant_counts_compare, x = 'varID', y = 'total', hue = 'study', order = order_vars)
# Add labels and title
plt.xlabel('Hotspot', fontsize = 11)
plt.ylabel('Total number of cases', fontsize = 11)
plt.legend(loc = 'upper right', fontsize = 10, frameon = False, handletextpad = 0.2)
plt.xticks(rotation = 90, fontsize = 9)
ax = plt.gca()  # Get current axes
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.savefig(f'/Users/barbarawalkowiak/Desktop/msc_thesis/results/reply/qc/5_cf_gu2023_total.pdf', bbox_inches='tight')
plt.close()

# %%
# FIGURE 3A OF RESPONSE 6: PLOT FOR ALL VARIANTS AGGREGATED TOGETHER 
    
bins = pd.qcut(netmhc1_df.age, q=4, labels=False)

# add column based on age bin
netmhc1_df['age_bin'] =  bins 

# count CH-positive and CH-negative in different age_bins
values_to_count = [0, 1, 2, 3] # 4 age bins, labelled 0, 1, 2, 3

counts_dt = pd.DataFrame(netmhc1_df.groupby(['ch_status', 'age_bin']).size()).reset_index()
counts_dt.columns.values[2] = 'count'
counts_dt = counts_dt.pivot(index=['age_bin'], columns=['ch_status'], values='count').reset_index()
counts_dt.columns.values[1] = 'ch_neg'
counts_dt.columns.values[2] = 'ch_pos'

counts_dt['fraction_CH'] = counts_dt['ch_pos'] / (counts_dt['ch_neg'] + counts_dt['ch_pos'])
counts_dt['std_err'] = np.sqrt(counts_dt['fraction_CH'] * (1 - counts_dt['fraction_CH'])) / np.sqrt(counts_dt['ch_neg'] + counts_dt['ch_pos']) 

counts_dt['fraction_CH'].fillna(0)
counts_dt['std_err'].fillna(0)    

# plot
ax = sns.stripplot(data = counts_dt, x = 'age_bin', y = 'fraction_CH', 
                size = 10, alpha = 1, dodge = False, jitter = False)

ax.set_xlabel(f'Age bin', fontsize = xaxis_font+3)
ax.set_ylabel('Fraction CH-positive individuals', fontsize = yaxis_font+3)
ax.set_title(f'All variants', fontsize = title_font+3)

ax.set_ylim(0, counts_dt['fraction_CH'].max() * 1.2)

# replace ticks with what the actual values for the age bins
x = [0, 1, 2, 3]
min0 = data[data['age_bin']==0].age.min()
max0 = data[data['age_bin']==0].age.max()
min1 = data[data['age_bin']==1].age.min()
max1 = data[data['age_bin']==1].age.max()
min2 = data[data['age_bin']==2].age.min()
max2 = data[data['age_bin']==2].age.max()
min3 = data[data['age_bin']==3].age.min()
max3 = data[data['age_bin']==3].age.max()

ax.set_xticks(x)
ax.set_xticklabels([f'{min0}-{max0}', f'{min1}-{max1}', f'{min2}-{max2}', f'{min3}-{max3}'], fontsize = xticks_font+3)

# add confidence intervals for each age bin 
for i, ab in enumerate(values_to_count):

    n = counts_dt[(counts_dt['age_bin']==ab)]['fraction_CH'].iloc[0]
    std_error = counts_dt[(counts_dt['age_bin']==ab)]['std_err'].iloc[0]

    # add error bars 
    ax.errorbar(x=i, y=n, yerr=[[std_error], [std_error]], fmt='none', capsize = 0.5, color='black', capthick=0)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

for j in range(40, len(axes)):
    fig.delaxes(axes[j])

plt.tight_layout()
plt.savefig('/Users/barbarawalkowiak/Desktop/msc_thesis/results/reply/qc/QC_age_dependence_allVariants.pdf', bbox_inches='tight')

# %%
# FIGURE 3B OF RESPONSE 6: PLOT FOR RARE VARIANTS AGGREGATED TOGETHER 
# For 'rare', I selected variants with fewer than 50 individuals 

bins = pd.qcut(netmhc1_df.age, q=4, labels=False)

# add column based on age bin
netmhc1_df['age_bin'] =  bins 

# count CH-positive and CH-negative in different age_bins
values_to_count = [0, 1, 2, 3] # 4 age bins, labelled 0, 1, 2, 3

# add CH-status for rare variants
rare_variants = ['IDH1_R132H', 'DNMT3A_R882P', 'DNMT3A_R882L', 'DNMT3A_R326G', 'DNMT3A_R736G', 'DNMT3A_R736L', 'DNMT3A_Y735F', 'MPL_W515L', 'TP53_R273H', 'TP53_R175H', 
'DNMT3A_P904R', 'NRAS_G12D', 'KRAS_G12D', 'IDH2_R172K', 'KRAS_G12S']
netmhc1_df['ch_status_rare'] = np.where(netmhc1_df['gene_var'].isin(rare_variants), 1, 0)

counts_dt = pd.DataFrame(netmhc1_df.groupby(['ch_status_rare', 'age_bin']).size()).reset_index()
counts_dt.columns.values[2] = 'count'
counts_dt = counts_dt.pivot(index=['age_bin'], columns=['ch_status_rare'], values='count').reset_index()
counts_dt.columns.values[1] = 'ch_neg'
counts_dt.columns.values[2] = 'ch_pos'

counts_dt['fraction_CH'] = counts_dt['ch_pos'] / (counts_dt['ch_neg'] + counts_dt['ch_pos'])
counts_dt['std_err'] = np.sqrt(counts_dt['fraction_CH'] * (1 - counts_dt['fraction_CH'])) / np.sqrt(counts_dt['ch_neg'] + counts_dt['ch_pos']) 

counts_dt['fraction_CH'].fillna(0)
counts_dt['std_err'].fillna(0)    

# plot 
ax = sns.stripplot(data = counts_dt, x = 'age_bin', y = 'fraction_CH', 
                size = 10, alpha = 1, dodge = False, jitter = False)

ax.set_xlabel(f'Age bin', fontsize = xaxis_font+3)
ax.set_ylabel('Fraction CH-positive individuals', fontsize = yaxis_font+3)
ax.set_title(f'Rare variants', fontsize = title_font+3)

ax.set_ylim(0, counts_dt['fraction_CH'].max() * 1.2)

# replace ticks with what the actual values for the age bins
x = [0, 1, 2, 3]
min0 = data[data['age_bin']==0].age.min()
max0 = data[data['age_bin']==0].age.max()
min1 = data[data['age_bin']==1].age.min()
max1 = data[data['age_bin']==1].age.max()
min2 = data[data['age_bin']==2].age.min()
max2 = data[data['age_bin']==2].age.max()
min3 = data[data['age_bin']==3].age.min()
max3 = data[data['age_bin']==3].age.max()

ax.set_xticks(x)
ax.set_xticklabels([f'{min0}-{max0}', f'{min1}-{max1}', f'{min2}-{max2}', f'{min3}-{max3}'], fontsize = xticks_font+3)

# add confidence intervals for each age bin 
for i, ab in enumerate(values_to_count):

    n = counts_dt[(counts_dt['age_bin']==ab)]['fraction_CH'].iloc[0]
    std_error = counts_dt[(counts_dt['age_bin']==ab)]['std_err'].iloc[0]

    # add error bars 
    ax.errorbar(x=i, y=n, yerr=[[std_error], [std_error]], fmt='none', capsize = 0.5, color='black', capthick=0)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

for j in range(40, len(axes)):
    fig.delaxes(axes[j])

plt.tight_layout()
plt.savefig('/Users/barbarawalkowiak/Desktop/msc_thesis/results/reply/qc/QC_age_dependence_rareVariants.pdf', bbox_inches='tight')

# %%
# FIGURE 3C OF RESPONSE 6: PLOT AGE DEPENDENCE FOR EACH VARIANT
# PLOT THE REL B/N VARIANT AND AGE W/O SPLITTING INTO GROUPS

rows = 8
cols = 5

with PdfPages('/Users/barbarawalkowiak/Desktop/msc_thesis/results/reply/qc/QC_age_dependence.pdf') as pdf:

    fig, axes = plt.subplots(rows, cols, figsize=(25, 40)) 
    axes = axes.flatten()

    for a, ax in enumerate(axes[:len(order_by_total)]):
        
        var = order_by_total[a]
        var_name = var.replace('_', ' ')
        data = netmhc1_df_labels_age[ (netmhc1_df_labels_age['variable']==f'score_{var}')]
        bins = pd.qcut(data.age, q=4, labels=False)

        # add column based on age bin
        data['age_bin'] =  bins 

        # count CH-positive and CH-negative in different age_bins
        values_to_count = [0, 1, 2, 3] # 4 age bins, labelled 0, 1, 2, 3

        counts_dt = pd.DataFrame(data.groupby(['CH_status', 'age_bin']).size()).reset_index()
        counts_dt.columns.values[2] = 'count'
        counts_dt = counts_dt.pivot(index=['age_bin'], columns=['CH_status'], values='count').reset_index()
        counts_dt.columns.values[1] = 'ch_neg'
        counts_dt.columns.values[2] = 'ch_pos'

        counts_dt['fraction_CH'] = counts_dt['ch_pos'] / (counts_dt['ch_neg'] + counts_dt['ch_pos'])
        counts_dt['std_err'] = np.sqrt(counts_dt['fraction_CH'] * (1 - counts_dt['fraction_CH'])) / np.sqrt(counts_dt['ch_neg'] + counts_dt['ch_pos']) 

        counts_dt['fraction_CH'].fillna(0)
        counts_dt['std_err'].fillna(0)

        sns.stripplot(data = counts_dt, x = 'age_bin', y = 'fraction_CH', 
                    size = 10, alpha = 1, dodge = False, jitter = False, ax = ax)

        ax.set_xlabel(f'Age bin', fontsize = xaxis_font+3)
        ax.set_ylabel('Fraction CH-positive individuals', fontsize = yaxis_font+3)
        ax.set_title(f'{var_name}', fontsize = title_font +4)

        ax.set_ylim(0, counts_dt['fraction_CH'].max() * 1.2)

        # replace ticks with what the actual values for the age bins
        x = [0, 1, 2, 3]
        min0 = data[data['age_bin']==0].age.min()
        max0 = data[data['age_bin']==0].age.max()
        min1 = data[data['age_bin']==1].age.min()
        max1 = data[data['age_bin']==1].age.max()
        min2 = data[data['age_bin']==2].age.min()
        max2 = data[data['age_bin']==2].age.max()
        min3 = data[data['age_bin']==3].age.min()
        max3 = data[data['age_bin']==3].age.max()

        ax.set_xticks(x)
        ax.set_xticklabels([f'{min0}-{max0}', f'{min1}-{max1}', f'{min2}-{max2}', f'{min3}-{max3}'], fontsize = xticks_font+3)
        
        # add confidence intervals for each age bin 
        for i, ab in enumerate(values_to_count):

            n = counts_dt[(counts_dt['age_bin']==ab)]['fraction_CH'].iloc[0]
            std_error = counts_dt[(counts_dt['age_bin']==ab)]['std_err'].iloc[0]

            # add error bars 
            ax.errorbar(x=i, y=n, yerr=[[std_error], [std_error]], fmt='none', capsize = 0.5, color='black', capthick=0)
        
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    for j in range(40, len(axes)):
        fig.delaxes(axes[j])

    plt.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)

# %%
# FIGURE 2 OF RESPONSE 6: PLOT FOR AGE DEPENDENCE FOR DIFFERENT NR OF READS

values_to_count = [0, 1, 2, 3] # 4 age bins, labelled 0, 1, 2, 3
dt_2reads = netmhc1_df[(netmhc1_df['var_depth']>=2)][['gene_var', 'Person_ID', 'VAF', 'ch_status', 'age', 'age_bin']]
ages_2reads = dt_2reads.groupby(['ch_status', 'age_bin']).size().reset_index(name='count')
ages_2reads['cat'] = '>= 2 reads'
dt_3reads = netmhc1_df[(netmhc1_df['var_depth']>=3)][['gene_var', 'Person_ID', 'VAF', 'ch_status', 'age', 'age_bin']]
ages_3reads = dt_3reads.groupby(['ch_status', 'age_bin']).size().reset_index(name='count')
ages_3reads['cat'] = '>= 3 reads'
dt_4reads = netmhc1_df[(netmhc1_df['var_depth']>=4)][['gene_var', 'Person_ID', 'VAF', 'ch_status', 'age', 'age_bin']]
ages_4reads = dt_4reads.groupby(['ch_status', 'age_bin']).size().reset_index(name='count')
ages_4reads['cat'] = '>= 4 reads'
ages_nrreads = pd.concat([ages_2reads, ages_3reads, ages_4reads])
ages_nrreads = ages_nrreads.reset_index()

# calculate fraction (total in the bin)
age_bins_counts = pd.DataFrame(netmhc1_df.age_bin.value_counts()).reset_index()
age_bins_counts = age_bins_counts.rename(columns={'count': 'total'}) 
ages_nrreads = pd.merge(ages_nrreads, age_bins_counts)

ages_nrreads['fraction_CH'] = ages_nrreads['count'] / ages_nrreads['total']
ages_nrreads['std_err'] = np.sqrt(ages_nrreads['fraction_CH'] * (1 - ages_nrreads['fraction_CH'])) / np.sqrt(ages_nrreads['total']) 

# plot 
ax = sns.stripplot(data = ages_nrreads, x = 'age_bin', hue = 'cat', y = 'fraction_CH', size = 7, alpha = 1, 
                    legend = True, jitter = False)

plt.xlabel(f'Age bin', fontsize = xaxis_font)
plt.ylabel('Fraction od CH-positive individuals', fontsize = yaxis_font)
plt.title(f'All variants', fontsize = title_font)

# add error bars 
for i, ab in enumerate(values_to_count):

    n1 = ages_nrreads[(ages_nrreads['age_bin']==ab) & (ages_nrreads['cat']=='>= 2 reads')]['fraction_CH'].iloc[0]
    n2 = ages_nrreads[(ages_nrreads['age_bin']==ab) & (ages_nrreads['cat']=='>= 3 reads')]['fraction_CH'].iloc[0]
    n3 = ages_nrreads[(ages_nrreads['age_bin']==ab) & (ages_nrreads['cat']=='>= 4 reads')]['fraction_CH'].iloc[0]
    
    std_error_1 = ages_nrreads[(ages_nrreads['age_bin']==ab) & (ages_nrreads['cat']=='>= 2 reads')]['std_err'].iloc[0]
    std_error_2 = ages_nrreads[(ages_nrreads['age_bin']==ab) & (ages_nrreads['cat']=='>= 3 reads')]['std_err'].iloc[0]
    std_error_3 = ages_nrreads[(ages_nrreads['age_bin']==ab) & (ages_nrreads['cat']=='>= 4 reads')]['std_err'].iloc[0]
    
    # add error bars 
    plt.errorbar(x=i, y=n1, yerr=[[std_error_1], [std_error_1]], fmt='none', capsize = 0.5, color='black', capthick=0)
    plt.errorbar(x=i, y=n2, yerr=[[std_error_2], [std_error_2]], fmt='none', capsize = 0.5, color='black', capthick=0)
    plt.errorbar(x=i, y=n3, yerr=[[std_error_3], [std_error_3]], fmt='none', capsize = 0.5, color='black', capthick=0)
    

# replace ticks with what the actual values for the age bins
x = [0, 1, 2, 3]
min0 = data[data['age_bin']==0].age.min()
max0 = data[data['age_bin']==0].age.max()
min1 = data[data['age_bin']==1].age.min()
max1 = data[data['age_bin']==1].age.max()
min2 = data[data['age_bin']==2].age.min()
max2 = data[data['age_bin']==2].age.max()
min3 = data[data['age_bin']==3].age.min()
max3 = data[data['age_bin']==3].age.max()

ax.set_xticks(x)
ax.set_xticklabels([f'{min0}-{max0}', f'{min1}-{max1}', f'{min2}-{max2}', f'{min3}-{max3}'])

plt.xticks(fontsize = xticks_font)
plt.yticks(fontsize = yticks_font)
plt.legend(loc = 'upper left', fontsize = legend_font, frameon = False, handletextpad = 0.2)
plt.ylim(bottom = 0)
ax = plt.gca()  # Get current axes
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.savefig(f'/Users/barbarawalkowiak/Desktop/msc_thesis/results/reply/qc/3_QC_age_dependence_allVars_fraction.pdf', bbox_inches='tight')

# %%
# ANALYSIS ADDRESSING RESPONSE 17 - SEX AND ANCESTRY DISTRIBUTION

# add sex and ancestry data 
pheno_df = pd.read_csv('/Users/barbarawalkowiak/Desktop/msc_thesis/data/pheno_info20230504.tsv', sep = '\t')
ancestry_df = pd.read_csv('/Users/barbarawalkowiak/Desktop/msc_thesis/data/ukb_iadmix_ancestry_props.csv')

# make sure we have consistent naming for Person_ID column 
pheno_df = pheno_df.rename(columns={'ID_v0': 'Person_ID'}) 
ancestry_df = ancestry_df.rename(columns={'ukbid': 'Person_ID'}) 

# clean up datasets: IDs < 0 participants resigned from the study 
pheno_df = pheno_df[pheno_df.Person_ID > 0] # only retain participants with IDs greater than 0
ancestry_df = ancestry_df[ancestry_df.Person_ID > 0] # only retain participants with IDs greater than 0

# check how many people you have the data for
ids_pheno_data = pheno_df.Person_ID.unique()
ids_ancestry_data = ancestry_df.Person_ID.unique()

print('Number of people with phenotypic data available:', len(ids_pheno_data))
print('Number of people with ancestry data available:', len(ids_ancestry_data))

# identify people for whom both are available
ids_pheno_ancestry = set(ids_pheno_data).intersection(set(ids_ancestry_data)) # intersection 
print('Number of people with phenotypic and ancestry data available:', len(ids_pheno_ancestry))

# %%    
# merge with pheno data 
netmhc1_df_labels_age_carriers_pheno = pd.merge(netmhc1_df_labels_age_carriers, pheno_df, on = 'Person_ID')
netmhc1_df_labels_age_carriers_ancestry = pd.merge(netmhc1_df_labels_age_carriers, ancestry_df, on = 'Person_ID')

# add sex 
counts_sex = netmhc1_df_labels_age_carriers_pheno.value_counts(['gene_var', 'Sex_v0']).reset_index(name = 'count')
counts_sex['gene_var2'] = counts_sex['gene_var'].str.replace('_', ' ')

# %%
# REPONSE 17: distirbution of scores (males and females - carriers of each variant)

# Compare core distribution in carriers (better vs worse binding)
colors = [col0r, col2r]

plt.figure(figsize = (16, 4))
ax = sns.stripplot(data = netmhc1_df_labels_age_carriers_pheno, 
                   dodge = True, jitter = False, x = 'gene_var', y = 'log_score', hue = 'Sex_v0', 
                size = 2.5, edgecolor = 'black', order = order_by_total, alpha = 0.6)

# GREY BACKGROUND (every other variant)
for i in range(1, len(order_by_total2), 2):
    ax.axvspan(i-0.5, i+0.5, color=col_background, alpha = 0.8)

# ADJUST HOW DODGED THE HUE IS (ie how separate are ppl who are better vs worse in MHC peptide binding)
for i, artist in enumerate(ax.collections):
    # Get the current positions of the points
    offsets = artist.get_offsets()
    dodge_extent = 0
    offsets[:, 0] += (i % 2) * dodge_extent - dodge_extent / 2
    # Update the positions
    artist.set_offsets(offsets)

# find median score and add onto the plot 
for i, category in enumerate(order_by_total):
            
    median_f = netmhc1_df_labels_age_carriers_pheno[(netmhc1_df_labels_age_carriers_pheno['gene_var'] == f'{category}') & (netmhc1_df_labels_age_carriers_pheno['Sex_v0'] == 'Female')].log_score.median()
    median_m = netmhc1_df_labels_age_carriers_pheno[(netmhc1_df_labels_age_carriers_pheno['gene_var'] == f'{category}') & (netmhc1_df_labels_age_carriers_pheno['Sex_v0'] == 'Male')].log_score.median()

    # Plot text for each hue group
    plt.text(i, median_f, '-', ha='right', va='center', fontsize=30, fontweight='bold', color = '#1f77b4')
    plt.text(i, median_m, '-', ha='left', va='center', fontsize=30, fontweight='bold', color = '#FF7F0E')

# add Mann-Whitney U test between groups (non-parametric t-test alternative)
p_values = []    

for i, category in enumerate(order_by_total):
    
    category_data = netmhc1_df_labels_age_carriers_pheno[netmhc1_df_labels_age_carriers_pheno['gene_var'] == f'{category}']
    max_value = category_data['log_score'].max()
    
    score_f = netmhc1_df_labels_age_carriers_pheno[(netmhc1_df_labels_age_carriers_pheno['gene_var'] == f'{category}') & (netmhc1_df_labels_age_carriers_pheno['Sex_v0'] == 'Female')].log_score.tolist()
    score_m = netmhc1_df_labels_age_carriers_pheno[(netmhc1_df_labels_age_carriers_pheno['gene_var'] == f'{category}') & (netmhc1_df_labels_age_carriers_pheno['Sex_v0'] == 'Male')].log_score.tolist()
    
    statistic, p_value = mannwhitneyu(score_f, score_m)
    p_values.append(p_value)

    # adjust to the number of tests performed 
    significance = ''
    if p_value > (0.05 / len(order_by_total2)):
        significance = 'ns'
    elif p_value < (0.01 / len(order_by_total2)):
        significance = '**'
    else:
        significance = '*'
    plt.text(i, 0.5 + max_value, significance, ha='center', va='center', fontsize=text_font)

# axes axes labels 
plt.xlabel(f'CH hotspot variant', fontsize = xaxis_font)
plt.ylabel('MHC-variant binding score', fontsize = yaxis_font)

plt.xlim(-0.5, len(order_by_total2)-0.5)
plt.ylim(-2, 3)

ax = plt.gca()  # Get current axes
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.yticks(fontsize = yticks_font)
plt.xticks(fontsize = xticks_font-2, rotation = 90)

# specify legend 
legend = plt.legend(loc = 'lower right', fontsize = legend_font, markerscale = 3, frameon = False, handletextpad = 0.2)
for legend_handle in legend.legendHandles:
    legend_handle.set_alpha(1)

# save the main figure 
plt.savefig(f'/Users/barbarawalkowiak/Desktop/msc_thesis/results/reply/qc/5_score_dist_by_sex.pdf', bbox_inches='tight')

# %%
# is there a significant difference to the distribution of sex in the population in the UKB overall (female-dominated generally)?
print('Nr of carriers with EUR ancestry:', len(netmhc1_df_labels_age_carriers_ancestry[netmhc1_df_labels_age_carriers_ancestry['jpt_main_ancestry']=='EUR'].Person_ID.unique()))
print('Nr of carriers with Mixed ancestry:', len(netmhc1_df_labels_age_carriers_ancestry[netmhc1_df_labels_age_carriers_ancestry['jpt_main_ancestry']=='Mixed'].Person_ID.unique()))
print('Nr of carriers with SAS ancestry:', len(netmhc1_df_labels_age_carriers_ancestry[netmhc1_df_labels_age_carriers_ancestry['jpt_main_ancestry']=='SAS'].Person_ID.unique()))
print('Nr of carriers with AFR ancestry:', len(netmhc1_df_labels_age_carriers_ancestry[netmhc1_df_labels_age_carriers_ancestry['jpt_main_ancestry']=='AFR'].Person_ID.unique()))
print('Nr of carriers with EAS ancestry:', len(netmhc1_df_labels_age_carriers_ancestry[netmhc1_df_labels_age_carriers_ancestry['jpt_main_ancestry']=='EAS'].Person_ID.unique()))

# %%
# ANALYSIS ADDRESSING COMMENT 17: distribution of coverage across variants

# plot distribution of coverage for each variant 

with PdfPages('/Users/barbarawalkowiak/Desktop/msc_thesis/results/reply/qc/QC_distirbution_of_depths.pdf') as pdf:

    fig, axes = plt.subplots(rows, cols, figsize=(25, 30)) 
    axes = axes.flatten()

    for a, ax in enumerate(axes[:len(order_by_total)]):

        var = order_by_total[a]
        parts = var.split('_')[0:2]
        part1, part2 = parts[0], parts[1]
        part1 = part1.replace('_', r'\_')
        part2 = part2.replace('_', r'\_')
        if var == "ALL_VARIANTS":
            var_name = part1 + '\n' + part2
        else:
            var_name = rf'$\mathit{{{part1}}}$' + '\n' + part2
        data = netmhc1_df[netmhc1_df['gene_var']==var]

        sns.histplot(data = data, x = 'depth', ax = ax)
        median = data.depth.median()
        ax.axvline(median, color='red', linewidth=2)

        ax.set_xlabel(f'Depth', fontsize = xaxis_font+3)
        ax.set_ylabel('Frequecy', fontsize = yaxis_font+3)
        ax.set_title(f'{var_name}', fontsize = title_font+4)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_xlim(0, 300)
        #ax.set_ylim(0, 180)

    for j in range(40, len(axes)):
        fig.delaxes(axes[j])

    plt.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)

# %%
