#!/usr/bin/env python

import HTSeq as hs
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import pysam
import scipy.stats as stats
import seaborn as sns
import sys


def read_counting(gene, bam):
    '''Takes a gene ID and counts the number of reads in the bam file'''

    reads = 0
    for interval in gene['exon_intervals']:
        reads += bam.count(gene['chr'], interval[0], interval[1])

    return reads


def rpkm(genes_df, bam_file, sample_name):
    '''Takes a datafrem with the relevant gene information, a sample bam file, and a sample_name. Calculates the rpkm values 
    for each gene ID in the input gene_df and adds them to the rpkm_df in which the columns are the sample names.'''

    total_reads = bam_file.count()

    for gene in genes_df.index:
        reads_nr = read_counting(genes_df.loc[gene], bam_file)
        rpkm = (10e9 * reads_nr) / (total_reads * genes_df.loc[gene, 'gene_length'])
        rpkm_df.loc[gene, sample_name] = round(rpkm, 1)


def quant_norm(input_df):
    '''Takes a dataframe and performs quantile normalization'''

    # Giving each value a ranking (handling same ranking by maintaining the same order)from highest to lowest in each column. 
    # Reshaping the dataframe (with the rankings as values) by making the previous columns the inner-most indexes.
    values_ranking = (input_df.rank(method="first").stack())

    # Making a new multiindex daframe in which both the inner and outer indexes are the same as 'value_ranking' and grouping them by
    # their rank. Then the mean of the values with the same ranking is calculated. This creates a Series in which the index 
    # is the rankings and the values are the mean of the initial values of the input dataframe with the same rank.
    mean_by_rank = (input_df.stack().groupby(values_ranking).mean())

    # Creating a new index for 'rank_mean' with values between each of the ranks to account for rank ties when ranking
    # by average in the last step
    new_ranks = ((mean_by_rank.index+0.5).to_list() + mean_by_rank.index.to_list())

    # Re-indexing and sorting 'rank_mean'. The resulting NA values (new index is longer than the Series) are replaced by the mean
    # of the previous and next values. These will become the values for the tied ranks in the last step.
    mapping_values = mean_by_rank.reindex(new_ranks).sort_index().interpolate()

    # Ranking the columns of the input dataframe (handling same ranking by averaging). Reshaping the dataframe to map the mean values
    # to the correct rankings and reshape it back.
    normalized = input_df.rank(method='average').stack().map(mapping_values).unstack()

    return normalized


def assign_category(dataframe):
    '''Takes a dataframe and creates a new column intended to catagorize each row of the dataframe by significance and regulation.'''
    if dataframe['fold_change'] < 0 and dataframe['p-val'] <= 0.05:
        val = 'Down Regulated'
    elif dataframe['fold_change'] > 0 and dataframe['p-val'] <= 0.05:
        val = 'Up Regulated'
    else:
        val = 'Non-significant'
    return val


organism = 'saccharomyces_cerevisiae'

path_to_output_files = f'/home/francisco/pipeline_project/pipeline_data/{organism}/output_files'

gtf_file = hs.GFF_Reader(f'/home/francisco/pipeline_project/pipeline_data/{organism}/input_files/gtf_{organism}.sorted.gtf.gz')


# Uses the information in the gtf_file to create a dictionary in which the keys are the gene IDs and the values are dictionaries
# containing the chromosome which contains the gene, the start and end position of the gene, number of exons in the gene, 
# start and end position of each exon, and the length of the gene (only ).
genes_info = {}
for feature in gtf_file:
    if (feature.type == 'gene') or (feature.type == 'exon'):

        # If the key already exists in the dictionary, add relevant information.
        try:
            genes_info[feature.attr['gene_id']]
            if feature.type == 'exon':
                interval = [eval(i) for i in (str(feature.iv).split('[')[1].split(')')[0].split(','))]

                genes_info[feature.attr['gene_id']]['exon_intervals'].append(interval)

                genes_info[feature.attr['gene_id']]['gene_length'] += interval[1]-interval[0]

                if int(feature.attr['exon_number']) > int(genes_info[feature.attr['gene_id']]['exon_number']):
                    genes_info[feature.attr['gene_id']]['exon_number'] = int(feature.attr['exon_number'])
            else:
                genes_info[feature.attr['gene_id']]['gene_interval'] = [eval(i) for i in (str(feature.iv).split('[')[1].split(')')[0].split(','))]

        # If the key does not exist, create the key and the corresponding value (dictionary).
        except KeyError:
            if feature.type == 'gene':
                genes_info[feature.attr['gene_id']] = {'chr': str(feature.iv).split(':')[0],
                                                       'gene_interval': [eval(i) for i in (str(feature.iv).split('[')[1].split(')')[0].split(','))],
                                                       'exon_number': 1,
                                                       'exon_intervals': [],
                                                       'gene_length': 0}
            else:
                interval = [eval(i) for i in (str(feature.iv).split('[')[1].split(')')[0].split(','))]
                genes_info[feature.attr['gene_id']] = {'chr': str(feature.iv).split(':')[0],
                                                       'gene_interval': None,
                                                       'exon_number': int(feature.attr['exon_number']),
                                                       'exon_intervals': [interval],
                                                       'gene_length': interval[1]-interval[0]}


genes_df = pd.DataFrame(genes_info).transpose()

# Dataframe in which the rpkm values will be stored.
rpkm_df = pd.DataFrame(index=genes_df.index)

samples = ['SRR28797006', 'SRR28797007', 'SRR28797008', 'SRR28797012', 'SRR28797013', 'SRR28797014']

for sample in samples:
    if os.path.exists(f'{path_to_output_files}/{sample}/{sample}_aligned_rd_rg_rc.bam'):
        rpkm(genes_df, pysam.AlignmentFile(f'{path_to_output_files}/{sample}/{sample}_aligned_rd_rg_rc.bam', 'rb'), sample)
    else:
        print(f'The sample {sample} is not in this directory.')


np.seterr(divide='ignore')

# Normalizing the RPKM values
normalized_df = quant_norm(rpkm_df)

# Applying log2 to the read_counts
sample_log_df = np.log2(normalized_df)

# Making negative values equal to zero
sample_log_df[sample_log_df < 0] = 0

# Change names of dataframe columns
sample_log_df.rename(columns={'SRR28797006': 'High1',
                              'SRR28797007': 'High2',
                              'SRR28797008': 'High3',
                              'SRR28797012': 'Low1',
                              'SRR28797013': 'Low2',
                              'SRR28797014': 'Low3'},
                              inplace=True)


# Making a new column with the p-value of a t test between the control and SPRRC samples
sample_log_df['p-val'] = stats.ttest_ind(sample_log_df[['High1', 'High2', 'High3']], sample_log_df[['Low1', 'Low2', 'Low3']], axis=1)[1]

# Making a new column with the log2 values the change in average gene expression between the control and SPRRC
sample_log_df['fold_change'] = np.log2(sample_log_df[['High1', 'High2', 'High3']].mean(axis=1) / sample_log_df[['Low1', 'Low2', 'Low3']].mean(axis=1))

# Creating a new column to better distinguish the significan values and up and down regulation
sample_log_df['category'] = sample_log_df.apply(assign_category, axis=1)

# Replacing all the -inf values with 0
sample_log_df.replace(-np.inf, 0, inplace=True)


# Violin Plot
sns.set_theme(font='Arial')
sns.set_style('white')
plt.figure(figsize=(8, 5))

vio_norm = sns.violinplot(data=sample_log_df[['High1', 'High2', 'High3', 'Low1', 'Low2', 'Low3']],
                     inner='box',
                     palette={sample: '#D81B60' if sample in sample_log_df.columns[:3]  
                              else '#1E88E5' for sample in sample_log_df.columns[:6]})

vio_norm.set_title("Normalized Read Counts per Sample", fontsize=18)
vio_norm.set_ylabel('$log_{10}$' + ' (Read Counts)', fontsize=10)
vio_norm.tick_params(labelsize=8)

plt.legend(handles=[mpatches.Patch(color='#1E88E5', label='Low'), 
                    mpatches.Patch(color='#D81B60', label='High')])
sns.despine()
plt.savefig('vio_norm', dpi=400)
plt.close()


# Vulcano Plot
sns.set_theme(font='Arial')
sns.set_style('darkgrid')
plt.figure(figsize=(8, 5))

vulcano_plot = sns.scatterplot(data=sample_log_df, 
                     x='fold_change', 
                     y=-np.log10(sample_log_df['p-val']),
                     hue='category',
                     s=10,
                     palette={'Up Regulated': 'red', 'Down Regulated': 'gold', 'Non-significant': 'black'},
                     hue_order=['Up Regulated', 'Down Regulated', 'Non-significant'])

vulcano_plot.axhline(y=-np.log10(0.05), 
           linewidth=2, 
           color='orange', ls=':')

vulcano_plot.set_title('Changes in Gene Expression', fontsize=18)
vulcano_plot.set_xlabel('$log_{2}$' + ' (Fold Change)', fontsize=10)
vulcano_plot.set_ylabel('$-log_{10}$' + ' (p-value)', fontsize=10)
vulcano_plot.set(xlim=(-6, 6), ylim=(0, 10))
vulcano_plot.tick_params(labelsize=8)
vulcano_plot.text(s='p-value = 0.05',
        x=-5.5, y=1.5,
        fontsize = 8,
        fontstyle = 'oblique')

plt.legend(loc='upper right',
           title='Gene Expression', 
           title_fontsize=12, 
           fontsize=10)
plt.savefig('vulcano_plot.png', dpi=400)
plt.close()


# Cluster Map
clustermap_data = sample_log_df[((sample_log_df['fold_change'] >= 0.5) | 
                              (sample_log_df['fold_change'] <= -0.5)) & (sample_log_df['p-val'] <= 0.05)]

clustermap = sns.clustermap(clustermap_data.iloc[:, : 6], z_score=0, cmap='RdYlBu_r', yticklabels=True)
clustermap.ax_heatmap.set_yticklabels(clustermap.ax_heatmap.get_ymajorticklabels(), fontsize=3)

plt.savefig('clustermap.png', dpi=400)
plt.close()


# Creating an overview table
overview_table = pd.DataFrame()
overview_table['Low_mean'] = round(sample_log_df[['Low1', 'Low2', 'Low3']].mean(axis=1), 2)
overview_table['Low_SD'] = round(np.std(sample_log_df[['Low1', 'Low2', 'Low3']], axis=1), 2)

overview_table['High_mean'] = round(sample_log_df[['High1', 'High2', 'High3']].mean(axis=1), 2)
overview_table['High_SD'] = round(np.std(sample_log_df[['High1', 'High2', 'High3']], axis=1), 2)

overview_table['fold_change'] = round(sample_log_df['fold_change'], 2)
overview_table['p-val'] = round(sample_log_df['p-val'], 2)

overview_table.to_csv('overview_table.csv')
