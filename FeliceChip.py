# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=1>

# Felice ChipSeq Results

# <markdowncell>

# Using the ChipSeq Results provided by Adam I was able to use the MACs v1.4 to find the peaks resulting from ChIP of AbcamA, AbcamB, OpbioA, OpbioB. For all comparisons I used the union of the RNAPol-IIA and RNAPol-IIB as the control. I used the option for MACs which controls for the different dataset sizes. All other arguements were the defaults. I'm assuming the Lanes listed as InputB and InputC are controls of some sort.

# <codecell>

import os
import os.path
from subprocess import check_call
import shlex
import tempfile
import glob
from pandas import *

os.chdir('/home/will/Tip60Analysis/Data/DerivedData/')

# <headingcell level=2>

# Mapping Summary

# <markdowncell>

# The MACS algorithm was able to run properly on all samples and did not produce any warnings.

# <headingcell level=3>

# Data Extraction

# <codecell>

macs_result_files = glob.glob('*/NA_peaks.bed')
macs_results = []
col_names = ['Chrom', 'Start', 'End', 'PeakName', 'Score']
for res in macs_result_files:
    anal = res.split('/', 1)[0]
    tdata = read_csv(res, sep='\t', names=col_names)
    tdata['Analysis'] = anal
    macs_results.append(tdata.copy())

macs_res = concat(macs_results, axis = 0, ignore_index=True)

# <headingcell level=3>

# Results

# <headingcell level=3>

# Number of Peaks

# <codecell>

print macs_res['Analysis'].value_counts()

# <markdowncell>

# This pretty consistent with my previous experience of 2000-10000 peaks.

# <headingcell level=3>

# Peak Distribution on Chromosomes

# <codecell>

chrom_dist = crosstab(cols = macs_res['Analysis'], 
                        rows = macs_res['Chrom'])
print chrom_dist.drop('dmel_mitochondrion_genome')

# <headingcell level=3>

# Data Extraction

# <codecell>

prom_files = glob.glob('*/promoters.bed')
prom_results = []
col_names = ['Chrom', 'Start', 'End', 'Genbank', 
                'JunkA', 'Strand', 'JunkB', 'JunkC',
                'JunkD','JunkE','JunkF','JunkG']
drop_cols = [col for col in col_names if col.startswith('Junk')]
for res in prom_files:
    anal = res.split('/', 1)[0]
    tdata = read_csv(res, sep='\t', names=col_names)
    tdata['Analysis'] = anal
    prom_results.append(tdata.drop(drop_cols, axis = 1))

prom_res = concat(prom_results, axis = 0, ignore_index=True)

# <headingcell level=3>

# Number of 'Controlled' Genes Found

# <codecell>

num_genes = prom_res.pivot_table(rows = 'Analysis', 
                                values = 'Genbank', 
                                aggfunc = lambda x: len(x.unique()))
print num_genes

# <markdowncell>

# Again, pretty consistent with my previous results. In this case I used the 10Kb upstream of a gene as the 'promoter region'.

# <headingcell level=3>

# Overlapping Genes

# <codecell>

from itertools import product
overlaps = DataFrame(index = num_genes.index,
                    columns = num_genes.index)
for a, b in product(num_genes.index, repeat = 2):
    maskA = prom_res['Analysis'] == a
    maskB = prom_res['Analysis'] == b
    genesA = prom_res['Genbank'][maskA]
    genesB = prom_res['Genbank'][maskB]
    
    overlaps[a][b] = len(set(genesA) & set(genesB))

print overlaps

# <markdowncell>

# You can see that there is quite a bit of overlap amongst all of the proteins.

# <headingcell level=3>

# Data Extraction

# <codecell>

charts = glob.glob('*/chart_*.txt')
group_data = []
for chart in charts:
    anal = chart.split('/')[0]
    tdata = read_csv(chart, sep = '\t')
    tdata['Analysis'] = anal
    group_data.append(tdata.copy())
    
all_data = concat(group_data, axis = 0, ignore_index=True)
    

# <headingcell level=3>

# Significant Terms and Pathways

# <codecell>

wanted_cols = ['Analysis', 'Category', 'Term', 
                'Count', '%', 'List Total', 
                'Pop Hits', 'Pop Total', 
                'PValue', 'Benjamini']
all_data[wanted_cols].to_excel('../Results/enrichment_results.xlsx', index = False)

# <codecell>

from matplotlib import pyplot as plt
import numpy as np
from pylab import get_cmap

sig_results = all_data['PValue'] < 0.01
sig_data = all_data[sig_results]
row_frac = 8.0/20.0

for cat in all_data['Category'].unique():
    wcat = sig_data['Category'] == cat
    if wcat.sum() == 0:
        continue
    res = crosstab(rows = sig_data['Term'][wcat], 
                    cols = sig_data['Analysis'][wcat], 
                    values = sig_data['PValue'][wcat], 
                    aggfunc = min)
    nrows = len(res.index)
    plt.figure(figsize = (4, int(row_frac*nrows)+1))
    plt.imshow(np.log10(res.values), aspect = 'auto', 
                interpolation = 'nearest', cmap = get_cmap('gray'))
    cbar = plt.colorbar()
    cbar.set_label('-log10(p-value)')
    plt.clim([-12,0])
    plt.xticks(range(len(res.columns)), 
                res.columns, rotation = 90);
    plt.yticks(range(len(res.index)), 
                res.index);
    plt.title(cat)
    
bigres = crosstab(rows = sig_data['Term'], 
                    cols = sig_data['Analysis'], 
                    values = sig_data['PValue'], 
                    aggfunc = min)
bigres.to_excel('../Results/enrichment_table.xlsx')

# <markdowncell>

# These results show that there is a good deal of overlap in the 'functional space' between the different TFs. For example, in the INTERPRO group you can see a strong signal in the 'Histone Binding'. There is also quiet a bit of chromosomal organization and chromatin binding/assembly/etc in the GO BP group. The full results are in the enrichment_table.tsv results.

# <headingcell level=3>

# Annotation Information

# <codecell>

import csv
conv_ids = []
with open('/home/will/Tip60Analysis/Data/Drosophila_melanogaster.gene_info') as handle:
    junk = handle.next()
    for row in csv.reader(handle, delimiter='\t'):
        conv_ids.append((row[2], row[1]))
        
gene_names = DataFrame(conv_ids, columns = ['FlyBase', 'Entrez'])
prom_res['Flybase'] = prom_res['Genbank'].map(lambda x: x.split('-')[0])

ndata = merge(gene_names, prom_res, 
                left_on = 'FlyBase', right_on = 'Flybase')
wanted_cols = ['Analysis', 'Entrez', 'Flybase' ,'Chrom', 'Start', 'End', 'Strand']
ndata[wanted_cols].to_excel('../Results/FoundPromoters.xlsx')

# <codecell>


