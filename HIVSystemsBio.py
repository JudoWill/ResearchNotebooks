# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

from pandas import *
import os, os.path
import csv

os.chdir('/home/will/HIVSystemsBio/')

# <codecell>

cocaine_genes = read_csv('CocaineGeneList.csv')
hiv_genes = read_csv('HIVGeneList.csv', sep = '\t')
biomart_conv = read_csv('mart_export.txt', sep = '\t')

# <codecell>

hiv_genes = merge(hiv_genes, biomart_conv,
                    left_on = 'Gene identifier',
                    right_on = 'Ensembl Gene ID',
                    how = 'inner')

# <codecell>

cocaine_genes

# <codecell>

printed = set()
with open('out_gene_list.tsv', 'w') as handle:
    writer = csv.writer(handle, delimiter = '\t')
    for gene, direc in hiv_genes[['EntrezGene ID', 'Expression']].dropna().values:
        geneid = int(gene)
        group = 'HIV-' + direc
        tup = (geneid, group)
        if tup not in printed:
            writer.writerow(tup)
            printed.add(tup)
            
    for gene in cocaine_genes['ID'].values:
        geneid = int(gene)
        writer.writerow((geneid, 'Cocaine'))
    

# <codecell>

both_genes = merge(hiv_genes, cocaine_genes,
                    left_on = 'Gene identifier',
                    right_on = 'ID',
                    how = 'inner')

# <codecell>

both_genes

# <codecell>


