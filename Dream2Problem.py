# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

from __future__ import division
import os, os.path
from pandas import read_csv, DataFrame, Series, Index, MultiIndex, concat
import csv
from collections import defaultdict
from itertools import islice
from StringIO import StringIO

# <codecell>

import sys
sys.path.append('/home/will/tdlData/will/IpythonNotebook/')
import DreamMicroUtils

os.chdir('/home/will/Dropbox/DREAMLargeData/WillStuff/')
datafile = '/home/will/Dropbox/DREAMProject/DREAM7_DrugSensitivity2/ACalifano_DLBCL_Ly3_14Comp_treatment.txt'

# <codecell>

def get_headers(handle):
    reader = csv.reader(handle, delimiter = '\t')
    ids = reader.next() #junk-line
    drugnames = reader.next()[2:]
    timepoints = reader.next()[2:]
    concs = reader.next()[2:]
    headers = [('probeid',  '', '', ''),
               ('genename', '', '', '')]
    found = set()
    for tup in zip(drugnames, [int(x) for x in timepoints], concs):
        num = 0
        while tup + (num,) in found:
            num += 1
        found.add(tup + (num,))
        headers.append(tup + (num,))
        
    headerindex = MultiIndex.from_tuples(headers, 
                                         names = ['drug', 'timepoint', 
                                                  'concentration', 'replicate'])
    return ids, headerindex

def iterate_lines(handle):
    reader = csv.reader(handle, delimiter = '\t')
    for row in reader:
        if ' /// ' not in row[1]:
            yield '\t'.join(row)
    

with open(datafile) as handle:
    ids, headers = get_headers(handle)
    odf = read_csv(StringIO('\n'.join(iterate_lines(handle))), sep = '\t', names = ids)
#odf.columns = headers

# <codecell>

gene_level_data = odf.drop(['AffyID'], axis = 1).groupby('Genename').agg('median')

# <codecell>

gene_level_data.columns = headers[2:]

# <codecell>

gene_level_data = gene_level_data.reorder_levels([0,2,1,3], axis=1)

# <codecell>

drugs = sorted(set(gene_level_data.columns.get_level_values(0))- set(['DMSO', 'Media']))
timepoints = sorted(set(gene_level_data.columns.get_level_values(2)))

# <codecell>

from itertools import product
from pandas import Panel

mi = MultiIndex.from_tuples(list(product(drugs, timepoints)), 
                            names = ['Drug', 'TimePoint'])
fold_change_panel = Panel(items = ['IC20', '1/10 of IC20'], 
                          major_axis = gene_level_data.index, 
                          minor_axis = mi)
pval_panel = Panel(items = ['IC20', '1/10 of IC20'], 
                          major_axis = gene_level_data.index, 
                          minor_axis = mi)

# <codecell>

import numpy as np
from scipy.stats import ttest_ind

untreated = concat([gene_level_data['DMSO'], 
                    gene_level_data['Media']],
                   axis = 1)

for drug, tp, c in product(drugs, timepoints, ['IC20', '1/10 of IC20']):
    
    m1 = gene_level_data[drug][c][tp]
    
    
    fold_change_panel[c][drug][tp] = np.log2(m1.mean(axis=1)/untreated.mean(axis=1))
    _, pvals = ttest_ind(m1, untreated, axis=1)
    #print pvals.shape, l.values.shape
    #raise KeyError
    pval_panel[c][drug][tp] = Series(pvals, index = untreated.index)

# <codecell>

pval_panel['IC20'].to_csv('/home/will/Downloads/Dream2pvals.csv')

# <codecell>


