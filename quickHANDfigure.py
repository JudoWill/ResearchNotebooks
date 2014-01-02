# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

from __future__ import division
from pandas import *
import os, os.path
import sys
import numpy as np

sys.path.append('/home/will/HIVReportGen/AnalysisCode/')
sys.path.append('/home/will/PySeqUtils/')
os.chdir('/home/will/HIVVariation/')

# <codecell>

from GeneralSeqTools import call_muscle

# <codecell>

store = HDFStore('/home/will/HIVReportGen/Data/SplitRedcap/2013-01-16/EntireCohort.hdf')
redcap_data = store['redcap']
seq_data = store['seq_data']

t = redcap_data['Event Name'].dropna().apply(lambda x: x.split(' - ')[0])
t.unique()
redcap_data['Patient visit number'] = redcap_data['Patient visit number'].combine_first(t)


wanted_cols = ['Patient ID', 'Patient visit number', 'Total Modified Hopkins Dementia Score']
wanted_redcap = redcap_data[wanted_cols]
data = merge(wanted_redcap, seq_data[['LTR']],
            left_on = ['Patient ID', 'Patient visit number'],
            right_index = True, how = 'inner')
data = data.rename(columns= {
                                'Patient visit number':'VisitNum',
                                'Date of visit':'Date',
                                'Total Modified Hopkins Dementia Score':'HIVD'\
                            }).dropna()
data.sort(['Patient ID', 'VisitNum'], inplace=True)

# <codecell>

hxb2_ltr = """TGGAAGGGCTAATTTACTCCCAAAAAAGACAAGATATCCTTGATCTGTGGGTC
TACCACACACAAGGCTACTTCCCTGATTGGCAGAACTACACACCAGGGCCAGG
GATCAGATATCCACTGACCTTTGGATGGTGCTTCAAGCTAGTACCAGTTGAGC
CAGAGAAGGTAGAAGAGGCCAATGAAGGAGAGAACAACAGCTTGTTACACCCT
ATGAGCCTGCATGGGATGGAGGACCCGGAGAAAGAAGTGTTAGTGTGGAAGTT
TGACAGCCGCCTAGCATTTCATCACATGGCCCGAGAGCTGCATCCGGAGTACT
ACAAGGACTGCTGACATCGAGCTTTCTACAAGGGACTTTCCGCTGGGGACTTT
CCAGGGAGGCGTGGCCTGGGCGGGACTGGGGAGTGGCGAGCCCTCAGATGCTG
CATATAAGCAGCTGCTTTTTGCCTGTACTGGGTCTCTCTGGTTAGACCAGATC
TGAGCCTGGGAGCTCTCTGGCTAACTAGGGAACCCACTGCTTAAGCCTCAATA
AAGCTTGCCTTGAGTGCTTCAAGTAGTGTGTGCCCGTCTGTTGTGTGACTCTG
GTAACTAGAGATCCCTCAGACCCTTTTAGTCAGTGTGGAAAATCTCTAGCA""".replace('\n', '')


def align_seq_ser(seq_series):
    
    nseqs = [(i, s) for i, s in zip(seq_series.index, seq_series.values)]
    nseqs += [('hxb2', hxb2_ltr)]
    daln = dict(call_muscle(nseqs))
    aln = [daln[str(s)] for s, _ in nseqs]
    aln_ser = Series(aln[:-1], seq_series.index)
    return aln_ser

data['LTR-align'] = data.groupby('Patient ID')['LTR'].apply(align_seq_ser)

# <codecell>

wanted_seq_cols = [108, 115, 120, 160, 168, 181, 251, 293] #1-based!!


def get_wanted_seq_cols(seq_series):
    
    nseqs = [(i, s) for i, s in zip(seq_series.index, seq_series.values)]
    nseqs += [('hxb2', hxb2_ltr)]
    daln = dict(call_muscle(nseqs))
    aln = [daln[str(s)] for s, _ in nseqs]
    outs = [[] for _ in range(len(aln)-1)]
    hxb2pos = 0
    for tup in zip(*aln):
        if tup[-1] != '-':
            hxb2pos += 1 #1-based!
        if hxb2pos in wanted_seq_cols:
            for out, let in zip(outs, tup):
                out.append(let)
    
    out_ser = Series(outs, seq_series.index)
    return out_ser
data['SNPCols'] = data.groupby('Patient ID')['LTR'].apply(get_wanted_seq_cols)

# <codecell>

wanted_seq_cols = [108, 115, 120, 160, 168, 181, 251, 293]
ref_val = ['A', 'A', 'C', 'C', 'G', 'A', 'G', 'G']

names = ['%i-%s' % (c, v) for c,v in zip(wanted_seq_cols, ref_val)]
def check_data(series):
    
    out = []
    for let, wlet in zip(series['SNPCols'], ref_val):
        if let == '-':
            out.append(np.nan)
        else:
            out.append(float(let == wlet))
    #print out
    return Series(out, index=names)

snp_data = data.apply(check_data, axis = 1)
ndata = concat([data, snp_data], axis = 1)
print ndata

# <codecell>

ndata['HIVD-I'] = ndata['HIVD']<10
nnames = [name + '-C' for name in names]
for nname, name in zip(nnames, names):
    ndata[nname] = ndata[name].dropna()==0

counts = ndata.groupby('HIVD-I')[names].sum()
ncounts = (~ndata[names].applymap(np.isnan)).sum(axis = 0) - counts
freqs = ndata.groupby('HIVD-I')[names].mean()
print ncounts
print counts

# <codecell>

from scipy.stats import fisher_exact

pvals = []
for name in names:
    table = [[counts[name][True], ncounts[name][True]], [counts[name][False], ncounts[name][False]]]
    _, pval = fisher_exact(table)
    pvals.append(pval)
    
pval_series = Series(pvals, names)
pval_df = DataFrame({'Pval':pval_series}).T
print pval_df

# <codecell>

concat([pval_df.reset_index(), counts.reset_index(), ncounts.reset_index()], axis =0 , ignore_index = True).to_excel('quick_hand_results.xlsx')

# <codecell>


