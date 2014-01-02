# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import os
import sys
sys.path.append('/home/will/PySeqUtils/')
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# <codecell>

os.chdir('/home/will/NewCoEvo/')
import GeneralSeqTools
import TFSeqTools

# <codecell>

import glob
from itertools import islice, imap
import os, os.path
import csv

def get_gi_acc(fname):
    gb = fname.split('/')[-1].split('.')[0]
    with open(fname) as handle:
        for line in handle:
            if line.startswith('ACCESSION'):
                acc = line.strip().split()[-1]
                return gb, acc
    raise AssertionError



gi_to_acc_dict = {}

fname = '/home/will/WLAHDB_data/gi_to_acc.csv'
if os.path.exists(fname):
    with open(fname) as handle:
        for row in csv.reader(handle):
            gi_to_acc_dict[row[0]] = row[1]
else:
    gb_files = glob.glob('/home/will/WLAHDB_data/GenbankDL/*.gb')
    with open(fname, 'w') as handle:
        writer = csv.writer(handle)
        for num, (gbm, acc) in enumerate(imap(get_gi_acc, gb_files)):
            if (num == 100) or (num % 50000 == 0):
                print num
            gi_to_acc_dict[gbm] = acc
            writer.writerow((gbm, acc))

# <codecell>

files = [('B', sorted(glob.glob('/home/will/WLAHDB_data/SeqDump/B_*')))]
seqs = []
for sub, sfiles in files:
    for f in sfiles:
        with open(f) as handle:
            base_name = f.rsplit(os.sep,1)[1].rsplit('.',1)[0]
            prot = base_name.split('_')[1]
            for name, seq in GeneralSeqTools.fasta_reader(handle):
                seqs.append({
                             'Seq':seq,
                             'ID':gi_to_acc_dict[name],
                             'Prot':prot,
                             })
            
seqdf = pd.DataFrame(seqs)

# <codecell>

pseqdf = pd.pivot_table(seqdf, 
                        rows = 'ID', 
                        cols = 'Prot', 
                        values = 'Seq', 
                        aggfunc = 'first')

# <codecell>

tmpenv = pseqdf[['gp120', 'gp41']].fillna('-').apply(lambda x: ''.join(x), axis=1)
tmpenv[tmpenv == '--'] = np.nan
pseqdf['env'] = pseqdf['env'].combine_first(tmpenv)

# <codecell>

ld = []
for f in glob.glob('/home/will/SubCData/LANLRes/*.xlsx'):
    tmp = pd.read_excel(f, 'Sheet1')
    print f
    if 'Coreceptor' in tmp:
        print tmp['Coreceptor'].unique()
    else:
        print 'not there!'
    
    ld.append(tmp.copy())
    
lanl_data = pd.concat(ld, axis = 0, ignore_index=True)

# <codecell>

from collections import Counter
uniq_items = Counter()
for key, val in lanl_data['Coreceptor'].dropna().to_dict().items():
    uniq_items += Counter(val.split())
print uniq_items

# <codecell>

from functools import partial

def is_trop(receptor, val):
    if val == val:
        return 1.0 if receptor in val else 0.0
    else:
        return np.nan
        
    
    
lanl_data['IsR5'] = lanl_data['Coreceptor'].map(partial(is_trop, 'CCR5'))
lanl_data['IsX4'] = lanl_data['Coreceptor'].map(partial(is_trop, 'CXCR4'))
lanl_data['HasTrop'] = lanl_data['IsR5'].notnull() | lanl_data['IsX4'].notnull()

# <codecell>

agg_lanl_data = lanl_data.groupby('Accession').first()

# <codecell>

tmp_lanl = agg_lanl_data[agg_lanl_data['HasTrop']]
wlanl, wseqs = tmp_lanl.align(pseqdf, join='inner', axis=0)

# <codecell>

wseqs.apply(lambda x: x.notnull(), axis=0).astype(float).sum()

# <codecell>

ref_func = TFSeqTools.align_to_ref

# <codecell>

import ConSeqs
from functools import partial

def get_ref_align(ref_func, conbseq, tseq):
    if tseq == tseq:
        query, ref = ref_func(tseq, conbseq)
        return ''.join(q for q,r in zip(query, ref) if r != '-')
    return np.nan



align_funcs = [('env', partial(get_ref_align, ref_func, ConSeqs.GetConSeq('env', alphabet='pro'))),
               ('ltr', partial(get_ref_align, ref_func, ConSeqs.GetConSeq('ltr', alphabet='dna'))),
               ('nef', partial(get_ref_align, ref_func, ConSeqs.GetConSeq('nef', alphabet='pro'))),
               ('pol', partial(get_ref_align, ref_func, ConSeqs.GetConSeq('pol', alphabet='pro'))),
               ('vif', partial(get_ref_align, ref_func, ConSeqs.GetConSeq('vif', alphabet='pro'))),
               ('vpr', partial(get_ref_align, ref_func, ConSeqs.GetConSeq('vpr', alphabet='pro'))),
               ('vpu', partial(get_ref_align, ref_func, ConSeqs.GetConSeq('vpu', alphabet='pro'))),
               ('gag', partial(get_ref_align, ref_func, ConSeqs.GetConSeq('gag', alphabet='pro'))),
               ('v3', partial(get_ref_align, ref_func, 'CTRPNNNTRKSIHIGPGRAFYTTGEIIGDIRQAHC')),
               ]

aseqs = {}
for col, afunc in align_funcs:
    
    aseqs[col] = wseqs[col].map(afunc)


aseqs_df = pd.DataFrame(aseqs)

# <codecell>

def decide_trop(ser):
    if (ser['IsX4'] == 1) and (ser['IsR5'] == 0):
        return 1.0
    elif (ser['IsX4'] == 0) and (ser['IsR5'] == 1):
        return 0.0
    return np.nan

naseqs_df, alanl = aseqs_df.align(wlanl, axis=0, join='inner')
wmask = alanl[['IsX4', 'IsR5']].any(axis=1) & ~alanl[['IsX4', 'IsR5']].all(axis=1)
wtrop = alanl[['IsX4', 'IsR5']][wmask].apply(decide_trop, axis=1)

# <codecell>

hwseqs = naseqs_df.apply(lambda x: x.notnull(), axis=0).astype(float)
hwseqs['Trop'] = wtrop
pd.pivot_table(hwseqs, 
               rows = 'Trop',
               aggfunc = 'sum')

# <codecell>

from operator import itemgetter
from sklearn.metrics import adjusted_mutual_info_score, normalized_mutual_info_score

results = []
for col in naseqs_df.columns:
    
    wanted_seqs, wanted_trops = naseqs_df[col].dropna().align(wtrop, join='inner')
    
    for pos in range(len(wanted_seqs.values[0])):
        getter = itemgetter(pos)
        wanted_col = wanted_seqs.map(getter)
        adj_score = adjusted_mutual_info_score(wanted_col.values, wanted_trops.values)
        norm_score = normalized_mutual_info_score(wanted_col.values, wanted_trops.values)
        results.append({
                        'Prot':col,
                        'Pos':pos+1,
                        'adj_score':adj_score,
                        'norm_score':norm_score
                        })
    
    
    

# <codecell>

results_df = pd.DataFrame(results).groupby(['Prot', 'Pos']).first()

# <codecell>

fig, ax = plt.subplots(1,1, figsize = (15, 5))
results_df.ix['env'].plot(ax = ax, alpha=0.5)
ax.set_xlim(0, results_df.ix['env'].index[-1])
ax.set_xlabel('Env Pos')
ax.set_ylabel('Score')

# <codecell>

prots = list(aseqs_df.columns)
fig, axs = plt.subplots(len(prots), 1, figsize = (10,15))
for ax, prot in zip(axs.flatten(), prots):
    
    results_df.ix[prot].plot(ax = ax, alpha=0.8, legend = False, linewidth=2)
    ax.set_xlim(0, results_df.ix[prot].index[-1])
    if ax.is_last_row():
        ax.set_xlabel('ConB Pos')
    ax.set_title(prot)
    ax.set_ylabel('Score')
    ax.set_ylim(0, ax.get_ylim()[1])
fig.tight_layout()

# <codecell>

genome_dict = {}
with open('HIV1_ALL_2012_genome_DNA.fasta') as handle:
    for name, seq in GeneralSeqTools.fasta_reader(handle):
        genome_dict[name.rsplit('.', 1)[-1]] = seq.replace('-', '')
genome_ser = pd.Series(genome_dict)

# <codecell>

trops = agg_lanl_data[['IsX4', 'IsR5']].dropna().apply(decide_trop, axis=1).dropna()
trop_omes, wgenomes = trops.align(genome_ser, join='inner')

# <codecell>

with open('/home/will/PySeqUtils/HIVDBFiles/HXB2Sequence.fasta') as handle:
    _, conb_genome = GeneralSeqTools.fasta_reader(handle).next()

# <codecell>

from Bio.Seq import Seq

def rev_trans(align_func, offset, seq):
    tseq = Seq(seq[:-(rframe+1)]).reverse_complement().translate(stop_symbol='X')
    return align_func(tseq.tostring())

rgenomes = {}
for rframe in range(3):
    
    cb_seq = Seq(conb_genome[:-(rframe+1)]).reverse_complement().translate(stop_symbol='X')
    afunc = partial(get_ref_align, ref_func, cb_seq)
    rev_ref_func = partial(rev_trans, afunc, rframe)
    rgenomes[rframe] = wgenomes.map(rev_ref_func)
    
rgenomes = pd.DataFrame(rgenomes)

# <codecell>

null_results = []
for col in rgenomes.columns:
    
    wanted_seqs, wanted_trops = rgenomes[col].dropna().align(trop_omes, join='inner')
    
    for pos in range(len(wanted_seqs.values[0])):
        getter = itemgetter(pos)
        wanted_col = wanted_seqs.map(getter)
        adj_score = adjusted_mutual_info_score(wanted_col.values, wanted_trops.values)
        norm_score = normalized_mutual_info_score(wanted_col.values, wanted_trops.values)
        null_results.append({
                             'Frame':col+1,
                             'Pos':pos+1,
                             'adj_score':adj_score,
                             'norm_score':norm_score
                             })

# <codecell>

null_res_df = pd.DataFrame(null_results)
null_res_df['adj_score'] = null_res_df['adj_score'].clip(0)
null_res_df['norm_score'] = null_res_df['norm_score'].clip(0)

# <codecell>

fig, axs = plt.subplots(2,1, sharex=True, figsize = (10, 10))
plot_cols = ['adj_score', 'norm_score']
edges = np.linspace(0, 0.2, 50)
for ax, col in zip(axs.flatten(), plot_cols):
    null_res_df[col].hist(bins = edges, ax = ax)
    ax.set_ylabel('#sites')
    ax.set_xlabel(col)
    

# <codecell>

all_rgenomes = {}
for rframe in range(3):
    print 'frame', rframe
    cb_seq = Seq(conb_genome[:-(rframe+1)]).reverse_complement().translate(stop_symbol='X')
    afunc = partial(get_ref_align, ref_func, cb_seq)
    rev_ref_func = partial(rev_trans, afunc, rframe)
    all_rgenomes[rframe] = genome_ser.map(rev_ref_func)
    
all_rgenomes = pd.DataFrame(all_rgenomes)

# <codecell>

from itertools import combinations
con_len = len(conb_genome)


pairwise_null = []
for p1, p2 in combinations(range(con_len), 2):
    
    for frame in all_rgenomes.columns:
        col1 = all_rgenomes[frame].map(itemgetter(p1))
        col2 = all_rgenomes[frame].map(itemgetter(p2))
        pairwise_null.append({
                              'Frame':frame+1,
                              'Pos1': p1,
                              'Pos2': p2,
                              'adj_score': adjusted_mutual_info_score(col1.values, col2.values),
                              'norm_score': normalized_mutual_info_score(col1.values, col2.values),
                          
                              })
    
    

# <codecell>




