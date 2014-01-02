# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

from __future__ import division
import os, os.path
import sys
sys.path.append('/home/will/PySeqUtils/')
from collections import Counter
os.chdir('/home/will/HIVCpG/')
from GeneralSeqTools import fasta_reader

# <codecell>

from itertools import izip, tee, imap, product

def yield_lets(infile):
    with open(infile) as handle:
        for name, seq in fasta_reader(handle):
            for l in imap(lambda x: x.upper(), seq):
                if l != '-':
                    yield l.upper()
                    
def yield_pairs(seq_iter):
    
    prev_iter, next_iter = tee(seq_iter, 2)
    _ = next_iter.next()
    
    for tup in izip(prev_iter, next_iter):
        yield ''.join(tup)
        
base_iter = yield_lets('lanlgenomes/hiv-db.fasta')
base_counts = Counter(yield_pairs(base_iter))
npairs = sum(base_counts.values())
keys = [''.join(tup) for tup in product('ACGT', repeat=2)]
base_freqs = dict([(key, base_counts[key]/npairs) for key in keys])
print base_freqs

# <codecell>

import pandas as pd
store = pd.HDFStore('/home/will/HIVReportGen/Data/SplitRedcap/2013-01-16/EntireCohort.hdf')
seq_data = store['seq_data']
redcap_data = store['redcap']
store.close()

# <codecell>

def yield_bin(row):
    
    for l, val in zip(row['LTR-bin-align'], row['LTR-seq-align']):
        if val not in wlets:
            yield np.nan
        else:
            yield l


tmp_seqs = seq_data[['LTR-seq-align', 'LTR-bin-align']].dropna()
nseqs = []
nvars = []
wlets = set('ACGT')
for key, row in tmp_seqs.iterrows():
    seq = row['LTR-seq-align']
    nlets = sum(l in wlets for l in seq)
    if nlets > 200:
        pat_info = [('Patient ID', key[0]), ('VisitNum', key[1])]
        tups = list(enumerate(yield_pairs(iter(seq))))+pat_info
        vtups = list(enumerate(yield_bin(row)))+pat_info
        nseqs.append(dict(tups))
        nvars.append(dict(vtups))
    
dinuc_seq_df = pd.DataFrame(nseqs).set_index(['Patient ID', 'VisitNum']).T
nuc_var_df = pd.DataFrame(nvars).set_index(['Patient ID', 'VisitNum']).T

# <codecell>

from pylru import lrudecorator

@lrudecorator(500)
def calc_fishers(intup):
    
    ge_row = [base_counts['CG'], npairs - base_counts['CG']]
    
    return fisher_exact([ge_row, list(intup)])[1]

# <codecell>

from scipy.stats import fisher_exact, chisquare
import numpy as np

def rolling_sig(inser):
    
    win_pos = inser.sum()
    win_neg = len(inser) - win_pos
    
    return calc_fishers((win_pos, win_neg))
    

sig_vals = pd.rolling_apply(dinuc_seq_df == 'CG', 100, rolling_sig, center=True)


# <codecell>

import matplotlib.pyplot as plt
(-np.log10(sig_vals)).mean(axis=1).plot(legend=False)
plt.ylabel('p-value')

# <codecell>

has_seq = pd.DataFrame({'has_seq':pd.Series(True, index=dinuc_seq_df.columns)})

t = redcap_data['Event Name'].dropna().apply(lambda x: x.split(' - ')[0])
redcap_data['Patient visit number'] = redcap_data['Patient visit number'].combine_first(t)
redcap_data["Drugs used (choice='Other')"] = redcap_data["Drugs used (choice='Other')"].dropna() == 'Checked'


drug_cols = ['Amphetamines',
             'Barbiturates',
            'Benzodiazepines',
            'Cannabinoid',
            'Cocaine + metabolite',
            'Opiates',
            'Phencyclidine'
            ]

admit_cols = ["Drugs used (choice='Marijuana')",
 "Drugs used (choice='Cocaine (crack, nasal, smoke, inject)')",
 "Drugs used (choice='Heroin (nasal, inject)')",
 "Drugs used (choice='Methamphetamine (smoke, nasal, inject)')",
 "Drugs used (choice='Benzodiazapine (i.e. valium, ativan, xanax, klonipin, etc)')",
 "Drugs used (choice='Narcotics')",
 "Drugs used (choice='Ecstasy')",
 "Drugs used (choice='PCP')",
 "Drugs used (choice='Ritalin')",
 "Drugs used (choice='Other')"]

wcols = ['Patient ID', 'Patient visit number'] + drug_cols + admit_cols
wanted_redcap = redcap_data[wcols].rename(columns = {
                                                     'Patient visit number':'VisitNum',
                                                     })
seq_redcap = pd.merge(wanted_redcap, has_seq,
                      left_on = ['Patient ID', 'VisitNum'],
                      right_index = True, how = 'outer')
seq_redcap

# <codecell>

PC = set()
PN = set()

for key, group in seq_redcap.groupby('Patient ID'):
    
    nseqs = group['has_seq'].sum()
    ever_drug = group[drug_cols].any(axis=0)
    admit_drug = group[admit_cols].any()
    
    always_coc = group['Cocaine + metabolite'].dropna().all()
    if nseqs < 3:
        pass
    elif ~ever_drug.any() and ~admit_drug.any():
        PN.add(key)
    elif (ever_drug.sum() == 1) and always_coc:
        PC.add(key)
        
print PN, PC

# <codecell>

pn_cols = [col for col in dinuc_seq_df.columns if col[0] in PN]
pc_cols = [col for col in dinuc_seq_df.columns if col[0] in PC]
all_vals = pd.rolling_apply(dinuc_seq_df == 'CG', 100, rolling_sig, center=True)
pc_vals = pd.rolling_apply(dinuc_seq_df[pc_cols] == 'CG', 100, rolling_sig, center=True)
pn_vals = pd.rolling_apply(dinuc_seq_df[pn_cols] == 'CG', 100, rolling_sig, center=True)

# <codecell>

fig, axs = plt.subplots(3, 1, sharex=True, figsize = (10,10))

groups = [('Conservation', 'Frequency', 1-nuc_var_df.dropna(thresh=100, axis = 0)),
          ('Non-Users', '-log(p-value)',-np.log10(pn_vals)),
          ('Pure-Cocaine', '-log(p-value)', -np.log10(pc_vals))]
for (name, label, group), ax in zip(groups, axs.flatten()):
    group.mean(axis = 1).plot(ax = ax)
    ax.set_ylabel(label)
    ax.set_title(name)
    
ax.set_xlabel('LTR-Position')
ax.set_xlim([1, 630])
plt.savefig('cpg_islands_by_cocaine.png', dpi=300)

# <codecell>

base_freqs['CG']*100

# <codecell>

len(PN), len(PC)

# <codecell>


