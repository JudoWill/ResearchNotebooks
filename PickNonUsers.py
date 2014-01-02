# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import pandas as pd
import sys
import os, os.path

sys.path.append('/home/will/PatientPicker/')

# <codecell>

import LoadingTools

# <codecell>

redcap_data = LoadingTools.load_redcap_data().set_index(['Patient ID', 'VisitNum'])
cyto_data = pd.read_csv('/home/will/HIVSystemsBio/NewCytokineAnalysis/CytoRawData.csv', sep='\t')
cyto_data['HasCyto'] = True
has_cyto = cyto_data.groupby(['Patient ID', 'VisitNum'])[['HasCyto']].all()

# <codecell>

cols = ['Psychomotor Speed Score',
 'Memory Recall Score',
 'Constructional Score',
 'TMHDS']
redcap_data['Psychomotor Speed Score'].unique()

# <codecell>

import glob
files = glob.glob('/home/will/HIVReportGen/Data/PatientFasta/*.fasta')
seqs = []
for f in files:
    fname = f.split('/')[-1]
    try:
        pid, vn, prot = fname.split('.')[0].split('-', 2)
    except ValueError:
        print fname
    seqs.append((pid, vn, prot, 1))
    
df = pd.DataFrame(seqs, columns = ['Patient ID', 'VisitNum', 'Prot', 'HasSeq'])
has_seq = pd.pivot_table(df, rows = ['Patient ID', 'VisitNum'], cols = 'Prot', values='HasSeq')

# <codecell>

import sys
sys.path.append('/home/will/PySeqUtils/')
import GeneralSeqTools

with open('/home/will/DrugStuff/pat_data.fasta') as handle:
    seqs = list(GeneralSeqTools.fasta_reader(handle))
    out = GeneralSeqTools.WebPSSM_V3_fasta(seqs)
    

# <codecell>

tmp = []
for row in out:
    parts = row[0].split('-')
    if len(parts) == 2:
        pat, vnum = parts
    else:
        pat, vnum, _ = parts
    tmp.append({'Patient ID':pat,
                'VisitNum':vnum,
                'IsR5':row[2]=='0',
                'IsX4':row[2]=='1',
                })
tropism = pd.DataFrame(tmp).groupby(['Patient ID', 'VisitNum']).first()

# <codecell>

redcap_data = pd.merge(redcap_data, has_cyto,
                       left_index=True, right_index=True,
                       how='outer')
redcap_data = pd.merge(redcap_data, has_seq,
                       left_index=True, right_index=True,
                       how='outer')
redcap_data = pd.merge(redcap_data, tropism,
                       left_index=True, right_index=True,
                       how='outer')
redcap_data = redcap_data.drop(['VisitNum', 'Patient ID'], axis=1)

# <codecell>

import numpy
def safe_sum(col):
    ncol = col.dropna()
    if len(ncol) == 0:
        return np.nan
    return ncol.sum()

def safe_col_apply_mean(col, func, indf):
    return indf[col].dropna().map(func).mean()

# <codecell>

from functools import partial
seq_cols = list(has_seq.columns)
test_drug_cols = [col for col in redcap_data.columns if col.startswith('Test-')]
admit_drug_cols = [col for col in redcap_data.columns if col.startswith('Admit-')]
race_cols = [col for col in redcap_data.columns if col.startswith('Race-')]
gender_cols = ['Male', 'Female']
mean_cols = test_drug_cols+admit_drug_cols+race_cols+seq_cols+['HasCyto', 'IsR5', 'IsX4', 'Hepatitis C status (HCV)']
agg_dict = dict([(col, safe_sum) for col in mean_cols])
agg_dict['VisitNum'] = 'count'

cut_list = [('%LowVL', 'Latest viral load', 50),
              ('%LowCD4', 'Latest CD4 count (cells/uL)', 200)]


# <codecell>

redcap_data['Male'] = redcap_data['Gender'] == 'Male'
redcap_data['Female'] = redcap_data['Gender'] == 'Female'
redcap_data[mean_cols] = redcap_data[mean_cols].applymap(float)
pat_sum = redcap_data.reset_index().groupby('Patient ID').agg(agg_dict)
for ncol, tcol, cut in cut_list:
    pat_sum[ncol] = (redcap_data[tcol]<cut).groupby(level='Patient ID').agg(mean)
pat_sum

# <codecell>

pat_sum[sorted(pat_sum.columns)].to_excel('/home/will/DrugStuff/large_pat_group.xlsx')

# <codecell>

non_users = pat_sum[(pat_sum[test_drug_cols]==0).all(axis=1)]

PN_pats = non_users.drop(test_drug_cols, axis=1)

# <codecell>

#PN_pats.to_excel('/home/will/HIVTropism/PN_pats.xlsx')

# <codecell>

import pickle

with open('/home/will/HIVTropism/trop_dict.pkl') as handle:
    trop_data = pickle.load(handle)
tmp = list(trop_data.items())

# <codecell>


for key, val in tmp:
    if key.startswith('A0'):
        print key, val

# <codecell>

(PN_pats['VisitNum']>=3).sum()

# <codecell>


                                             

# <codecell>


# <codecell>


