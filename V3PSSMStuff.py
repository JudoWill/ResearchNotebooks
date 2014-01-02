# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import os, os.path
from pandas import DataFrame, Series, MultiIndex
from collections import defaultdict
from glob import glob
from itertools import groupby, product

# <codecell>

def fasta_reader(handle):
    
    for key, lines in groupby(handle, lambda x: x.startswith('>')):
        if key:
            name = next(lines)[1:].strip()
        else:
            seq = ''.join(line.strip() for line in lines)
            yield name, seq
    

# <codecell>

infiles = glob('/home/will/HIVReportGen/Data/TrainingSequences/*')
seq_data = defaultdict(dict)

for f in infiles:
    fname = f.split(os.sep)[-1].split('.')[0]
    prot, name = fname.split('-',1)
    if 'multi' in name:
        with open(f) as handle:
            for name, seq in fasta_reader(handle):
                seq_data[prot][name] = ''.join(s for s in seq if s.isalpha())
    else:
        with open(f) as handle:
            _, seq = next(fasta_reader(handle))
            seq_data[prot][name] = seq
            
    

# <codecell>

rmkeys = set()
for key, seq in seq_data['V3'].items():
    if len(seq) != 35:
        rmkeys.add(key)
print(len(rmkeys))

# <codecell>

for key in rmkeys:
    seq_data['V3'].pop(key);

# <codecell>

for key in seq_data.keys():
    with open('/home/will/HIVReportGen/Data/TrainingSequences/' + key + '.fasta', 'w') as handle:
        for sname, seq in seq_data[key].items():
            handle.write('>%s\n%s\n' % (sname, seq))

# <codecell>

import re
reg_obj = re.compile('(A\d+)_(R\d+)-(Primer\d+)')
nseq_file = '/home/will/HIVReportGen/Data/PatientFasta/pbmc_new_all[1].fasta'
path_func = lambda x,y,z: '/home/will/HIVReportGen/Data/PatientFasta/%s-%s-%s.fasta' % (x,y,z)
with open(nseq_file) as handle:
    for name, seq in fasta_reader(handle):
        pat, visit, primer = reg_obj.findall(name)[0]
        nfname = path_func(pat, visit, primer)
        with open(nfname, 'w') as handle:
            handle.write('>%s-%s-%s\n%s' % (pat, visit, primer, seq))
        

# <codecell>

import pickle
import StringIO

with open('/home/will/HIVReportGen/Data/TrainingSequences/training_seqs.pkl', 'r') as handle:
    training_seqs = pickle.load(handle)
with open('/home/will/HIVReportGen/Data/TrainingSequences/training_pssm.pkl') as handle:
    training_scores = pickle.load(handle)

# <codecell>

training_scores.to_csv('/home/will/Downloads/PSSMscores.csv')

# <codecell>

from pandas import *

fields = ['name','score','pred','x4.pct','r5.pct','geno','pos.chg','net.chg','percentile']
PSSMscores = read_csv('/home/will/Dropbox/HIVseqs/output_data.tsv', names = fields, sep = '\t')

# <codecell>

r5_cut = -6.95
x4_cut = -2.88

def BenMeth(val):
    
    if val < r5_cut:
        return 'R5'
    elif val > x4_cut:
        return 'X4'
    else:
        return 'None'

# <codecell>

PSSMscores['n-name'] = 'name'
scores = PSSMscores.groupby('name', as_index = False).agg({'score':np.mean, 'n-name':'count'})

# <codecell>

scores['Tropism'] = scores['score'].map(BenMeth)

# <codecell>

bscores = scores[scores['n-name'] == 1].drop(['n-name'], axis = 1)
bscores.to_csv('/home/will/Downloads/BiggerPSSMscores.csv')

# <codecell>

bscores['Tropism'].value_counts()

# <codecell>


