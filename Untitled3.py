# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import os
import sys
sys.path.append('/home/will/PySeqUtils/')

# <codecell>

import GeneralSeqTools
import glob

# <codecell>

import pandas as pd
files = sorted(glob.glob('/home/will/HIVTropism/LANLdata/SubB*.fasta'))

seqs = []
for f in files:
    prot_name = f.split('/')[-1].split('.')[0].split('-')[1]
    print prot_name
    with open(f) as handle:
        for name, seq in GeneralSeqTools.fasta_reader(handle):
            seqs.append({
                         'GI':name,
                         'Seq':seq.replace('-', '').upper(),
                         'Prot':prot_name
                         })
            

# <codecell>

seq_df = pd.pivot_table(pd.DataFrame(seqs),
                        rows = 'GI',
                        cols = 'Prot',
                        values = 'Seq',
                        aggfunc = 'first')

# <codecell>

from Bio import Seq
from Bio.Alphabet import generic_dna
res = Seq.Seq('ATG', alphabet=generic_dna).translate()
res.tostring()

# <codecell>

def translate(inseq):
    return Seq.Seq(inseq, alphabet=generic_dna).translate().tostring()
benj_seqs = seq_df[['LTR', 'Tat_1', 'Tat_2', 'Vpr', 'V3']].dropna()['Tat_2'].map(translate)

# <codecell>

with open('/home/will/Downloads/tat2_for_benj.fasta', 'w') as handle:
    GeneralSeqTools.fasta_writer(handle, benj_seqs.to_dict().items())

# <codecell>


