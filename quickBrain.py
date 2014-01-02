# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

from pandas import *
from subprocess import check_call
from tempfile import NamedTemporaryFile as NTF
import os, os.path
import numpy as np
from scipy.stats import ttest_ind
from itertools import groupby,combinations, islice
from operator import itemgetter
import sys

from random import shuffle
import csv, shlex, shutil

os.chdir('/home/will/HIVTropism/')
sys.path.append('/home/will/HIVReportGen/AnalysisCode/')
sys.path.append('/home/will/PySeqUtils/')


# <codecell>

from GeneralSeqTools import fasta_reader, 

# <codecell>

tat_ex1_pos = (5830, 6044) #0 based

pos_data = read_csv('simple_results.txt', sep = '\t')

with open('Tat1-AB1_passed-cleaned.fasta') as handle:
    seq_data = list(fasta_reader(handle))

# <codecell>

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_protein

def translate_to_tat(inseq, start_pos, end_pos, rep = 0):
    
    if rep == 2:
        print 'ESCAPED!!!'
        
        return 
    nstart = tat_ex1_pos[0] - start_pos
    nend = tat_ex1_pos[1] - end_pos
    
    if (nstart >= 0):
        nseq = inseq[nstart:nend]
        tseq = Seq(nseq, generic_dna)
        aa = tseq.translate()
        return aa.tostring()
    else:
        #print 'starts late'
        nseq = 'N'*abs(nstart)+inseq
        #print 'fixing'
        #print inseq
        #print nseq
        #print len(inseq), len(nseq), tat_ex1_pos[1]-tat_ex1_pos[0]
        return translate_to_tat(nseq, tat_ex1_pos[0], end_pos, rep = rep+1)
        
    
    

        
    

# <codecell>

prots = []
names = []
for (name, seq), (idx, row) in zip(seq_data, pos_data.iterrows()):
    aa = translate_to_tat(seq, row['GenomeStart']-1, row['GenomeEnd']-1)
    if sum(1 for l in aa if l == 'X') < 10:
        
        prots.append(aa)
        names.append(name)

# <codecell>

seq_df = DataFrame({
                        'Tat':Series(prots, index = names)
                    })

# <codecell>

from SeqProcessTools import align_seq_data_frame

out_align = align_seq_data_frame(seq_df, '/home/will/HIVReportGen/Data/BlastDB/ConBseqs.txt')

# <codecell>

cohort_data = read_csv('Texas_Cohort_Data.txt', sep = '\t')
nout_align = out_align.reset_index()
nout_align['short_ind'] = nout_align['index'].map(lambda x: x.split('-')[0])


nout = merge(nout_align, cohort_data, 
            left_index = 'short_ind', right_index = 'Patient ID')

nout = nout.drop(nout['NeuroCog'] == 'Not Tested', axis = 0)
print nout

# <codecell>

from SeqAnalysisTools import calculate_fisher_exact

impaired = nout['NeuroCog'] != 'Normal'

res = calculate_fisher_exact(nout['Tat-bin-align'][impaired], nout['Tat-bin-align'][~impaired])

# <codecell>

#import numpy as np
plt.figure(figsize = (10,10))
res.plot()
plt.ylabel('p-value')
plt.xlabel('Tat EX-1')
plt.ylim([0,1.1])

# <codecell>

res.to_csv('p_value_res.csv')

# <codecell>

nout.dropna()['NeuroCog'].value_counts()

# <codecell>

nout['NeuroCog']

# <codecell>


