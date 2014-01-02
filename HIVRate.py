# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import sys
import os, os.path
import pandas
sys.path.append('/home/will/PySeqUtils/')
from GeneralSeqTools import fasta_reader, seq_align_to_ref
os.chdir('/home/will/HIVRate/')

# <codecell>

ltr_seqs = []
with open('LTRsearch.fasta') as handle:
    for name, seq in fasta_reader(handle):
        gi = name.split('|')[1]
        ltr_seqs.append((gi, seq.replace('-', '')))

# <codecell>

from Bio import SeqIO
from dateutil import parser


gi2date = {}
with open('ltrsequence.gbx.xml') as handle:
    for rec in SeqIO.parse(handle, 'genbank'):
        tdate = parser.parse(rec.annotations['date'])
        gi = rec.annotations['gi']
        gi2date[gi] = tdate

# <codecell>

hxb2 = """TGGAAGGGCTAATTCACTCCCAACGAAGACAAGATATCCTTGATCTGTGGATCTACCACACACAAGGCTACTTCCCT
GATTAGCAGAACTACACACCAGGGCCAGGGATCAGATATCCACTGACCTTTGGATGGTGCTACAAGCTAGTACCAGT
TGAGCCAGAGAAGTTAGAAGAAGCCAACAAAGGAGAGAACACCAGCTTGTTACACCCTGTGAGCCTGCATGGAATGGATG
ACCCGGAGAGAGAAGTGTTAGAGTGGAGGTTTGACAGCCGCCTAGCATTT--CATCACATGGCCCGAG------------
------------AGCTGCATCCGGAGTACTTC---------AAGAACTGCT-----------------------------
-------GACATCGA------------------------GCTTG---CT----------------------ACAA---GG
GACTTTCCGCTGGGGACTTTCCAG-------------GGAGGCGTGGCCTGGGCGGGACT---GGGGAGTGGCGA---GC
CCTCAGATCCTGCATATAAGCAGC---TGCTTTTTGCCTGTACTGGGTCTCTCTGGTTAGACCAGATCTGAGCCTGGGAG
CTCTCTGGCT---AACTAGGGAACCCACTGCTTAAGCCTCAATAAAGCTTGCCT---TGAGTGCTTC---AAGTAGTGTG
TGC---CCGTCTG---TTGTGTGACTCTGGT---AACTAGAGATCCC---TCAGAC---CCT---TTTAGTCAGTGTGG-
--AAAATCTCT""".replace('\n', '').replace('-', '')

# <codecell>

aln_seqs = []
for gi, aln_seq in seq_align_to_ref(ltr_seqs, hxb2, max_workers = 50):
    try:
        aln_seqs.append((gi2date[gi], aln_seq))
    except KeyError:
        pass

# <codecell>

sorted_alns = sorted(aln_seqs, key = lambda x: x[0])

# <codecell>

from itertools import groupby
from collections import defaultdict

found_seqs = set()
nums = []


for date, seq in sorted_alns:
    isNew = False
    if seq not in found_seqs:
        new_seq_num += 1
        found_seqs.add(seq)
        isNew = True
    nums.append((date, 1, isNew))
    
ltr_df = pandas.DataFrame(nums, columns = ['Date', 'isSeq', 'isNew'])

# <codecell>

nltr = ltr_df.groupby('Date').sum()
nltr.cumsum().plot(legend = False)

# <codecell>

hxb2

# <codecell>

sorted_alns[-1]

# <codecell>


