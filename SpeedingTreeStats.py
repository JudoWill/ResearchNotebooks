# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import sys
sys.path.append('/home/will/PySeqUtils/')

# <codecell>

import TreeingTools
import GeneralSeqTools
import dendropy

# <codecell>

with open('/home/will/SubCData/mafft_ep.fasta') as handle:
    seqs = list(GeneralSeqTools.fasta_reader(handle))

# <codecell>

import os, os.path
import csv
from itertools import product
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from operator import methodcaller
from itertools import groupby
from Bio.Seq import Seq
from Bio import Motif
from Bio.Alphabet import IUPAC
from StringIO import StringIO

from subprocess import check_output, check_call
from tempfile import NamedTemporaryFile as NTF
import shlex

# <codecell>

tmp = Motif.Thresholds.ScoreDistribution(mot, precision = 50)
tmp.threshold_fnr?

# <codecell>

def yield_motifs():
    with open('/home/will/LTRtfAnalysis/Jaspar_PWMs.txt') as handle:
        for key, lines in groupby(handle, methodcaller('startswith', '>')):
            if key:
                name = lines.next().strip().split()[-1].lower()
            else:
                tmp = ''.join(lines)
                mot = Motif.read(StringIO(tmp), 'jaspar-pfm')
                yield name, mot
                yield name+'-R', mot.reverse_complement()
    tmp = u"""A 0  0 6 1 0 0 0 4 2 2 0 0 3 
C  1 1 1 0 5 6 4 1 0 0 0 3 5 5 4 0  
G  0 6 0 1 1 0 0 0 0 7 1 1 0 0 1 0 
T  6 0 0 0 1 1 3 5 7 0 0 0 0 2 2 4"""
    mot = Motif.read(StringIO(tmp), 'jaspar-pfm')
    yield 'coup2', mot
    yield 'coup2-R', mot.reverse_complement()

            
pwm_dict = {}
for num, (name, mot) in enumerate(yield_motifs()):
    if num % 100 == 0:
        print num
    
    #low_thresh = Motif.Thresholds.ScoreDistribution(mot, precision = 50).threshold_fpr(0.01)
    #high_thresh = Motif.Thresholds.ScoreDistribution(mot, precision = 50).threshold_fnr(0.3)
    pwm_dict[name] = mot
    


# <codecell>

wanted_mots = ['ap1', 'ap1-R',
               'cebpa', 'cebpa-R',
               'creb1', 'creb1-R',
               'coup2', 'coup2-R',
               'ets1','ets1-R',
               #'fev', 'fev-R',
               'foxc1',	'foxc1-R',
               #'gata2',	'gata2-R',
               #'gata3',	'gata3-R',
               #'hnf4a',	'hnf4a-R',
               #'hoxa5',	'hoxa5-R',
               'nf-kappab','nf-kappab-R',
               'nfatc2', 'nfatc2-R',
               'nr2f1','nr2f1-R',
               #'tfap2a', 'tfap2a-R',    
               #'znf354c','znf354c-R',
               'sp1', 'sp1-R']

# <codecell>

with open('/home/will/SubCData/C_ltr.fasta') as handle:
    raw_seqs = list(GeneralSeqTools.fasta_reader(handle))

# <codecell>

from operator import itemgetter


def run_mafft(inseqs):
    
    orig_order = [name for name, _ in inseqs]
    with NTF(suffix = '.fasta') as handle:
        GeneralSeqTools.fasta_writer(handle, inseqs)
        handle.flush()
        os.fsync(handle)
        
        cmd = 'mafft --quiet --op 10 --ep 0.123 %s' % handle.name
        out = check_output(shlex.split(cmd))
        
    out_dict = dict(GeneralSeqTools.fasta_reader(StringIO(out)))
        
    return [(name, out_dict[name]) for name in orig_order]


def align_blocks(inseqs, start = 0, winsize = 50):
    
    if start != 0:
        slicer = itemgetter(slice(0, start))
        yield [(name, slicer(seq)) for name, seq in inseqs]
    
    for num in range(start, len(inseqs[0][1]), winsize):
        
        slicer = itemgetter(slice(num, num+winsize))
        yield [(name, slicer(seq)) for name, seq in inseqs]
        
    slicer = itemgetter(slice(num, len(inseqs[0][1])))
    yield [(name, slicer(seq)) for name, seq in inseqs]
    
    
def join_blocks(blocks):
    
    final_seqs = ['']*len(blocks[0])
    
    for tup in blocks:
        seqs = [s for _, s in tup]
        for lets in zip(*seqs):
            if any(l != '-' for l in lets):
                
                for pos, l in enumerate(lets):
                    final_seqs[pos] += l
            else:
                print 'dropped column!'
    names = [n for n, _ in tup]
    return [(n, s) for n, s in zip(names, final_seqs)]
                

def refine_alignment(inseqs):
    
    winsizes = [100, 75, 50]
    starts = [0, 25, 50, 75]
    for win, start in product(winsizes, starts):
        blocks = align_blocks(inseqs, winsize=win, start = start)
        ablocks = []
        print win, start
        for num, block in enumerate(blocks):
            ablocks.append(run_mafft(block))
        inseqs = join_blocks(ablocks)
    return inseqs

# <codecell>

aligned_seqs = run_mafft(raw_seqs)

# <codecell>

refined = refine_alignment(aligned_seqs)

# <codecell>

with open('/home/will/SubCData/refined.fasta', 'w') as handle:
    GeneralSeqTools.fasta_writer(handle, refined)
    
    
    
    

# <codecell>

refined = join_blocks(aligned_blocks)

# <codecell>

aligned_seqs[0]

# <codecell>

refined[0]

# <codecell>

from itertools import chain
def score_tf_align(mots, seq):
    
    nseq = Seq(seq.replace('-', ''), 
               alphabet=IUPAC.unambiguous_dna)
    all_scores = -np.inf*np.ones((len(seq), 1)).flatten()
    for mot in mots:    
        scores = mot.scanPWM(nseq)
        score_iter = chain(iter(scores.flatten()), [np.nan]*len(mot))
        oscore = []
        for num, l in enumerate(seq):
            if l != '-':
                oscore.append(score_iter.next())
            else:
                oscore.append(-np.inf)
        all_scores = np.maximum(all_scores.flatten(), np.array(oscore).flatten())
                
    return all_scores
    

# <codecell>

from itertools import islice
fig, axs = plt.subplots(100, 1, sharex=True, figsize=(10,30))

for ax, (name, seq) in zip(axs.flatten(), seqs):
    out_scores = score_tf_align([pwm_dict[mot] for mot in wanted_mots], seq)
    ax.plot(out_scores)
    ax.set_xticklabels([])
    ax.set_yticklabels([])

# <codecell>

plt.plot(out_scores)

# <codecell>

fig, axs = plt.subplots(len(wanted_mots), 1, figsize = (10, 20), sharex=True)
for mot_name , ax in zip(wanted_mots, axs.flatten()):
    score_mat = []
    for name, seq in seqs:
        score_mat.append(score_tf_align(pwm_dict[mot_name], seq))
    score_np = np.array(score_mat)
    ax.plot(np.mean(score_np>0, axis = 0))
    ax.set_title(mot_name)
fig.tight_layout()

# <codecell>

plt.plot()

# <codecell>


