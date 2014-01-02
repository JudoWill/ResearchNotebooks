# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

from pandas import DataFrame, Series, merge, read_csv, MultiIndex, Index, concat
from subprocess import check_call
from tempfile import NamedTemporaryFile as NTF
import os, os.path
import numpy as np
from scipy.stats import ttest_ind
from itertools import groupby,combinations, islice
from operator import itemgetter
from Bio import Phylo
import networkx
import sys
import pickle

from random import shuffle
import csv, shlex, shutil

os.chdir('/home/will/HIVTropism//')
sys.path.append('/home/will/HIVReportGen/AnalysisCode/')
sys.path.append('/home/will/PySeqUtils/')

# <codecell>

from SeqProcessTools import read_pat_seq_data, load_training_seq_data, align_seq_data_frame
from TreeingTools import make_mrbayes_trees, run_bats, get_pairwise_distances, check_distance_pvals
from GeneralSeqTools import fasta_reader, WebPSSM_V3_fasta, yield_chunks
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import glob
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
from itertools import chain

# <codecell>


with open('trop_dict.pkl') as handle:
    trop_dict = pickle.load(handle)

with open('wanted_data.pkl') as handle:
    wanted_data = pickle.load(handle)

 
from itertools import product
def yield_regions(trop_dict):
    
    regions = ['gp41-seq-align',
               'gp120-seq-align',
                'LTR-seq-align',
                'Nef-seq-align']
    win_sizes = [5,35,10]#,15,20,40,45]
    
    for region in regions:
        prot = region.split('-')[0]
        seq_ser = wanted_data[region].dropna()
        seqs = [(name, ''.join(list(seq))) for name, seq in zip(seq_ser.index, seq_ser.values)]
        seq_len = len(seqs[0][1])
        
        for win, start in product(win_sizes, range(seq_len)):
            stop = start+win
            if stop < seq_len:
                nseqs = [(name, seq[start:stop]) for name, seq in seqs]
                yield prot, start, win, nseqs, trop_dict
                

import dendropy
import TreeingTools
from Bio.Alphabet import generic_dna, generic_protein
def calculate_region(arg):
    prot, start, win, nseqs, trop_dict = arg
    
    fname = 'phyliptrees/%s-%i-%i.tree' % (prot, start, win)
    
    if os.path.exists(fname):
        contree = dendropy.Tree.get_from_path(fname, 'nexus')
        treeset = dendropy.TreeList.get_from_path(fname + 'set', 'nexus')
    else:
        
        alphabet = generic_protein if prot != 'LTR' else generic_dna
        contree = TreeingTools.phylip_tree(nseqs, alphabet=alphabet)
        treeset = dendropy.TreeList([contree])
        contree.write_to_path(fname, 'nexus')
        treeset.write_to_path(fname + 'set', 'nexus')
    
    
    try:
        bats_res = TreeingTools.run_bats(treeset, trop_dict, nreps = 1000)
    except:
        bats_res = None
    
    try:
        dmat = TreeingTools.get_pairwise_distances(contree)
        benj_res = TreeingTools.check_distance_pvals(dmat, trop_dict, nreps = 50)
    except:
        benj_res = None
    
    return prot, win, start, bats_res, benj_res
    

# <codecell>

from itertools import groupby, imap
from operator import itemgetter
from types import StringType
from concurrent.futures import ThreadPoolExecutor

bats_fields = ['Prot',
                 'Start',
                 'WinSize',
                'null mean',
                 'significance',
                 'upper 95% CI',
                 'observed mean',
                 'upper 95% CU',
                 'Statistic',
                 'lower 95% CI',
                 ]
benj_fields = ['Prot',
                'Start',
                'WinSize',
                'Group2Mean',
                'Group2Std',
                'Group2Name',
                'Group1Mean',
                'Group1Std',
                'RawPval',
                'AdjPval',
                'Group1Name']
benj_writer = csv.DictWriter(open('phylip_BenjRes.tsv', 'w'), benj_fields, delimiter = '\t')
bats_writer = csv.DictWriter(open('phylip_BatsRes.tsv', 'w'), bats_fields, delimiter = '\t')

benj_writer.writeheader()
bats_writer.writeheader()

multi = False
print 'Starting multiprocessing!'
if multi:
    pool = ThreadPoolExecutor(max_workers = 30)
    results = pool.map(calculate_region, yield_regions(trop_dict))
else:
    results = imap(calculate_region, yield_regions(trop_dict))

for prot, win, start, bats_res, benj_res in results:
    if type(bats_res) == StringType:
        print 'eror making tree: ', prot, win, start, bats_res
        continue
    
    tdict = {
             'Prot':prot,
             'Start':start,
             'WinSize':win
             }
    if benj_res is None:
        print 'Error making benj_res at', prot, start, win
    else:
        benj_res.update(tdict)
        benj_writer.writerow(benj_res)
    
    if bats_res is None:
        print 'Error making BATS_res at', prot, start, win
        
    else:
        for row in bats_res:
            if None in row:
                row.pop(None)
            row.update(tdict)
            bats_writer.writerow(row)
        
if multi:
    pool.shutdown()


