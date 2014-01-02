# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

from __future__ import division
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

from types import ListType
from itertools import combinations, groupby, islice, imap
from collections import Counter
from subprocess import check_output
from operator import itemgetter
from StringIO import StringIO
import csv
import shlex
import sys
import glob
sys.path.append('/home/will/PySeqUtils/')
sys.path.append('/home/will/DeepPipeline/AnalysisCode/')
from GeneralSeqTools import fasta_reader, fasta_writer
import HapReconTools

# <codecell>

base_path = '/home/will/DeepPipeline/Data/ShoRAHruns/%s/tmp.bam'
files = ['DrexelMed.A0017.lastz', 'DrexelMed.A0017.R02.lastz',
         'DrexelMed.A0107.lastz', 'DrexelMed.A0107.R02.lastz']
hap_dict = {}

for name in files:
    print 'working on %s' % name
    fname = base_path % name
    
    print 'getting reads'
    reads = sorted(HapReconTools.read_from_bam(fname), key = itemgetter(1))
    ir_reads = list(HapReconTools.yield_IR_reads(reads))
    
    print 'generating graph'
    graph = HapReconTools.generate_hap_graph(ir_reads)
    
    print 'generating haplotypes'
    paths = HapReconTools.simple_paths(graph)
    out_haps = HapReconTools.assemble_haps(paths, graph)
    
    print 'estimating frequencies'
    hap_dict[name] = HapReconTools.estimate_freqs(reads, out_haps, 0.001, tol=0.1)
    

# <codecell>

reload(HapReconTools)
nhap_dict = {}

for name in files:
    print 'working on %s' % name
    fname = base_path % name
    
    print 'getting reads'
    reads = sorted(HapReconTools.read_from_bam(fname), key = itemgetter(1))
    ir_reads = list(HapReconTools.yield_IR_reads(reads))
    
    print 'generating graph'
    graph = HapReconTools.generate_hap_graph(ir_reads)
    
    print 'generating haplotypes'
    paths = HapReconTools.simple_paths(graph)
    out_haps = HapReconTools.assemble_haps(paths, graph)
    
    print 'estimating frequencies'
    nhap_dict[name] = HapReconTools.estimate_freqs(reads, out_haps, 0.01, tol=0.1)

# <codecell>

thap_dict = {}

for name in files:
    print 'working on %s' % name
    fname = base_path % name
    
    print 'getting reads'
    reads = sorted(HapReconTools.read_from_bam(fname), key = itemgetter(1))
    ir_reads = list(HapReconTools.yield_IR_reads(reads))
    
    print 'generating graph'
    graph = HapReconTools.generate_hap_graph(ir_reads)
    
    print 'generating haplotypes'
    paths = HapReconTools.simple_paths(graph)
    out_haps = HapReconTools.assemble_haps(paths, graph)
    
    print 'estimating frequencies'
    thap_dict[name] = HapReconTools.estimate_freqs(reads, out_haps, 0, quick_pass_reps=0, tol=0.1)

# <codecell>

from tempfile import NamedTemporaryFile as NTF
from subprocess import check_output, check_call
import shlex
import os
from concurrent.futures import ProcessPoolExecutor
import csv
from StringIO import StringIO
from itertools import islice
from functools import partial

def check_seqs(db_path, seqs):
    cmd = "blastn -db %(db)s -query %(q)s -outfmt '10 qseqid sseqid pident nident length' -num_threads 20 -max_target_seqs 1"
    fields = ['SeqA', 'SeqB', 'pident', 'nident', 'length']
    dpath =  '/home/will/tmpstuf/haptest/tmpseqs/'
    
    with NTF(suffix='.fa', dir=dpath, delete=False) as check_handle:
        
        fasta_writer(check_handle, seqs)
        check_handle.flush()
        os.fsync(check_handle.fileno())
        
        tdict = {
                 'db':db_path,
                 'q':check_handle.name
                 }
        cmd_list = shlex.split(cmd % tdict)
        out = check_output(cmd_list)
        reader = csv.DictReader(StringIO(out), fieldnames=fields)
        return list(reader)
    

def yield_blocks(iterable, block_size):
    
    block = list(islice(iterable, block_size))
    while block:
        yield block
        block = list(islice(iterable, block_size))
    

def blast_all_v_all(seqsA, seqsB, block_size=20):
        
    dpath = '/home/will/tmpstuf/haptest/tmpseqs/'
    with NTF(suffix='.fa', dir=dpath, delete=False) as db_handle:
        fasta_writer(db_handle, seqsA)
        db_handle.flush()
        os.fsync(db_handle.fileno())
        
        cmd = 'makeblastdb -in %s -dbtype nucl' % db_handle.name
        cmd_list = shlex.split(cmd)
        check_call(cmd_list)
        
        align_func = partial(check_seqs, db_handle.name)
        check_iterable = islice(yield_blocks(iter(seqsB), 50), 20)
        blocks = []
        with ProcessPoolExecutor(max_workers=5) as pool:
            res_iter = pool.map(align_func, check_iterable)
            for num, block in enumerate(res_iter):
                blocks += block
        return blocks

# <codecell>

def add_names(seqs, base):
    tbase = base + '%03i'
    return [(tbase % num, seq) for num, ((seq, pos), _) in enumerate(seqs)]

paired_visits = [('DrexelMed.A0107.lastz', 'DrexelMed.A0107.R02.lastz'),
                 ('DrexelMed.A0017.lastz', 'DrexelMed.A0017.R02.lastz')]
scatter_data = []
for vA, vB in paired_visits:

    vA_names = add_names(hap_dict[vA], 'R00-')
    vB_names = add_names(hap_dict[vB], 'R02-')
    scores = blast_all_v_all(vA_names, vB_names)
    for row in scores:
        sA_pos = int(row['SeqB'].split('-')[-1])
        sB_pos = int(row['SeqA'].split('-')[-1])
        freqA = hap_dict[vA][sA_pos][-1]
        cons = float(row['nident'])/float(row['length'])
        scatter_data.append((freqA, cons, vA.split('.')[1]))
    

# <codecell>

import pandas as pd
tdata = pd.DataFrame(scatter_data, columns = ['SourceFreq', 'Cons', 'Pat'])
cdict = {'A0107':'g', 'A0017':'r'}
tdata['Color'] = tdata['Pat'].map(lambda x: cdict[x])

# <codecell>

plt.figure(figsize = (10, 10))
ax = plt.subplot(111)
ax.scatter(tdata['SourceFreq'], tdata['Cons'], c = list(tdata['Color']), alpha=0.5)
ax.set_xscale('log')

# <codecell>

tdata['Pat']

# <codecell>


