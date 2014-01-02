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

def simple_translate(inseq):
    seq = Seq(inseq, alphabet=generic_dna)
    return seq.translate().tostring()


seq_files = glob.glob('LANLdata/*.fasta')
seq_data = []
for f in seq_files:
    sub, prot = f.split(os.sep)[-1].split('-')[:2]
    with open(f) as handle:
        for name, seq in fasta_reader(handle):
            nseq = ''.join(l for l in seq if l.isalpha())
            if prot != 'LTR':
                nseq = simple_translate(nseq)
            seq_data.append((name, sub, prot, nseq))
            
seq_df = DataFrame(seq_data, columns=['Name', 'Sub', 'Prot', 'Seq'])


# <codecell>

test_seqs = seq_df[seq_df['Prot'] == 'Nef']['Seq'].map(lambda x: x[:30])

# <codecell>

from tempfile import mkdtemp
from TreeingTools import tmp_directory
import shutil
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna, generic_protein
import contextlib
from StringIO import StringIO
import dendropy
from Bio.Alphabet import IUPAC


@contextlib.contextmanager
def push_dir(path):
    
    cur_path = os.path.abspath(os.curdir)
    try:
        os.chdir(path)
        yield
    finally:
        os.chdir(cur_path)
    
def clean_seqs(inseq, alphabet):
    if alphabet == generic_protein:
        valid_lets = set(IUPAC.IUPACProtein.letters)
    elif alphabet == generic_dna:
        valid_lets = set(IUPAC.IUPACUnambiguousDNA.letters)
    else:
        raise KeyError, 'Unknown alphabet.'
    out = ''
    for l in inseq:
        out += l if l.upper() in valid_lets else '-'
        
    return out
    
    
def write_phylip_seqs(seqs, outhandle, alphabet = generic_protein):
    
    trans_names = {}
    out_seqs = []
    for num, (name, seq) in enumerate(seqs):
        new_name = 'Seq-%i' % num
        trans_names[name] = new_name
        out_seqs.append(SeqRecord(Seq(clean_seqs(seq, alphabet), 
                                      alphabet = alphabet), 
                                  id = str(new_name)))
    SeqIO.write(out_seqs, outhandle, 'phylip')
    return trans_names

def make_phylip_seq_dist_mat(inseqs, alphabet, tmp_path = None):
    
    commands = "2\ny\n"
    
    with tmp_directory(dir=tmp_path, rm_dir=False) as r_path:
        with push_dir(r_path):
            with open('cmd_file', 'w') as cmd_handle:
                cmd_handle.write(commands)
            cmd_handle = open('cmd_file')
            with open('infile', 'w') as handle:
                trans_names = write_phylip_seqs(inseqs, handle, alphabet = alphabet)
            if alphabet == generic_protein:
                cmd = 'phylip protdist'
            elif alphabet == generic_dna:
                cmd = 'phylip dnadist'
            else:
                raise KeyError, 'Unknown alphabet.'
            cmd = shlex.split(cmd)
            try:
                check_call(cmd, stdin = cmd_handle)
            except:
                print r_path
                raise KeyError
            with open('outfile') as handle:
                dist_data = handle.read()
            
    return trans_names, dist_data
                
def make_phylip_tree(dist_data, tmp_path = None):
    
    commands = "2\n3\ny\n"
    with tmp_directory(dir=tmp_path, rm_dir=True) as r_path:
        with push_dir(r_path):
            with open('cmd_file', 'w') as cmd_handle:
                cmd_handle.write(commands)
            cmd_handle = open('cmd_file')
            with open('infile', 'w') as handle:
                handle.write(dist_data)
            cmd = shlex.split('phylip neighbor')
            check_call(cmd, stdin = cmd_handle)
            with open('outtree') as handle:
                tree = dendropy.Tree.get_from_stream(handle, schema='newick')
    return tree

def phylip_tree(seqs, alphabet = generic_protein, tmp_path = '/home/will/tmpstuf/ntest/'):
    
    trans_names, dist_data = make_phylip_seq_dist_mat(seqs, alphabet, tmp_path=tmp_path)
    out_tree = make_phylip_tree(dist_data)
    for orig_name, new_name in trans_names.items():
        print orig_name, new_name
        node = out_tree.find_node_with_taxon_label(new_name)
        if node:
            node.taxon = orig_name
    return out_tree

# <codecell>

names, dmat = make_phylip_seq_dist_mat(test_seqs[0:200].to_dict().items(), generic_protein)

# <codecell>

names.items()[:5]
rev_dict = dict((val, key) for key, val in names.items())
handle = StringIO(dmat)
nseqs = int(handle.next())
dist_groups = []

for line in handle:
    if line.startswith('Seq-'):
        dist_groups.append('')
    dist_groups[-1] += line

omat = {}
tmpl = 'Seq-%i'
for seq_num, group in enumerate(dist_groups[:5]):
    parts = group.split()[1:]
    for onum, val in enumerate(parts):
        nkey = (rev_dict[tmpl % seq_num], rev_dict[tmpl % onum])
        nval = float(val)
        if nval >= 0:
            omat[nkey] = nval


# <codecell>

omat.items()[:5]

# <codecell>

from concurrent.futures import ProcessPoolExecutor
from itertools import imap

def take(iterable, N):
    return list(islice(iterable, N))
    
tmp_path = '/home/will/tmpstuf/ntest/'

#with ProcessPoolExecutor(max_workers = 20) as e:
    
iterable = iter(list(test_seqs.to_dict().items()))
items = []
block = take(iterable, 100)
while block:
    items.append(block)
    block = take(iterable, 100)
    
    

# <codecell>


