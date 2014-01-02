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
from GeneralSeqTools import fasta_reader, WebPSSM_V3_fasta, yield_chunks
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import glob
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
from itertools import chain, product

# <codecell>

with open('trop_dict.pkl') as handle:
    trop_dict = pickle.load(handle)

with open('wanted_data.pkl') as handle:
    wanted_data = pickle.load(handle)

trans_dict = wanted_data['Name'].to_dict()
ntrop_dict = dict((trans_dict[key], val) for key, val in trop_dict.items())
trop_dict = ntrop_dict

wanted_data = wanted_data.set_index('Name')

# <codecell>

wanted_data['Tropism'][wanted_data['gp120-seq-align'].notnull()].value_counts()

# <codecell>

from GeneralSeqTools import fasta_writer
fourkb_cols = ['gp120-seq-align', 'Nef-seq-align', 'Vpr-seq-align', 
                 'Tat-1-seq-align', 'Tat-2-seq-align', 'LTR-seq-align']
four = wanted_data[fourkb_cols].dropna()
wseqs = set()
with open('/home/will/Dropbox/HIVseqs/BensTropismLabels.csv') as handle:
    for row in csv.DictReader(handle, delimiter=','):
        wseqs.add(row['Patient ID'])

        
for col in four.columns:
    found = set()
    prot = col.rsplit('-', 2)[0]
    fname = 'AlignForBenj/fourKB_%s.fasta' % prot
    with open(fname, 'w') as handle:
        for seq, name in zip(four[col], four.index):
            if name in wseqs and name not in found:
                fasta_writer(handle, [(name+'-'+trop_dict[name], ''.join(seq))])
                found.add(name)
    print prot, len(found)

# <codecell>

foukb_lanl = ['AB078005', 'AB221126', 'AB253432', 'AB286955', 
              'AB287365', 'AB287367', 'AB287368', 'AB287369', 
              'AB480695', 'AB485642', 'AB565479', 'AB565496', 
              'AB565497', 'AB565499', 'AB565500', 'AB565502', 
              'AB604946', 'AB604948', 'AB604950', 'AB604951', 
              'AB641836', 'AF003887', 'AF003888', 'AF004394', 
              'AF042100', 'AF042101', 'AF538302', 'AF538303', 
              'AF538307', 'AJ271445', 'AY173953', 'AY352275', 
              'AY835748', 'AY835754', 'AY835759', 'AY835762', 
              'AY835766', 'AY835769', 'AY835770', 'AY835774', 
              'AY835777', 'AY835779', 'AY970950', 'DQ007901', 
              'DQ007902', 'DQ007903', 'DQ295192', 'DQ295193', 
              'DQ295194', 'DQ295195', 'DQ358809', 'DQ837381', 
              'DQ990880', 'EF057102', 'EF363123', 'EF363124', 
              'EF363126', 'EF363127', 'GU647196', 'GU733713', 
              'JN944928', 'JN944936', 'JN944939', 'JN944940', 
              'JN944942', 'JN944943', 'JN944944', 'JN944946', 
              'JN944948', 'JQ316126', 'JQ316128', 'JQ316131', 
              'JQ316132', 'JQ316134', 'JQ316135', 'JQ341411', 
              'JQ429433', 'M17449', 'U34604']

# <codecell>

from collections import Counter
trops = []
for p in wanted_data['gp120-seq-align'].dropna().index:
    trops.append(trop_dict.get(p, None))
    
Counter(trops)

# <codecell>

wseqs = set(wanted_data['gp120'].dropna().index)
cols = ['gp120-seq-align', 'Nef-seq-align', 'Vpr-seq-align', 
                 'Tat-1-seq-align', 'Tat-2-seq-align', 'LTR-seq-align']
for col in cols:
    found = set()
    prot = col.rsplit('-', 2)[0]
    fname = 'AlignForBenj/has_env_%s.fasta' % prot
    df = wanted_data[col].dropna()
    with open(fname, 'w') as handle:
        for seq, name in zip(df, df.index):
            if name in wseqs and name not in found:
                fasta_writer(handle, [(name+'-'+trop_dict[name], ''.join(seq))])
                found.add(name)

# <codecell>

def yield_regions(trop_dict):
    
    regions = ['LTR-seq-align',
               'gp41-seq-align',
               'gp120-seq-align',
               'Nef-seq-align',
               'Vpr-seq-align',
               'Tat-1-seq-align',
               'Tat-2-seq-align',
                ]
    tail_cols = ['gp120', 'gp41', 'Nef', 'Vpr', 
                 'Tat-1', 'Tat-2', 'LTR']
    fourkb_cols = ['gp120', 'Nef', 'Vpr', 
                 'Tat-1', 'Tat-2', 'LTR']
    groups = [('fourkb', wanted_data[fourkb_cols].dropna().index),
              ('full_env', wanted_data[['gp120', 'gp41']].dropna().index),
              ('full_tail', wanted_data[tail_cols].dropna().index),
              ]
    subs = ['SubB']
    win_sizes = [5, 10, 15, 20, 30, 35]
    
    for region, (gname, ind), sub in product(regions, groups, subs):
        prot = region.split('-')[0]
        gwanted = wanted_data.ix[ind]
        mask = gwanted['Sub'] == sub
        seq_ser = gwanted[mask][region].dropna()
        print prot, gname, sub, len(seq_ser)
        seqs = [(name, ''.join(list(seq))) for name, seq in zip(seq_ser.index, seq_ser.values)]
        seq_len = len(seqs[0][1])
        
        for win, start in product(win_sizes, range(seq_len)):
            stop = start+win
            if stop < seq_len:
                nseqs = [(name, seq[start:stop]) for name, seq in seqs]
                yield gname, sub, prot, start, win, nseqs, trop_dict

# <codecell>

import dendropy

# <codecell>

from Bio.Alphabet import generic_dna, generic_protein
import TreeingTools
def calculate_region(arg):
    gname, sub, prot, start, win, nseqs, trop_dict = arg
    
    treename = 'quicktrees/%s-%s-%s-%i-%i.tree' % (gname, sub, prot, start, win)
    matfname = 'quicktrees/%s-%s-%s-%i-%i.pkl' % (gname, sub, prot, start, win)
    
    if os.path.exists(treename):
        #benj_res = 'Already Processed'
        #return gname, sub, prot, win, start, benj_res
        
        with open(matfname) as handle:
            dmat = pickle.load(handle)
            
        with open(treename) as handle:
            tree = dendropy.Tree.get_from_stream(handle, 'newick')
        
    else:
        
        is_aa = prot != 'LTR'
        alphabet = generic_protein if is_aa else generic_dna
        
        try:
            tree, dmat = TreeingTools.phylip_tree_collapse_unique(nseqs, alphabet=alphabet)
        except ValueError:
            benj_res = 'Too few unique sequences to process'
            return gname, sub, prot, win, start, benj_res
        except:
            benj_res = 'uncaught exception in dist-mat'
            return gname, sub, prot, win, start, benj_res
        print 'writing'
        with open(matfname, 'w') as handle:
            pickle.dump(dmat, handle)
        with open(treename, 'w') as handle:
            tree.write_to_stream(handle, 'newick')
    
    try:
        benj_res = TreeingTools.check_distance_pvals(dmat, trop_dict, nreps = 50)
    except AssertionError:
        benj_res = 'too few groups'
        return  gname, sub, prot, win, start, benj_res
    except:
        benj_res = 'uncaught exception'
        return  gname, sub, prot, win, start, benj_res
    
    
    try:
        out = TreeingTools.evaluate_association_index(tree, trop_dict)
        benj_res['AI'], benj_res['AI-pval'], benj_res['AI-null'] = out
    except:
        benj_res['AI'], benj_res['AI-pval'], benj_res['AI-null'] = ('error', 'error', 'error')
    
    return gname, sub, prot, win, start, benj_res

# <codecell>

def quick_yield_regions(trop_dict):
    
    nseqs = wanted_data.ix[foukb_lanl]['gp120-seq-align']
    aseqs = wanted_data['gp120-seq-align'].dropna()
    
    regions = [('C1', 0, 101),
               ('V1', 101, 127),
               ('V2', 127, 166),
               ('C2', 166, 266),
               ('V3', 266, 301),
               ('C3', 301, 355),
               ('V4', 355, 388),
               ('C4', 388, 430),
               ('V5', 430, 439),
               ('C5', 439, 462)]
    
    for name, start, stop in regions:
        
        seqs = [(n, ''.join(s[start:stop])) for n, s in zip(nseqs.index, nseqs.values)]
        yield '44kb', 'gp120', name, start, stop, seqs, trop_dict
        
        seqs = [(n, ''.join(s[start:stop])) for n, s in zip(aseqs.index, aseqs.values)]
        yield 'All', 'gp120', name, start, stop, seqs, trop_dict
    
    
    

# <codecell>

from itertools import groupby, imap
from types import StringType

benj_fields = ['GroupName',
               'Subtype',
               'Prot',
                'Start',
                'WinSize',
                'Group2Mean',
                'Group2Std',
                'Group2Name',
                'Group1Mean',
                'Group1Std',
                'RawPval',
                'AdjPval',
                'Group1Name',
                'AI',
                'AI-pval',
                'AI-null']
fname = 'gp120_new_BenjRes.tsv'
handle = open(fname, 'w')
benj_writer = csv.DictWriter(handle, benj_fields, delimiter = '\t')
benj_writer.writeheader()
results = imap(calculate_region, quick_yield_regions(trop_dict))

for gname, sub, prot, win, start, benj_res in results:
    
    #print prot, start, win
    tdict = {
             'Prot':prot,
             'Start':start,
             'WinSize':win,
             'GroupName':gname,
             'Subtype':sub,
             }
    if type(benj_res) is StringType:
        if (benj_res == 'Already Processed') or benj_res.startswith('Too few unique sequences'):
            continue
        print benj_res, prot, start, win
    else:
        benj_res.update(tdict)
        benj_writer.writerow(benj_res)
handle.close()

# <codecell>

from itertools import groupby, imap
from operator import itemgetter

from concurrent.futures import ThreadPoolExecutor


benj_fields = ['GroupName',
               'Subtype',
               'Prot',
                'Start',
                'WinSize',
                'Group2Mean',
                'Group2Std',
                'Group2Name',
                'Group1Mean',
                'Group1Std',
                'RawPval',
                'AdjPval',
                'Group1Name',
                'AI',
                'AI-pval',
                'AI-null']
fname = 'more_phylip_BenjRes.tsv'
benj_writer = csv.DictWriter(open(fname, 'w'), benj_fields, delimiter = '\t')
   

benj_writer.writeheader()

multi = True
print 'Starting multiprocessing!'
if multi:
    pool = ProcessPoolExecutor(max_workers = 30)
    results = pool.map(calculate_region, yield_regions(trop_dict))
else:
    results = imap(calculate_region, islice(yield_regions(trop_dict), 0,35))

for gname, sub, prot, win, start, benj_res in results:
    
    #print prot, start, win
    tdict = {
             'Prot':prot,
             'Start':start,
             'WinSize':win,
             'GroupName':gname,
             'Subtype':sub,
             }
    if type(benj_res) is StringType:
        if (benj_res == 'Already Processed') or benj_res.startswith('Too few unique sequences'):
            continue
        print benj_res, prot, start, win
    else:
        benj_res.update(tdict)
        benj_writer.writerow(benj_res)
    
        
if multi:
    pool.shutdown()

# <codecell>


# <codecell>

#with open('allgp120.fasta', 'w') as handle:
tres = []
for key, row in wanted_data[['gp120-seq-align', 'Tropism']].dropna().iterrows():
    oname = key+'-'+row['Tropism']
    tres.append((oname, ''.join(row['gp120-seq-align'])))
    
    

# <codecell>

tree, dmat = TreeingTools.phylip_tree_collapse_unique(tres, alphabet=generic_protein)

# <codecell>

with open('gp120tree.nexus', 'w') as handle:
    tree.write_to_stream(handle, 'nexus')

# <codecell>

import networkx
with open('gp120tree.dot') as handle:
    new_tree = networkx.read_dot(handle)

# <codecell>

pos = networkx.spring_layout(new_tree, dim=100)
#networkx.draw_spring(new_tree, 
#                     with_labels = False,
#                     dim = 10)

# <codecell>

pos.items()[-10:]

# <codecell>


