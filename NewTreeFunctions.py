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

v3_seqs = [(name, seq) for name, sub, prot, seq in seq_data if prot == 'V3']

have_data = []
need_data = set([name for name, _ in v3_seqs])
count = 0
print 'Getting WebPSSM scores'
while need_data and count < 5:
    count += 1
    print len(need_data)
    gen = ((name, seq) for name, seq in v3_seqs if name in need_data)
    chunks = yield_chunks(gen, 50)
    with ThreadPoolExecutor(max_workers = 10) as e:
        res = e.map(WebPSSM_V3_fasta, chunks)
        have_data += list(chain.from_iterable(res))
    need_data -= set(row[0] for row in have_data)
            
            

# <codecell>

seq_df['Prot'].unique()

# <codecell>

import numpy as np
def safe_float(inval):
    try:
        return float(inval)
    except ValueError:
        return np.nan

fields = ['name','score','pred','x4.pct','r5.pct','geno','pos.chg','net.chg','percentile']
pssm_data = DataFrame(have_data, columns=fields)
pssm_data['pred'] = pssm_data['pred'] == '1'

float_cols = [1, 3, 4, 6, 7, 8]
for col in float_cols:
    pssm_data[fields[col]] = pssm_data[fields[col]].map(safe_float)
    
valid_pssm = pssm_data[pssm_data['percentile']<0.95]

# <codecell>

trop_scores = valid_pssm.groupby('name')[['score']].mean()
#trop_scores.to_excel('NewPSSMScores.xls')

# <codecell>

grouped_seq_df = seq_df.pivot(index = 'Name', columns='Prot', values='Seq')
grouped_seq_df = merge(grouped_seq_df, trop_scores, 
                        left_index = True, right_index = True)
print grouped_seq_df

# <codecell>

from collections import defaultdict
def decide_tropism(inval):
    if inval < -6.95:
        return 'R5'
    elif inval > -2.88:
        return 'X4'
    return np.nan

grouped_seq_df['Tropism'] = grouped_seq_df['score'].map(decide_tropism)
trops = grouped_seq_df[['Tropism']].dropna()

trop_dict = defaultdict(lambda :'R5')
for ind, trop in zip(trops.index, trops['Tropism'].values):
    trop_dict[ind] = trop

# <codecell>


grouped_seq_df.dropna(subset = ['gp120'])[['score']].to_excel('NewPSSMScores.xlsx')

# <codecell>

wanted_seq_data = grouped_seq_df.dropna(subset = ['gp120'])

# <codecell>

print 'aligning'
align_data = align_seq_data_frame(wanted_seq_data,  '/home/will/HIVReportGen/Data/BlastDB/ConBseqs.txt')

# <codecell>

align_data['Tropism'] = align_data['score'].map(decide_tropism)

# <codecell>

wanted_data = align_data.dropna(subset = ['Tropism'])

# <codecell>

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
                

# <codecell>

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
    
        
    
#prot, start, win, nseqs = regions.next()
#tmp = calculate_region(prot, start, win, nseqs, trop_dict)
#print tmp

# <codecell>


#grouper = itemgetter(0,1)
#grouper(tmp)
#reload(TreeingTools)

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

multi = True
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

# <codecell>

raise KeyError

# <codecell>

from dendropy.treecalc import PatristicDistanceMatrix
from sklearn.cross_validation import StratifiedShuffleSplit
from random import shuffle

def make_trees(all_seqs, col, start, stop):
    
    test_seqs = []
    for (pid, vn), seq in all_seqs.iterrows():
        try:
            tseq = ''.join(list(seq[col][start:stop]))
        except IndexError:
            continue
        key = pid+'-'+vn
        if (key in trop_dict):
            test_seqs.append((key, tseq))
    print len(test_seqs)
    return make_mrbayes_trees(test_seqs)


def get_pairwise_distances(con_tree):
   
    taxons = con_tree.taxon_set
    pdm = PatristicDistanceMatrix(con_tree)
    pdm.calc()
    dmat = {}
    for p1, p2 in combinations(taxons, 2):
        d = pdm(p1, p2)
        dmat[(p1.label, p2.label)] = d
        dmat[(p2.label, p1.label)] = d
        
    return dmat

def check_distance_pvals(mat_data, trop_dict):
    nreps = 500
    frac = 0.5
    g1dist = []
    g2dist = []
    for (key1, key2), dist in mat_data.items():
        if (trop_dict[key1]=='R5') and (trop_dict[key2] == 'R5'):
            g1dist.append(dist)
        elif (trop_dict[key1] == 'X4') and (trop_dict[key2] == 'X4'):
            g2dist.append(dist)
    nitems = int(min(frac*len(g1dist), frac*len(g2dist)))
    print len(g1dist), len(g2dist)
    _, raw_pval = ttest_ind(g1dist, g2dist)
    cor_pvals = []
    for _ in range(nreps):
        shuffle(g1dist)
        shuffle(g2dist)
        _, pval = ttest_ind(g1dist[:nitems], g2dist[:nitems])
        cor_pvals.append(pval)
    return raw_pval, np.mean(cor_pvals), np.mean(g1dist), np.mean(g2dist), np.std(g1dist), np.std(g2dist)

# <codecell>

from concurrent.futures import ProcessPoolExecutor
from itertools import imap
import csv

def do_analysis(tup):
    start, window = tup
    stop = start + window
    print 'starting process'
    con_tree, multi_trees = make_trees(all_seqs[['gp120-seq-align']].dropna(), 'gp120-seq-align', start, stop)
    #print con_tree
    #print multi_trees
    #bats_res = run_bats(multi_trees[:100], trop_dict)
    
    dmat = get_pairwise_distances(con_tree)
    mgd_res = check_distance_pvals(dmat, trop_dict)
    
    return start, window, mgd_res



window = 35
final_res = []
win = 35
inputs = [(start, win) for win in [5,10,15,20] for start in range(0, 462-win)]
if not os.path.exists('/home/will/tmpstuf/results.csv'):
    with open('/home/will/tmpstuf/results.csv', 'a') as ohandle:
        writer = csv.writer(ohandle, delimiter = '\t')
        with ProcessPoolExecutor(max_workers = 20) as ex:
            res_list = ex.map(do_analysis, inputs)
            #res_list = imap(do_analysis, inputs)
            for start, win, mgd_res in res_list:
                print start, win, 'finished'
                writer.writerow((start, win)+mgd_res)
raise KeyError

# <codecell>

headers = ['Start', 'WinSize', 'RawP', 'AdjP', 'R5Mean', 'X4Mean', 'R5std', 'X4std']
final_data = read_csv('/home/will/tmpstuf/results.csv', sep = '\t', names=headers)

final_data['Prot'] = 'gp120'

nheaders = ['Prot', 'Start', 'WinSize', 'RawP', 'AdjP', 'R5Mean', 'X4Mean', 'R5std', 'X4std']

# <codecell>

#one_based
gp120_features = [('V1', 100, 127),
                  ('V2', 127, 166),
                  ('V3', 266, 301),
                  ('V4', 355, 388),
                  ('V5', 430, 439)]
from matplotlib.patches import Rectangle

# <codecell>

final_data.head()
nfinal_data = final_data.groupby(['WinSize', 'Start']).agg('mean')

# <codecell>

plt.figure(figsize = (10,10))
winsizes = sorted(final_data['WinSize'].unique())
plt.hold(True)
for wsize in winsizes:
    df = nfinal_data.ix[wsize].dropna()
    yvals = (-np.log10(df['AdjP'])).values
    xvals = (df.index + wsize/2).values
    plt.plot(xvals, yvals, label = str(wsize))
plt.xlabel('Gp120-pos')
plt.ylabel('-log10(p-val)')

for name, start, stop in gp120_features:
    rect = Rectangle([start, 0], stop-start, 350, facecolor = 'r', alpha = 0.2)
    plt.gca().add_patch(rect)
    plt.text((start+stop)/2, 330, name)
    #plt.vlines([start, stop], 0, 300)

plt.legend(loc='upper left')
plt.ylim([0,350])
plt.xlim([0, 460])
plt.hold(False)
plt.savefig('gp120-multi-win.png')

# <codecell>

plt.figure(figsize = (10,10))

plt.hold(True)
for wsize in winsizes:
    if wsize < 11:
        df = nfinal_data.ix[wsize].dropna()
        yvals = (-np.log10(df['AdjP'])).values
        xvals = (df.index + wsize/2).values
        plt.plot(xvals, yvals, label = str(wsize))
plt.xlabel('Gp120-pos')
plt.ylabel('-log10(p-val)')

plt.legend(loc='upper left')
plt.ylim([0,350])
plt.xlim([266, 301])
plt.hold(False)

# <codecell>

smoothed_pvals = rolling_mean(nfinal_data['AdjP'], 10)
print smoothed_pvals

# <codecell>

from pandas.stats.moments import rolling_mean

plt.figure(figsize = (10,10))

plt.hold(True)
smoothed_pvals = rolling_mean(nfinal_data['AdjP'], 10)

for wsize in winsizes:
    df = smoothed_pvals.ix[wsize].dropna()
    #print wsize, df
    yvals = (-np.log10(df)).values
    xvals = (df.index + wsize/2).values
    plt.plot(xvals, yvals, label = str(wsize))
plt.xlabel('Gp120-pos')
plt.ylabel('-log10(p-val)')

for name, start, stop in gp120_features:
    rect = Rectangle([start, 0], stop-start, 25, facecolor = 'r', alpha = 0.2)
    plt.gca().add_patch(rect)
    #plt.text((start+stop)/2, 330, name)
    #plt.vlines([start, stop], 0, 300)

plt.legend(loc='upper left')
plt.ylim([0,25])
plt.xlim([0, 460])
plt.hold(False)
plt.savefig('gp120-multi-smoothed.png')

# <codecell>

import pickle

# <codecell>

with open('wanted_data.pkl') as handle:
    wanted_data = pickle.load(handle)

# <codecell>

import TreeingTools
seq_data = wanted_data['Nef-seq-align'].dropna().map(lambda x: ''.join(x[:30])).to_dict().items()
with open('test_nef_seq.phylip', 'w') as handle:
    TreeingTools.write_phylip_seqs(seq_data, handle)

# <codecell>


