# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

from __future__ import division
from pandas import *
import os, os.path
import csv
from itertools import product
import numpy as np
import matplotlib.pyplot as plt

os.chdir('/home/will/Tip60Data/')

# <codecell>

microdata = read_csv('microdata.tsv', sep='\t', index_col = 0)
microdata.columns = MultiIndex.from_tuples([col.split('_') for col in microdata.columns], names = ['Experiment', 'Run'])
agg_microdata = microdata.groupby(level = 'Experiment', axis = 1).mean()

# <codecell>

check_cols = ['dTIP60E431Q', 'dTIP60WT']
fold_changes = agg_microdata.copy()
for col in check_cols:
    fold_changes[col] = agg_microdata[col]/agg_microdata['Control']
    
fold_changes = fold_changes.drop(['Control'], axis = 1).applymap(np.log2)

# <codecell>

microdata

# <codecell>

from scipy.stats import norm, ttest_1samp
fold_changes.plot(kind = 'kde')
plt.title('KS Density of fold-changes')

# <codecell>

pvals = fold_changes.copy()
for col in check_cols:
    mval = fold_changes[col].mean()
    stdval = fold_changes[col].std()
    pvals[col] = 1 - norm.cdf(abs(fold_changes[col]), loc=mval, scale = stdval)

# <codecell>

cut = 0.05/len(pvals.index)
good_mask = (pvals<cut).any(axis = 1)
print pvals[good_mask].to_string()

# <codecell>

from operator import methodcaller
from itertools import groupby
from Bio.Seq import Seq
from Bio import Motif
from StringIO import StringIO

def yield_motifs():
    motifdir = '/home/will/Tip60Data/TFdata/'
    with open(motifdir + 'matrix_only.txt') as handle:
        for key, lines in groupby(handle, methodcaller('startswith', '>')):
            if key:
                name = lines.next().strip().split()[-1].lower()
            else:
                tmp = ''.join(lines)
                mot = Motif.read(StringIO(tmp), 'jaspar-pfm')
                yield name, mot
                yield name+'-R', mot.reverse_complement()

            
pwm_dict = {}
for num, (name, mot) in enumerate(yield_motifs()):
    if num % 100 == 0:
        print num
    thresh = Motif.Thresholds.ScoreDistribution(mot, precision = 50).threshold_fpr(0.0001)
    pwm_dict[name] = (mot, thresh)

# <codecell>

from itertools import imap
from operator import itemgetter
def unique_justseen(iterable, key=None):
    "List unique elements, preserving order. Remember only the element just seen."
    # unique_justseen('AAAABBBCCDAABBB') --> A B C D A B
    # unique_justseen('ABBCcAD', str.lower) --> A B C A D
    return imap(next, imap(itemgetter(1), groupby(iterable, key)))

def take(n, iterable):
    "Return first n items of the iterable as a list"
    return list(islice(iterable, n))

# <codecell>

def scan_seqs(tup):
    pwm_tup, row = tup
    seq = Seq(row['Seq'])
    start = int(row['Start'])
    ch = row['Chromosome']
    strand = row['Strand']
    name, (mot, thresh) = pwm_tup
    results = []
    for loc, m in mot.search_pwm(seq, threshold=thresh):
        results.append((name, ch, start+loc, start+loc+len(mot), strand))
    
    return results
    

# <codecell>

from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import partial
from itertools import islice

interval_fields = ['Chromosome','Start','Stop','Refseq','Junk1','Strand','Junk2','Junk3','Junk4','Junk5','Junk6','Junk7','Seq']

with open('seqdata/PromoterSeqdata.interval') as handle:
    with open('seqdata/TFBindingPos.interval', 'w') as ohandle:
        reader = csv.DictReader(handle, delimiter = '\t', fieldnames=interval_fields)
        nreader = unique_justseen(reader, key = lambda x: x['Seq'])
        
        writer = csv.writer(ohandle, delimiter = '\t')
        writer.writerow(['TFname', 'Chr', 'Start', 'Stop', 'strand'])
        
        check_items = product(sorted(pwm_dict.items()), nreader)
        blocksize = 50000
        with ProcessPoolExecutor(max_workers = 30) as executor:
            block = take(blocksize, check_items)
            start = 0
            while block:
                res = executor.map(scan_seqs, block)
                for r in res:
                    if r:
                        writer.writerows(r)
                block = take(blocksize, check_items)    
                start += blocksize
                num = int(start/len(pwm_dict))
                print num
            

# <codecell>

def make_filename(inp):
    return inp.replace(' ', '-').replace(':', '-')


with open('seqdata/TFBindingPos.interval') as handle:
    reader = csv.reader(handle, delimiter = '\t')
    junk = reader.next()
    grouper = lambda x: x[0]
    for key, rows in groupby(reader, key = grouper):
        fname = make_filename(key)
        with open('seqdata/TFhits/'+fname + '.bed', 'a') as ohandle:
            writer = csv.writer(ohandle, delimiter = '\t')
            for row in rows:
                if (int(row[2]) > 0) and (int(row[3]) > 0):
                    writer.writerow(row[1:])

# <codecell>

wanted_genes = set([n.split('.')[0] for n in pvals[good_mask].index])

with open('seqdata/PromoterInervals') as handle:
    with open('seqdata/SigGenePromoters.bed', 'w') as ohandle:
        reader = csv.reader(handle, delimiter = '\t')
        writer = csv.writer(ohandle, delimiter = '\t')
        for row in reader:
            if row[3] in wanted_genes:
                print row[3]
                wanted = [row[0], row[1], row[2], row[5]]
                writer.writerow(wanted)

# <codecell>

wanted_genes

# <codecell>


