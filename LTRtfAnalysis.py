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

from operator import methodcaller
from itertools import groupby
from Bio.Seq import Seq
from Bio import Motif
from StringIO import StringIO

from subprocess import check_output, check_call
from tempfile import NamedTemporaryFile as NTF
import shlex

os.chdir('/home/will/LTRtfAnalysis/')

# <codecell>

def yield_motifs():
    with open('Jaspar_PWMs.txt') as handle:
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

import pickle

with open('/home/will/HIVReportGen/Data/TrainingSequences/training_seqs.pkl') as handle:
    training_seqs = pickle.load(handle)
    
with open('/home/will/HIVReportGen/Data/TrainingSequences/training_pssm.pkl') as handle:
    pssm_data = pickle.load(handle)
    
with open('/home/will/HIVReportGen/Data/PatientFasta/seq_data.pkl') as handle:
    pat_seq_data = pickle.load(handle)

# <codecell>

def fasta_reader(handle):

    name = None
    for key, lines in groupby(handle, lambda x:x.startswith('>')):
        if key:
            name = next(lines)[1:].strip()
        else:
            seq = ''.join(line.strip() for line in lines)
            yield name, seq

def fasta_writer(seq_tups, handle, alpha_only = True):

    for name, seq in seq_tups:
        if alpha_only:
            seq = ''.join(s for s in seq if s.isalpha())
        handle.write('>%s\n%s\n' % (name, seq))

def seq_number_gen(inseq):
    c = -1
    for let in inseq:
        if let != '-':
            c+=1
        yield c
            
            


def map_to_conb(conb_seq, other_seq):
    with NTF(mode = 'w') as handle:
        fasta_writer([('ConB', conb_seq), ('Other', other_seq)], handle)
        handle.flush()
        os.fsync(handle)
        with NTF(mode = 'rt') as ohandle:
            cmd = shlex.split('muscle -quiet -nocore -in %s -out %s' % (handle.name, ohandle.name))
            check_call(cmd)
            seq_dict = dict(fasta_reader(ohandle))
    
    mapping = dict(zip(seq_number_gen(seq_dict['Other']),
                       seq_number_gen(seq_dict['ConB'])))
    return mapping

# <codecell>

ltr_data = training_seqs['LTR'].combine_first(pat_seq_data['LTR']).dropna()



# <codecell>

with open('ltr_seqs.fasta', 'w') as handle:
    for (p,v), seq in zip(ltr_data.index, ltr_data.values):
        handle.write('>%s-%s\n%s\n' % (p,v,seq))

        
        
        

# <codecell>

pssm_data = read_csv('/home/will/HIVReportGen/Data/TrainingSequences/pssm_data.csv', index_col = [0,1])
def decide_tropism(inval):
    if inval < -6.95:
        return True
    elif inval > -2.88:
        return False
    return np.nan
tropism_data = pssm_data['score'].map(decide_tropism).dropna()
trop_dict = {}
for (pat, visit), val in zip(tropism_data.index, tropism_data.values):
    trop_dict[pat+'-'+visit] = val
    
    
with open('/home/will/Dropbox/HIVseqs/BensTropismLabels.csv') as handle:
    reader = csv.DictReader(handle)
    for row in reader:
        trop_dict['%s-%s' % (row['Patient ID'], row['Visit'])] = row['Prediction'] == 'TRUE'

# <codecell>

hxb2_ltr = """TGGAAGGGCTAATTTACTCCCAAAAAAGACAAGATATCCTTGATCTGTGGGTC
TACCACACACAAGGCTACTTCCCTGATTGGCAGAACTACACACCAGGGCCAGG
GATCAGATATCCACTGACCTTTGGATGGTGCTTCAAGCTAGTACCAGTTGAGC
CAGAGAAGGTAGAAGAGGCCAATGAAGGAGAGAACAACAGCTTGTTACACCCT
ATGAGCCTGCATGGGATGGAGGACCCGGAGAAAGAAGTGTTAGTGTGGAAGTT
TGACAGCCGCCTAGCATTTCATCACATGGCCCGAGAGCTGCATCCGGAGTACT
ACAAGGACTGCTGACATCGAGCTTTCTACAAGGGACTTTCCGCTGGGGACTTT
CCAGGGAGGCGTGGCCTGGGCGGGACTGGGGAGTGGCGAGCCCTCAGATGCTG
CATATAAGCAGCTGCTTTTTGCCTGTACTGGGTCTCTCTGGTTAGACCAGATC
TGAGCCTGGGAGCTCTCTGGCTAACTAGGGAACCCACTGCTTAAGCCTCAATA
AAGCTTGCCTTGAGTGCTTCAAGTAGTGTGTGCCCGTCTGTTGTGTGACTCTG
GTAACTAGAGATCCCTCAGACCCTTTTAGTCAGTGTGGAAAATCTCTAGCA""".replace('\n', '')

# <codecell>

hxb2_ltr[329:337]

# <codecell>



known_binding_pos = [('AP-1 IV', 'ap1', 104),
     ('AP-1 III','ap1', 119),
     ('AP-1 II','ap1', 154),
     #('GRE', pwm_dict['GRE'][0], 191),
     ('AP-1 I','ap1', 213),
     ('C/EBP II', 'cebpa', 280),
     #('USF-1', pwm_dict['USF-1'][0], 221),
     ('ETS-1', 'ets1', 304),
     #('Lef-1', pwm_dict['Lef-1'][0], 317),
     ('ATF/Creb', 'creb1', 329),
     ('C/EBP I', 'cebpa', 337),
     ('NFkB II', 'nf-kappab', 349),
     ('NFkB I', 'nf-kappab', 362),
     ('Sp III', 'sp1', 376),
     ('Sp II', 'sp1', 387),
     ('Sp I', 'sp1', 398),
     ('AP-1','ap1', 539),     
     ('AP-1','ap1', 571),
     #('Oct-1', pwm_dict['OCT1'][0], 440),
]

wanted_pwms = [('ap1', pwm_dict['ap1'][0]),
                ('cebpa', pwm_dict['cebpa'][0]), 
                ('ets1', pwm_dict['ets1'][0]),
                ('creb1', pwm_dict['creb1'][0]),
                ('nf-kappab', pwm_dict['nf-kappab'][0]),
                ('sp1', pwm_dict['sp1-R'][0]),]

# <codecell>

seq = 'TGACATCG'
pos = hxb2_ltr.find(seq)

print pos, pos+len(seq)

# <codecell>

from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import partial
from itertools import islice
from itertools import imap, chain
from operator import itemgetter
from collections import defaultdict

from scipy.optimize import minimize_scalar

def unique_justseen(iterable, key=None):
    "List unique elements, preserving order. Remember only the element just seen."
    # unique_justseen('AAAABBBCCDAABBB') --> A B C D A B
    # unique_justseen('ABBCcAD', str.lower) --> A B C A D
    return imap(next, imap(itemgetter(1), groupby(iterable, key)))

def take(n, iterable):
    "Return first n items of the iterable as a list"
    return list(islice(iterable, n))


def scan_seqs(seq, pwm_tup):
    
    seq = Seq(seq)
    name, mot, thresh_fpr = pwm_tup
    thresh = Motif.Thresholds.ScoreDistribution(mot, precision = 50).threshold_fpr(thresh_fpr)
    results = []
    for loc, m in mot.search_pwm(seq, threshold=thresh):
        if loc > 0:
            results.append((name, loc, loc+len(mot)))
    
    return results

def check_seq(seq, mapping, pwms, thresh_fpr, executor):
    
    tups = [(name, pwm, thresh_fpr) for name, pwm in pwms]
    anal_fun = partial(scan_seqs, seq)
    res = executor.map(anal_fun, tups)
    for tf, start, stop in chain.from_iterable(res):
        yield tf, mapping[start], mapping[stop-1]
    
def check_all_seqs(seqs, mappings, wanted_pwms, thresh_fpr, executor):
    
    all_res = []
    for n, (seq, mapping) in enumerate(zip(seqs, mappings)):
        all_res.append((n, check_seq(seq, mapping, wanted_pwms, thresh_fpr, executor)))
    for n, res in all_res:
        for tf, start, stop in res:
            yield n, tf, start, stop
    
    
    
def obj_fun(thresh_fpr, seqs, mappings, pwms, allowed_binding_pos, correct_binding_pos, executor):
    
    correct_found = 0
    extra_found = 0
    missing = 0
    
    
    for pat, rows in groupby(check_all_seqs(seqs, mappings, pwms, thresh_fpr, executor), key = lambda x: x[0]):
        found_tfs = defaultdict(int)
        for _, tf, start, stop in rows:
            found_tfs[allowed_binding_pos.get((tf, start), None)] += 1
            #if thresh_fpr < 0.001:
            #    print tf, start, stop
        
        for binding_pos in correct_binding_pos:
            if found_tfs[binding_pos] == 0:
                missing += 1
            elif found_tfs[binding_pos] == 1:
                correct_found += 1
            else:
                correct_found += 1
                extra_found += found_tfs[binding_pos]-1
        extra_found += found_tfs[None]
    #if thresh_fpr < 0.001:
    #    raise KeyError
  
    print thresh_fpr, correct_found, missing, extra_found
    return -correct_found + extra_found
    

final_results = []
ltrseqs = ltr_data.dropna().head(n=50)


        
with ProcessPoolExecutor(max_workers = 20) as executor:
    mapping_fun = partial(map_to_conb, hxb2_ltr)
    seq_mappings = list(executor.map(mapping_fun, ltrseqs.values))
    
    for row in wanted_pwms:
        print row
        allowed_binding_pos = dict()
        for _, tf, pos in known_binding_pos:
            if tf == row[0]:
                for nudge in range(-5,6):
                    allowed_binding_pos[(tf, pos+nudge)] = '%s-%i' % (tf, pos)

        correct_binding_pos = set(allowed_binding_pos.values())
    
        res = minimize_scalar(obj_fun, bounds = [0,0.1], method = 'bounded', 
                args = (ltrseqs.values, seq_mappings, [row], allowed_binding_pos, correct_binding_pos, executor))
    
print res
        
            
        

# <codecell>

correct_binding_pos

# <codecell>

tfdata = DataFrame(final_results, columns = ['Patient ID', 'Visit Number', 'TFName', 'Start', 'Stop'])
tfdata

# <codecell>

tf_counts = tfdata[['Patient ID', 'TFName']].groupby('TFName').count()['Patient ID']
print tf_counts

# <codecell>

tf_grouped = tfdata.groupby(['TFName', 'Patient ID', 'Visit Number', 'Start']).first()
print tf_grouped

# <codecell>

def crazy_iterable():
    tindex = list(ltrseqs.index.copy())
    print tindex[:10]
    tindex.sort(key = lambda x: trop_dict.get('%s-%s'%x, 'Unknown'))
    for key, inds in groupby(tindex, key = lambda x: trop_dict.get('%s-%s'%x, 'Unknown')):
        if key == True:
            key = 'R5'
        elif key == False:
            key = 'X4'
            
        for n, (p, v) in enumerate(inds):
            yield (key, n, p, v)
pat_inds = DataFrame(list(crazy_iterable()), columns = ['Tropism', 'Row', 'Patient ID', 'Visit Number'])
map_sizes = pat_inds[['Tropism', 'Row']].groupby('Tropism').max()['Row']
   

# <codecell>

def make_filename(inp):
    return inp.replace(' ', '-').replace(':', '-')

# <codecell>

from pylab import get_cmap


order = ['R5', 'X4', 'Unknown']
cmap = get_cmap('Greys')
for tf, num in zip(tf_counts.index, tf_counts.values):
    if num > 40:
        map_dict = {
            'Unknown':np.zeros((map_sizes['Unknown']+1,700)),
            'R5':np.zeros((map_sizes['R5']+1,700)),
            'X4':np.zeros((map_sizes['X4']+1,700)),
        }
        tmp = tf_grouped.ix[tf].reset_index()
        merged = merge(tmp, pat_inds, left_on = ['Patient ID', 'Visit Number'], right_on = ['Patient ID', 'Visit Number'])
        for _, row in merged.iterrows():
            map_dict[row['Tropism']][row['Row'], row['Start']:row['Stop']] += 1
    
        fig, axes = plt.subplots(3,1, sharex = True, figsize = (10,10))
        plt.title(tf)
        for key, ax in zip(order, axes.flatten()):
            ax.imshow(map_dict[key], cmap = cmap, aspect = 'auto')
            ax.set_ylim(0, map_sizes[key]+1)
            if key != 'Unknown':
                ax.set_ylabel(key + ' TF:' + tf)
            else:
                ax.set_ylabel(key)
        plt.xlabel('LTR Position')
        
        fname = make_filename(tf)
        #plt.savefig('figures/%s.png' % fname)
        

# <codecell>

map_sizes

# <codecell>

from sklearn.cluster import MeanShift


cluster_data = DataFrame(columns = ['Patient ID', 'Visit Number', 'TFName', 'Start', 'Cluster'])

for tf, num in zip(tf_counts.index, tf_counts.values):
    
    data = tf_grouped.ix[tf].reset_index()
    data['TFName'] = tf
    clust = MeanShift(bandwidth = 10)
    res = clust.fit_predict(data[['Start']].values)
    data['Cluster'] = res
    cluster_data = concat([cluster_data, data], axis = 0, ignore_index = True)



# <codecell>

res = crosstab(rows = [cluster_data['Patient ID'], cluster_data['Visit Number']], cols = [cluster_data['TFName'], cluster_data['Cluster']])

# <codecell>

from sklearn.cluster import k_means, mean_shift

centroids, labels = mean_shift(res.values)

labels = Series(labels, index = res.index)
labels.sort()

plt.figure(figsize = (20,20))

plt.imshow(res.ix[labels.index].values)

# <codecell>

labels

# <codecell>


