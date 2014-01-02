# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

from __future__ import division
from pandas import *
import os, os.path
import sys
import numpy as np

sys.path.append('/home/will/HIVReportGen/AnalysisCode/')
sys.path.append('/home/will/PySeqUtils/')
os.chdir('/home/will/HIVVariation/')
from GeneralSeqTools import call_muscle

# <codecell>

from GeneralSeqTools import fasta_reader

seq_data = []
with open('PBMC_analyzed.clean.fasta') as handle:
    for name, seq in fasta_reader(handle):
        try:
            pid, vn = name.split('-')[0:2]
        except ValueError:
            print name
            raise ValueError
        seq_data.append((pid, vn, seq))
        
seq_df = DataFrame(seq_data, columns=['Patient ID', 'VisitNum', 'Seq'])

# <codecell>

wanted_seq_cols = [340, 381] #1-based!!
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

sp3_pos = hxb2_ltr.find('GAGGCGTGGC')
sp2_pos = hxb2_ltr.find('TGGGCGGGAC')
sp1_pos = hxb2_ltr.find('GGGGAGTGGC')

sp3_locs = [sp3_pos+1+i for i in range(10)]
sp2_locs = [sp2_pos+1+i for i in range(10)]
sp1_locs = [sp1_pos+1+i for i in range(10)]

wanted_seq_cols += sp3_locs+sp2_locs+sp1_locs


def get_wanted_seq_cols(seq_series):
    
    nseqs = [('test', seq_series.values[0]),
             ('hxb2', hxb2_ltr)]
    daln = dict(call_muscle(nseqs))
    out = dict([(str(x), None) for x in wanted_seq_cols])
    scols = set(wanted_seq_cols)
    hxb2pos = 0
    for hxb2_l, test_l in zip(daln['hxb2'], daln['test']):
        if hxb2_l != '-':
            hxb2pos += 1 #1-based!
            
        if hxb2pos in scols:
            out[str(hxb2pos)] = test_l
    
    out_ser = Series(out)
    return out_ser
out = seq_df[['Seq']].apply(get_wanted_seq_cols, axis = 1)

# <codecell>

col_groups = [('SP1', [str(x) for x in sp1_locs]),
              ('SP2', [str(x) for x in sp2_locs]),
              ('SP3', [str(x) for x in sp3_locs])]
keep_cols = set(['340', '381'])
for tf, group in col_groups:
    
    tmp = out[group].apply(lambda x: ''.join(x.values), axis=1)
    out[tf] = tmp
    wdrop = [col for col in group if col not in keep_cols]
    out = out.drop(wdrop, axis=1)
    

# <codecell>

nseq_df = concat([seq_df, out], axis=1)
nseq_df['340-381'] = nseq_df[['340', '381']].apply(lambda x: ''.join(x.values), axis = 1)

# <codecell>

from itertools import product
counts = []
for col in ['340', '381']:
    hxb2let = hxb2_ltr[int(col)-1]
    ncol = str(col) + '-' + hxb2let
    for let in 'ACGT':
        count = (nseq_df[col] == let).sum()
        counts.append((ncol, let, count))

single_done = set()
for l1, l2 in product('ACGT', 'ACGT'):
    if (l1 == l2):
        if l1 in single_done:
            continue
        single_done.add(l1)
        
    tmp = l1+l2
    count = nseq_df['340-381'] == tmp
    counts.append(('340-381', tmp, count.sum()))

        
df = DataFrame(counts, columns=['Col', 'Let', 'Count'])

# <codecell>

pos_counts = pivot_table(df, rows = ['Col', 'Let'], values='Count', aggfunc = 'sum')
print pos_counts
DataFrame({'Counts':pos_counts}).to_excel('CountsForSonya.xlsx')

# <codecell>

pos_counts.groupby(level='Col').sum()

# <codecell>

wanted_pats = set(nseq_df['Patient ID'][nseq_df['340-381'] == 'TT'])
wanted_seqs = nseq_df[nseq_df['Patient ID'].map(lambda x: x in wanted_pats)]
print wanted_seqs.to_string()

# <codecell>

store_loc = os.path.join('home', 'will', 'HIVReportGen', 'Data',
                        'SplitRedcap', '2013-01-16',
                        'EntireCohort.hdf')
store = HDFStore('/'+store_loc)
redcap_data = store['redcap']

#fix the Visit Number issue
t = redcap_data['Event Name'].dropna().apply(lambda x: x.split(' - ')[0])
vnum_field = 'Patient visit number'
redcap_data[vnum_field] = redcap_data[vnum_field].combine_first(t)

#Fix the Other-drug issue
ofield = "Drugs used (choice='Other')"
redcap_data[ofield] = redcap_data[ofield].dropna() == 'Checked'

# <codecell>

wanted_cols = {'Patient ID':'Patient ID', 
                'Patient visit number':'VisitNum', 
                'Date of visit':'VisitDate',
                'HIV seropositive date':'HIV seropositive date',
                'Total Modified Hopkins Dementia Score':'HIVD'
                }
wanted_redcap = redcap_data[wanted_cols.keys()].rename(columns = wanted_cols)
has_date_mask = wanted_redcap['HIV seropositive date'].notnull()

wanted_redcap['YearsSeropositive'] = np.nan
years_sero_func = lambda x: x['VisitDate'].year - x['HIV seropositive date'].year
years_sero_data = wanted_redcap[['VisitDate', 'HIV seropositive date']].reset_index().dropna().apply(years_sero_func, axis = 1)
years_sero_data.index = has_date_mask[has_date_mask].index
wanted_redcap['YearsSeropositive'][has_date_mask] = years_sero_data

wanted_redcap

# <codecell>

wanted_red_pats = wanted_redcap[wanted_redcap['Patient ID'].map(lambda x: x in wanted_pats)]
out_data = merge(wanted_red_pats, wanted_seqs.drop(['Seq'], axis = 1),
                 left_on = ['Patient ID', 'VisitNum'],
                 right_on = ['Patient ID', 'VisitNum'],
                 how = 'outer')
out_data.sort(['Patient ID', 'VisitNum'])
print out_data.to_string()
out_data.to_excel('ThreeFiveTPatInfo.xlsx')
                 

# <codecell>

from Bio import Motif
from Bio.Seq import Seq, Alphabet
from StringIO import StringIO

def make_seq(seq, comp = False):
    if comp:
        return Seq(seq,Alphabet.IUPAC.unambiguous_dna).reverse_complement()
    else:
        return Seq(seq,Alphabet.IUPAC.unambiguous_dna) 

def score_seq(mot, seq, comp = False):
    return mot.scanPWM(make_seq(seq, comp = comp))[0]


tmp = u"""A 1 2 0 0 0 2 0 0 1 2 
C 1 1 0 0 5 0 1 0 1 0
G 4 4 8 8 2 4 5 6 6 0
T 2 1 0 0 1 2 2 2 0 6"""

sp1_mot = Motif.read(StringIO(tmp), 'jaspar-pfm')

# <codecell>

test_seqs = [('sp3', 'GAGGCGTGGC'),
             ('sp2', 'TGGGCGGGAC'),
             ('sp1', 'GGGGAGTGGC')]
res = []
for name, base_seq in test_seqs:
    bs = list(base_seq)
    
    mat = np.zeros((6, len(bs)+2))
    for n in range(len(bs)):
        olet = bs[n]
        for ln, let in enumerate('ACTG'):
            bs[n] = let
            mat[ln+1, n+1] = score_seq(sp1_mot, ''.join(bs))
        bs[n] = olet
    res.append((name, base_seq, mat.copy(), score_seq(sp1_mot, base_seq)))

# <codecell>

def add_letters(lets):
    yvals = {'A':1, 'C':2, 'T':3,'G':4}
    for n, l in enumerate(lets):
        yval = yvals[l]
        plt.annotate(l, xy=(n+1, yval), fontsize=60)

def add_sig(sig_mask):
    nrows, ncols = sig_mask.shape
    rows, cols = np.where(sig_mask)
    for row, col in zip(rows, cols):
        plt.annotate('*', xy=(col, row), fontsize=60)
        

# <codecell>

from pylab import get_cmap

#plt.figure(figsize = (10, 90))
#fig, subs = plt.subplots(3,1, figsize = (5,15))
width_per_let = 0.75
for name, seq, r, conB_score in res:
    if True:
        plt.figure(figsize = (width_per_let*len(seq),3.4))
        plt.title(name)
        plt.pcolor(-(r-conB_score), cmap = get_cmap('RdBu'))
        plt.yticks([1.5,2.5,3.5,4.5], ('A', 'C', 'T', 'G'))
        plt.xticks(np.arange(1.5, len(seq)+1), range(1,len(seq)+1))
        plt.ylim([1,5])
        plt.xlim([1,len(seq)+1])
        plt.clim([-10, 10])
        add_letters(seq)
        add_sig(abs(r[0:-1, 0:-1]-conB_score)>3.84)
        plt.savefig('vSNP-%s.png' % name)
        plt.close()

# <codecell>

from functools import partial

calc_sp_val = partial(score_seq, sp1_mot)
for col in ['SP1', 'SP2', 'SP3']:
    wmask = nseq_df[col].map(lambda x: '-' not in x)
    nseq_df[col+'-bind'] = np.nan
    nseq_df[col+'-bind'][wmask] = nseq_df[col][wmask].map(calc_sp_val)
    
    
    

# <codecell>

cbseq = 'GAGGCGTGGC'
cb_bind = calc_sp_val(cbseq)

variants = []
for key, df in nseq_df.groupby('SP3'):
    if '-' not in key:
        d = ''
        for num, (c, k) in enumerate(zip(cbseq, key),1):
            if c != k:
                d += str(num)+k
        variants.append((key, d, cb_bind, df['SP3-bind'].mean(), len(df)))
var_df = DataFrame(variants, columns = ['SP3-Seq', 
                                        'SP3-Variant',
                                        'ConB-Binding',
                                        'Variant-Binding',
                                        'NumSamples'])
var_df.to_excel('sp3_variants.xlsx', index = False)

# <codecell>

cols = ['SP1-bind', 'SP2-bind', 'SP3-bind']
lets = 'ACGT'
bins = np.linspace(-5,10, 20)
tit = '%s-%s %s N:%i'
for pos in ['340', '381']:
    fig, axs = plt.subplots(4,3, sharex=True, sharey=True, figsize=(10,10))
    for ax, (let, col) in zip(axs.flatten(), product(lets, cols)):

        mask = nseq_df[pos] == let
        counts, _ = np.histogram(nseq_df[col][mask].dropna().values, bins=bins)
        fracs = counts / counts.sum()
        ax.bar(bins[1:], fracs, bins[1]-bins[0])
        ntit = tit % (pos, let, col, len(nseq_df[col][mask].dropna()))
        ax.set_title(ntit)
        if col == 'SP1-bind':
            ax.set_ylabel('Frequency')
        if let == 'T':
            ax.set_xlabel('Binding score')
    fig.tight_layout()
    plt.savefig('binding-' + pos + '-fig.png')

# <codecell>


