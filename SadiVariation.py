# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import os, os.path
import sys
import pandas as pd
import numpy as np

os.chdir('/home/will/SadiVariation/')
sys.path.append('/home/will/PySeqUtils/')

# <codecell>

from GeneralSeqTools import fasta_reader, fasta_writer, WebPSSM_V3_series
import glob

# <codecell>

files = [('x4_seqs.fasta.old', 'x4_seqs.fasta'),
         ('r5_seqs.fasta.old', 'r5_seqs.fasta')]
for ifile, ofile in files:
    with open(ifile) as handle:
        with open(ofile, 'w') as ohandle:
            for name, seq in fasta_reader(handle):
                fasta_writer(ohandle, [(name, seq[1:-1])])

# <codecell>

subtype_files = glob.glob('/home/will/WLAHDB_data/SubtypeGuess/*.gb')
subtypes = []
for f in subtype_files:
    gb = f.rsplit(os.sep, 1)[-1].split('.')[0]
    with open(f) as handle:
        subtype = handle.next().strip()
        if subtype != 'Unk':
            subtypes.append((int(gb), subtype))
subtype_df = pd.DataFrame(subtypes, columns = ['GI', 'Subtype'])

subtype_ser = subtype_df.groupby('GI')['Subtype'].first()

# <codecell>

with open('hxb2.needle') as handle:
    aligned_seqs = list(fasta_reader(handle))

# <codecell>

from scipy.stats import linregress

hxb2_inds = [350, 364, 375, 388, 398, 456]
our_inds = [-103, -99, -85, -74, -63, 0]
m, b, _, _, _ = linregress(hxb2_inds, our_inds)
new_our_inds = np.arange(0, 634)
new_hxb2_inds = np.ceil(m*new_our_inds+b)

# <codecell>

starts = range(0, len(aligned_seqs), 2)
aligned = []
for s in starts:
    _, hxb2_seq = aligned_seqs[s]
    gi, gi_seq = aligned_seqs[s+1]
    aseq = ''.join(q for q, r in zip(gi_seq, hxb2_seq) if r.isalpha())
    aligned.append((int(gi), np.array(list(aseq))))

# <codecell>

aligned_ser = pd.DataFrame(aligned, columns=['GI', 'alignment']).groupby('GI').first()['alignment']
subs, _ = subtype_ser.align(aligned_ser, join='right')

# <codecell>

subs.value_counts()

# <codecell>

aset = set(gi for gi, _ in aligned)
with open('/home/will/WLAHDB_data/SeqDump/B_v3.fasta') as handle:
    v3_seqs = []
    for gi, seq in fasta_reader(handle):
        if int(gi) in aset:
            v3_seqs.append((int(gi), seq))
    print len(v3_seqs)

# <codecell>

trop_dict = dict((int(gi), 'X4' if trop > 0.5 else 'R5') for gi, trop in WebPSSM_V3_series(v3_seqs))
v3_ser = pd.Series(trop_dict)
trop_data, _ = v3_ser.align(aligned_ser, join='right')

# <codecell>

trop_data.value_counts()

# <codecell>

lanl_data = pd.read_csv('/home/will/HIVTropism/R5Cluster/LANLResults.tsv', sep = '\t')

# <codecell>

wtissue, _ = lanl_data['NewSimpleTissue'].dropna().align(aligned_ser, join = 'right')
wcor_receptor, _ = lanl_data['Coreceptor'].dropna().align(aligned_ser, join = 'right')

# <codecell>

check_seqs = [('CEBP-US2', 'ATTTCATCA', -170, -162),
              ('ATF-CREB', 'CTGACATCG', -123, -115),
              ('CEBP-US1', 'AGCTTTCTACAA', -114, -103),
              ('NFKB-II', 'AGGGACTTTCC', -103, -93),
              ('NFKB-I', 'GGGGACTTTCC', -99, -89),
              ('SP-III', 'GAGGCGTGG', -85, -77),
              ('SP-II', 'TGGGCGGGA', -74, -66),
              ('SP-I', 'GGGGAGTGG', -63, -55),
              ('AP1-I', 'TTGAGTGCT', 85, 93),
              ('AP1-II', 'TGTTGTGTGAC', 121, 138),
              ('AP1-III', 'TTTAGTCAG', 153, 161),
              ('DS3-B', 'TCAGTGTGGAAAATC', 158, 175),
              ('DS3-C', 'GTAGTGTGGAAAATC', 158, 175),
              ('DS3-D', 'TCAGTGTGGAAAATC', 158, 175),
              ('DS3-A', 'ACTGTGTAAAAATC', 158, 175)
              ]
slop = 20
check_dict = {'Subs':subs, 'tissue': wtissue, 'Coreceptor': trop_data}
for name, seq, start, stop in check_seqs:
    mask = (new_hxb2_inds>(start-slop)) & (new_hxb2_inds<(stop+slop))
    extracted = aligned_ser.map(lambda x: ''.join(x[mask]).replace('-', ''))
    
    check_dict[name] = extracted.map(lambda x: seq in x).map(float)
    check_dict[name][extracted.map(len)==0] = np.nan
df = pd.DataFrame(check_dict)

# <codecell>


# <codecell>

print pd.pivot_table(df, rows = 'Subs', aggfunc='mean').T*100

# <codecell>

print 9.0/290

# <codecell>

print pd.pivot_table(df, rows = 'Coreceptor', aggfunc='mean').T*100

# <codecell>

def find_best(tf_seq, ltr_seq):
    width = len(tf_seq)
    scores = []
    for start in range(0, len(ltr_seq)-width):
        scores.append(sum(s == t for s, t in zip(tf_seq, ltr_seq[start:(start+width)])))
    if scores:
        return max(scores)/float(width)
    else:
        return np.nan

best_dict = {'Subs':subs, 'tissue': wtissue, 'Coreceptor': trop_data}
for name, seq, start, stop in check_seqs:
    print name
    mask = (new_hxb2_inds>(start-slop)) & (new_hxb2_inds<(stop+slop))
    extracted = aligned_ser.map(lambda x: ''.join(x[mask]).replace('-', ''))
    best_dict[name] = extracted.map(lambda x: find_best(seq, x))
bdf = pd.DataFrame(best_dict)

# <codecell>

print (1-pd.pivot_table(bdf, rows = 'Subs')).T*100

# <codecell>

print (1-pd.pivot_table(bdf, rows = 'Coreceptor')).T*100

# <codecell>

slop=5
mask = (new_hxb2_inds>(159-slop)) & (new_hxb2_inds<(174+slop))
ds3_ser = aligned_ser.map(lambda x: ''.join(x[mask]).replace('-', ''))
ds3_ser[ds3_ser.map(len)==0] = np.nan
with open('r5_seqs.fasta', 'w') as handle:
    fasta_writer(handle, ds3_ser[trop_data == 'R5'].dropna().to_dict().items())

with open('x4_seqs.fasta', 'w') as handle:
    fasta_writer(handle, ds3_ser[trop_data == 'X4'].dropna().to_dict().items())
    
with open('subC_seqs.fasta', 'w') as handle:
    fasta_writer(handle, ds3_ser[subs == 'C'].dropna().to_dict().items())
          
#print ds3_ser[trop_data == 'X4'].dropna()

# <codecell>

quick_seqs = []
with open('r5_seqs.fasta') as handle:
    for name, seq in fasta_reader(handle):
        quick_seqs.append({
                           'Name':name,
                           'Trop':'R5',
                           'Seq':seq
                           })
with open('x4_seqs.fasta') as handle:
    for name, seq in fasta_reader(handle):
        quick_seqs.append({
                           'Name':name,
                           'Trop':'X4',
                           'Seq':seq
                           })

# <codecell>

from Bio.Seq import Seq
from Bio import Motif
from StringIO import StringIO
from itertools import groupby
from operator import methodcaller
from Bio.Alphabet import IUPAC

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
pwm_dict = {}
for num, (name, mot) in enumerate(yield_motifs()):
    if num % 100 == 0:
        print num
    pwm_dict[name] = mot

# <codecell>

from functools import partial
from scipy.stats import ttest_ind, gaussian_kde, chi2_contingency

def score_seq(mot, seq):
    bseq = Seq(seq, alphabet=IUPAC.unambiguous_dna)
    scores = mot.scanPWM(bseq)
    return np.max(scores)

def make_cdfs(kde, points):
    
    cdf = []
    for point in points:
        cdf.append(kde.integrate_box_1d(-np.inf, point))
    return 1-np.array(cdf)
    
    
wanted_mots = ['cebpa-R',
               #'nfatc2',
               'nfatc2-R']
fig, axs = plt.subplots(2,1, sharex=True, figsize = (10, 5))
quick_seqs_df = pd.DataFrame(quick_seqs)

r5_mask = quick_seqs_df['Trop'] == 'R5'
x4_mask = quick_seqs_df['Trop'] == 'X4'
for ax, mot in zip(axs.flatten(), wanted_mots):
    quick_seqs_df[mot] = quick_seqs_df['Seq'].map(partial(score_seq, pwm_dict[mot]))
    
    r5_vals = quick_seqs_df[mot][r5_mask].dropna().values
    x4_vals = quick_seqs_df[mot][x4_mask].dropna().values
    
    r5_kde = gaussian_kde(r5_vals)
    x4_kde = gaussian_kde(x4_vals)
    
    points = np.linspace(0, 15)
    ax.plot(points, make_cdfs(r5_kde, points), 'b', label = 'R5')
    ax.plot(points, make_cdfs(x4_kde, points), 'g', label = 'X4')
    ax.set_title(mot)
    if ax.is_last_row():
        ax.set_xlabel('TF Score')
    else:
        ax.legend()
    
    ax.set_ylabel('Frac Sequences')
    thresh = Motif.Thresholds.ScoreDistribution(pwm_dict[mot], precision = 100).threshold_fpr(0.005)
    ax.vlines(thresh, 0, 1)
    
    ch2table = [[(r5_vals>thresh).sum(), (r5_vals<=thresh).sum()],
                [(x4_vals>thresh).sum(), (x4_vals<=thresh).sum()],]
    
    _, pval, _, _ = chi2_contingency(ch2table)
    print mot, np.mean(r5_vals), np.mean(x4_vals), pval
plt.tight_layout()
#plt.savefig('/home/will/Dropbox/Wigdahl HIV Lab/SadiTFFigure/TFscores.png', dpi = 300)

# <codecell>

ax.vlines?

# <codecell>

from sklearn.svm import OneClassSVM

data = np.array([2,5,30,4,8,3,5,4,2,5,3,4,5,4]).reshape((-1, 1))
tmp = OneClassSVM().fit(data).predict(data)
data[tmp>0]

# <codecell>

count = 0
for num, f in enumerate(glob.glob('/home/will/WLAHDB_data/RegionSplit/ltr/*.fasta')):
    if num % 5000 == 0:
        print num, count
    with open(f) as handle:
        if len(handle.read()) > 10:
            count += 1

# <codecell>


