# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import os, os.path
from itertools import groupby
import numpy as np
import matplotlib.pyplot as plt
from Bio import Motif
from Bio.Seq import Seq, Alphabet

os.chdir('/home/will/Dropbox/IpythonShare/HIVTfFig/')

# <codecell>

from collections import Counter
seqfiles = [('pwms/AP1.fasta', 'AP1'),
             ('pwms/coup1.fasta', 'COUP1'),
             ('pwms/coup2.fasta', 'COUP2')]
count_dict = {}
order = 'ACGT'
for f, name in seqfiles:
    seqs = [line.strip().upper() for line in open(f) if not line.startswith('>')]
    cols = []
    for col in zip(*seqs):
        cols.append(Counter(col))
    #print(cols)
    ostr = u''
    for l in order:
        ostr += l
        for count in cols:
            ostr += ' ' + str(count[l])
        ostr += '\n'
    count_dict[name] = ostr

# <codecell>

print(count_dict['COUP2'])

# <codecell>

from io import StringIO

motif_dict = {}
for key, tmp in count_dict.items():
    print key, type(tmp)
    motif_dict[key] = Motif.read(StringIO(tmp), 'jaspar-pfm')

tmp = u"""A  0  0 16  5  3  0 16
C  1  0  2 12  0 15  0
G  0 15  0  1  1  3  1
T  17  3  0  0 14  0  1"""

motif_dict['AP1'] = Motif.read(StringIO(tmp), 'jaspar-pfm')

tmp = u"""A  3  1  4  2  4  2 18 18  0
C  0  1  1  9  2 15  0  0  6
G  0  4  6  2 10  0  0  0  2
T  15 12  7  5  2  1  0  0 10"""

motif_dict['CEBP'] = Motif.read(StringIO(tmp), 'jaspar-pfm')

tmp = u"""A  0  0  0  4  2  0  1  0  6  3
C  32 30 35 27  5 28 31 24 25 26
G  1  1  0  0 15  1  0  3  0  3
T  2  4  0  4 13  6  3  8  4  3"""

motif_dict['SP1'] = Motif.read(StringIO(tmp), 'jaspar-pfm')

tmp = u"""A  0  1 12  6  0  0  0  1  2  6  6  1  3  0
C  0  0  0  7 13  3  2  0  0  4  5 10  6  3
G  2 12  1  0  0  0  0  0 11  3  1  1  0  3
T 11  0  0  0  0 10 11 12  0  0  1  1  4  7"""

motif_dict['NR2F1'] = Motif.read(StringIO(tmp), 'jaspar-pfm')

tmp = u""" A 7  4  1  1  0  0 12 11 10  0  0  0  0 11
C 3  0  1  1  3 11  0  0  0  0  0  1 12  2
G 3  6 10  5  4  0  0  2  3 13  7  0  0  0
T 0  3  1  6  6  2  1  0  0  0  6 12  1  0"""


motif_dict['HNF4A'] = Motif.read(StringIO(tmp), 'jaspar-pfm')

tmp = u"""A 0  0 6 1 0 0 0 4 2 2 0 0 3 
C  1 1 1 0 5 6 4 1 0 0 0 3 5 5 4 0  
G  0 6 0 1 1 0 0 0 0 7 1 1 0 0 1 0 
T  6 0 0 0 1 1 3 5 7 0 0 0 0 2 2 4"""

motif_dict['COUP-T2'] = Motif.read(StringIO(tmp), 'jaspar-pfm')


tmp = u"""A  8 13  0  3  2  0 14  3
C  1  0  0  0  2 13  0  8
G  3  1 13 11  0  0  0  2
T  1  0  1  0 10  1  0  0"""

motif_dict['NR4A2'] = Motif.read(StringIO(tmp), 'jaspar-pfm')

tmo = """A  4  2  4  6  7  0  0 16
C  4  5  4  4  3  0  0  0
G  5  6  3  5  2 16  0  0
T  3  3  5  1  4  0 16  0"""

motif_dict['FoxC1'] = Motif.read(StringIO(tmp), 'jaspar-pfm')

tmp = u"""A  2  9  0  0 13 12  7  0
C  9  3  0  0  0  0  0  0
G  1  1 13 13  0  0  6  0
T  1  0  0  0  0  1  0 13"""


motif_dict['FEV'] = Motif.read(StringIO(tmp), 'jaspar-pfm')

tmp = u"""A 33  4  4 42 41  5  3
C  3  0  0  0  0  3  5
G  6 37 38  0  1 33  2
T  0  1  0  0  0  1 32"""

motif_dict['SPI1'] = Motif.read(StringIO(tmp), 'jaspar-pfm')

tmp = u"""A  0   0   0  22  19  55  53  19   9
C  0 185 185  71  57  44  30  16  78
G 185   0   0  46  61  67  91 137  79
T   0   0   0  46  48  19  11  13  19"""

motif_dict['TFAP2A'] = Motif.read(StringIO(tmp), 'jaspar-pfm')

tmp = u"""A 4 17  0  0  0  5
C 16  0  1 39 39  3
G  4  0  0  1  0 17
T 16 23 39  0  1 15"""

motif_dict['ETS1'] = Motif.read(StringIO(tmp), 'jaspar-pfm')

tmp = u"""A 13  0 52  0 25
C 13  5  0  0  7
G 18 48  1  0 15
T  9  0  0 53  6"""


motif_dict['GATA2'] = Motif.read(StringIO(tmp), 'jaspar-pfm')

tmp = u"""A  7 10  6 13  4 21  0 22
C  1  4  3  4 10  0  2  1
G  4  2  6  4  2  2  0  0
T 11  7  8  2  7  0 21  0"""


motif_dict['FOXC1'] = Motif.read(StringIO(tmp), 'jaspar-pfm')

tmp = u"""A  2  7  0  6 14 14  0  0
C 13  0  7  0  0  0  0  1
G  0  5  5  1  0  2  0  6
T  0  4  4  9  2  0 16  9"""

motif_dict['HOXA5'] = Motif.read(StringIO(tmp), 'jaspar-pfm')

tmp = u"""A 25  0 61  0 39 15
C 14  1  0  0  1  3
G  4 62  1  5  4 37
T 20  0  1 58 19  8"""

motif_dict['GATA3'] = Motif.read(StringIO(tmp), 'jaspar-pfm')

tmp = u"""A  4 17  0  0  0  5
C 16  0  1 39 39  3
G  4  0  0  1  0 17
T 16 23 39  0  1 15"""

motif_dict['ETS1'] = Motif.read(StringIO(tmp), 'jaspar-pfm')

tmp = u"""A  7  3  0  0 16  0
C  6  2 16 16  0 15
G  3  0  0  0  0  1
T  0 11  0  0  0  0"""


motif_dict['ZNF354C'] = Motif.read(StringIO(tmp), 'jaspar-pfm')

tmp = u"""A  358   81   81   91 1226 3298
C 1832   67   67   88 5364  981
G  176  300 6713 6643  160 1186
T 4546 6464   51   90  162 1447"""


motif_dict['NFIC'] = Motif.read(StringIO(tmp), 'jaspar-pfm')

tmp = u"""A  0  8  0  0  0  0
C 19  2  1  0  3  0
G  0  2  4  0 19  1
T  3 10 17 22  0 21"""


motif_dict['SOX10'] = Motif.read(StringIO(tmp), 'jaspar-pfm')

tmp = u"""A  0  0 11  0  1  0  2  8
C  1  1  0  9  0  3  7  0
G  1 10  0  2 10  0  1  1
T  9  0  0  0  0  8  1  2"""

motif_dict['CREB'] = Motif.read(StringIO(tmp), 'jaspar-pfm')

tfs = list(motif_dict.keys())

for m in tfs:
    if not m.endswith('-R'):
        motif_dict[m+'-R'] = motif_dict[m].reverse_complement()
    


# <codecell>

test_seqs = [('SP-1', 1, 'GAGGCGTGGC')
             #('AP1', 105, 'GGGATCA'),
             #('AP1', 120, 'CTGACCT'),
             #('AP1', 155, 'TGAGCCA'),
             #('CEBP', 281, 'ATTTCATCA'),
             #('CREB', 330, 'TGACATCG'),
             #('CEBP', 340, 'CTTTCTACA'),
             #('SP1', 377, 'GAGGCGTGGC'),
             #('NR2F1', 95, 'ACCAGGGCCAGGGAT'),
             #('NR2F1', 118, 'CACTGACCTTTGGAT'),
             #('NR2F1-R', 98, 'AGGGCCAGGGATCAG'),
             #('NR2F1', 121, 'TGACCTTTGGATGGT'),
             #('HNF4A', 99, 'GGGCCAGGGATCAG'),
             #('COUP-T2', 95, 'ACCAGGGCCAGGGA'),
             #('COUP-T2', 118, 'CACTGACCTTTGGA'),
             #('COUP-T2', 121, 'TGACCTTTGGATGG'),
             #('COUP-T2-R', 98, 'AGGGCCAGGGATCA'),
             #('NR4A2', 64, 'AAGGCTAC'),
             #('FOXC1-R', 69, 'TACTTCCC'),
             #('FEV-R', 70, 'ACTTCCCT'),
             #('SPI1-R', 70, 'ACTTCCC'),
             ('TFAP2A', 101, 'GCCAGGGATC'),
             #('ETS1-R', 105, 'GGGATC'),
             #('GATA2', 106, 'GGATC'),
             #('GATA2-R', 107, 'GATCA'),
             #('FOXC1', 108, 'ATCAGATA'),
             #('HOXA5', 110, 'CAGATATC'),
             #('GATA2', 111, 'AGATA'),
             #('GATA3', 111, 'AGATAT'),
             #('GATA3-R', 113, 'ATATCC'),
             #('ETS1', 114, 'TATCCA'),
             #('GATA2-R', 114, 'TATCC'),
             #('ZNF354C', 115, 'ATCCAC'),
             #('NR4A2', 120, 'CTGACCTT'),
             #('AP1-R', 155, 'TGAGCCA'),
             #('NFIC-R', 157, 'AGCCAG'),
             #('GATA3', 161, 'AGAGAA'),
             #('NR4A2', 174, 'GAGGCCAA'),
             #('NFIC', 176, 'GGCCAA'),
             #('SOX10-R', 178, 'CCAATG')
             
]


# <codecell>

def make_seq(seq, comp = False):
    if comp:
        return Seq(seq,Alphabet.IUPAC.unambiguous_dna).reverse_complement()
    else:
        return Seq(seq,Alphabet.IUPAC.unambiguous_dna) 

def score_seq(mot, seq, comp = False):
    return mot.scanPWM(make_seq(seq, comp = comp))[0]



# <codecell>

res = []
lets = 'ACTG'
for name, startpos, base_seq in test_seqs:
    bs = list(base_seq)
    TM = motif_dict[name]
    mat = np.zeros((6, len(bs)+2))
    for n in range(len(bs)):
        olet = bs[n]
        for ln, let in enumerate(lets):
            bs[n] = let
            mat[ln+1, n+1] = score_seq(TM, ''.join(bs))
        bs[n] = olet
    res.append((name, startpos, base_seq, mat.copy(), score_seq(TM, base_seq)))

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
for name, startpos, seq, r, conB_score in res:
    if True:
        plt.figure(figsize = (width_per_let*len(seq),3.4))
        plt.title(name + '-' + str(startpos))
        plt.pcolor(-(r-conB_score), cmap = get_cmap('RdBu'))
        plt.yticks([1.5,2.5,3.5,4.5], ('A', 'C', 'T', 'G'))
        plt.xticks(np.arange(1.5, len(seq)+1), range(startpos,startpos+len(seq)+1))
        plt.ylim([1,5])
        plt.xlim([1,len(seq)+1])
        plt.clim([-10, 10])
        add_letters(seq)
        add_sig(abs(r[0:-1, 0:-1]-conB_score)>3.84)
        #plt.colorbar()
        plt.savefig('%s-%i-color.png' % (name, startpos))

# <codecell>

for name, startpos, seq, r, _ in res:
    if True:
        bnum = score_seq(motif_dict[name], seq)
        print(bnum)
        plt.figure(figsize = (15,5))
        plt.title(name + '-' + str(startpos))
        mask = 0.5*((-(r-bnum))>3.84) + -0.5*((-(r-bnum))<-3.84)
        plt.pcolor(mask, cmap = get_cmap('RdBu'))
        plt.yticks([1.5,2.5,3.5,4.5], ('A', 'C', 'T', 'G'))
        plt.xticks(np.arange(1.5, len(seq)+1), range(startpos,startpos+len(seq)+1))
        plt.ylim([1,5])
        plt.xlim([1,len(seq)+1])
        mval = abs(r).max()
        plt.clim([-0.5, 0.5])
        add_letters(seq)
        plt.colorbar()
        plt.savefig('%s-%i-stat.png' % (name, startpos))

# <codecell>


# <codecell>

len(check_seq)

# <codecell>


