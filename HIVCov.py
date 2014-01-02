# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import os, os.path
from itertools import groupby
import numpy as np
import csv
os.chdir('/home/will/Dropbox/IpythonShare/HIVCov/')

# <codecell>

def fasta_iter(handle):
    
    for key, lines in groupby(handle, lambda x: x.startswith('>')):
        if key:
            name = lines.next().strip()[1:]
        else:
            seq = ''.join(line.strip() for line in lines)
            yield name, seq
            

# <codecell>

with open('gp120.aln') as handle:
    aligns = list(csv.reader(handle, delimiter = '\t'))
names = [n for n, s in aligns]
seqs = [s for n, s in aligns]
align_mat = np.array([list(s) for s in seqs])

# <codecell>

from itertools import product
allowed_lets = list('ACDEFGHIKLMNPQRSTVWY')

align_mask = np.array([align_mat == l for l in allowed_lets]).astype(np.float32)

# <codecell>

#(lets, seqs, columns)
check_seq = align_mat[:,5]
check_mask = align_mask[:,:,5]

# <codecell>

pAbase = align_mask.mean(axis = 1)
pAcheck = pAbase[:,5]

# <codecell>

MI = np.zeros((754,))
for x, y in product(range(len(allowed_lets)), repeat = 2):
    if pAcheck[x]==0:
        continue
    
    pXY = (np.tile(check_mask[x,:], (754,1)).T & align_mask[y,:,:]).mean(axis = 0)
    mask = (pAbase[y,:] == 0) | (pXY == 0)
    t = pXY[~mask]*np.log10(pXY[~mask]/(pAcheck[x]*pAbase[y,~mask]))
    if (t != t).any():
        print('breaking')
        break
    MI[~mask] += pXY[~mask]*np.log10(pXY[~mask]/(pAcheck[x]*pAbase[y,~mask]))

# <codecell>

def calculate_MI_double(check_mask, pAcheck, align_mask, pAbase, allowed_lets = list('ACDEFGHIKLMNPQRSTVWY')):
    MI = np.zeros((align_mask.shape[2],))
    for x, y in product(range(len(allowed_lets)), repeat = 2):
        if pAcheck[x]==0:
            continue
    
        pXY = (np.tile(check_mask[x,:], (pAbase.shape[1],1)).T * align_mask[y,:,:]).mean(axis = 0)
        mask = (pAbase[y,:] == 0) | (pXY == 0)
        MI[~mask] += pXY[~mask]*np.log10(pXY[~mask]/(pAcheck[x]*pAbase[y,~mask]))
    return MI

# <codecell>

dcheck_mask = check_mask.astype(np.float32)
dalign_mask = align_mask.astype(np.float32)

inds = np.arange(check_mask.shape[1])
Tres = calculate_MI_double(dcheck_mask[:,inds], pAcheck, dalign_mask, pAbase)
Rres = np.zeros((100,align_mat.shape[1]))
for rep in range(100):
    np.random.shuffle(inds)
    Rres[rep,:] = calculate_MI_double(dcheck_mask[:,inds], pAcheck, dalign_mask, pAbase)
    
    

# <codecell>

def calulcate_nullMI(check_mask, pAcheck, align_mask, pAbase, allowed_lets = list('ACDEFGHIKLMNPQRSTVWY'), nreps = 100):
    inds = np.arange(check_mask.shape[1])
    Rres = np.zeros((nreps,align_mat.shape[1]))
    for rep in range(nreps):
        Rres[rep,:] = calculate_MI_double(check_mask[:,inds], pAcheck, align_mask, pAbase)
    return Rres.mean(axis = 0)

# <codecell>

#nullInfo = np.zeros((align_mat.shape[1], align_mat.shape[1]))

for col in range(54, align_mat.shape[1]):
    if (col % 100 == 0) | (col == 1) | (col == 5) | (col == 10):
        print(col)
    check_mask = align_mask[:,:,col]
    nullInfo[col,:] = calulcate_nullMI(check_mask, pAcheck, align_mask, pAbase, nreps = 50)


# <codecell>

check_cols = [(i, align_mask[:,:,i], pAcheck, align_mask, pAbase) for i in range(align_mat.shape[1])]

# <codecell>

from IPython.parallel import Client

rc = Client()
lview = rc.load_balanced_view()

# <codecell>

@lview.parallel()
def multi_linker(tup):
    col, check_mask, pAcheck, align_mask, pAbase = tup
    return col, calulcate_nullMI(check_mask, pAcheck, align_mask, pAbase, nreps = 50)

# <codecell>

larg_res = multi_linker.map(check_cols[:50])

# <codecell>

for r in larg_res:
    print(r)

# <codecell>


