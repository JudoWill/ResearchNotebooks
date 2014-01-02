# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import fisher_exact, chi2

# <codecell>

x = np.arange(0.01, 1.0, 0.05)
y = np.arange(0.01, 1.0, 0.05)
XX, YY = np.meshgrid(x, y)
npats = np.arange(10, 5000, 10)

# <codecell>

def check_fisher(x,y,npat):
    tpat = npat/2
    #print(x,y,tpat)
    xpos = int(tpat*x)
    xneg = int(tpat*(1-x))
    ypos = int(tpat*y)
    yneg = int(tpat*(1-y))
    _, pval = fisher_exact([[xpos, ypos], [xneg, yneg]])
    return pval < 0.01

# <codecell>

reqpats = np.ones_like(XX)*np.nan
for ind in range(XX.shape[0]*XX.shape[1]):
    minpats = 20
    maxpats = 1000
    while (maxpats-minpats)>10:
        npat = int(minpats + (maxpats-minpats)/2)
        #print(npat)
        if check_fisher(XX.flat[ind], YY.flat[ind], npat):
            maxpats = npat
        else:
            minpats = npat
    reqpats.flat[ind] = npat

# <codecell>

from scipy.interpolate import interp2d

func = interp2d(x,y,reqpats)

# <codecell>

nX = np.arange(0.01, 1.0, 0.01)
nY = np.arange(0.01, 1.0, 0.01)
plt.figure(figsize = (15,15))
nreq = func(nX, nY)
#print(nreq)
#raise KeyError
plt.imshow(reqpats, extent = (0,100,0,100), interpolation = 'bilinear')
plt.colorbar()
plt.ylabel('% G1 Mutant')
plt.xlabel('% G2 Mutant')
#plt.xticks(np.arange(0.01, 1.0, 0.05))

# <codecell>

from collections import defaultdict
from scipy.stats import chi2, randint

def log_factorial(n):
    return sum(np.log10(x) for x in range(1,n+1))

def multi_nomial_dist(observed_count, total_count = None):

    #print observed_count
   
    if total_count is None:
        total_count = dict(observed_count.items())
    
    for key in observed_count:
        total_count[key] = max(observed_count[key], total_count[key])
        
    tp = count2prob(dict(total_count.items()), want_dec = False)
    N = int(sum(list(observed_count.values())))
    nf_log = log_factorial(N)

    d_log = 0
    for n in observed_count.values():
        d_log += log_factorial(n)
        
    p = nf_log-d_log
    for k, nnp in tp.items():
        p += observed_count[k]*np.log10(nnp)#.log10()
        
    return 10**float(p)

def countdict(intup):
    r = defaultdict(int)
    for n in intup:
        if n.isalpha() or len(n)>1:
            r[n] += 1
    return r

def count2prob(d, want_dec = False):
    n = sum(list(d.values()))
    for key in d.keys():
        d[key] = d[key]/n
    return d
 
    
def likelihood_ratio(g1align, g2align):
    
    g1count = countdict(g1align)
    g2count = countdict(g2align)
    if (sum(list(g1count.values())) < 5) | (sum(list(g2count.values())) < 5):
        return None, None
    
    self_p = multi_nomial_dist(g1count)
    g2_p = multi_nomial_dist(g1count, total_count = g2count)
    
    ratio = -2*(np.log(g2_p)-np.log(self_p))
    df = len(g1count)
    pval = 1-chi2.cdf(ratio, df)
    #print self_p, r5_p, ratio, df, pval
    return ratio, pval

def gen_seqs(base_seq):
    rv = randint(0, len(base_seq))
    while True:
        yield base_seq[rv.rvs()]

# <codecell>

def predict_npats(check_fun, s1, s2, s1iter, s2iter, max_pats = 500):
    
    pval = 1
    npats = 0
    while (pval < 0.05) & (npats < max_pats):
        s1 += next(s1iter)
        s2 += next(s2iter)
        
    

# <codecell>

from scipy.stats import scoreatpercentile
g1 = 'ATTTTCTTTATCT'
g2 = 'ATTCCCGCTGCTT'




types = [('optimistic', gen_seqs(g1), gen_seqs(g2)),
        ('neutral', gen_seqs(g1+g2), gen_seqs(g1+g2))]

for name, g1iter, g2iter in types:
    numpats = []
    for nrep in range(100):
        #print(nrep)
        ng1 = g1
        ng2 = g2
        pval = 1
        while (pval > 0.05) and (len(ng1) < 500): 
            ng1 += next(g1iter)
            ng2 += next(g2iter)
            _, pval = likelihood_ratio(ng1, ng2)
        numpats.append(len(ng1))
    bmean = np.mean(numpats)
    print(name, bmean)
        

# <codecell>

from scipy.stats import norm
kde = norm(np.mean(numpats), np.std(numpats))
print(name,kde.mean(), kde.interval(0.95))

# <codecell>

tmp.rvs()

# <codecell>


