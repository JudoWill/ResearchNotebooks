# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import pandas as pd
import matplotlib.pyplot as plt
import os, os.path
import numpy as np
os.chdir('/home/will/HIVTropism/')

# <codecell>

tdata = pd.read_csv('BenjRes.tsv', sep='\t')

# <codecell>

wanted = tdata[tdata['Prot'] == 'gp120']

# <codecell>

plt.figure(figsize = (10,10))
for rsize in [5,10,15]:
    mask = wanted['WinSize'] == rsize
    tmp = wanted[mask]
    out = pd.rolling_mean(-np.log10(tmp['AdjPval']), rsize)
    plt.plot(tmp['Start'].values, out.values)
    

# <codecell>


