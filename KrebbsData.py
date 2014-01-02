# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import sys
import numpy as np
import pandas as pd

sys.path.append('/home/will/PySeqUtils/')
import PlottingTools

fname = '/home/will/KrebbsData/Results.csv'

# <codecell>

from sklearn.covariance import EllipticEnvelope

data = pd.read_csv(fname, sep = '\t')
def safe_float(val):
    try:
        return float(val)
    except ValueError:
        return np.nan

cytos = ['VEGF','IL-1beta','G-CSF','EGF','IL-10','HGF','FGF-basic',
'IFN-alpha','IL-6','IL-12','Rantes','Eotaxin','IL-13','IL-15',
'IL-17','MIP-1alpha','GM-CSF','MIP-1beta','MCP-1','IL-5',
'IFN-gamma','TNF-alpha','IL-RA','IL-2','IL-7',
'IP-10','IL-2R','MIG','IL-4','IL-8']

for col in cytos:
    data[col] = data[col].map(safe_float)
    try:
        env = EllipticEnvelope().fit(data[col].dropna().values.reshape(-1,1))
        mask = env.predict(data[col].values.reshape(-1,1))
        data[col][mask == -1] = np.nan
    except:
        pass
        
    
    #print mask
    #break
    

# <codecell>

pos = dict(zip('ABCDEFGH', range(8)))
def xpos(val):
    _, p = val.split('(')
    return pos[p.split(',')[1][0]]

def ypos(val):
    _, p = val.split('(')
    _, r = p.split(',')
    return int(''.join(l for l in r[:-1] if not l.isalpha()))
    
data['xpos'] = data['Location'].map(xpos)
data['ypos'] = data['Location'].map(ypos)

# <codecell>

data['NumBad'] = data[cytos].applymap(np.isnan).sum(axis=1)

# <codecell>

counts = pd.pivot_table(data, rows = 'xpos', cols = 'ypos', values = 'NumBad', aggfunc = 'sum').fillna(len(cytos))

# <codecell>

import matplotlib.pyplot as plt

fig = PlottingTools.make_heatmap_df(counts, figsize = (10,4), colormap = 'copper_r')
plt.colorbar()

# <codecell>

data.head().T

# <codecell>

grouped = data.groupby(['CellLine', 'Age', 'Patient', 'Exposed', 'WellSide', 'TimePoint'], as_index=False).mean()

# <codecell>

from itertools import combinations
from scipy.stats import ttest_ind
sf_data = data[(data['CellLine'] == 'Seminal Fluid')] 
old_sf =  sf_data['Age'] == 'Aged'
young_sf = sf_data['Age'] == 'Young'
pooled_sf = sf_data['Age'] == 'Pooled'

comp = [('Old', old_sf),
        ('Young', young_sf),
        ('Pooled', pooled_sf)]

results = []
for (g1_name, g1), (g2_name, g2) in combinations(comp, 2):
    for cyto in cytos:
        g1vals = sf_data[g1][cyto].dropna()
        g2vals = sf_data[g2][cyto].dropna()
        if (len(g1vals) > 2) and (len(g2vals) > 2):
            tstat, pval = ttest_ind(g1vals, g2vals)
            results.append({
                            'Group1':g1_name,
                            'Group2':g2_name,
                            'Tstat': tstat,
                            'Pval': pval,
                            'G1N':len(g1vals),
                            'G2N':len(g2vals),
                            'G1mean':g1vals.mean(),
                            'G2mean':g2vals.mean(),
                            'Cyto': cyto
                            })
    
res = pd.DataFrame(results)

# <codecell>

from statsmodels.graphics.boxplots import beanplot
labels = ['Old', 'Young', 'Pooled']
for cyto in cytos:
    bins = []
    labs = []
    for gname, gmask in comp:
        td = sf_data[cyto][gmask].dropna()
    
        if len(td) > 2:
            bins.append(td.copy())
            labs.append(gname)
    if len(bins) > 1:
        fig = plt.figure()
        ax = plt.subplot(111)
        beanplot(bins, labels=labs, ax=ax)
        ax.set_ylim([0, ax.get_ylim()[1]])
        ax.set_ylabel(cyto)
        ax.set_title(cyto)
        plt.savefig('/home/will/KrebbsData/draft_figures/violin_plot_%s.png' % cyto)
        plt.close()
    

# <codecell>

print res[(res['Pval'] < 0.1)][['Cyto','Group1', 'G1mean','Group2','G2mean','Pval']]

# <codecell>

wanted_cl = set(['Ect1', 'End1', 'VK2'])
cl_data = data[data['CellLine'].map(lambda x: x in wanted_cl)]
exposed = cl_data['Exposed'].unique()
exposed_colors = 'kbrgm'
wellsides = ['Top', 'Bottom']

# <codecell>

for cyto in cytos:
    fig, axs = plt.subplots(3,2, figsize = (10,10), sharey = True)
    for row, cl in enumerate(wanted_cl):
        for col, side in enumerate(wellsides):
            ax = axs[row, col]
            if ax.is_first_row():
                ax.set_title('%s %s' % (cyto, side))
            if ax.is_first_col():
                ax.set_ylabel('log(%s exp)' % cl)
            tmask = data['CellLine'] == cl
            tmask &= data['WellSide'] == side
            data['TimePoint'][tmask] = data['TimePoint'][tmask].map(safe_float)
            tmp = pd.pivot_table(data[tmask], cols = ['Exposed', 'TimePoint'], values = cyto, rows = 'Location')
            tmp.applymap(np.log10).boxplot(ax=ax, rot = 90)
            #data[tmask].pivot('Exposed', 'TimePoint', cyto).boxplot(ax=ax)
            #data[tmask].boxplot(column = cyto, by = ['Exposed', 'TimePoint'], ax = ax)
    fig.tight_layout()
    plt.savefig('/home/will/KrebbsData/draft_figures/boxplot_%s.png' % cyto)
    plt.close()
    

# <codecell>

cl_data['nTimePoint'] = cl_data['TimePoint'].map(str)
gdata = cl_data.groupby(['nTimePoint', 'Exposed', 'CellLine', 'WellSide']).mean()

# <codecell>

from scipy.stats import norm
effect = np.log2(gdata.ix['4']/gdata.ix['24'])[cytos]
mock_effects = effect.ix['Mock'].values.reshape(-1,1)
mock_effects[np.isnan(mock_effects)] = []
rvs = norm(*norm.fit(np.abs(mock_effects)))

pval_func = lambda x: rvs.sf(abs(x))
pvals = effect.applymap(pval_func)

# <codecell>

from statsmodels.stats.multitest import multipletests
tpvals = pvals.values.flatten()
tpvals[np.isnan(tpvals)] = []
_, bh_pvals, _, _ = multipletests(tpvals, alpha = 0.05, method = 'bonferroni')
tdict = dict(zip(tpvals, bh_pvals))
tdict[np.nan] = np.nan

bon_pvals = pvals.applymap(tdict.get)
wmask = bon_pvals < 0.1

print effect.where(wmask)['IFN-gamma'].dropna()
print bon_pvals.where(wmask)['IFN-gamma'].dropna()


# <codecell>

fname = '/home/will/KrebbsData/MouseResults.csv'
mouse_data = pd.read_csv(fname, sep = '\t')

mouse_cytos = ['IL-1beta','IL-10','IL-6','IL-12','GM-CSF','IL-5', 'IFN-gamma','TNF-alpha','IL-2','IL-4']

for col in mouse_cytos:
    mouse_data[col] = mouse_data[col].map(safe_float)
    try:
        env = EllipticEnvelope().fit(mouse_data[col].dropna().values.reshape(-1,1))
        mask = env.predict(mouse_data[col].values.reshape(-1,1))
        mouse_data[col][mask == -1] = np.nan
    except:
        pass

# <codecell>

check_cytos = ['IL-1beta','IL-6','IL-5','TNF-alpha','IL-2','IL-4']
wanted_reag = ['CpG', 'N-9', 'Imiquimod'] #PBS on all
wanted_time = ['2', '6', '12', '24']
control_mask = mouse_data['Reagent'] == 'PBS'
for cyto in check_cytos:
    fig, axs = plt.subplots(1,3, figsize = (10,5), sharey = True)
    for ax, reag in zip(axs.flatten(), wanted_reag):
        boxes = [mouse_data[cyto][control_mask].dropna().values]
        reag_mask = mouse_data['Reagent'] == reag
        for tm in wanted_time:
            boxes.append(mouse_data[cyto][reag_mask & (mouse_data['TimeGroup']==tm)].values)
        ax.boxplot(boxes)
        if ax.is_first_col():
            ax.set_ylabel(cyto)
        ax.set_title(reag)
        ax.set_xticklabels(['PBS']+wanted_time)
        
    fig.tight_layout()
    plt.savefig('/home/will/KrebbsData/draft_figures/mouse_boxplot_%s.png' % cyto)
    plt.close()

# <codecell>


# <codecell>


