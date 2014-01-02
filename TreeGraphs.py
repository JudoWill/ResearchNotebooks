# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

from __future__ import division
import os, os.path
import pandas as pd
import matplotlib.pyplot as plt
import sys

os.chdir('/home/will/HIVTropism/')

# <codecell>

def fix_tat(indata, key):
    fix_mask = indata['Prot'] == 'Tat'
    tmp_data = indata[fix_mask].sort(key)
    add_inds = set()
    for key, group in tmp_data.groupby(key):
        #print group
        #print (group['Start'] == 0).map(float).cumsum()
        pos_fix = (group['Start'] == 0).map(float).cumsum()==2
        add_inds |= set(pos_fix[pos_fix].index)
    return add_inds

# <codecell>

tat_ex1_width = 67

tdf = pd.read_csv('bigger_phylip_BenjRes.tsv', sep = '\t')
tdf['Start'] += (tdf['WinSize']/2).map(int)
tdf['Start'][tdf['Prot'] == 'gp120'] += 30
finds = fix_tat(tdf, ['GroupName', 'Subtype', 'Prot', 'WinSize'])
tdf['Start'].ix[list(finds)] += tat_ex1_width
df = tdf.groupby(['GroupName', 'Subtype', 'Prot', 'Start', 'WinSize'], as_index=False).first()

#wmask = (df['GroupName'] == 'full_tail') & (df['Subtype'] == 'SubB')
#df = df[wmask]

# <codecell>

#odf = pd.read_csv('phylip_BenjRes.tsv', sep = '\t')
#finds = fix_tat(odf, ['Prot', 'WinSize'])
#odf['Start'].ix[list(finds)] += tat_ex1_width
#ndf = odf.groupby(['Prot', 'Start', 'WinSize']).first()
#wdf = df.groupby(['Prot', 'Start', 'WinSize']).first()

# <codecell>

#tmp = ndf.combine_first(wdf)
#print tmp['AdjPval'].head()
#print ndf['AdjPval'].head()
#print wdf['AdjPval'].head()
#df = tmp.reset_index()

# <codecell>

def make_prot_plot(ax, xpos, mgd_logpvals, ai_logpvals, xmax = None,
                   mgd_yticks = None, mgd_maxL = None,
                   ai_yticks = None, ai_maxL = None,
                   title_prepend = '', annot_func = None):
    
    mgd_color = 'b'
    ai_color = 'r'
    mgd_ax = ax
    ai_ax = mgd_ax.twinx()
    
    mgd_ax.plot(xpos, mgd_logpvals, color=mgd_color)
    ai_ax.plot(xpos, ai_logpvals, color=ai_color)
    
    if mgd_maxL:
        mgd_ax.set_ylim([0, mgd_maxL])
    if mgd_yticks:
        mgd_ax.set_yticks(mgd_yticks)
    mgd_ax.tick_params(axis='y', colors=mgd_color)
    
    if ai_maxL:
        ai_ax.set_ylim([0, ai_maxL])
    if ai_yticks:
        ai_ax.set_yticks(ai_yticks)
    ai_ax.tick_params(axis='y', colors=ai_color)
    
    mgd_ax.set_ylabel(title_prepend + '-log10(MGD-Pval)', 
                      fontdict = {'color':mgd_color})
    ai_ax.set_ylabel(title_prepend + '-log10(AI-Pval)', 
                     fontdict = {'color':ai_color})
    
    if xmax:
        mgd_ax.set_xlim([0, xmax])
    
    if annot_func:
        annot_func(mgd_ax)

def annot_from_features(feature_list, ax,
                        **kwargs):
    
    max_val = ax.get_ylim()[1]
    for name, start, stop in feature_list:
        rect = Rectangle([start, 0], stop-start, max_val, **kwargs)
        mgd_ax.add_patch(rect)

# <codecell>

fname_tmp = 'NewerFigures/%s_%s.png'
winsizes = sorted(df['WinSize'].unique())
for gname, group_data in df.groupby('GroupName'):
    for prot, prot_data in group_data.groupby('Prot'):
        #print gname, prot
        fig, axs = plt.subplots(6,1, sharey=True, sharex=True, figsize=(10,10))
        for cut, ax in zip(winsizes, axs.flatten()):
            group = prot_data[prot_data['WinSize'] == cut]
            x = group['Start']
    
            y = pd.rolling_max(-group['AdjPval'].map(np.log10), 10)
            y2 = pd.rolling_max(-group['AI-pval'].map(np.log10), 10)
    
            make_prot_plot(ax, x, y, y2, title_prepend = '%i-mer\n' % cut)
            if ax.is_first_row():
                ax.set_title(', '.join([gname, prot]))
        fig.savefig(fname_tmp % (gname, prot))
        plt.close(fig)
            

# <codecell>


# <codecell>

fig, axs = plt.subplots(2,2, figsize = (10,10))

for ax, (prot, group) in zip(axs.flatten(), df.groupby('Prot')):
    for wsize, wgroup in group.groupby('WinSize'):
        ax.plot(wgroup['PlotPos'], 
                -pd.rolling_mean(wgroup['AdjPval'].map(np.log10), 10), 
                label = str(wsize))
    ax.set_title(prot)

# <codecell>

gp120_data = df[(df['Prot'] == 'gp120') & (df['WinSize'] == 35)]
tmp = gp120_data[['Start', 'AdjPval']]
tmp['Start'] += 17
tmp.to_excel('/home/will/Dropbox/R21Grants/2013-NINDS/gp120scanning.xlsx', index=False)

# <codecell>

from matplotlib.patches import Rectangle
gp120_features = [('V1', 131, 157),
                  ('V2', 157, 196),
                  ('V3', 296, 331),
                  ('V4', 385, 418),
                  ('V5', 460, 469)]

gp120_data = df[(df['Prot'] == 'gp120') & (df['WinSize'] == 35)]

plt.figure(figsize = (10,3))
mgd_ax = plt.gca()
ai_ax = mgd_ax.twinx()

x = gp120_data['PlotPos']
y = pd.rolling_max(-gp120_data['AdjPval'].map(np.log10), 10)
y2 = pd.rolling_max(-gp120_data['AI-pval'].map(np.log10), 10)
mgd_ax.plot(x, y)
ai_ax.plot(x, y2, 'r')
mgd_ax.set_xlim([1, 511])
ai_ax.set_xlim([1, 460])
mgd_ax.set_xticks(range(0,460,50))

mgd_ax.set_ylim([0, y.max()*1.1])

for name, start, stop in gp120_features:
    rect = Rectangle([start, 0], stop-start, y.max()*1.1, facecolor = 'r', alpha = 0.2)
    mgd_ax.add_patch(rect)
    #plt.text((start+stop)/2, y.max()*1.05, name)

plt.title('gp120 35-mer Scanning')
mgd_ax.set_ylabel('-log10(MGD-Pval)', fontdict = {'color':'b'})
ai_ax.set_ylabel('-log10(AI-Pval)', fontdict = {'color':'r'})
plt.xlabel('gp120 position')
#plt.savefig('Figures/gp120-35mer.png', dpi=1500)

# <codecell>

from functools import partial

gp120_anot = partial(annot_from_features, gp120_features, facecolor='r', alpha=0.2)
gp120_plot_func = partial(make_prot_plot, xmax=511, annot_func = gp120_anot,
                          mgd_yticks = [5, 10, 20, 30, 40], mgd_maxL = 45,
                          ai_yticks = [2, 4, 6, 8], ai_maxL = 9)

gp120_multidata = df[(df['Prot'] == 'gp120')]
fig, axs = plt.subplots(5,1,figsize = (10,10), sharex = True, sharey=True)

for mgd_ax, winsize in zip(axs.flatten(), [5, 10, 15, 20, 35]):
    group = gp120_multidata[gp120_multidata['WinSize'] == winsize]
    x = group['PlotPos']
    
    y = pd.rolling_max(-group['AdjPval'].map(np.log10), 10)
    y2 = pd.rolling_max(-group['AI-pval'].map(np.log10), 10)
    
    gp120_plot_func(mgd_ax, x, y, y2,                  
                   title_prepend = '%i-mer\n' % winsize)
mgd_ax.set_xlabel('gp120 position')
fig.savefig('Figures/gp120-multi-mer.png', dpi=1500)

# <codecell>

gp41_plot_func = partial(make_prot_plot, xmax=304, 
                          mgd_yticks = [20, 40, 60], mgd_maxL = 65,
                          ai_yticks = [1, 2, 3, 4], ai_maxL = 4.2)
gp41_multidata = df[(df['Prot'] == 'gp41')]
fig, axs = plt.subplots(5,1,figsize = (10,10), sharex = True, sharey=True)

for mgd_ax, winsize in zip(axs.flatten(), [5, 10, 15, 20, 35]):
    group = gp41_multidata[gp41_multidata['WinSize'] == winsize]
    x = group['PlotPos']
    
    y = pd.rolling_max(-group['AdjPval'].map(np.log10), 10)
    y2 = pd.rolling_max(-group['AI-pval'].map(np.log10), 10)
    
    gp41_plot_func(mgd_ax, x, y, y2,                  
                   title_prepend = '%i-mer\n' % winsize)

mgd_ax.set_xlabel('gp41 position')
fig.savefig('Figures/gp41-multi-mer.png', dpi=1500)

# <codecell>

nef_plot_func = partial(make_prot_plot, xmax=207, 
                          mgd_yticks = [5, 10, 15, 20], mgd_maxL = 21,
                          ai_yticks = range(2, 11, 2), ai_maxL = 11)
nef_multidata = df[(df['Prot'] == 'Nef')]
fig, axs = plt.subplots(5,1,figsize = (10,10), sharex = True, sharey=True)

for mgd_ax, winsize in zip(axs.flatten(), [5, 10, 15, 20, 35]):
    group = nef_multidata[nef_multidata['WinSize'] == winsize]
    x = group['PlotPos']
    
    y = pd.rolling_max(-group['AdjPval'].map(np.log10), 10)
    y2 = pd.rolling_max(-group['AI-pval'].map(np.log10), 10)
    
    nef_plot_func(mgd_ax, x, y, y2,                  
                   title_prepend = '%i-mer\n' % winsize)

mgd_ax.set_xlabel('Nef position')
fig.savefig('Figures/Nef-multi-mer.png', dpi=1500)

# <codecell>

vpr_plot_func = partial(make_prot_plot, xmax=97, 
                          mgd_yticks = None, mgd_maxL = 5,
                          ai_yticks = None, ai_maxL = 2)
vpr_multidata = df[(df['Prot'] == 'Vpr')]
fig, axs = plt.subplots(5,1,figsize = (10,10), sharex = True, sharey=True)

for mgd_ax, winsize in zip(axs.flatten(), [5, 10, 15, 20, 35]):
    group = vpr_multidata[vpr_multidata['WinSize'] == winsize]
    x = group['PlotPos']
    
    y = pd.rolling_max(-group['AdjPval'].map(np.log10), 10)
    y2 = pd.rolling_max(-group['AI-pval'].map(np.log10), 10)
    
    vpr_plot_func(mgd_ax, x, y, y2,                  
                   title_prepend = '%i-mer\n' % winsize)

mgd_ax.set_xlabel('Vpr position')
fig.savefig('Figures/Vpr-multi-mer.png', dpi=1500)

# <codecell>

vpr_plot_func = partial(make_prot_plot, xmax=97, 
                          mgd_yticks = None, mgd_maxL = 5,
                          ai_yticks = None, ai_maxL = 2)
vpr_multidata = df[(df['Prot'] == 'Vpr')]
fig, axs = plt.subplots(5,1,figsize = (10,10), sharex = True, sharey=True)

for mgd_ax, winsize in zip(axs.flatten(), [5, 10, 15, 20, 35]):
    group = vpr_multidata[vpr_multidata['WinSize'] == winsize]
    x = group['PlotPos']
    
    y = pd.rolling_max(-group['AdjPval'].map(np.log10), 10)
    y2 = pd.rolling_max(-group['AI-pval'].map(np.log10), 10)
    
    vpr_plot_func(mgd_ax, x, y, y2,                  
                   title_prepend = '%i-mer\n' % winsize)

mgd_ax.set_xlabel('Vpr position')

