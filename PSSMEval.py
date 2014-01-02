# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import sys
sys.path.append('/home/will/PySeqUtils/')
import os, os.path
os.chdir('/home/will/HIVTropism/R5Cluster/')
from itertools import islice
import pandas as pd
from GeneralSeqTools import fasta_reader, call_muscle

# <codecell>

import autoreload
import SeqSklearn
%load_ext autoreload
%autoreload 2

# <codecell>

from sklearn.preprocessing import label_binarize
def decide_tropism(inval):
    if inval < -6.95:
        return 'R5'
    elif inval > -2.88:
        return 'X4'
    elif inval < -4.92:
        return 'R5-P'
    elif inval >= -4.92:
        return 'X4-P'
    
    return np.nan

# <codecell>

lanl_data = pd.read_csv('LANLResults.tsv', sep='\t')
def safe_float(num):
    
    try:
        return float(num)
    except ValueError:
        return np.nan
    
def convert_group_series(inser, cluster_only=False):
    
    groups = inser.dropna().unique()
    c_col = inser.name+'-cluster'
    if cluster_only:
        odf = pd.DataFrame(np.nan, 
                           columns=[c_col],
                           index=inser.index)
        for cnum, ocol in enumerate(groups,1):
            odf[c_col][inser==ocol] = cnum
        
    else:
        new_labels = [inser.name+'-'+col for col in groups]
        odf = pd.DataFrame(np.nan, 
                           columns=new_labels+[c_col],
                           index=inser.index)
    
        for cnum, (ocol, ncol) in enumerate(zip(groups, new_labels),1):
            odf[ncol] = (inser==ocol).astype(float)
            odf[ncol][inser.isnull()] = np.nan
            odf[c_col][inser==ocol] = cnum
        
    return odf
    
    
    

treated_lanl_dict = {'CD8':lanl_data['CD8 count'],
                     'CD4':lanl_data['CD4 count'],
                     'Log-VL':lanl_data['Viral load'].map(np.log10),
                     'Years-Infected':lanl_data['Days from Infection']/365,
                     'Years-Seropositve':lanl_data['Days from Seroconversion'].map(safe_float)/365,
                     'Accession':lanl_data['Accession'],
                     'Gender':lanl_data['Patient Sex'],
                     'Tissue':lanl_data['Sample Tissue'],
                     'STissue':lanl_data['NewSimpleTissue'],
                     'DrugNaive':lanl_data['Drug Naive'],
                     'GeoRegion':lanl_data['Georegion'],
                     
                }

treated_lanl = pd.DataFrame(treated_lanl_dict)
treated_lanl = pd.concat([treated_lanl, 
                          convert_group_series(treated_lanl['Gender']),
                          convert_group_series(treated_lanl['Tissue'], cluster_only=True),
                          convert_group_series(treated_lanl['STissue'], cluster_only=True),
                          convert_group_series(treated_lanl['DrugNaive']),
                          convert_group_series(treated_lanl['GeoRegion'], cluster_only=True),
                          ], 
                         axis=1)
treated_lanl

# <codecell>

def safe_mean(inser):
    return inser.mean()

def most_common(inser):
    try:
        return inser.value_counts().index[0]
    except IndexError:
        return np.nan

V3_tropism = pd.read_csv('LANLTropism.tsv', sep='\t')
V3_tropism = V3_tropism[V3_tropism['percentile']<0.95]
trop_data = V3_tropism.groupby('name')['score'].mean()
pssm_bins = [-15.0, -13.0, -11.0, -9.0, -6.96, -4.92, -2.88, 1]
trop_data = pd.DataFrame({
                           'PSSMScore':trop_data,
                           'PSSMTrop':trop_data.map(decide_tropism),
                           'PSSMSubs':pd.Series(np.digitize(trop_data.values, pssm_bins), 
                                                index=trop_data.index),
                           })
ntrop_data = pd.concat([trop_data, convert_group_series(trop_data['PSSMTrop'])], axis=1)
print ntrop_data.head()

# <codecell>

nlanl_data = pd.merge(treated_lanl, ntrop_data,
                      left_on='Accession',
                      right_index=True, 
                      how='outer').groupby('Accession').first()
nlanl_data

# <codecell>

tmp_data = pd.read_csv('cluster_dump.tsv', sep=',')

# <codecell>

pivot_data = pd.pivot_table(tmp_data, rows='Accession', 
                            cols=['Prot','Start'], 
                            values='Cluster')

# <codecell>

from sklearn.metrics import adjusted_rand_score
from scipy.stats import norm
from scipy.stats import mode

ns_score = adjusted_rand_score
pvals = []
nreps = 250
nclusters = []
prots = ['LTR', 'gp120', 'Int', 'gp41']

for prot in prots:
    for col_start in pivot_data[prot].columns:
        if (col_start % 50 == 0):
            print prot, col_start
    
        col, ps_sub = nlanl_data['PSSMSubs'].dropna().align(pivot_data[prot][col_start].dropna(), join='inner')
        ts = ns_score(col.values, ps_sub.values)
        runs = [ns_score(col.values, np.random.permutation(ps_sub.values)) for _ in range(nreps)]
        args = norm.fit(runs)
        fit_model = norm(*args)
        pval = fit_model.sf(ts)
        pvals.append({
                      'Start':col_start,
                      'Prot':prot,
                      'TrueScore':ts,
                      'MeanScore':args[0],
                      'Pval':pval
                      })
        
        tmp = pd.DataFrame({'V3':col, 'other':ps_sub})
        out_vals = tmp.groupby('other')['V3'].transform(lambda x: mode(x)[0])
        nclusters += [(prot, col_start+17, acc, ncl) for acc, ncl in  out_vals.to_dict().items()]

        
    

# <codecell>

from statsmodels.stats.multitest import multipletests

pval_df = pd.DataFrame(pvals)



_, bh_pvals, _, _ = multipletests(pval_df['Pval'], alpha = 0.01, method = 'fdr_bh')
pval_df['BHPval'] = bh_pvals

# <codecell>

mask = (pval_df['Prot']=='gp120') & (pval_df['Start']>325)& (pval_df['Start']<375)
pval_df[mask]['BHPval'].head(n=50).map(np.log10)

# <codecell>

wprots = set(['LTR', 'gp120', 'Int', 'gp41'])
sig_dict = {}
for prot, rows in pval_df.groupby('Prot'):
    if prot in wprots:
        tot = pd.Series(np.nan, index=range(rows['Start'].max()+17))
        lpvals = -np.log10(rows['BHPval']).replace({-np.inf:-500})
        lpvals.index = rows['Start']+17
        smoothed_p = pd.rolling_mean(lpvals, 5)
        sig_dict[prot], _ = smoothed_p.align(tot)
        plt.figure(figsize=(8,4))
        plt.plot(smoothed_p.index, smoothed_p)
        plt.ylabel('-log10(pval)')
        plt.xlabel('ConB Position')
        plt.title(prot)
        plt.savefig('final_figures/%s_clustering.png' % prot, dpi=1000)

# <codecell>

cnames = ['Prot', 'FStart', 'Accession', 'NewCluster']
new_clust_df = pd.pivot_table(pd.DataFrame(nclusters, columns = cnames),
                              rows='Accession',
                              cols=['Prot', 'FStart'],
                              values='NewCluster',
                              aggfunc='first')


# <codecell>

from PlottingTools import make_heatmap
from matplotlib import colors
from sklearn.cluster import KMeans

cmap = colors.ListedColormap(['white', 'red', 'blue', 'green', 'cyan', 'purple', 'magenta'])

for prot, sig_vals in sig_dict.items():

    tdf, tmask = new_clust_df[prot].dropna(axis=0).align(sig_vals, axis=1, join='outer')
    tdf, scores = tdf.align(nlanl_data['PSSMScore'].dropna(), axis=0, join='inner')
    tvals = tdf.values

    order = np.argsort(scores.values)
    tvals[:, tmask.values<1.5] = np.nan
    tvals = tvals[order,:]
    
    xpos = np.array(tdf.columns)
    if prot == 'gp120':
        xpos += 30
        pass
    
    xx, yy = np.meshgrid(xpos, range(tvals.shape[0]))

    #print xx
    fig, (ax1, ax2) = plt.subplots(2,1, figsize=(10,5), sharex=True)

    #plt.figure(figsize=(10,5))
    ax1.set_title(prot)
    ax1.plot(xpos, sig_vals)
    ax1.set_ylabel('-log10(p-value)')
    ax1.set_xlim([0, max(xpos)+17])

    #plt.figure()
    ax2.pcolormesh(xx,yy, tvals, vmin=0, vmax=8, cmap=cmap)
    #ax2.imshow(tvals+1, aspect='auto', cmap=cmap)
    ax2.set_yticks([])
    ax2.set_ylabel('ordered by PSSM with R5 at top')
    ax1.set_xlim([0, max(xpos)+17])
    plt.savefig('final_figures/%s_scanning_group.png' % prot, dpi=1000)

# <codecell>

10**(-sig_dict['Int'].max())

# <codecell>

best_pos = sig_dict['Int'].idxmax()
print best_pos
int_clusts = pivot_data['Int'][best_pos-17]
wanted_lanl, int_pos_clusts = nlanl_data.align(int_clusts.dropna(), axis=0, join='inner')
wanted_lanl['IntClusts'] = int_pos_clusts

med_pssm = wanted_lanl.groupby('IntClusts')['PSSMScore'].median()
ranks = med_pssm.rank(method='first')
rep_dict = ranks.to_dict()
wanted_lanl['IntClusts'] = wanted_lanl['IntClusts'].replace(rep_dict)


keys = [('Log-VL', 'boxplot'),
        ('CD4', 'boxplot'),
        ('Years-Seropositve', 'boxplot'),
        ('CD8', 'boxplot'),
        ('PSSMScore', 'boxplot'),
        ('Count', 'count'),]

fig, axs = plt.subplots(3,2, figsize=(10,10), sharex=True)

for ax, (col,typ) in zip(axs.flatten(), keys):
    
    if col == 'Count':
        boxes = [len(vals) for _, vals in wanted_lanl.groupby('IntClusts')]
    else:
        boxes = [vals.dropna() for _, vals in wanted_lanl.groupby('IntClusts')[col]]
        nums = [len(vals.dropna()) for _, vals in wanted_lanl.groupby('IntClusts')[col]]
        print col, nums
    if typ == 'boxplot':
        ax.boxplot(boxes)
    elif typ == 'frac':
        ax.bar(np.arange(0.5,len(boxes)), [b.mean() for b in boxes])
        ax.set_ylim([0, 1])
    elif typ == 'count':
        ax.bar(np.arange(0.5,len(boxes)), boxes)
        
        
    ax.set_title('IntClusts-' + col)
    #if ax.is_last_row():
    #    num_clust = int(wanted_lanl['IntClusts'].max())
    #    ax.set_xticklabels(['CL-%i' % i for i in range(1, len(boxes)+1)])

# <codecell>


