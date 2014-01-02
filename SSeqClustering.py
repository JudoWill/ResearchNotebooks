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

lanl_data = pd.read_csv('LANLResults.tsv', sep='\t')

# <codecell>

lanl_data

# <codecell>

from sklearn.preprocessing import label_binarize
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
    
    
    

treated_lanl_dict = {'CD8':lanl_data['CD8 count'].clip(0, 3000),
                     'CD4':lanl_data['CD4 count'].clip(0, 2000),
                     'Log-VL':lanl_data['Viral load'].map(np.log10),
                     'Years-Infected':lanl_data['Days from Infection']/365,
                     'Years-Seropositve':lanl_data['Days from Seroconversion'].map(safe_float)/365,
                     'Accession':lanl_data['Accession'],
                     'Gender':lanl_data['Patient Sex'],
                     'Tissue':lanl_data['Sample Tissue'],
                     'STissue':lanl_data['NewSimpleTissue'],
                     'DrugNaive':lanl_data['Drug Naive'],
                     'GeoRegion':lanl_data['Georegion'],
                     'PatientID':lanl_data['Patient Id']
                     
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
pssm_bins = [-13.0, -11.0, -9.0, -6.96, -4.92, -2.88]

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
                      how='outer')
nlanl_data

# <markdowncell>

# Just to make sure. There should be an increased Log-VL in X4 as compared to R5.

# <codecell>

nlanl_data[['Log-VL', 'PSSMTrop']].dropna().boxplot(column=['Log-VL'], 
                                                            by='PSSMTrop')

# <codecell>

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from collections import Counter, defaultdict
from GeneralSeqTools import fasta_writer

def filter_seq(handle, trans):
    for name, seq in fasta_reader(handle):
        tseq = ''.join(l for l in seq if l.isalpha())
        l = len(tseq)
        if (l == 105):
            if trans:
                rseq = Seq(tseq, generic_dna).translate()
                yield name, rseq.tostring()
            else:
                yield name, tseq



with open('V3filter.nt.fasta.raln') as handle:
    seq_list = list(fasta_reader(handle))
                

with open('V3filter.aa.fasta.raln') as handle:
    aa_seq_list = list(fasta_reader(handle))
    

# <codecell>

import numpy as np
from sklearn.base import BaseEstimator, TransformerMixin, ClassifierMixin
from itertools import product
from scipy.sparse import csr_matrix, eye

class BioTransformer(BaseEstimator, TransformerMixin):
    
    def __init__(self, typ = 'nuc'):
        self.typ = typ
        
    def fit(self, *args):
        return self
    
    def transform(self, X):
        if self.typ == 'nuc':
            letters = 'ACGT-'
        else:
            letters = 'ARNDCEQGHILKMFPSTWYV-'
        
        nrows, ncols = X.shape
        #out = eye(nrows, ncols*len(letters), format='csr')
        data = []
        rows = []
        cols = []
        for row in range(nrows):
            for num, (col,l) in enumerate(product(range(ncols), letters)):
                if X[row, col].upper()==l:
                    data.append(1)
                    rows.append(row)
                    cols.append(num)
                
        return csr_matrix((np.array(data), (np.array(rows), np.array(cols))), 
                          shape=(nrows, ncols*len(letters)), dtype=float).todense()

# <codecell>

tmp_seqs = np.array([list(seq) for name, seq in seq_list])
names = [name for name, _ in seq_list]
nuc_df = pd.DataFrame(BioTransformer(typ='nuc').transform(tmp_seqs), index=names)

tmp_aa_seqs = np.array([list(seq) for name, seq in aa_seq_list])
aa_names = [name for name, _ in aa_seq_list]
aa_df = pd.DataFrame(BioTransformer(typ='aa').transform(tmp_aa_seqs), index=aa_names)


# <codecell>

from sklearn.pipeline import Pipeline
from sklearn.cluster import KMeans, MiniBatchKMeans, AffinityPropagation
from sklearn.grid_search import GridSearchCV
from sklearn.metrics import silhouette_score, adjusted_rand_score, normalized_mutual_info_score
from sklearn.metrics.pairwise import euclidean_distances
from sklearn.cross_validation import Bootstrap


def sil_linker(predictor, X):
    clusters = predictor.predict(X)
    if len(set(clusters)) == 1:
        clusters[-1] += 1
    return silhouette_score(X, clusters)

def rand_linker(predictor, X, y):
    
    clusters = predictor.predict(X)
    return normalized_mutual_info_score(y, clusters)


#nuc_dist = euclidean_distances(tmp_bin)
#aa_dist = euclidean_distances(tmp_aa_bin)

# <markdowncell>

# Start with some simple analysis using only North American sequences from PBMC/serum and sampled in North America

# <codecell>

pat_lanl_data = nlanl_data.groupby('Accession').first()

# <codecell>

wanted_mask = (pat_lanl_data['GeoRegion'] == 'North America') & \
((pat_lanl_data['STissue'] == 'PBMC') | (pat_lanl_data['STissue'] == 'plasma/serum/blood'))
NA_blood_df, NA_wanted_lanl = aa_df.align(pat_lanl_data[wanted_mask], axis=0, join='inner')

# <codecell>

from sklearn.cross_validation import permutation_test_score

def check_num_clusts(X_data, n_clusts=range(2, 30), n_iter=10):
    
    cv = Bootstrap(X_data.shape[0], n_iter=n_iter, train_size=0.7)
    
    pred = GridSearchCV(KMeans(), {'n_clusters':n_clusts}, 
                        n_jobs=-1, pre_dispatch=30, cv = cv,
                        scoring=sil_linker,
                        refit=True)
    pred.fit(X_data)
    res = []
    for d in pred.grid_scores_:
        for r in d.cv_validation_scores:
            res.append((d.parameters['n_clusters'], r))
    return pred, pd.DataFrame(res, columns=['n_clusters','score'])

def check_trop_score(X_data, trop_clusters):
    
    cv = Bootstrap(X_data.shape[0], n_iter=3, train_size=0.7)
    pred = KMeans(n_clusters=len(set(trop_clusters)))
    t_score, scores, pval = permutation_test_score(pred, X_data, 
                                                   n_permutations=100,
                                                   y = trop_clusters,
                                                   n_jobs=20,
                                                   scoring=rand_linker,
                                                   cv=cv)
    return t_score, scores, pval

# <codecell>

_, trop_data = check_num_clusts(NA_blood_df.values, 
                                n_iter=10, 
                                n_clusts=range(2, 60))
t = trop_data.groupby('n_clusters')['score'].mean()
e = trop_data.groupby('n_clusters')['score'].std()
plt.errorbar(t.index, t.values, yerr=e.values)
plt.title('Clustering of North American V3 sequnces')
plt.xlabel('Cluster Size')
plt.xlim([1.5, 60])
plt.ylim([0, 1])
plt.ylabel('Silhouette Score')
plt.savefig('final_figures/long_NA_v3_clustering.png', dpi = 1000)

# <codecell>

from SeqSklearn import BinBasedCluster

bin_clust = BinBasedCluster(bins = pssm_bins)

pca_trans, biny, xx, yy, Z = bin_clust.make_vern_points(NA_blood_df.values, NA_wanted_lanl['PSSMScore'])


# <codecell>

from pylab import get_cmap
plt.figure(figsize=(10,10))
jitter = 0.1*np.random.randn(*pca_trans.shape)+pca_trans
plt.scatter(jitter[:,0], jitter[:,1], vmax = 0, c=NA_wanted_lanl['PSSMScore'], cmap=get_cmap('copper_r'), alpha=0.5)
cbar = plt.colorbar()
cbar.set_label('PSSMScore')
plt.ylabel('PC-1')
plt.xlabel('PC-2')

# <codecell>

v3_clust_pred = KMeans(n_clusters=7, n_jobs=-1)
v3_cluster_vals = pd.Series(v3_clust_pred.fit_predict(NA_blood_df.values), NA_blood_df.index)

NA_wanted_lanl['V3Clusters'] = v3_cluster_vals

# <codecell>

vals = NA_wanted_lanl['PSSMScore'].dropna().values
counts, edges = np.histogram(vals, bins = 500, normed=True)

# <codecell>

print edges[x4][:-1].shape, counts[x4[:-1]].shape

# <codecell>

plt.figure(figsize=(10,10))
r5 = edges<-6.96
r5p = (edges>-6.96) & (edges<-4.92)
x4p = (edges>-4.92) & (edges<-2.88)
x4 = (edges>-2.88)
ax = plt.gca()
ax.bar(edges[r5], counts[r5[:-1]], width=0.001, edgecolor = 'green')
ax.bar(edges[r5p], counts[r5p[:-1]], width=0.001, edgecolor = 'orange')
ax.bar(edges[x4p], counts[x4p[:-1]], width=0.001, edgecolor = 'blue')
ax.bar(edges[x4][:-1], counts[x4[:-1]], width=0.001, edgecolor = 'red')
ax.set_ylabel('Frac Sequences')
ax.set_xlabel('PSSM Score')
plt.savefig('final_figures/pssm_scores.png', dpi=1000)

# <codecell>


# <codecell>

pssm_ranks = NA_wanted_lanl.groupby('V3Clusters')['PSSMScore'].median().rank(method='first').to_dict()
NA_wanted_lanl['V3Clusters'] = NA_wanted_lanl['V3Clusters'].replace(pssm_ranks)

# <codecell>

keys = [('Log-VL', 'boxplot'),
        ('CD4', 'boxplot'),
        ('Years-Seropositve', 'boxplot'),
        ('CD8', 'boxplot'),
        ('PSSMScore', 'boxplot'),
        ('Count', 'count'),]

fig, axs = plt.subplots(3,2, figsize=(10,10), sharex=True)

for ax, (col,typ) in zip(axs.flatten(), keys):
    
    if col == 'Count':
        boxes = [len(vals) for _, vals in NA_wanted_lanl.groupby('V3Clusters')]
    else:
        boxes = [vals.dropna() for _, vals in NA_wanted_lanl.groupby('V3Clusters')[col]]
        nums = [len(vals.dropna()) for _, vals in NA_wanted_lanl.groupby('V3Clusters')[col]]
        print col, nums
    if typ == 'boxplot':
        ax.boxplot(boxes)
    elif typ == 'frac':
        ax.bar(np.arange(0.5,len(boxes)), [b.mean() for b in boxes])
        ax.set_ylim([0, 1])
    elif typ == 'count':
        ax.bar(np.arange(0.5,len(boxes)), boxes)
        
        
    ax.set_title('clustering-' + col)
    if ax.is_last_row():
        num_clust = int(NA_wanted_lanl['V3Clusters'].max())
        ax.set_xticklabels(['CL-%i' % i for i in range(1, len(boxes)+1)])
        
plt.savefig('final_figures/clustering_pat_results.png', dpi=1000)

# <codecell>

def num_not_null(ser):
    return ser.notnull().sum()

agg_dict = {'Log-VL': num_not_null,
            'CD4':num_not_null,
            'Years-Seropositve':num_not_null,
            'CD8':num_not_null,
            'PSSMScore':num_not_null}
print NA_wanted_lanl.groupby('PSSMSubs').agg(agg_dict).sum()
print len(NA_wanted_lanl.index)

# <codecell>

keys = [('Log-VL', 'boxplot'),
        ('CD4', 'boxplot'),
        ('Years-Seropositve', 'boxplot'),
        ('CD8', 'boxplot'),
        ('PSSMScore', 'boxplot'),
        ('Count', 'count'),]

fig, axs = plt.subplots(3,2, figsize=(10,10), sharex=True)

for ax, (col,typ) in zip(axs.flatten(), keys):
    
    if col == 'Count':
        boxes = [len(vals) for _, vals in NA_wanted_lanl.groupby('PSSMSubs')]
    else:
        boxes = [vals.dropna() for _, vals in NA_wanted_lanl.groupby('PSSMSubs')[col]]
        nums = [len(vals.dropna()) for _, vals in NA_wanted_lanl.groupby('PSSMSubs')[col]]
        print col, nums
    if typ == 'boxplot':
        ax.boxplot(boxes)
    elif typ == 'frac':
        ax.bar(np.arange(0.5,len(boxes)), [b.mean() for b in boxes])
        ax.set_ylim([0, 1])
    elif typ == 'count':
        ax.bar(np.arange(0.5,len(boxes)), boxes)
        
        
    ax.set_title('binning-' + col)
    if ax.is_last_row():
        num_clust = int(NA_wanted_lanl['PSSMSubs'].max())
        ax.set_xticklabels(['CL-%i' % i for i in range(1, len(boxes)+1)])
        
plt.savefig('final_figures/binning_pat_results.png', dpi=1000)

# <markdowncell>

# If we do the clustering analysis on this sequence data we would expect to see two main clusters (X4 and R5). So when scanning with the Silhouette Score we expect to see a peak at 2 and then trail off. Then, to confirm that this is X4/R5 clustering we would expect the Adjusted Rand Index to have a positive value when comparing these two clusters to known X4/R5 tropism.

# <codecell>


# <codecell>

from scipy.stats import hypergeom

def fisher_pval():
    
    pass

def fisher_odds(val):
    pass

#order=['X4', 'X4-P', 'R5-P']+['R5-%i' %i for i in range(1,7)]
larger_wanted_mask = (pat_lanl_data['GeoRegion'] == 'North America') & pat_lanl_data['PSSMSubs'].notnull()
NA_lanl = pat_lanl_data[larger_wanted_mask]
tissue_num = NA_lanl['STissue'].value_counts()
cluster_num = NA_lanl['PSSMSubs'].value_counts()
#cluster_num.index = pd.Index(order)

counts = pd.crosstab(NA_lanl['STissue'], NA_lanl['PSSMSubs'])
#counts.columns = pd.Index(order)

pvals = []
for tissue, row in counts.iterrows():
    for cluster, val in row.to_dict().items():
        
        rv = hypergeom(larger_wanted_mask.sum(), tissue_num[tissue], cluster_num[cluster])
        pvals.append({
                      'Tissue':tissue,
                      'Cluster':cluster,
                      'Pval': rv.sf(val),
                      'Expected': rv.median(),
                      'Observed': val
                      })


pval_df = pd.DataFrame(pvals)



#(counts.T/counts.sum(axis=1)).T.to_excel('seq_counts.xlsx')

# <codecell>

pdata = pd.pivot_table(pval_df, rows='Tissue', cols='Cluster', 
                       values=['Pval', 'Observed', 'Expected'])
pdata['Pval'] = -pdata['Pval'].applymap(np.log10)

#pdata['Pval'][order].to_excel('new_seq_counts.xlsx', sheet_name='Pvals')
#pdata['Observed'][order].to_excel(writer, sheet_name='Observed')
#pdata['Expected'][order].to_excel(writer, sheet_name='Expected')

# <codecell>

tmp_data.index

# <codecell>

from PlottingTools import make_heatmap_df
order = ['R5-1', 'R5-2', 'R5-3', 'R5-4', 'R5-P', 'X4-P', 'X4']
tmp_data = pdata['Pval'].copy()
tmp_data.columns= pd.Index(order)
tmp_data[tmp_data<2] = np.nan
order = ['CSF', 'brain', 'meninges', 
         'PBMC', 'monocyte', 'plasma/serum/blood', 'resting CD4+ T cells',
         'lymph node', 'spleen',
         'GALT', 'colon',
         'BAL', 'lung', 'liver', 'cervix/vagina',
         'semen/testis', 'urethra']
fig = make_heatmap_df(tmp_data.ix[order], colormap='copper_r', figsize=(10,10))
cbar = plt.colorbar()
cbar.set_clim([0, 40])
cbar.set_ticks(range(0, 40, 5))
plt.savefig('final_figures/bining_tissue_figure_shortened.png', dpi=1000)

# <codecell>

mask = ((trim_lanl['PSSMScore'] > -2.98) & (trim_lanl['STissue'] == 'brain'))
wanted_acc = set(trim_lanl['Accession'][mask])
tdata = trim_lanl[['Accession','PSSMScore']][mask].groupby('Accession').first()
out_seqs = []
for name, seq in aa_seq_list:
    if name in wanted_acc:
        out_seqs.append({
                         'Accession':name,
                         'PSSM':tdata.ix[name]['PSSMScore'],
                         'Seq':seq
                         })
out_df = pd.DataFrame(out_seqs)
        

# <codecell>

with open('extra_brain.fasta', 'w') as handle:
    found = set()
    for name, seq in aa_seq_list:
        if name in wanted_acc:
            if seq not in found:
                fasta_writer(handle, [(name+'NEWSEQS!!!!!!', seq)])
                found.add(seq)

# <codecell>

out_df.to_csv('brain_x4.tsv', sep='\t')

# <codecell>

ax = trim_lanl.boxplot(column='PSSMScore', by = 'STissue', vert=False, figsize=(10,10))
ax.set_ylim([-1, ax.get_ylim()[1]])
order = ['R5-E', 'R5-1', 'R5-2', 'R5-3', 'R5-P', 'X4-P', 'X4']
for line, name in zip(pssm_bins, order):
    ax.annotate(name, (line-1.5, -0.5), fontsize=10)
ax.annotate('X4', (-1.5, -0.5), fontsize=10)

# <codecell>

print pssm_bins
trim_lanl = nlanl_data[['Accession', 'STissue', 'PSSMScore', 'DrugNaive-yes']]
trim_lanl['PSSMScore'] = trim_lanl['PSSMScore'].clip(-20, 0)
fig, axs = plt.subplots(1,2, figsize=(10,10), sharey=True)
tissues = trim_lanl['STissue'].unique()
for (is_naive, tl), ax in zip(trim_lanl.groupby('DrugNaive-yes'), axs.flatten()):
    
    boxes = dict([(tiss, tgroup['PSSMScore'].dropna()) for tiss, tgroup in tl.groupby('STissue')])
    ax.boxplot([boxes.get(tiss, []) for tiss in tissues], vert=False)
    #tl.boxplot(by=['STissue'], column='PSSMScore', vert=False, ax=ax)
    if is_naive:
        ax.set_title('Drug Naive')
    else:
        ax.set_title('Seen Therapy')
    ax.set_ylim([-1, ax.get_ylim()[1]])
    ax.set_yticklabels(tissues)
    ax.vlines(pssm_bins, *ax.get_ylim())
    order = ['R5-E', 'R5-1', 'R5-2', 'R5-3', 'R5-P', 'X4-P', 'X4']
    for line, name in zip(pssm_bins, order):
        ax.annotate(name, (line-1.5, -0.5), fontsize=10)
    ax.annotate('X4', (-1.5, -0.5), fontsize=10)

# <headingcell level=2>

# Other Seq Stuff

# <codecell>

NC_data = pd.read_csv('HIVSeqDB-MetaData.tsv', sep = '\t')

# <codecell>

NC_wanted = pd.DataFrame(
                         {
                          'Accession':NC_data['Genbank Locus'],
                          'HAART':NC_data['Antiretroviral Treatment'],
                          'Neuropath':NC_data['Neuropathological diagnosis'],
                          'Neurocog':NC_data['Neurocognitive diagnosis'],
                          }
                         )
NC_merge = pd.merge(NC_wanted, nlanl_data,
                    left_on = 'Accession',
                    right_on = 'Accession')

# <codecell>

NC_merge['Accession'][NC_merge['PSSMScore']>=-2.88]

# <codecell>

NC_merge.boxplot(column='PSSMScore', by='Neurocog', vert=False, figsize=(10,4))

# <codecell>

NC_merge.boxplot(column='PSSMScore', by='Neuropath', vert=False, figsize=(10,6))

# <codecell>

pd.crosstab(NC_merge['PSSMSubs'], NC_merge['Neurocog']).T

# <codecell>

orer

# <codecell>

order

# <codecell>


