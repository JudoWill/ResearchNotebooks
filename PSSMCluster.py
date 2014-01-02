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

autoreload?

# <codecell>

import autoreload
import SeqSklearn
%load_ext autoreload
%autoreload 2

# <headingcell level=2>

# Pulling out LANL Data

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

lanl_data = pd.read_csv('LANLResults.tsv', sep='\t')

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
    
    
    

treated_lanl_dict = {'Log-CD8':lanl_data['CD8 count'].map(np.log10),
                     'Log-CD4':lanl_data['CD4 count'].map(np.log10),
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

# <headingcell level=2>

# Bring in Seqs

# <codecell>

import glob
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from GeneralSeqTools import fasta_writer, seq_align_to_ref, WebPSSM_V3_fasta
from concurrent.futures import ThreadPoolExecutor
import csv
from itertools import imap
import pickle



def trans_seq(handle, wanted_seqs, trans=True):
    for name, seq in fasta_reader(handle):
        if name not in wanted_seqs:
            continue
            
        tseq = ''.join(l for l in seq if l.isalpha())
        if trans:
            rseq = Seq(tseq, generic_dna).translate()
            yield name, ''.join(l for l in rseq.tostring() if l.isalpha())
        else:
            yield name, tseq

with open('/home/will/HIVReportGen/Data/BlastDB/ConBseqs.txt') as handle:
    ref_seqs = dict()
    for line in handle:
        parts = line.strip().split('\t')
        ref_seqs[parts[0]] = parts[1]
            

files = glob.glob('/home/will/HIVTropism/LANLdata/SubB-*.fasta')
oseqs = []
if os.path.exists('aligned_seq.pkl'):
    with open('aligned_seq.pkl') as handle:
        oseqs = pickle.load(handle)
else:
    for f in files:
        fname = f.rsplit(os.sep,1)[-1].split('.')[0]
        parts = fname.split('-')
        prot_name = parts[1].replace('_','-')
        if prot_name == 'V3':
            continue
        with open(f) as handle:
        
            prot_seqs = list(trans_seq(handle, set(trop_data.index), trans=prot_name != 'LTR'))
            print prot_name, len(prot_seqs)
            aligned_seqs = seq_align_to_ref(prot_seqs, ref_seqs[prot_name], max_workers=20)
            for name, seq in aligned_seqs:
                oseqs.append({
                              'Accession':name,
                              'Seq':seq,
                              'Prot':prot_name
                              })
            

# <codecell>

aligned_seqs = pd.pivot_table(pd.DataFrame(oseqs),
                              rows='Accession',
                              cols='Prot',
                              values='Seq',
                              aggfunc='first')

# <headingcell level=3>

# Scanning Co-Linear

# <codecell>

def var_func(X, y):
    scores = np.squeeze(np.asarray(np.var(X, axis=0)))
    pvals = np.ones_like(scores)
    return scores, pvals

# <codecell>

from sklearn.pipeline import Pipeline
from sklearn.grid_search import GridSearchCV, RandomizedSearchCV
from sklearn.cluster import KMeans
from sklearn.feature_selection import SelectPercentile, SelectKBest



param_grid = {'n_clusters':range(2,30)}

region_transform = Pipeline(steps = [('SeqTransormed', SeqSklearn.BioTransformer(typ='aa')),
                                     ('VarFilter', SelectKBest(var_func, k=20))])

region_predictor = GridSearchCV(KMeans(), 
                                param_grid, n_jobs=5, pre_dispatch=10)

region_scorer = KMeans()

# <codecell>

from operator import itemgetter
import warnings
import csv

#handle = open('cluster_dump.tsv', 'w')
writer = csv.DictWriter(handle, ['Prot', 'Start', 'Accession', 'Cluster'])
writer.writeheader()
scan_width = 35
out_data = []
for col in ['Int', 'LTR', 'Nef', 'Tat-2', 'Tat-1', 'PR', 'RT', 'Vif', 'Vpr', 'gp120', 'gp41']:
    break
    seq_df, known_subs = aligned_seqs[col].dropna().align(nlanl_data['PSSMSubs'].dropna(), join='inner')
    seq_list_seqs = [list(l) for l in seq_df.values]
    seq_len = len(seq_list_seqs[0])
    for start in range(0,seq_len-scan_width):
        if start < 20:
            print col, start
        elif start > (seq_len-50):
            print col, start-seq_len
        getter = itemgetter(*range(start,start+scan_width))
        seq_region = np.array(map(getter, seq_list_seqs))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            region_transform.fit(seq_region, np.ones((seq_region.shape[0], 1)))
            seq_data = region_transform.transform(seq_region)
        region_predictor.fit(seq_data)
        clusters = region_predictor.predict(seq_data)
        for c, acc in zip(clusters, seq_df.index):
            out_data.append({
                             'Prot':col,
                             'Start':start,
                             'Accession':acc,
                             'Cluster':c
                             })
            writer.writerow({
                             'Prot':col,
                             'Start':start,
                             'Accession':acc,
                             'Cluster':c
                             })
        
    

# <codecell>

handle.close()

# <codecell>

tmp = [len(aligned_seqs[col].dropna().values[0]) for col in aligned_seqs.columns]

# <codecell>

aligned_seqs.columns

# <codecell>

#from sklearn.cross_validation import permutation_test_score, Bootstrap
##
#
#tregion = region_scorer.set_params(n_clusters=4)
#ts, scores, pval = permutation_test_score(tregion, seq_data, y=known_subs.values,
#                                          scoring = SeqSklearn.normalized_mutual_info_score_linker,                       
#                                          n_permutations=100,
#                                          cv=Bootstrap(seq_region.shape[0], train_size=0.7, n_iter=1))
#print ts, pval

# <codecell>

seq_df, known_subs = aligned_seqs['Int'].dropna().align(nlanl_data['PSSMScore'].dropna(), join='inner')

# <codecell>

from SeqSklearn import BinBasedCluster

clust = BinBasedCluster(bins=pssm_bins)
getter = itemgetter(*range(109-17,109+17))
seq_list_seqs = [list(l) for l in seq_df.values]
seq_region = np.array(map(getter, seq_list_seqs))
region_transform.fit(seq_region, np.ones((seq_region.shape[0], 1)))
seq_data = region_transform.transform(seq_region)

pca_trans, biny, xx, yy, Z = clust.make_vern_points(seq_data, known_subs.values)

# <codecell>

from pylab import get_cmap
plt.figure(figsize=(10,10))
jitter = 0.01*np.random.randn(*pca_trans.shape)+pca_trans
plt.scatter(jitter[:,0], jitter[:,1], vmax = 0, c=known_subs, cmap=get_cmap('copper_r'), alpha=0.5)
cbar = plt.colorbar()
cbar.set_label('PSSMScore')
plt.ylabel('PC-1')
plt.xlabel('PC-2')

# <codecell>

from sklearn.cross_validation import permutation_test_score, Bootstrap
from SeqSklearn import silhouette_score_linker

def check_num_clusts(X_data, n_clusts=range(2, 30), n_iter=10):
    
    cv = Bootstrap(X_data.shape[0], n_iter=n_iter, train_size=0.7)
    
    pred = GridSearchCV(KMeans(), {'n_clusters':n_clusts}, 
                        cv = cv,
                        scoring=silhouette_score_linker,
                        verbose=1,
                        refit=True)
    pred.fit(X_data)
    res = []
    for d in pred.grid_scores_:
        for r in d.cv_validation_scores:
            res.append((d.parameters['n_clusters'], r))
    return pred, pd.DataFrame(res, columns=['n_clusters','score'])


_, trop_data = check_num_clusts(seq_data, 
                                n_iter=5, 
                                n_clusts=range(2, 60))
t = trop_data.groupby('n_clusters')['score'].mean()
e = trop_data.groupby('n_clusters')['score'].std()
plt.errorbar(t.index, t.values, yerr=e.values)
#plt.title('Clustering of North American V3 sequnces')
plt.xlabel('Cluster Size')
plt.xlim([1.5, 60])
plt.ylabel('Silhouette Score')

# <codecell>

t = trop_data.groupby('n_clusters')['score'].mean()
e = trop_data.groupby('n_clusters')['score'].std()
plt.errorbar(t.index, t.values, yerr=e.values)
#plt.title('Clustering of North American V3 sequnces')
plt.xlabel('Cluster Size')
plt.xlim([1.5, 60])
plt.ylabel('Silhouette Score')

# <codecell>

pd.DataFrame({
              'val':t,
              'err':e
              }).to_csv('make_out.csv')

# <codecell>

t = pd.read_csv('make_out.csv', names = ['nc', 'err', 'val'])
plt.errorbar(t['nc'], t['val'].values, yerr=t['err'].values)
plt.xlabel('Cluster Size')
plt.xlim([1.5, 60])
plt.ylim([0,1])
plt.ylabel('Silhouette Score')
plt.savefig('final_figures/tall_int_109_clustering.png', dpi=1000)

# <codecell>


