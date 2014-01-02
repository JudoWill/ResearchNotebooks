# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

from __future__ import division
import pandas as pd
import numpy as np
import os
import sys
import gspread
from StringIO import StringIO
import csv
sys.path.append('/home/will/PySeqUtils/')
sys.path.append('/home/will/PatientPicker/')
os.chdir('/home/will/LTRGraphs/')

# <codecell>

import LoadingTools
import GeneralSeqTools
import TFSeqTools
import ConSeqs

# <codecell>

redcap_data = LoadingTools.load_redcap_data()
redcap_data['SingleID'] = redcap_data[['Patient ID', 'VisitNum']].apply(lambda x: '-'.join(x), axis = 1)

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
    
    return 'Unk'

gc = gspread.login('judowill@gmail.com', 'gfzasfxjagrxdmqq')
spread = gc.open("V3 Project")
worksheet = spread.worksheet('PBMC Progress Report')
handle = StringIO()
writer = csv.writer(handle, delimiter = '\t')
rows = worksheet.get_all_values()
writer.writerow(rows[0])
for row in rows[1:]:
    if row[0].startswith('A'):
        try:
            writer.writerow(row)
        except UnicodeEncodeError:
            print row
handle.seek(0)

df = pd.read_csv(handle, sep = '\t', parse_dates = [5])
df['HasSeq'] = df['V3 Amino Acid Sequence'].notnull()
df['Date'] = df['Date'].map(pd.to_datetime)
df['Date'][df['Date'] == 'nan'] = np.nan
df['TropismPrediction'] = df['PSSM Score'].map(decide_tropism)
df['SingleID'] = df[['Patient', 'Visit']].apply(lambda x: '-'.join(x), axis = 1)

# <codecell>

pat_data = pd.merge(redcap_data, df,
                    left_on ='SingleID',
                    right_on = 'SingleID',
                    how = 'outer').groupby('SingleID').first()

# <codecell>

import glob
ltr_files = sorted(glob.glob('/home/will/HIVReportGen/Data/PatientFasta/*LTR.fasta'))
ltr_seqs = {}
for f in ltr_files:
    with open(f) as handle:
        _, seq = GeneralSeqTools.fasta_reader(handle).next()
        fname = os.path.basename(f).rsplit('-', 1)[0]
        ltr_seqs[fname] = seq

# <codecell>

ltr_df = pd.DataFrame({
                       'LTR':pd.Series(ltr_seqs)
                       })
ltr_df.head()

# <codecell>

conb_ltr = ConSeqs.GetConSeq('ltr')
conb_ltr

# <codecell>

from functools import partial
def scan_seq(mot, name, seq):
    if len(seq) < len(mot):
        return pd.Series([np.nan, np.nan], index = [name+'-Score', name+'-Seq'])
    score, _, seq = TFSeqTools.simple_score_pwm(seq, mot)
    return pd.Series([score, seq], index = [name+'-Score', name+'-Seq'])

def region_extractor(conb_ltr, start, stop, seq):
    oseq = TFSeqTools.slice_to_ref(seq, conb_ltr, start, stop)
    nseq = oseq.replace('-', '')
    if len(nseq):
        return nseq
    else:
        return np.nan

pwm_dict = TFSeqTools.Load_PWMS()
             
regions = [('AP1-IV', 104-5, 104+15, pwm_dict['ap1']),
           ('AP1-III', 119-5, 119+15, pwm_dict['ap1']),
           ('AP1-II', 154-5, 154+15, pwm_dict['ap1']),
           ('AP1-I', 213-5, 213+15, pwm_dict['ap1']),
           ('CEBP-II', 280-5, 280+15, pwm_dict['cebpa']),
           ('ETS', 304-5, 304+15, pwm_dict['ets1']),
           ('ATF/Creb', 329-5, 329+15, pwm_dict['creb1']),
           ('CEBP-I', 329-5, 329+15, pwm_dict['cebpa']),
           ('NFKB-II', 349-5, 349+15, pwm_dict['nf-kappab']),
           ('NFKB-I', 362-5, 362+15, pwm_dict['nf-kappab']),
           ('SpIII', 376-5, 376+15, pwm_dict['sp1']),
           ('SPII', 387-5, 387+15, pwm_dict['sp1']),
           ('SPI', 398-5, 398+15, pwm_dict['sp1'])]

ltr_df_cp = ltr_df.copy()
for tfname, start, stop, mot in regions:
    r_ex = partial(region_extractor, conb_ltr, start, stop)
    tf_checker = partial(scan_seq, mot, tfname)
    ltr_df_cp[tfname+'-Region'] = ltr_df_cp['LTR'].map(r_ex)
    out = ltr_df_cp[tfname+'-Region'].dropna().apply(tf_checker)
    ltr_df_cp = pd.merge(ltr_df_cp, out,
                         left_index=True,
                         right_index=True,
                         how='left')
    

# <codecell>

cols = dict([('Patient ID', 'PatientID'),
        ('VisitNum', 'VisitNum'),
        ('Age', 'Age'),
        ('Gender', 'Gender'),
        ('Latest viral load', 'VL'),
        ('Latest CD4 count (cells/uL)', 'CD4'),
        ('Latest CD8 count (cells/uL)', 'CD8'),
        ('Current ART status', 'ART'),
        ('Hepatitis C status (HCV)', 'HCV'),
        ('TropismPrediction', 'Tropism')
        ])

tfnames = [tf for tf, _, _, _ in regions]
score_cols = [tf+'-Score' for tf, _, _, _ in regions]
seq_cols = [tf+'-Seq' for tf, _, _, _ in regions]

wanted_pat = pat_data[cols.keys()].dropna()
wanted_scores = ltr_df_cp[score_cols+seq_cols].dropna()
wanted_scores['TFJoin'] = wanted_scores[seq_cols].apply(lambda x: ''.join(x), axis=1)
#wanted_scores = wanted_scores.drop(seq_cols, axis=1)

check_data = pd.concat(wanted_pat.align(wanted_scores, axis=0, join='inner'), axis=1).rename(columns = cols)
check_data = check_data.fillna(check_data[score_cols].min())

ncols = dict((col, col.replace('-', '_').replace('/', '_')) for col in check_data.columns)
check_data = check_data.rename(columns = ncols)

# <codecell>

import TreeingTools

tree = TreeingTools.run_FastTree(check_data['TFJoin'].to_dict().items(),
                                 alphabet=TreeingTools.generic_dna)

# <codecell>

import networkx as nx
from itertools import combinations
import csv
with open('ltr_tree.nwk', 'w') as handle:
    tree.write_to_stream(handle, schema = 'phylip', exclude_chars=True)

# <codecell>



# <codecell>

def get_node_names(node):
    if node.taxon is None:
        return node.oid
    else:
        return node.taxon.label
        
names = [get_node_names(node) for node in tree.nodes()]

tree_graph = nx.DiGraph()
tree_graph.add_nodes_from(names)

for parent_tree_node in tree.nodes():
    parent_graph_name = get_node_names(parent_tree_node)
    for child_node in parent_tree_node.child_nodes():
        edge_len = child_node.edge_length
        child_graph_name = get_node_names(child_node)
        tree_graph.add_edge(parent_graph_name, child_graph_name, weight=edge_len)

# <codecell>

check_data['LogVL'] = check_data['VL'].map(np.log10)
check_data.to_csv('ltr_node_attr.csv')
nx.write_gml(tree_graph, 'tree_graph.gml')

# <codecell>

check_data.std()

# <codecell>

check_data['AP1_I_Score'].describe()

# <codecell>



nscore_cols = [col for col in check_data.columns if 'Score' in col]
pat_graph = nx.Graph()
with open('node_attributes.csv', 'w') as handle:
    writer = csv.DictWriter(handle, check_data.columns)
    for key, row in check_data.iterrows():
        pat_graph.add_node(key, **row.dropna().to_dict())
        writer.writerow(row.dropna().to_dict())
    
with open('edge_table.csv', 'w') as handle:
    writer = csv.DictWriter(handle, ['Source', 'Target']+nscore_cols)
    for (key1, row1), (key2, row2) in combinations(check_data.iterrows(), 2):
    
        tdists = (row1[nscore_cols]-row2[nscore_cols])
        tdict = tdists.dropna().to_dict()
        tdict['Source'], tdict['Target'] = (key1, key2)
        writer.writerow(tdict)
        mean_dists = tdists.mean()
        if mean_dists == 0:
            pat_graph.add_edge(key1, key2)


# <codecell>

nx.write_gml(pat_graph, 'ltr_tf_graph.gml')

# <codecell>

dists = []
for (key1, row1), (key2, row2) in combinations(check_data.iterrows(), 2):
    
    dists.append(row1[score_cols]-row2[score_cols])
dists_df = pd.DataFrame(dists)
    

# <codecell>

print diffs_df.std()

# <codecell>

from scipy.stats import gaussian_kde
#fig, ax = plt.subplots(1,1, figsize=(20,5))
vals = diffs_df.fillna(diffs_df.min()).mean(axis=1)
#kde = gaussian_kde(vals)
rank_vals = vals.rank()



#ax.set_xticks(np.arange(0, 10, 0.5))

# <codecell>

print vals[rank_vals.idxmax()], vals[rank_vals.idxmin()]

# <codecell>

vals

