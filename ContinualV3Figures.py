# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=1>

# Continual V3 Progress Report

# <markdowncell>

# This notebook is intended to keep track of continual results of the various V3 tropism analyses that we plan to publish. This includes (but is not limited too): differences in clinical parameters, LTR SNPs, cytokine profiles, etc. This script auto-pulls from the Google-spreadsheet that Greg/Kyle/etc are putting data into so it should be able to auto-update as they complete sequencing.

# <codecell>

from __future__ import division
import numpy as np
from scipy.stats import linregress
import matplotlib.pyplot as plt
from matplotlib import dates
import pandas as pd
import gspread
from StringIO import StringIO
import csv
import sys
sys.path.append('/home/will/PatientPicker/')
sys.path.append('/home/will/PySeqUtils/')
import LoadingTools
from GeneralSeqTools import fasta_reader, fasta_writer, seq_align_to_ref
base_path = '/home/will/Dropbox/Wigdahl HIV Lab/V3Progress/'

# <headingcell level=2>

# Pull out Google-Docs Data

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

with open('/home/will/IpythonNotebook/secret.txt') as handle:
    line = handle.next()
    login, pword = line.split('\t')
gc = gspread.login(login, pword)
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
df['Prediction'] = df['PSSM Score'].map(decide_tropism)
df['Prediction'][df['PCR'].notnull()] = df['Prediction'][df['PCR'].notnull()].fillna('Attempted')

# <codecell>

df

# <headingcell level=3>

# Generate a Linear Regression of the data

# <codecell>

from datetime import timedelta
tdf = df.groupby(['Patient', 'Visit'], as_index=False).first()

date_based = pd.pivot_table(tdf[tdf['Date'].notnull()], rows = 'Date', 
                            cols = 'Prediction',
                            values = 'Patient',
                            aggfunc = 'count')
date_cum = pd.expanding_sum(date_based)['2013-10':]

date_cum['Total'] = date_cum.sum(axis=1)
td = date_cum[['Total']].reset_index()
td['dDate'] = (td['Date'] - pd.to_datetime('2013-7-1')).apply(lambda x: x / timedelta64(1, 'D'))

m, b, r, p, e = linregress(td['dDate'], td['Total'])
num_days = (len(tdf)-b)/m
nd = pd.DataFrame({
                       'Date':pd.date_range(start = '2013-7-1', 
                                            freq = 'M',
                                            periods = np.ceil(num_days/30))
                       })
nd['dDate'] = (nd['Date'] - pd.to_datetime('2013-7-1')).apply(lambda x: x / timedelta64(1, 'D'))
nd['GuessNum'] = m*nd['dDate'] + b
nd = nd.set_index('Date')
pdata = nd['GuessNum'][nd['GuessNum']<len(tdf)]

# <codecell>

from operator import itemgetter
pssm_bins = [-15.0, -13.0, -11.0,   -9.0,    -6.96, -4.92, -2.88, 1]
pssm_names = ['R5-4', 'R5-3', 'R5-2', 'R5-1', 'R5-P', 'X4-P', 'X4-R']
tdf['PSSMCluster'] = pd.Series(np.digitize(tdf['PSSM Score'].values.astype(float), pssm_bins).astype(float),
                               index=tdf.index)
tdf['PSSMCluster'][tdf['PSSM Score'].isnull()] = np.nan
rdict = dict(enumerate(pssm_names,1))
tdf['PSSMCluster'] = tdf['PSSMCluster'].map(lambda x: rdict.get(x, np.nan))
pd.crosstab(tdf['PSSMCluster'], tdf['Prediction'])

# <codecell>

bar_cols = ['X4', 'X4-P', 'R5', 'R5-P', 'Attempted']
colors = 'rmbcg'
fig, axs = plt.subplots(2, 1, figsize = (8,8))
for ax in axs.flatten():

    bottoms = np.zeros_like(date_cum['X4'].values)
    for col, c in zip(bar_cols, list(colors)):
        ax.bar(list(date_cum.index), date_cum[col].values, color = c, width=5, bottom = bottoms)
        bottoms += date_cum[col].values
    
    if ax.is_first_row():
        ax.legend(bar_cols, 'lower right')
    
    if ax.is_first_row():
        ldate = date_cum.index[-1]+timedelta(days=60)
        tmp = pdata[pd.to_datetime('2013-9'):ldate]
        ax.plot(tmp.index, tmp.values)
    else:
        ax.plot(pdata.index, pdata.values)
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=45)
    ax.set_ylabel('Samples')
fig.tight_layout()
plt.savefig(base_path+'SeqProgress.png')

# <markdowncell>

# This figure shows the cumulative number of sequences that we've typed so far. The blue trendline shows the projected rate of completion.

# <headingcell level=2>

# Pull out Redcap Data

# <codecell>

pat_data = LoadingTools.load_redcap_data().groupby(['Patient ID', 'VisitNum']).first()
tdf.set_index(['Patient', 'Visit'], inplace=True)

# <codecell>

pat_data['Tropism'] = np.nan
_, trops = pat_data['Tropism'].align(tdf['Prediction'], join = 'left')
pat_data['Tropism'] = trops.replace('Attempted', np.nan)
pat_data['PSSMCluster'] = np.nan
_, pat_data['PSSMCluster'] = pat_data['PSSMCluster'].align(tdf['PSSMCluster'], join = 'left')

# <codecell>

def safe_time_min(inser):
    
    if inser.isnull().all():
        return pd.NaT
    else:
        return inser.dropna().min()
    
qc_cols = [('Race-Asian','min'), 
           ('Race-Indian','min'),
           ('Race-Black','min'),
           ('Race-Hawaiian','min'),
           ('Race-White','min'),
           ('Race-Multiple','min'),
           ('Race-Unknown','min'),
           ('HIV Seropositive Date', 'min'),
           ('Nadir CD4 count (cells/uL)', 'min'),
           ('Nadir CD8 count (cells/uL)', 'min'),
           ('Peak viral load (copies/mL)', 'max')]

date_cols = ['HIV Seropositive Date',
             'Date Of Visit',
             'Date of initial CD4 count',
             'Date of nadir CD4 count',
             'Date of latest CD4 count',
             'Date of initial CD8 count',
             'Date of nadir CD8 count',
             'Date of latest CD8 count',
             'Date of initial viral load',
             'Date of peak viral load',
             'Date of latest viral load']

for col in date_cols:
    pat_data[col] = pd.to_datetime(pat_data[col], coerce = True)

for col, func in qc_cols:
    pat_data[col] = pat_data[col].groupby(level = 'Patient ID').transform(func)
    

# <codecell>

type(pat_data['Date Of Visit'].values[0])

# <codecell>

def to_days(ival):
    return ival/np.timedelta64(1, 'D')

date_check_cols = [('Latest viral load', 'Date of latest viral load'),
                   ('Latest CD4 count (cells/uL)', 'Date of latest CD4 count'),
                   ('Latest CD8 count (cells/uL)', 'Date of latest CD8 count')]
cutoff = 90
for col, dcol in date_check_cols:
    date_deltas = (pat_data['Date Of Visit'] - pat_data[dcol]).map(to_days).abs()
    mask = date_deltas<cutoff
    pat_data['Close ' + col] = np.nan
    pat_data['Close ' + col][mask] = pat_data[col][mask]

# <codecell>

def to_years(ival):
    return ival/np.timedelta64(1, 'D')/365

pat_data['YearsSero'] = (pat_data['Date Of Visit'] - pat_data['HIV Seropositive Date']).apply(to_years)
log_cols = [('Latest viral load', 'Log-Latest-VL'),
            ('Peak viral load (copies/mL)', 'Log-Peak-VL'),
            ('Close Latest viral load', 'Close-Log-Latest-VL'),
            ]
for orig, new in log_cols:
    pat_data[new] = pat_data[orig].map(np.log10)

# <codecell>

pat_data['Tropism'].value_counts()

# <codecell>

from statsmodels.graphics.api import violinplot
from scipy.stats import ttest_ind, kruskal
from statsmodels.stats.power import tt_ind_solve_power

def generate_violion_plots(plot_col, group_col, group_order, ax):
    
    boxes = []
    mus = []
    stds = []
    g_order = []
    for group in group_order:
        mask = group_col == group
        tmp = plot_col[mask].dropna()
        if len(tmp) > 2:
            g_order.append(group)
            boxes.append(tmp.copy().values)
            mus.append(plot_col[mask].mean())
            stds.append(plot_col[mask].std())
        
    if len(boxes) == 2:
        ef = abs(np.diff(mus))/(np.sum(stds))
        ratio = len(boxes[1])/len(boxes[0])
        n0 = tt_ind_solve_power(effect_size=ef, alpha = alpha, power = power, ratio = ratio)
        sizes = [str(int(n0)), str(int(n0*ratio))]
        _, pval = ttest_ind(*boxes)
    else:
        sizes = ['']*len(boxes)
        _, pval = kruskal(*boxes)
    
    labels = ['%s n=%i/%s' % (t, len(b), n) for t, b, n in zip(g_order, boxes, sizes)]
    violinplot(boxes, ax = ax, labels = labels)
    return pval, ax
    
        
    

# <codecell>

checks = [('VL', 'Log-Latest-VL'),
          ('Close-VL', 'Close-Log-Latest-VL'),
          ('Peak-VL', 'Log-Peak-VL'),
          ('CD4', 'Latest CD4 count (cells/uL)'),
          ('Close-CD4', 'Close Latest CD4 count (cells/uL)'),
          ('Nadir-CD4', 'Nadir CD4 count (cells/uL)'),
          ('CD8', 'Latest CD8 count (cells/uL)'),
          ('Close-CD8', 'Close Latest CD8 count (cells/uL)'),
          ('Nadir-CD8', 'Nadir CD8 count (cells/uL)'),
          ('TMHDS', 'TMHDS'),
          ('Years-Sero', 'YearsSero')]
alpha = 0.05
power = 0.8

trops = ['X4', 'R5']
fig, axs = plt.subplots(4,3, figsize = (10,10))
for (name, col), ax in zip(checks, axs.flatten()):
    
    pval, ax = generate_violion_plots(pat_data[col], 
                                      pat_data['Tropism'], 
                                      trops, ax)
    ax.set_title(name + ' pval:%f' % pval)
    ax.set_ylim([0, ax.get_ylim()[1]])

plt.tight_layout()
plt.savefig(base_path + 'uncorrected_clinical_params.png')

# <markdowncell>

# This figure shows the difference in the X4 vs R5 population for the entire sequenced cohort. This was done at the 'sample level'. We can see a clear difference in Nadir-CD4, Close-CD4 (cd4 within 90 days of visit) and Latest CD4. However, there are numerous confounders in this analysis so I choose a subset of R5 patients such that the two groups are matched for Age, gender, race, etc.

# <codecell>

fig, axs = plt.subplots(4,3, figsize = (15,15))
for (name, col), ax in zip(checks, axs.flatten()):
    
    pval, ax = generate_violion_plots(pat_data[col], 
                                      pat_data['PSSMCluster'], 
                                      pssm_names, ax)
    ax.set_title(name + ' pval:%f' % pval)
    ax.set_ylim([0, ax.get_ylim()[1]])
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=70 )
plt.tight_layout()
plt.savefig(base_path + 'uncorrected_clinical_params_MULTI_R5.png')

# <codecell>

from sklearn.metrics.pairwise import euclidean_distances
from sklearn.preprocessing import normalize
pat_data['IsMale'] = pat_data['Gender'] == 'Male'

control_cols = ['Age', 
                'IsMale', 
                'Race-Asian',
                'Race-Indian',
                'Race-Black',
                'Race-Hawaiian',
                'Race-White',
                'Race-Multiple',
                'HAART-Naive',
                'HAART-Non-Adherent',
                'HAART-Off',
                'HAART-On',
                'YearsSero']

r5_mask = pat_data['Tropism'] == 'R5'
x4_mask = pat_data['Tropism'] == 'X4'

r5_data = pat_data[r5_mask][control_cols].dropna()
x4_data = pat_data[x4_mask][control_cols].dropna()

dists = euclidean_distances(normalize(r5_data.values.astype(float)), normalize(x4_data.values.astype(float)))

# <codecell>

def assign_best(dists):
    
    valid_dists = dists.copy()
    out = []
    for col in range(valid_dists.shape[1]):
        pos = valid_dists[:, col].argmin()
        out.append(pos)
        valid_dists[pos, :] = np.inf
    return np.array(out)
    
def match_items(small_set, large_set):
    
    small_data = normalize(small_set.values.astype(float))
    large_data = normalize(large_set.values.astype(float))
    
    dists = euclidean_distances(large_data, small_data)
    large_locs = set(assign_best(dists))
    
    
    mask = pd.Series([i in large_locs for i in range(len(large_set.index))],
                     index = large_set.index)
    
    return small_set.index, large_set[mask].index
        
x4_inds, r5_inds = match_items(x4_data, r5_data)

# <codecell>

fig, axs = plt.subplots(4,3, figsize = (10,10))
pat_data['Keep'] = False
pat_data['Keep'][x4_inds] = True
pat_data['Keep'][r5_inds] = True
    
keep_mask = pat_data['Keep']
for (name, col), ax in zip(checks, axs.flatten()):
    
    pval, ax = generate_violion_plots(pat_data[col][keep_mask], 
                                      pat_data['Tropism'][keep_mask], 
                                      trops, ax)
    ax.set_title(name + ' pval:%f' % pval)
    ax.set_ylim([0, ax.get_ylim()[1]])

plt.tight_layout()
plt.savefig(base_path + 'matched_clinical_params.png')

# <markdowncell>

# Here is the data for the matched cohorts. We still see a strong effect on the CD4.

# <codecell>

fig, axs = plt.subplots(4,3, figsize = (15,15))
for (name, col), ax in zip(checks, axs.flatten()):
    
    pval, ax = generate_violion_plots(pat_data[col][keep_mask], 
                                      pat_data['PSSMCluster'][keep_mask], 
                                      pssm_names, ax)
    ax.set_title(name + ' pval:%f' % pval)
    ax.set_ylim([0, ax.get_ylim()[1]])
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=70 )
plt.tight_layout()
plt.savefig(base_path + 'matched_clinical_params_MULTI_R5.png')

# <headingcell level=2>

# Pull out the cytokine data

# <codecell>

from sklearn.covariance import EllipticEnvelope
cytos = sorted(['IL.8','VEGF','IL.1beta',
        'G.CSF','EGF','IL.10','HGF',
        'FGF.basic','IFN.alpha','IL.6',
        'IL.12','Rantes','Eotaxin',
        'GM.CSF','MIP.1beta',
        'MCP.1','IL.5','IL.13', 'IFN.gamma','TNF.alpha',
        'IL.RA','IL.2','IL.7','IP.10',
        'IL.2R','MIG','IL.4','IL.15',
        'IL.17','MIP.1alpha']) + ['Th1', 'Th2']

cyto_data_raw = pd.read_csv('/home/will/HIVSystemsBio/NewCytokineAnalysis/CytoRawData.csv', sep = '\t')
cyto_data_raw['Th1'] = cyto_data_raw['IFN.gamma'] + \
                            cyto_data_raw['IL.2']+cyto_data_raw['TNF.alpha']
cyto_data_raw['Th2'] = cyto_data_raw['IL.4'] + \
                            cyto_data_raw['IL.5']+cyto_data_raw['IL.10']

# <codecell>

cyto_data = cyto_data_raw.groupby(['Patient ID', 'VisitNum']).mean()
tranfer_cols = ['Log-Latest-VL', 
                'Tropism',
                'Keep',
                'IsMale',
                'Race-Black',
                'Age',
                'HAART-Naive',
                'HAART-Non-Adherent',
                'HAART-Off',
                'HAART-On',
                'Hepatitis C status (HCV)']
for col in tranfer_cols:
    _, cyto_data[col] = cyto_data.align(pat_data[col], join='left', axis = 0)
cyto_data['HCV'] = cyto_data['Hepatitis C status (HCV)']

# <codecell>

for col in cytos:
    env = EllipticEnvelope(contamination=0.05)
    env.fit(cyto_data[col].dropna().values.reshape(-1, 1))
    mask = env.predict(cyto_data[col].values.reshape(-1,1))
    cyto_data[col][mask==-1] = np.nan

# <codecell>


fig, axs = plt.subplots(11,3, figsize = (10,20))

for ax, col in zip(axs.flatten(), cytos):
    
    boxes = []
    mus = []
    stds = []
    for trop in trops:
        mask = cyto_data['Tropism'] == trop
        #mask &= cyto_data['Keep']
        boxes.append(cyto_data[col][mask].dropna().values)
        mus.append(cyto_data[col][mask].mean())
        stds.append(cyto_data[col][mask].std())
    
    ef = abs(np.diff(mus))/(np.sum(stds))
    ratio = len(boxes[1])/len(boxes[0])
    n0 = tt_ind_solve_power(effect_size=ef, alpha = alpha, power = power, ratio = ratio)
    sizes = [n0, n0*ratio]
    _, pval = ttest_ind(*boxes)
    labels = ['%s n=%i/%i' % (t, len(b), n) for t, b, n in zip(trops, boxes, sizes)]
    violinplot(boxes, ax = ax, labels = labels)
    # 
    ax.set_title(col + ' pval:%f' % pval)
    ax.set_ylim([0, ax.get_ylim()[1]])
plt.tight_layout()
plt.savefig(base_path+'uncorrected_cyto_data.png')

# <codecell>

fig, axs = plt.subplots(11,3, figsize = (10,20))

for ax, col in zip(axs.flatten(), cytos):
    
    boxes = []
    mus = []
    stds = []
    for trop in trops:
        mask = cyto_data['Tropism'] == trop
        mask &= cyto_data['Keep']
        boxes.append(cyto_data[col][mask].dropna().values)
        mus.append(cyto_data[col][mask].mean())
        stds.append(cyto_data[col][mask].std())
    
    ef = abs(np.diff(mus))/(np.sum(stds))
    ratio = len(boxes[1])/len(boxes[0])
    n0 = tt_ind_solve_power(effect_size=ef, alpha = alpha, power = power, ratio = ratio)
    sizes = [n0, n0*ratio]
    _, pval = ttest_ind(*boxes)
    labels = ['%s n=%i/%i' % (t, len(b), n) for t, b, n in zip(trops, boxes, sizes)]
    violinplot(boxes, ax = ax, labels = labels)
    # 
    ax.set_title(col + ' pval:%f' % pval)
    ax.set_ylim([0, ax.get_ylim()[1]])
plt.tight_layout()
plt.savefig(base_path+'matched_cyto_data.png')

# <codecell>

import statsmodels.api as sm

con_cols = ['Log-Latest-VL', 
            'IsR5',
            'IsMale',
            'Age',
            'HAART-Naive',
            'HAART-Off',
            'HAART-On',
            'HCV']
cyto_data['IsR5'] = 1.0
cyto_data['IsR5'][cyto_data['Tropism'].isnull()] = np.nan
cyto_data['IsR5'][cyto_data['Tropism']=='X4'] = 0.0
aa_mask = cyto_data['Race-Black'] == True

fig, axs = plt.subplots(11,3, figsize = (10,20))
for ax, col in zip(axs.flatten(), cytos):
    
    tdata = cyto_data[con_cols + [col]][aa_mask].dropna()
    res = sm.GLM(tdata[col],tdata[con_cols].astype(float)).fit()
    pval = res.pvalues['IsR5']

    boxes = []
    mus = []
    stds = []
    for trop in trops:
        mask = cyto_data['Tropism'] == trop
        mask &= aa_mask
        boxes.append(cyto_data[col][mask].dropna().values)
        mus.append(cyto_data[col][mask].mean())
        stds.append(cyto_data[col][mask].std())
    
    ef = abs(np.diff(mus))/(np.sum(stds))
    ratio = len(boxes[1])/len(boxes[0])
    n0 = tt_ind_solve_power(effect_size=ef, alpha = alpha, power = power, ratio = ratio)
    sizes = [n0, n0*ratio]
    
    labels = ['%s n=%i/%i' % (t, len(b), n) for t, b, n in zip(trops, boxes, sizes)]
    violinplot(boxes, ax = ax, labels = labels)
    # 
    ax.set_title(col + ' pval:%f' % pval)
    ax.set_ylim([0, ax.get_ylim()[1]])
plt.tight_layout()
plt.savefig(base_path + 'corrected_cyto_data.png')


# <codecell>

import glob

ltr_files = sorted(glob.glob('/home/will/HIVReportGen/Data/PatientFasta/*LTR.fasta'))
ltr_seqs = []
for f in ltr_files:
    with open(f) as handle:
        ltr_seqs += list(fasta_reader(handle))
print len(ltr_seqs)

# <codecell>

conb_ltr = """TGGAAGGGCTAATTCACTCCCAACGAAGACAAGATATCCTTGATCTGTGGATCTACCACACACAA
GGCTACTTCCCTGATTAGCAGAACTACACACCAGGGCCAGGGATCAGATATCCACTGACCTTTGGATGGTGCTACAAGC
TAGTACCAGTTGAGCCAGAGAAGTTAGAAGAAGCCAACAAAGGAGAGAACACCAGCTTGTTACACCCTGTGAGCCTGCA
TGGAATGGATGACCCGGAGAGAGAAGTGTTAGAGTGGAGGTTTGACAGCCGCCTAGCATTTCATCACATGGCCCGAGAG
CTGCATCCGGAGTACTTCAAGAACTGCTGACATCGAGCTTGCTACAAGGGACTTTCCGCTGGGGACTTTCCAGGGAGGC
GTGGCCTGGGCGGGACTGGGGAGTGGCGAGCCCTCAGATCCTGCATATAAGCAGCTGCTTTTTGCCTGTACTGGGTCTC
TCTGGTTAGACCAGATCTGAGCCTGGGAGCTCTCTGGCTAACTAGGGAACCCACTGCTTAAGCCTCAATAAAGCTTGCC
TTGAGTGCTTCAAGTAGTGTGTGCCCGTCTGTTGTGTGACTCTGGTAACTAGAGATCCCTCAGACCCTTTTAGTCAGTG
TGGAAAATCTCT""".replace('\n', '')

ltr_align = list(seq_align_to_ref(ltr_seqs, conb_ltr, max_workers = 20))

# <codecell>

ltr_bin_align = []
align_inds = []
for key, seq in ltr_align:
    align_inds.append(tuple(key.split('-')))
    hit_first = False
    bins = []
    for g, cb in zip(seq, conb_ltr):
        if ~hit_first and (g == '-'):
            bins.append(np.nan)
            continue
        hit_first = True
        bins.append(1.0 if g==cb else 0.0)
    ltr_bin_align.append(np.array(bins))
ltr_bin_array = np.array(ltr_bin_align)

# <codecell>

columns = ['ConB-%03i' % i for i in range(1, len(conb_ltr)+1)]
seq_df = pd.DataFrame(ltr_bin_array, 
                      index = pd.MultiIndex.from_tuples(align_inds, 
                                                        names = ['Patient', 'Visit']), 
                      columns = columns)

# <codecell>

from scipy.stats import fisher_exact

_, seq_df['Tropism'] = seq_df.align(tdf['Prediction'], 
                                    join = 'left', 
                                    axis = 0)
r5_mask = seq_df['Tropism'] == 'R5'
x4_mask = seq_df['Tropism'] == 'X4'

# <codecell>

pvals = []
for col in columns:
    r5_col = seq_df[col][r5_mask]
    x4_col = seq_df[col][x4_mask]
    if (r5_col.isnull().sum()>5) and (x4_col.isnull().sum()>5):
        f_table = [[(r5_col == 1).sum(), (x4_col == 1).sum()],
                   [(r5_col == 0).sum(), (x4_col == 0).sum()]]
        _, pval = fisher_exact(f_table)
        pvals.append(pval)
    else:
        pvals.append(np.nan)

# <codecell>

ltr_features = pd.read_excel(base_path+'LTR_features.xlsx', 'Sheet1')
ltr_features

# <codecell>

bottoms = 3.0+(ltr_features['Column']*0.5)
heights = 0.5*np.ones_like(ltr_features['EndPos'].values)
widths = ltr_features['EndPos'] - ltr_features['StartPos']
lefts = ltr_features['StartPos']
colors= ltr_features['Color'].values

# <codecell>

from statsmodels.sandbox.stats.multicomp import multipletests
apvals = np.array(pvals)
tmask = ~np.isnan(apvals) & (apvals < 1)
reject, adj_pvals, _, _ = multipletests(apvals[tmask], 0.1, 'fdr_bh')
areject = np.zeros_like(apvals)
areject[tmask] = reject
fig, ax = plt.subplots(figsize = (10,5))
ax.plot(-np.log10(pvals))
for num, (col, pval, mc) in enumerate(zip(columns, pvals, areject.flatten())):
    if pval < (0.05):
        label = '%s p=%f' % (col, pval)
        if mc:
            label += '*'
        ax.annotate(label, (num, -np.log10(pval)))
rects = ax.bar(lefts.values, heights, width=widths.values, bottom=bottoms, color = list(colors))
for left, bottom, label in zip(lefts, bottoms, ltr_features['Name']):
    ax.annotate(label, (left+5, bottom+0.25), 
                rotation=70, 
                xytext = (left, 7.5), 
                ha = 'left',
                arrowprops = {'width':1, 'headwidth':1})
    
plt.savefig(base_path + 'snp_figure.png')

# <codecell>

from Bio.Seq import Seq
from Bio import Motif
from StringIO import StringIO
from itertools import groupby
from operator import methodcaller

def yield_motifs():
    with open('/home/will/LTRtfAnalysis/Jaspar_PWMs.txt') as handle:
        for key, lines in groupby(handle, methodcaller('startswith', '>')):
            if key:
                name = lines.next().strip().split()[-1].lower()
            else:
                tmp = ''.join(lines)
                mot = Motif.read(StringIO(tmp), 'jaspar-pfm')
                yield name, mot
                yield name+'-R', mot.reverse_complement()
pwm_dict = {}
for num, (name, mot) in enumerate(yield_motifs()):
    if num % 100 == 0:
        print num
    pwm_dict[name] = mot
    
    
tmp = u"""A 0  0 6 1 0 0 0 4 2 2 0 0 3 
C  1 1 1 0 5 6 4 1 0 0 0 3 5 5 4 0  
G  0 6 0 1 1 0 0 0 0 7 1 1 0 0 1 0 
T  6 0 0 0 1 1 3 5 7 0 0 0 0 2 2 4"""

pwm_dict['coup2'] = Motif.read(StringIO(tmp), 'jaspar-pfm')
pwm_dict['coup2-R'] = Motif.read(StringIO(tmp), 'jaspar-pfm').reverse_complement()

# <codecell>

from Bio.Alphabet import IUPAC

def score_seq(seq, mot):
    bseq = Seq(seq, alphabet=IUPAC.unambiguous_dna)
    scores = mot.scanPWM(bseq)
    for pos, score in enumerate(scores.flatten(),1):
        if ~np.isnan(score):
            tseq = seq[pos:pos+len(mot)]
            yield pos, tseq, score
    

wanted_mots = ['ap1', 'ap1-R',
               'cebpa', 'cebpa-R',
               'creb1', 'creb1-R',
               'coup2', 'coup2-R',
               'ets1','ets1-R',
               #'fev', 'fev-R',
               'foxc1',	'foxc1-R',
               #'gata2',	'gata2-R',
               #'gata3',	'gata3-R',
               #'hnf4a',	'hnf4a-R',
               #'hoxa5',	'hoxa5-R',
               'nf-kappab','nf-kappab-R',
               'nfatc2', 'nfatc2-R',
               'nr2f1','nr2f1-R',
               #'tfap2a', 'tfap2a-R',    
               #'znf354c','znf354c-R',
               'nf-kappab', 'nf-kappab-R',
               'sp1', 'sp1-R']



big_res = []
for mot_name in wanted_mots:
    print mot_name
    mot = pwm_dict[mot_name]
    thresh = Motif.Thresholds.ScoreDistribution(mot, precision = 50).threshold_fpr(0.01)
    for name, seq in ltr_align:
        pid, vnum = name.split('-')
        for pos, tf_seq, score in score_seq(seq, mot):
            big_res.append({
                            'Patient':pid,
                            'Visit':vnum,
                            'TF':mot_name,
                            'PosSeq':tf_seq,
                            'Score':score,
                            'ConBPos':pos
                            })
tf_df = pd.DataFrame(big_res)

# <codecell>

def process_vals(inser):
    return pd.rolling_max(inser, 8, min_periods=1)

gkey = ['Patient', 'Visit', 'TF']

tf_df['MaxScore'] = tf_df.groupby(gkey)['Score'].transform(process_vals)

# <codecell>

tf_pivot = pd.pivot_table(tf_df, 
                          rows = ['Patient', 'Visit'],
                          cols = ['TF', 'ConBPos'],
                          values = ['MaxScore', 'PosSeq'],
                          aggfunc = 'first')

# <codecell>

def calc_adjust(col, bind_data):
    confounder_cols = ['Age', 
                       'IsMale',
                       'Race-Black',
                       'Race-White',
                       'HAART-Naive',
                       'HAART-Non-Adherent',
                       'HAART-Off',
                       'HAART-On',
                       'YearsSero']
    condata = pat_data[sorted(set(confounder_cols+[col]))].dropna()
    condata, bondata = condata.align(bind_data, join = 'inner', axis = 0)
    sumres = sm.GLM(bondata, condata.astype(float)).fit()
    return sumres.pvalues[col]
    

pos_scores = tf_pivot['MaxScore']

results = []
for mot_name in wanted_mots:
    print mot_name
    score_data = pos_scores[mot_name]
    for pat_col_name, pat_col in checks:
        clin_data = pat_data[pat_col].dropna()
        for col in score_data.columns:
            bind_data = score_data[col].dropna().astype(float)
            if len(bind_data) < 10:
                continue
            bdata, cdata = bind_data.align(clin_data, join = 'inner')
            if (len(bdata) > 100):
                m, b, r, p, _ = linregress(bdata.values, cdata.values)
                results.append({
                                'Motif':mot_name,
                                'ClinicalVal':pat_col_name,
                                'Slope':m,
                                'Intercept':b,
                                'R2': r**2,
                                'pval':p,
                                'ConBCol':col,
                                'N':len(bdata),
                                'AdjPval':calc_adjust(pat_col, bind_data)
                                })
    

# <codecell>

cor_results = pd.DataFrame(results)

nres = pd.pivot_table(cor_results,
                     rows = ['Motif', 'ConBCol'],
                     cols = 'ClinicalVal',
                     values = ['pval', 'Slope', 'R2', 'N', 'AdjPval'],
                     aggfunc = 'first')

# <codecell>

all_pvals = np.array([d['pval'] for d in results if d['N'] > 200])
reject, adj_pvals, _, _ = multipletests(all_pvals, 0.2, 'fdr_bh')
cutoff = all_pvals[reject].max()
print cutoff
lcutoff = -np.log10(cutoff)

# <codecell>

cols = [col for _, col in checks]
ncols = pd.Index([col for col, _ in checks])
small_pat_data = pat_data[cols]
small_pat_data.columns = ncols
tf_pat_data = pd.merge(tf_df,
                       small_pat_data,
                       left_on = ['Patient', 'Visit'],
                       right_index = True, 
                       how = 'inner')

# <codecell>

from types import TupleType

crazy_big = pd.merge(tf_pat_data,
                     nres,
                     left_on = ['TF', 'ConBPos'],
                     right_index = True,
                     how = 'inner')
ncols = []
vcols = set(c for c, _ in checks)
for col in crazy_big.columns:
    if (type(col) is TupleType) and (col[1] in vcols):
        ncols.append((col[-1], col[0]))
    elif (type(col) is TupleType):
        print col
    else:
        ncols.append(('APatient', col))
        
crazy_big.columns = pd.MultiIndex.from_tuples(ncols, names = ['Value', 'Analysis'])
crazy_big.sortlevel(axis=1, inplace=True)

# <codecell>

print crazy_big.head().T.to_string()

# <codecell>

pmask = crazy_big.xs('pval', level = 'Analysis', axis=1) <= cutoff*10
pmask &= crazy_big.xs('AdjPval', level = 'Analysis', axis=1) <= cutoff*10
rmask = (crazy_big.xs('R2', level = 'Analysis', axis=1) > 0.05)
nmask = (crazy_big.xs('N', level = 'Analysis', axis=1) >= 200)
wanted_mask = (rmask & pmask & nmask).any(axis = 1)
wdata = crazy_big[wanted_mask]
wdata

# <codecell>

def coerce_cols(tup):
    try:
        return '--'.join(tup)
    except TypeError:
        return coerce_cols(tup[0])

gkey = [('APatient', 'TF'), ('APatient', 'ConBPos'), ('APatient', 'PosSeq')]
grouped_data = wdata.groupby(gkey).agg(['mean', 'count'])

# <codecell>

drop_cols = [tup for tup in grouped_data.columns if (tup[-1] == 'count') and (tup[0] != 'APatient')]
out_data = grouped_data.drop(drop_cols, axis=1)
for cname, _ in checks:
    mask = (out_data[cname]['pval']['mean'] > 0.1) & (out_data[cname]['AdjPval']['mean'] > 0.1)
    tmask = np.tile(mask.values, (len(out_data[cname].columns),1))
    out_data[cname] = out_data[cname].where(np.transpose(tmask))
out_data.reset_index(inplace=True)
out_data.columns = [coerce_cols(tup) for tup in out_data.columns]
out_data.to_excel(base_path+'TFbinding.xlsx', index=False)

# <codecell>

from itertools import product
clin_vals = sorted(cor_results['ClinicalVal'].unique())
motifs = sorted(cor_results['Motif'].unique())
pos = range(1, 630)
nindex = pd.MultiIndex.from_tuples(list(product(clin_vals, motifs, pos)), names = ['ClinicalVal', 'Motif', 'ConBCol'])

# <codecell>

cor_results.head()

# <codecell>

cor_results['LogP-signed'] = (-np.log10(cor_results['pval']))*np.sign(cor_results['Slope'])
plot_vals = cor_results.groupby(['ClinicalVal', 'Motif', 'ConBCol'])['LogP-signed'].first()

# <codecell>

plot_vals.head()

# <codecell>

def make_annotation(ax, tf, positions, color):
    if len(positions) == 1:
        label = '%s-%i' % (tf, min(positions))
    else:
        label = '%s-[%i-%i]' % (tf, min(positions), max(positions))
    ax.annotate(label,
                (max(lpos), val), 
                textcoords = 'offset points',
                xytext = (0, np.random.normal()*4),
                color = color)
    
def simple_color(val):
    if val > 0:
        return 'b'
    return 'r'
    

fig, axs = plt.subplots(4,3, figsize = (15, 10), sharex=True, sharey=True)
for ax, (col, _) in zip(axs.flatten(), checks):
    vals = -nres['pval'][col].applymap(np.log10).fillna(0)
    
    xpos = np.array([int(v.split('-')[1]) for v in vals.columns])
    yvals = vals[(vals > 1).any(axis = 1)]
    
    direct = nres['Slope'][col][(vals > 1).any(axis = 1)]
    corr = nres['R2'][col][(vals > 1).any(axis = 1)]
    
    ax.bar(left = xpos, height = direct.T.values, c = 'g', width = 1/len(direct.index))
    raise KeyboardInterrupt
    ax.set_ylim([0, 10])
    
    for tf, row in yvals.iterrows():
        lpos = []
        for c, val in sorted(row[row>lcutoff].to_dict().items()):
            if corr[c][tf] < 0.5:
                continue
            pos = int(c.split('-')[1])
            color = 'b' if direct[c][tf] > 0 else 'r'
            #print col, tf, c, val, corr[c][tf]
            if (len(lpos)==0) or (lpos[-1]+5>pos):
                lpos.append(pos)
            else:
                make_annotation(ax, tf, lpos, color)
                lpos = []
        if len(lpos)>0:
            make_annotation(ax, tf, lpos, color)
            
    
    ax.set_title(col)
    if ax.is_first_col():
        ax.set_ylabel('Slope')
    if ax.is_last_row():
        ax.set_ylabel('HXB2 Pos')
plt.tight_layout()
#plt.savefig(base_path+'tf_figure.png', dpi=200)

# <codecell>


# <codecell>


