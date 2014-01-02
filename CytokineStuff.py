# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=1>

# Cytokine Analysis

# <markdowncell>

# This analysis uses data from Luminex plates measuring the amount of particular cytokines from patient samples. These patiens were chosen such that they are either Pure Cocaine Users (PC), Pure Benzo (PB), Cocaine + other drug (MDU) or pure NonUsers (PN). These drugs of abuse were tested using blood tests.

# <headingcell level=2>

# Data Extraction

# <codecell>

from __future__ import division
import os, os.path
from pandas import *
import numpy as np
import csv
from itertools import groupby
from collections import defaultdict
from copy import deepcopy
from types import TupleType
import matplotlib.pyplot as plt
from itertools import product, imap, islice
from patsy import dmatrices
from patsy.contrasts import Treatment
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
import shutil
import glob
import shlex
from subprocess import check_call, CalledProcessError
from tempfile import mkdtemp
from concurrent.futures import ProcessPoolExecutor
from types import TupleType
from rpy2.robjects import Formula
import rpy2.robjects as robjects
import pandas.rpy.common as com
from rpy2.rinterface import RRuntimeError
import sys


sys.path.append('/home/will/PySeqUtils/')
sys.path.append('/home/will/PatientPicker/')
os.chdir('/home/will/HIVSystemsBio/NewCytokineAnalysis/')

import Rtools

SAVE_FIGURES = False

# <codecell>

import LoadingTools
redcap_data = LoadingTools.load_redcap_data()
redcap_data = redcap_data.groupby(['Patient ID', 'VisitNum']).first()

# <codecell>

def count_with_skips(inser, nskips):
    skips = 0
    for row in inser.values:
        if row:
            skips = 0
        else:
            skips += 1
        if skips > nskips:
            return False
    return True

name_mappings = {'Test-Benzodiazepine':'SBe', 
                'Test-Cannabinoid':'SCa', 
                'Test-Cocaine':'PC',
                'Test-Opiates':'PO',
                'Test-Amphetamines':np.nan,
                'Test-Barbiturates':np.nan,
                'Test-Phencyclidine':np.nan
                }

def niz_groupings(indf):
    
    inds = [v for p, v in indf.index]
    
    if len(indf.index) < 3:
        return Series(np.nan, index=inds)
    indf = indf.dropna()
    ever_used = indf.any(axis = 0)
    if not ever_used.any():
        return Series('PN', index=inds)
    all_used = indf.all(axis = 0)
    if all_used.sum() == 1:
        return Series(name_mappings[all_used.idxmax()], index=inds)
    elif all_used.sum() > 1:
        return Series('MDU', index=inds)
    
    pure_cols = []
    non_pure_cols = []
    for col in indf.columns:
        if count_with_skips(indf[col], 1):
            pure_cols.append(col)
        else:
            non_pure_cols.append(col)
    if ever_used[non_pure_cols].any():
        return Series(np.nan, index=inds)
            
    if len(pure_cols) == 1:
        return Series(name_mappings[pure_cols[0]], index=inds)
    else:
        return Series('MDU', index=inds)

admit_cols = [col for col in redcap_data.columns if col.startswith('Admit')]
admit_data = redcap_data[admit_cols].any(axis = 1).groupby(level = 'Patient ID').transform('any')
drug_names = ['Test-Benzodiazepine', 'Test-Cannabinoid', 'Test-Cocaine', 'Test-Opiates']
niz_groupings = redcap_data[drug_names].groupby(level = 'Patient ID').apply(niz_groupings)
niz_groupings[(niz_groupings == 'PN') & (admit_data)] = np.nan

# <codecell>

def safe_days(val):
    try:
        return val/np.timedelta64(1, 'D')
    except:
        return np.nan

def safe_float(val):
    try:
        return float(val)
    except:
        return np.nan

def get_days_since_bl(ser):
    
    fdate = ser.dropna().min()
    diff_dates = (ser-fdate).apply(safe_days)
    return diff_dates


def guess_race(indf):
    race_cols = [col for col in indf.columns if col.startswith('Race-')]
    race = indf[race_cols].sum().idxmax()
    vnums = [v for _, v in indf.index]
    rser = Series([race]*len(indf), index = vnums)
    return rser

def convert_haart(inval):
    tdict = {'on': 'cH',
             'non-adherent': 'dH',
             'off': 'dH',
             'naive': 'nH'}
    return tdict.get(inval, np.nan)





cols = {'Age':redcap_data.groupby(level = 0)['Age'].transform('min'),
        'NumTotalVisits': redcap_data.groupby(level = 0)['Age'].transform(len).map(safe_float),
        'Latest CD4 count': redcap_data['Latest CD4 count (cells/uL)'].map(safe_float),
        'Current Alcohol use': redcap_data['Current Alcohol Use'],
        'Current Tobacco use': redcap_data['Current Tobacco Use'],
        #'Days since baseline': redcap_data.groupby(level = 0)['Date Of Visit'].transform(get_days_since_bl),
        'Gender': redcap_data['Gender'],
        'Grouping': niz_groupings,
        'Race': redcap_data.groupby(level=0).apply(guess_race),
        'HAART': redcap_data['Current ART status'].map(convert_haart),
        'Hepatitis C status (HCV)': redcap_data.groupby(level=0)['Hepatitis C status (HCV)'].transform(expanding_max),
        'HIVD score': redcap_data['TMHDS'].map(safe_float),
        'HIVD.I': redcap_data['TMHDS']<10,
        'Hepatitis B status (HBV)': redcap_data.groupby(level=0)['Hepatitis B status (HBV)'].transform(expanding_max).map(safe_float),
        'Latest CD8 count': redcap_data['Latest CD8 count (cells/uL)'].map(safe_float),
        'Nadir CD4 count': redcap_data.groupby(level=0)['Nadir CD4 count (cells/uL)'].transform(expanding_min).map(safe_float),
        'Nadir CD8 count': redcap_data.groupby(level=0)['Nadir CD8 count (cells/uL)'].transform(expanding_min).map(safe_float),
        'Peak viral load': redcap_data.groupby(level=0)['Peak viral load (copies/mL)'].transform(expanding_max).map(safe_float),
        'Latest Viral load': redcap_data['Latest viral load'].map(safe_float),
        'Years Seropositive': redcap_data['Years Seropositive'].map(safe_float),
        'TOSample.Benzodiazepines': redcap_data.groupby(level=0)['Test-Benzodiazepine'].transform(expanding_mean).map(safe_float),
        'TOSample.Cannabinoid': redcap_data.groupby(level=0)['Test-Cannabinoid'].transform(expanding_mean).map(safe_float),
        'TOSample.Cocaine': redcap_data.groupby(level=0)['Test-Cocaine'].transform(expanding_mean).map(safe_float),
        'TOSample.Opiates': redcap_data.groupby(level=0)['Test-Opiates'].transform(expanding_mean).map(safe_float),
        'ALL.Benzodiazepines': redcap_data.groupby(level=0)['Test-Benzodiazepine'].transform('mean').map(safe_float),
        'ALL.Cannabinoid': redcap_data.groupby(level=0)['Test-Cannabinoid'].transform('mean').map(safe_float),
        'ALL.Cocaine': redcap_data.groupby(level=0)['Test-Cocaine'].transform('mean').map(safe_float),
        'ALL.Opiates': redcap_data.groupby(level=0)['Test-Opiates'].transform('mean').map(safe_float),
        'ATSample.Benzodiazepines': redcap_data['Test-Benzodiazepine'].map(safe_float),
        'ATSample.Cannabinoid': redcap_data['Test-Cannabinoid'].map(safe_float),
        'ATSample.Cocaine': redcap_data['Test-Cocaine'].map(safe_float),
        'ATSample.Opiates': redcap_data['Test-Opiates'].map(safe_float)}

known_pat_data = DataFrame(cols)

# <codecell>

cytos = sorted(['IL.8','VEGF','IL.1beta',
        'G.CSF','EGF','IL.10','HGF',
        'FGF.basic','IFN.alpha','IL.6',
        'IL.12','Rantes','Eotaxin',
        'GM.CSF','MIP.1beta',
        'MCP.1','IL.5','IL.13', 'IFN.gamma','TNF.alpha',
        'IL.RA','IL.2','IL.7','IP.10',
        'IL.2R','MIG','IL.4','IL.15',
        'IL.17','MIP.1alpha']) + ['Th1', 'Th2']

# <codecell>



raw_cyto_data = read_csv('CytoRawData.csv', 
                            sep = '\t').groupby(['Patient ID', 'VisitNum', 'SampleNum']).first()

raw_cyto_data['Th1'] = raw_cyto_data['IFN.gamma'] + \
                            raw_cyto_data['IL.2']+raw_cyto_data['TNF.alpha']
raw_cyto_data['Th2'] = raw_cyto_data['IL.4'] + \
                            raw_cyto_data['IL.5']+raw_cyto_data['IL.10']
raw_cyto_data[cytos] = raw_cyto_data[cytos].applymap(safe_float)
#known_pat_data = read_csv('CytoPatData.csv', 
#                            index_col=[0,1], 
#                            sep = '\t')

norm_cyto_data = Rtools.quantile_norm_with_R(raw_cyto_data[cytos].applymap(safe_float))





# <markdowncell>

# The quantile normalization has placed all of the data into a normal form with the same distribution across all cytokines.

# <codecell>

             
pat_cyto_data = merge(known_pat_data, norm_cyto_data.reset_index(),
                        left_index = True, 
                        right_on = ['Patient ID', 'VisitNum'],
                        how = 'inner')

pat_cyto_data = pat_cyto_data.set_index(['Patient ID', 
                                        'VisitNum', 
                                        'SampleNum'])

# <codecell>

from functools import partial

pat_cyto_data['CD4:CD8'] = pat_cyto_data['Latest CD4 count'] / \
                            pat_cyto_data['Latest CD8 count']
known_pat_data['CD4:CD8'] = known_pat_data['Latest CD4 count'] / \
                            known_pat_data['Latest CD8 count']

pat_cyto_data['Log10 Viral Load'] = pat_cyto_data['Latest Viral load'].map(np.log10)
known_pat_data['Log10 Viral Load'] = known_pat_data['Latest Viral load'].map(np.log10)

pat_cyto_data['Th1:Th2'] = pat_cyto_data['Th1']/pat_cyto_data['Th2']
raw_cyto_data['Th1:Th2'] = raw_cyto_data['Th1']/raw_cyto_data['Th2']

pat_cyto_data['HIV-Healthy'] = (pat_cyto_data['Latest Viral load'] <= 100) & \
                                (pat_cyto_data['Latest CD4 count'] >= 250)
pat_cyto_data['Drug Status'] = pat_cyto_data['Grouping']
pat_cyto_data['OlderThan50'] = pat_cyto_data['Age']>50
pat_cyto_data['OlderThan60'] = pat_cyto_data['Age']>60
pat_cyto_data['AgePast50'] = (pat_cyto_data['Age'] - 50).map(partial(max, 0))
pat_cyto_data['AgePast60'] = (pat_cyto_data['Age'] - 60).map(partial(max, 0))

# <codecell>

tmp = pat_cyto_data[['Log10 Viral Load', 'Th2']].dropna()
plt.scatter(tmp['Th2'].values, tmp['Log10 Viral Load'].values)

# <codecell>

pat_cyto_data[['Log10 Viral Load', 'Th1']].plot

# <headingcell level=2>

# Quality Control

# <headingcell level=3>

# Normalize Data

# <codecell>

pat_cyto_data['Race'].unique()

# <codecell>

from collections import defaultdict

fig, axes = plt.subplots(6,6, figsize = (10,10), sharex = True)
fig.tight_layout()
masks = [('PN', pat_cyto_data['Grouping'] == 'PN'),
         ('PC', pat_cyto_data['Grouping'] == 'PC'),
         ('SCa', pat_cyto_data['Grouping'] == 'SCa'),
         ('SBe', pat_cyto_data['Grouping'] == 'SBe'),
         ('MDU', pat_cyto_data['Grouping'] == 'MDU'),
        ('HCV+', pat_cyto_data['Hepatitis C status (HCV)']==True),
        ('HCV-', pat_cyto_data['Hepatitis C status (HCV)']==False),
        ('AA', pat_cyto_data['Race'] == 'Race-Black'),
        ('N/AA', pat_cyto_data['Race'] != 'Race-Black')]
for ax, cyto in zip(axes.flatten(), cytos):
    
    items = []
    for key, mask in masks:
        items.append(pat_cyto_data[cyto][mask].dropna())
    
    plt.sca(ax)
    plt.boxplot(items, sym = '')
    plt.xticks(range(1, len(masks)+1), [name for name, _ in masks], rotation = 90)
    plt.title(cyto)
    plt.ylabel('NormUnits')
    
if SAVE_FIGURES:
    plt.savefig('figures/QCFigures/normalized_relative_cytokine_levels.png')


# <headingcell level=3>

# Raw Cytokine Values

# <codecell>

def drop_outliers(series, rep = 0):
    if rep == 3:
        return series
    mval = series.mean()
    stdval = series.std()
    lval = mval -2*stdval
    hval = mval +2*stdval
    mask = (series < hval) & (series > lval)
    #print mask.isnull().sum()
    return drop_outliers(series[mask], rep = rep+1)
    
    
    
    
fig, axes = plt.subplots(6,6, figsize = (10,10), sharex = True)
fig.tight_layout()
for ax, cyto in zip(axes.flatten(), cytos):
    
    items = []
    for key, mask in masks:
        items.append(drop_outliers(pat_cyto_data[cyto][mask].dropna()))
    
    plt.sca(ax)
    plt.boxplot(items, sym = '')
    plt.xticks(range(1, len(masks)+1), [name for name, _ in masks], rotation = 90)
    plt.title(cyto)
    plt.ylabel('RealUnits')
    
if SAVE_FIGURES:
    plt.savefig('figures/QCFigures/raw_relative_cytokine_levels.png')

# <headingcell level=3>

# Missing Data

# <codecell>

from pylab import get_cmap
from PlottingTools import make_heatmap_df

have_values = pat_cyto_data.apply(lambda x: x.notnull())
grouped_vals = have_values.groupby(level = 'Patient ID').mean()


#plt.figure(figsize = (10,20))
#plt.imshow(grouped_vals.T.values, 
#            interpolation = 'nearest', 
#            cmap = get_cmap('gray'), 
#            aspect='auto')

fig = make_heatmap_df(grouped_vals.T, colormap = 'gray', figsize = (10,20));

plt.xticks(range(len(grouped_vals.index)),[]);
plt.clim([0,1]);
plt.xlabel('Patient');

if SAVE_FIGURES:
    plt.savefig('figures/QCFigures/missing_values.png')

# <headingcell level=2>

# Statisical Analysis

# <headingcell level=3>

# Data Extraction

# <markdowncell>

# Adding some derived data into the mix.

# <codecell>

from scipy.stats import ks_2samp, ttest_ind, mannwhitneyu
from scipy.stats import kruskal
from types import ListType
def test_categoical(cyto_series, cat_series):
    
    data = DataFrame({
                        'cyto':cyto_series,
                        'cat':cat_series
                    })
    test_data = []
    for _, group in data.groupby('cat'):
        test_data.append(group['cyto'].values)
    if len(test_data) < 2:
        return 1
    elif len(test_data) == 2:
        _, pval = mannwhitneyu(*test_data)
    else:
        _, pval = kruskal(*test_data)
    return pval
    
def test_multibinary(cyto_series, binary_df):
    pvals = []
    for col in binary_df.columns:
        pvals.append(test_categoical(cyto_series, binary_df[col]))
    return min(pvals)

def test_continous(cyto_series, value_series):
    
    tmp = DataFrame({'cyto':cyto_series,
                    'value':value_series}).dropna().reset_index(drop = True)
    
    try:
        res = ols(y=tmp['cyto'], x = tmp['value'])
        return res.f_stat['p-value']
    except:
        return None
    


confounders = [('Age', test_continous),
('Latest CD4 count', test_continous),
('Current Alcohol use', test_categoical),
('Current Tobacco use', test_categoical),
('Days since baseline', test_continous),
('Gender', test_categoical),
('Grouping', test_categoical),
('Race', test_categoical),
('HAART', test_categoical),
('Hepatitis C status (HCV)', test_categoical),
('HIVD score', test_continous),
('HIVD.I', test_categoical),
('Hepatitis B status (HBV)', test_categoical),
('Latest CD8 count', test_continous),
('Nadir CD4 count', test_continous),
('Peak viral load', test_continous),
('Nadir CD8 count', test_continous),
('Latest Viral load', test_continous),
('Years Seropositive', test_continous),
('CD4:CD8', test_continous),
('Log10 Viral Load', test_continous),
('Drug Status', test_categoical),
('HIV-Healthy', test_categoical),
]

# <codecell>

from itertools import combinations

def mask2df(cyto_data, mask_tupA, mask_tupB):
    nameA, maskA = mask_tupA
    nameB, maskB = mask_tupB
    valid = maskA | maskB
    wdata = cyto_data.ix[valid]
    yield nameA, nameB, 'All Race',wdata, maskA.ix[valid]
        
    black_mask = cyto_data['Race'] == 'Race-Black'
    valid = (maskA | maskB) & black_mask
    wdata = cyto_data.ix[valid]
    yield nameA, nameB, 'AA', wdata, maskA.ix[valid]
    

gp_names = [(gp, pat_cyto_data['Grouping']==gp) for gp in pat_cyto_data['Grouping'].dropna().unique()]
drug_groups = list(combinations(gp_names, 2)) + [(('HCV-', pat_cyto_data['Hepatitis C status (HCV)']==False), ('HCV+', pat_cyto_data['Hepatitis C status (HCV)']==True))]

#nameA, nameB, race, cyto, confounder, p-val
test_res = []

for tupA, tupB in drug_groups:
    
    for nameA, nameB, race, data, mask in mask2df(pat_cyto_data, tupA, tupB):
        for cyto in cytos:
            for conf, test in confounders:
                pval = test(data[cyto], data[conf])
                test_res.append((nameA, nameB, race, cyto, conf, pval))
association_results_df = DataFrame(test_res, columns = ['GroupA', 'GroupB', 'Race', 'Cytokine', 'Column', 'p-value'])

# <codecell>

pivoted_data = pivot_table(association_results_df, 
                            rows = ['Race', 'GroupA', 'GroupB', 'Column'], 
                            cols = 'Cytokine', 
                            values = 'p-value')
filtered_data = pivoted_data[pivoted_data<0.05]

# <markdowncell>

# They want a specific order to the heatmaps.

# <codecell>


drug_order = [
'Age', 
'Gender',
'Race', 
'Current Alcohol use', 
'Current Tobacco use', 
'Latest CD4 count', 
'Nadir CD4 count', 
'Latest CD8 count',
'Nadir CD8 count', 
'CD4:CD8', 
'Latest Viral load', 
'Peak viral load', 
'Log10 Viral Load', 
'HIV-Healthy',
'HAART', 
'Hepatitis B status (HBV)', 
'Hepatitis C status (HCV)', 
'HIVD-I', 
'HIVD score', 
'Years Seropositive', 
'Days since baseline'] 

hcv_order = [
'Age', 
'Gender', 
'Race', 
'Current Alcohol use', 
'Current Tobacco use', 
'Drug Status', 
'Latest CD4 count', 
'Nadir CD4 count', 
'Latest CD8 count', 
'Nadir CD8 count', 
'CD4:CD8', 
'Latest Viral load', 
'Peak viral load', 
'Log10 Viral Load', 
'HIV-Healthy',
'HAART', 
'Hepatitis B status (HBV)', 
'HIVD-I', 
'HIVD score', 
'Years Seropositive', 
'Days since baseline']


presentation_order = {
 ('PN', 'MDU'):drug_order,
 ('PN', 'PC'):drug_order,
 ('PN', 'SBe'):drug_order,
 ('MDU', 'SCa'):drug_order,
 ('MDU', 'PC'):drug_order,
 ('MDU', 'SBe'):drug_order,
 ('SCa', 'PC'):drug_order,
 ('SCa', 'SBe'):drug_order,
 ('PC', 'SBe'):drug_order,
    ('HCV-','HCV+'):hcv_order,
}

# <codecell>

def convert_for_R(indf, cyto, items, factor_refs):
    """Converts a dataframe into an R dataframe"""
    
    
    fcols = {
                'Hepatitis C status (HCV)':'HCV',
                'HIVD-I':'HIVDI',
                'Log10 Viral Load':'LVL',
                'TOSample-Benzodiazepines':'TOSample.Benzodiazepines',
                'TOSample-Cannabinoid':'TOSample.Cannabinoid',
                'TOSample-Cocaine':'TOSample.Cocaine',
                'TOSample-Opiates':'TOSample.Opiates',
                'ALL-Benzodiazepines':'ALL.Benzodiazepines',
                'ALL-Cannabinoid':'ALL.Cannabinoid',
                'ALL-Cocaine':'ALL.Cocaine',
                'ALL-Opiates':'ALL.Opiates'
            }

    
    wanted_cols = [col for col, _ in items] + [cyto]
    if 'Age' not in set(wanted_cols):
        wanted_cols.append('Age')
    ndf = indf.rename(columns = fcols)[wanted_cols].dropna()
    pat_index = ndf.index
    ndf = ndf.reset_index()
    #since RPY2 chokes on boolean columns
    bool_cols = ['HCV', 'HIVDI', 'OlderThan50', 'OlderThan60']
                    
    for col in bool_cols:
        try:
            ndf[col] = ndf[col].apply(str)
        except KeyError:
            pass
    #print ndf
    rpy_ndf = com.convert_to_r_dataframe(ndf)
    factor_convert_cols = []
    neqn = []
    for col, term in items:
        if len(ndf[col].dropna().unique()) > 1:
            neqn.append(term)
        if (col in factor_refs) and (factor_refs[col] in set(ndf[col])):
            factor_convert_cols.append((col, factor_refs[col]))
    
    rpy_ndf = Rtools.convert_columns_to_factors(rpy_ndf, factor_convert_cols)
            
    res_eqn = ' + '.join(neqn)
    
    formula = Formula(cyto + ' ~ ' + res_eqn)
    return rpy_ndf, formula, pat_index

# <codecell>

def convert_output_rows(tTable):
    """Converts all the different representations of drug columns into a standard setup."""
       
    
    change_rules = [('PC', 'Cocaine'),
                    ('Cocaine', 'Cocaine'),
                    ('MDU', 'MDU'),
                    ('SCa', 'Cannabinoid'),
                    ('Cannabinoid', 'Cannabinoid'),
                    ('SBe', 'Benzodiazepines'),
                    ('Benzodiazepines', 'Benzodiazepines'),
                    ('Opiates', 'Opiates'),
                    ('HCV', 'HCV'),
                    ('Age', 'Age')]
    
    change_rows = {}
    for row, (check_term, res_row) in product(tTable.index, change_rules):
        if check_term.lower() in row.lower():
            change_rows[row] = res_row
    try:
        return tTable.rename(index = change_rows)
    except:
        print tTable
        print change_rows
        tTable.rename(index = change_rows)

def RuiRegressTest(indf, items, factor_refs, cyto, debug = False):
    """Runs the Linear Mixed Effects model."""
    
    rpy_ndf, formula, pat_index = convert_for_R(indf, cyto, items, factor_refs)
    rand_formula = Formula('~1|Patient.ID')
    
    class Blackhole(object):

        def write(self, string):
            pass
        
        def flush(self, *args, **kwargs):
            pass
    
    if debug:
        result_dict = Rtools.R_linear_mixed_effects_model(rpy_ndf, 
                                                            formula, 
                                                            rand_formula)
    else:
        stdout = sys.stdout
        try:
            sys.stdout = Blackhole()
            result_dict = Rtools.R_linear_mixed_effects_model(rpy_ndf, 
                                                                formula, 
                                                                rand_formula)
        except RRuntimeError:
            return None
        finally:
            sys.stdout = stdout
       
        
    tTable = convert_output_rows(result_dict['tTable'])
    return tTable.T, pat_index

# <markdowncell>

# Subcohorts and correction equations we want to look at:

# <codecell>

rui_groups = set(['PN', 'PC', 'MDU'])
subcohorts = [  ('AA', pat_cyto_data['Race'] == 'Race-Black'),
                ('AA-Restricted', (pat_cyto_data['Race'] == 'Race-Black') & pat_cyto_data['Grouping'].notnull()),
                #('All', pat_cyto_data['Age'] > 0),
                #('HIV-Healthy', pat_cyto_data['HIV-Healthy'] == True),
                #('HIV-NONHealthy', pat_cyto_data['HIV-Healthy'] == False),
                #('HIV-LowCD4', pat_cyto_data['Latest CD4 count'] < 250),
                #('HIV-HighCD4', pat_cyto_data['Latest CD4 count'] > 250),
                #('HIV-HighVL', pat_cyto_data['Latest Viral load'] > 100),
                #('HIV-LowVL', pat_cyto_data['Latest Viral load'] <= 100),
                #('High-CD4-CD8', pat_cyto_data['CD4:CD8'] > 1),
                #('Low-CD4-CD8', pat_cyto_data['CD4:CD8'] <= 1)
                ]

test_cols = [tup for tup in combinations(pat_cyto_data['Grouping'].dropna().unique(),2)]+[('HCV-', 'HCV+')]
test_cytos = [(tup[0], tup[1], cyto) for tup, cyto in product(test_cols, cytos)]
columns = MultiIndex.from_tuples(test_cytos)

final_test_res = []
final_resid = []
drug_groups = list(pat_cyto_data['Grouping'].dropna().unique())
eqns = [
('HCVOnlyCorrection', (('Gender', 'as.factor(Gender)'), 
                    ('Race', 'as.factor(Race)'), 
                    ('HAART', 'as.factor(HAART)'), 
                    ('HCV', 'as.factor(HCV)'),
                    ('LVL', 'LVL'))),
('GroupingCorrection', (('Grouping', 'as.factor(Grouping)'),
                    ('Gender', 'as.factor(Gender)'), 
                    ('Race', 'as.factor(Race)'), 
                    ('HAART', 'as.factor(HAART)'), 
                    ('HCV', 'as.factor(HCV)'),
                    ('LVL', 'LVL'))),
('HCV_AtSample', (('ATSample.Benzodiazepines', 'ATSample.Benzodiazepines'),
                  ('ATSample.Cannabinoid', 'ATSample.Cannabinoid'),
                  ('ATSample.Cocaine', 'ATSample.Cocaine'),
                  ('ATSample.Opiates', 'ATSample.Opiates'),
                ('Gender', 'as.factor(Gender)'), 
                ('Race', 'as.factor(Race)'), 
                ('HAART', 'as.factor(HAART)'), 
                ('HCV', 'as.factor(HCV)'),
                ('LVL', 'LVL'))),
('HCV_AllSample', (('ALL.Benzodiazepines', 'ALL.Benzodiazepines'),
                  ('ALL.Cannabinoid', 'ALL.Cannabinoid'),
                  ('ALL.Cocaine', 'ALL.Cocaine'),
                  ('ALL.Opiates', 'ALL.Opiates'),
                ('Gender', 'as.factor(Gender)'), 
                ('Race', 'as.factor(Race)'), 
                ('HAART', 'as.factor(HAART)'), 
                ('HCV', 'as.factor(HCV)'),
                ('LVL', 'LVL'))),
('HCV_ToSample', (('TOSample.Benzodiazepines', 'TOSample.Benzodiazepines'),
                  ('TOSample.Cannabinoid', 'TOSample.Cannabinoid'),
                  ('TOSample.Cocaine', 'TOSample.Cocaine'),
                  ('TOSample.Opiates', 'TOSample.Opiates'),
                ('Gender', 'as.factor(Gender)'), 
                ('Race', 'as.factor(Race)'), 
                ('HAART', 'as.factor(HAART)'), 
                ('HCV', 'as.factor(HCV)'),
                ('LVL', 'LVL'))),
]

eqns += [(name+'_NVL', vals[:-1]) for name, vals in eqns]

age_factors = [('Age', 'Age')]
group_references = ['PN']
    
    

def linker_fun(tup):
    
    (eqn_name, items), (cohort_name, cohort_mask), cyto, age_factor, group_ref = tup
    data = pat_cyto_data.ix[cohort_mask]
    #print len(data.index)
    #print eqn_name
    if ('Sample' in eqn_name) and (group_ref != 'PN'):
        return None
    #print 'after', eqn_name
    factor_refs = {
                    'HCV':'False',
                    'Grouping':group_ref,
                    'Race':'Race-Black',
                    'HIVDI':'False',
                  }
    #print 'testing'
    try:
        output = RuiRegressTest(data, items+(age_factor,), factor_refs, cyto, debug = False)
    except:
        print data
        raise TypeError
    if output is None:
        #print 'wrong'
        return None
    #print 'yes'
    ttables, pat_index = output
    ttables['Cytokine'] = cyto
    ttables['Correction'] = eqn_name
    ttables['CohortName'] = cohort_name
    ttables['DrugReference'] = group_ref
    ttables['AgeFactor'] = age_factor[0]
    pat_groups = pat_cyto_data['Grouping'].ix[pat_index].groupby(level = 'Patient ID').first()
    pat_groups = pat_groups.fillna('None')
    ttables['NumSamples'] = len(pat_index)
    ttables['NumPatients'] = len(pat_groups)
    counts = pat_groups.value_counts()
    for ind, count in zip(counts.index, counts.values):
        ttables['NumPatients-'+ind] = count
    return ttables
    
test_items = product(eqns, subcohorts, cytos, age_factors, group_references)
num_to_do = len(eqns)*len(subcohorts)*len(cytos)*len(age_factors)*len(group_references)

#with ProcessPoolExecutor(max_workers = 20) as E:
    #results = E.map(linker_fun, test_items)
results = imap(linker_fun, test_items)
for num, res in enumerate(results):
    if res is None:
        #print 'wrongness'
        continue
    #print 'rightness'
    #print res.T.to_string()
    #raise KeyError
    final_test_res.append(res)
    if num % 1000 == 0:
        print num, num_to_do


        
        

# <codecell>

print final_test_res[0].T.to_string()

# <codecell>

adj_data = concat([df.reset_index() for df in final_test_res], axis = 0, ignore_index=True)

# <codecell>

#adj_data = concat(final_test_res, axis = 0).reset_index().rename(columns = {'level_1':'ResultType'})
print adj_data.head(n=6).T.to_string()

# <markdowncell>

# Now we need to aggregate the data into something more useful.

# <codecell>


test_cols = ['Opiates',
'Cocaine',
'Benzodiazepines',
'Cannabinoid',
'MDU',
'HCV',
'LVL',
'Age']
sample_cols = [
'NumPatients',
'NumPatients-MDU',
'NumPatients-None',
'NumPatients-PC',
'NumPatients-PN',
'NumPatients-SBe',
'NumPatients-SCa',
'NumSamples']

wanted_cols = test_cols + sample_cols

grouping_key = ['CohortName', 'Correction','Cytokine']
masks = [('Pval', adj_data['index'] == 'p-value'),
        ('Effect', adj_data['index'] == 'Value')]
concat_data = []
for name, mask in masks:
    ndata = adj_data.ix[mask]
    odata = {}
    for idx, group in ndata.groupby(grouping_key):
        odata[idx] = group[test_cols].mean()
    grouped_data = DataFrame(odata).T
    grouped_data.index = MultiIndex.from_tuples(grouped_data.index, names = grouping_key)
    grouped_data.columns = MultiIndex.from_tuples([(name, col) for col in grouped_data.columns])
    concat_data.append(grouped_data.copy())

num_samples = adj_data.groupby(grouping_key).agg(dict([(col, 'max') for col in sample_cols])).fillna(0)
num_samples.columns = MultiIndex.from_tuples([('NumSamples', col) for col in num_samples.columns])
concat_data.append(num_samples)
big_table = concat(concat_data, axis = 1)

# <codecell>

from functools import partial

def do_col_analysis(typ, indf):
    odf = indf.copy()
    func = do_bh if typ == 'BH' else do_bon
    #print indf
    for col in indf.columns:
        #print indf[col]
        #print indf[col].notnull().mean()
        if indf[col].notnull().mean()>0.2:
            #print 'doing'
            odf[col] = func(indf[col])
        else:
            #print 'skipping'
            odf[col] = np.nan
    return odf

def do_bh(ser):
    
    nser = ser.copy()
    nval = nser.values
    mask = ~np.isnan(nval)
    _, bh_pvals, _, _ = multipletests(nval[mask], alpha = 0.05, method = 'fdr_bh')
    nser.values[mask] = bh_pvals
    return nser

def do_bon(ser):
    
    nser = ser.copy()
    nval = nser.values
    mask = ~np.isnan(nval)
    _, bh_pvals, _, _ = multipletests(nval[mask], alpha = 0.05, method = 'bonferroni')
    nser.values[mask] = bh_pvals
    return nser

# <codecell>

bh_qvals = big_table['Pval'].groupby(level = ['CohortName', 'Correction']).apply(partial(do_col_analysis, 'BH'))
bon_qvals = big_table['Pval'].groupby(level = ['CohortName', 'Correction']).apply(partial(do_col_analysis, 'Bon'))
 
bh_qvals.columns = MultiIndex.from_tuples([('BH-Qval', col) for col in bh_qvals.columns])
bon_qvals.columns = MultiIndex.from_tuples([('Bon-Qval', col) for col in bon_qvals.columns])

final_table = concat([big_table, bh_qvals, bon_qvals], axis = 1)
final_table.to_excel('tables/ModelBasedResults.xlsx')

# <codecell>

store = HDFStore('ModelOutputs/results.hdf')
store['final_table'] = final_table
store['pat_cyto_data'] = pat_cyto_data
store['raw_cyto_data'] = raw_cyto_data
store.close()

# <codecell>

from types import StringType, ListType
def drop_data(inp):
    if inp < 0.1:
        return -np.log10(inp)
    else:
        return np.nan
    

def make_pval_heatmap(input_data, pval_col, wanted_corrections, wanted_cohorts, wanted_reference, 
                        filename = None, colormap = 'copper_r', labels = 'log'):
    
    if type(pval_col) == StringType:
        pval_data = input_data[[pval_col]]
        nwanted_corrections = list(wanted_corrections)
    elif type(pval_col) == ListType:
        pval_data = input_data[pval_col]
        nwanted_corrections = list(product(pval_col, wanted_corrections))
    else:
        raise TypeError
        
    corrected_data = pval_data.ix[wanted_corrections]
    tmp_data = corrected_data.reset_index()
    #print tmp_data
    
    results = pivot_table(tmp_data, rows = 'Cytokine', cols = 'Correction', values = pval_col)
    #print results
    #plt.figure(figsize = (10,10))
    tdata = results.applymap(drop_data)[nwanted_corrections].ix[cytos].dropna(axis = 1, how='all')
    
    try:
        fig = make_heatmap_df(tdata, colormap = colormap, 
                        figsize = (0.75*len(nwanted_corrections),10), 
                        grid_kwargs=grid_kwargs)
        
    except ValueError:
        return
    #fig.tight_layout()
    #plt.imshow(tdata.values, interpolation = 'nearest', cmap = get_cmap('copper_r'))
    #plt.xticks(range(len(tdata.columns)), tdata.columns, rotation = 90, fontsize = 16)
    #plt.yticks(range(len(tdata.index)), tdata.index, fontsize = 16)
    #print -np.log10([1, 0.1, 0.01, 0.001, 0.0001])
    cbar = plt.colorbar(ticks = -np.log10([0.1, 0.05, 0.01, 0.001, 0.0001]))
    plt.clim([0, 4])
    if labels != 'log':
        cbar.ax.set_yticklabels(['0.1', '0.05', '0.01', '0.001', '0.0001'])
        cbar.set_label('P-value')
    else:
        cbar.set_label('-log10(P-value)')
    plt.title(pval_col)
    if filename:
        plt.savefig(filename, pad_inches = 2)

        
cohorts = ['AA']
wanted = [('HCV', 'HCV'),
          ('Cocaine', 'Cocaine'),
          ('Cannabinoid', 'Cannabinoid')]
color_maps = ['copper_r']
label_types = ['norm']
grid_kwargs = {'lw':2, 'color':'w'}
new_wanted_corrections = [('HCV', ['HCVOnlyCorrection', 'HCV_AllSample']),
                          ('Cocaine', ['HCVOnlyCorrection', 
                                        'HCV_AllSample', 
                                        'HCV_ToSample',
                                        'HCV_AtSample']),
                          ('Cannabinoid', ['HCVOnlyCorrection', 
                                           'HCV_AllSample', 
                                           'HCV_ToSample',
                                           'HCV_AtSample']),
                          ]
iterable = product(cohorts, wanted, color_maps, label_types, new_wanted_corrections)
for cohort, (name, cols), cmap, label_type, (cor_name, wanted_corrections) in iterable:
    fname = 'figures/HeatMaps/%s-%s-%s-PvalFigure-fixed-AA-%s-%s.png' % (cohort, name, cor_name, cmap, label_type)
    #fname = None
    make_pval_heatmap(big_table.ix[cohort]['Pval'], cols, wanted_corrections, 
    'AA', 'PN', filename = fname, colormap = cmap, labels = label_type)
    plt.close()

# <codecell>

def fix_nums(inp):
    if inp:
        return 0.05
    else:
        return 1
    
def nfix_nums(inp):
    if (inp == 1).all():
        return 1
    else:
        return -np.log10(np.mean(inp))

lit_review = read_csv('FormattedInput/NewLitReviewTable.csv', sep = ',')


lit_review['IsSig'] = (lit_review['WasSig']==1).map(fix_nums).combine_first(lit_review['Pval'])
lit_review['IsSig'][lit_review['WasSig']==0] = 1

print lit_review.columns
tmp = [('Cocaine', 'WithCocaine'),
        ('Cannabinoid', 'WithCannabinoid'),
        ('HCV', 'WithCoinfecHCV'),
        ('HCV', 'WithCoinfecHIV')]
big_lit = lit_review
for cohort, (col, gp) in product(['AA', 'AA-Restricted'], tmp):
    niz_res = big_table.ix[cohort]['Pval'][col].reset_index().rename(columns = {col:'IsSig'})
    niz_res['Author'] = 'Niz'
    niz_res['Year'] = 2013
    niz_res['WithHIV'] = 1.0
    niz_res[gp] = 1.0
    niz_res['Grouping'] = niz_res['Correction'].map(lambda x: cohort+'-'+x)
    niz_res['Sample'] = 'Plasma'
    niz_res['Cohort'] = cohort
    big_lit = concat([niz_res, big_lit], axis = 0, ignore_index=True)

skip = set(['GroupingCorrection', 'HCV_AllSample', 'HCV_ToSample'])
drop_mask = big_lit['Correction'].map(lambda x: x in skip) & (big_lit['Cohort'] == 'AA-Restricted')

big_lit = big_lit[~drop_mask]

big_lit['Label'] = big_lit[['Author', 'Grouping', 'Sample']].apply(lambda x: ', '.join(x), axis = 1)

col_order = {
             'Cocaine':['EGF', 'Eotaxin', 'FGF.basic', 'G.CSF', 'GM.CSF', 'HGF', 'IFN.alpha', 'IFN.gamma',
                           'IL.1', 'IL.10', 'IL.12', 'IL.13', 'IL.15', 'IL.17', 'IL.18', 'IL.1RA', 'IL.1beta',
                           'IL.2', 'IL.2R', 'IL.3', 'IL.4', 'IL.5', 'IL.6', 'IL.7', 'IL.8', 'IL.RA', 'IP.10',
                           'MCP.1', 'MIG', 'MIP.1alpha', 'MIP.1beta', 'MPA', 'NAP.2', 'Rantes', 'TGF.beta',
                           'TNF.alpha', 'TNF.beta', 'VEGF', 'sCD14', 'sCD40L', 'Th1', 'Th2'],
             'Cannabinoid':['EGF', 'Eotaxin', 'FGF.basic', 'G.CSF', 'GM.CSF', 'HGF', 'IFN.alpha', 'IFN.gamma',
                               'IL.10', 'IL.12', 'IL.13', 'IL.15', 'IL.17', 'IL.18', 'IL.1RA', 'IL.1beta',
                               'IL.2', 'IL.2R', 'IL.3', 'IL.4', 'IL.5', 'IL.6', 'IL.7', 'IL.8', 'IL.RA', 'IP.10',
                               'MCP.1', 'MIG', 'MIP.1alpha', 'MIP.1beta', 'Rantes', 'TGF.beta',
                               'TNF.alpha', 'TNF.beta', 'VEGF', 'sCD14', 'Th1', 'Th2'],
             'Coinfect':['EGF', 'Eotaxin', 'FGF.basic', 'G.CSF', 'GM.CSF', 'HGF', 'IFN.alpha', 'IFN.gamma',
                               'IL.10', 'IL.12', 'IL.13', 'IL.15', 'IL.17', 'IL.18', 'IL.1beta',
                               'IL.2', 'IL.2R', 'IL.3', 'IL.4', 'IL.5', 'IL.6', 'IL.7', 'IL.8', 'IL.RA', 'IP.10',
                               'MCP.1', 'MIG', 'MIP.1alpha', 'MIP.1beta', 'Rantes', 
                               'TNF.alpha',  'VEGF', 'sCD27', 'sTNFRII', 'Th1', 'Th2'],
             }


row_order = []
found = set()
for it in big_lit['Label'].values:
    if it not in found:
        row_order.append(it)
        found.add(it)

niz_order = ['Niz, AA-GroupingCorrection, Plasma', 
            'Niz, AA-Restricted-HCV_AtSample, Plasma', 
            'Niz, AA-HCV_AtSample, Plasma', 
            'Niz, AA-HCV_ToSample, Plasma',
            'Niz, AA-HCV_AllSample, Plasma', 
            ]
        
review_figs = [('Coinfect', ((big_lit['WithCoinfecHIV']==1 )| (big_lit['WithCoinfecHCV']==1))),
               ('Cannabinoid', ((big_lit['WithHIV']==1) | (big_lit['WithCannabinoid']==1)) & (big_lit['WithCocaine'].fillna(0)==0)),
                ('Cocaine', ((big_lit['WithHIV']==1) | (big_lit['WithCocaine']==1)) & (big_lit['WithCannabinoid'].fillna(0)==0))]
for name, mask in review_figs:
    #print big_lit.ix[mask].tail().T.to_string()
    lit_table = pivot_table(big_lit.ix[mask], 
                            rows = 'Label', 
                            cols = 'Cytokine', 
                            values = 'IsSig', 
                            aggfunc=nfix_nums).ix[row_order]
    row_order = [ind for ind in lit_table.index if not ind.startswith('Niz,')]
    row_order = niz_order + row_order

    
    for col in col_order:
        if col not in lit_table:
            lit_table[col] = np.nan
    tmpa = lit_table.ix[row_order]
    tmpb = col_order[name]
    out_lit = lit_table.ix[row_order][col_order[name]].abs().dropna(axis = 0, how = 'all')
    make_heatmap_df(out_lit, figsize = (30,10), colormap = 'copper_r', grid_kwargs=grid_kwargs)
    #print lit_table[col_order].ix['Katona, Cannabinoids, MS patients who used marijuana serum samples']
    #raise KeyError
    plt.xticks(range(len(col_order[name])), col_order[name], rotation = 90, fontsize = 16);
    plt.yticks(range(len(out_lit.index)), out_lit.index, fontsize = 16);
    plt.hlines([4.5], -0.5, len(col_order[name]), color = 'r', linewidth = 5);
    last = None
    for num, val in enumerate(out_lit.index):
        if num > 4:
            nval = tuple(val.split(', ')[1:])
            if nval != last:
                plt.hlines([num-0.5], -0.5, len(col_order[name]), color = 'g', linewidth = 5)
                last = nval
            
    #plt.hlines([4.5, 8.5, 9.5, 23.5], -0.5, len(col_order), color = 'g', linewidth = 2);
    cbar = plt.colorbar(ticks = -np.log10([1, 0.1, 0.05, 0.01, 0.001, 0.0001]))
    plt.clim([0, 4])
    cbar.ax.set_yticklabels(['1', '0.1', '0.05', '0.01', '0.001', '0.0001'])
    cbar.set_label('P-value')
    #plt.gcf().tight_layout()
    
    plt.savefig('figures/LitReviewHeatmaps/%s-litreview.png' % name)
    plt.close()


# <markdowncell>

# This plot shows the cytokines which were tested in our dataset or what we found in the literature. In this figure the Copper boxes we tested but are not significant and the black boxes are significant. White boxes were not tested.

# <codecell>

def drop_outliers(series):
    mu = series.mean()
    std = series.std()
    
    mask = (series>(mu+1.5*std)) | (series<(mu-1.5*std))
    series[mask] = np.nan
    return series

def safe_med(inser):
    try:
        return np.median(inser.dropna())
    except:
        return np.nan

nraw_cyto_data = raw_cyto_data.groupby(level = ['Patient ID', 'VisitNum']).median()

raw_cyto_data_no_outliers = nraw_cyto_data.apply(drop_outliers)

# <codecell>

have_values = raw_cyto_data_no_outliers.apply(lambda x: x.notnull())
grouped_vals = have_values[cytos].groupby(level = 'Patient ID').mean()
grouped_vals.mean()

# <codecell>

from functools import partial
from scipy.stats import linregress

def regress_pvals(x,y, nreps = 30000):
    
    _, _, r, _, _ = linregress(x, y = y)
    lx = list(x)
    r = r**2
    num = 0
    for n in range(nreps):
        shuffle(lx)
        _, _, nr, _, _ = linregress(lx, y = y)
        num += nr**2 > r
    #print num, nreps, float(num)/float(nreps)
    return float(num)/float(nreps)

def check_outliers(data, nstd=5):
    
    med = np.median(data.unique())
    mad = (data-med).abs().mean()
    #print med, mad
    
    high_mask = (data > med+nstd*mad)
    
    while high_mask.any() & (data[~high_mask].max()*1.3 >= data[high_mask].min()):
        #print 'reclaiming', data[~high_mask].max(), data[high_mask].min()
        high_mask[data[high_mask].idxmin()] = False
        
    return high_mask

def prep_broken_figure(tmp_cyto_data, main_min, main_max):
    
    cout = tmp_cyto_data['out']
    main_ax = None
    out_ax = None
    #main_min = tmp_cyto_data['cyto'][~cout].min()
    #main_max = tmp_cyto_data['cyto'][~cout].max()
    if cout.any():
        out_min = tmp_cyto_data['cyto'][cout].min()
        out_max = tmp_cyto_data['cyto'][cout].max()

        if (out_min > main_max) & (main_min < out_min):
            #outliers are above
            top_ax = plt.axes([0.1, 0.75, 0.8, 0.2])
            bot_ax = plt.axes([0.1, 0.05, 0.8, 0.65])
            
            main_ax = bot_ax
            out_ax = top_ax
            
            out_ax.hold(True)
            out_ax.set_ylim([out_min*0.9, out_max*1.1])
            out_ax.set_xlim([-0.05,1.05])
            out_ax.set_yticks(out_ax.get_yticks()[::3])
            
            top_ax.spines['bottom'].set_visible(False)
            top_ax.set_xticks([])
            bot_ax.spines['top'].set_visible(False)
            bot_ax.tick_params(labeltop='off') # don't put tick labels at the top
            bot_ax.xaxis.tick_bottom()
    
    if main_ax is None:
        main_ax = plt.axes()
    
    main_ax.hold(True)
    main_ax.set_ylim([main_min, main_max])
    main_ax.set_xlim([-0.05,1.05])
    main_ax.set_yticks(main_ax.get_yticks()[::2])
    
    return main_ax, out_ax

# <codecell>

import PlottingTools
from statsmodels.graphics.boxplots import beanplot
from sklearn.covariance import EllipticEnvelope

def prep_data(name, column, x_data_col, valid_pats):
    
    valid_pats = valid_pats.groupby(level = ['Patient ID', 'VisitNum']).any().dropna()
    tcol, _ = known_pat_data[x_data_col].align(valid_pats[valid_pats])
    tcolumn = column.groupby(level = ['Patient ID', 'VisitNum']).median()
    drug_frac, cyto_data = tcol.dropna().align(tcolumn[valid_pats].dropna())
    tmp = DataFrame({'drug':drug_frac, 'cyto':cyto_data}).dropna()
    vals = tmp['cyto'].values.reshape(-1,1)
    pred = EllipticEnvelope(contamination=0.05).fit(vals)
    mask = pred.predict(vals)<0
    return tmp, mask


def violin_plots(name, column, x_data_col, valid_pats):
    wanted_cuts = np.array([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    wanted_x = np.array([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1])
    
    results = []
    tmp, mask = prep_data(name, column, x_data_col, valid_pats)
    tmp['cyto'][mask] = np.nan
    tmp = tmp.dropna()
    boxes = []
    pos = []
    labels = []
    for num, (low, hi) in enumerate(zip(wanted_cuts, list(wanted_cuts[1:])+[1.1])):
        tmask = (tmp['drug'] >= low) & (tmp['drug'] < hi)
        if tmask.sum()>2:
            boxes.append(tmp['cyto'][tmask])
            pos.append(num)
            labels.append(low)
        elif tmask.sum()==1:
            tser = Series(np.linspace(0,0.01,3))+tmp['cyto'][tmask].values[0]
            boxes.append(tser.copy())
            pos.append(num)
            labels.append(low)
    #print boxes    
    plt.figure(figsize=(8,8))
    ax = plt.subplot(111)
    beanplot(boxes, positions=pos, labels=labels, ax=ax)
    plt.ylabel(name)
    
    m, b, r, sp, _ = linregress(10*tmp['drug'].values, y = tmp['cyto'].values)
    ny = m*wanted_x*10 + b
    ax.plot(wanted_x*10, ny, color = 'g', lw = 3)


def make_imagesc_figure(name, column, x_data_col, valid_pats, 
                        filename = None, draw_annotations = True, 
                        mdu_pats = None, nreps = 100, with_vl = True):
    wanted_cuts = np.array([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    wanted_x = np.array([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1])
    
    results = []
    tmp, cout = prep_data(name, column, x_data_col, valid_pats)
    
    tmp['out'] = cout
    cyto_bins = np.linspace(tmp['cyto'][~cout].min(), tmp['cyto'][~cout].max(), 15)
    H, x, y = np.histogram2d(tmp['drug'][~cout].values, tmp['cyto'][~cout].values,
                                bins = [wanted_x, cyto_bins])
    
    m, b, r, sp, _ = linregress(tmp['drug'].values, y = tmp['cyto'].values)
    p = regress_pvals(tmp['drug'].values, tmp['cyto'].values, nreps=nreps)
    
    ny = m*wanted_x + b
    
    wanted_cmap = get_cmap('copper_r')
    wanted_cmap.set_under('w')
    fig = plt.figure(figsize = (8,8), dpi = 200)
    main_ax, out_ax = prep_broken_figure(tmp, y[0], y[-1])
        
    main_ax.pcolormesh(x-0.05,y,H.transpose(), cmap = wanted_cmap, vmax = 10, vmin=1)
       
    PlottingTools.add_grid(main_ax, x-0.05,y, color = 'w', lw = 2)
    
    for low, hi in zip(wanted_cuts, list(wanted_cuts[1:])+[1.1]):
        tmask = (tmp['drug'] >= low) & (tmp['drug'] < hi)
        med_val = tmp['cyto'][tmask].median()
        main_ax.hlines(med_val, low-0.05, hi-0.05, color = 'b', lw = 3)
    main_ax.plot(wanted_x, ny, color = 'g', lw = 3) 
    main_ax.set_ylabel(name)
    if out_ax is not None:
        out_ax.scatter(tmp['drug'], tmp['cyto'])
    return m, b, r**2, p, sp, len(tmp)

tcols = [(cyto, pat_cyto_data[cyto]) for cyto in cytos]
wanted_cols = [('Log-Viral-Load', pat_cyto_data['Log10 Viral Load']),
               ('Th1', pat_cyto_data['Th1']),
               ('Th2', pat_cyto_data['Th2']),
               ('Th1:Th2', pat_cyto_data['Th1:Th2']),
               ('CD4', pat_cyto_data['Latest CD4 count']),
               ('CD8', pat_cyto_data['Latest CD8 count']),
               ('CD4:CD8', pat_cyto_data['CD4:CD8'])
                ] + tcols

fname = 'figures/TrendingFigures/%(drug)s_%(sample)s.png'

checks = ['ALL']
drugs = ['Cocaine', 'Cannabinoid']
two_cuts = set(['Cannabinoid', 'Opiates', 'Benzodiazepines'])
with_vl = [('VL', True)]
iterable = product(drugs, wanted_cols, checks, with_vl)
hname = 'figures/TrendingFigures/heatmap_%(drug)s_%(sample)s_%(VL)s.png'
vname = 'figures/TrendingFigures/violin_%(drug)s_%(sample)s_%(VL)s.png'
trend_results = []
for drug, (name, col), check, (vlabel, vl) in iterable:
    
    tcols = set(['ALL.Benzodiazepines', 'ALL.Cannabinoid', 'ALL.Opiates', 'ALL.Cocaine'])
    tcols.discard('ALL.'+drug)
    check_col = check + '.' + drug
    mask = (pat_cyto_data['Race'] == 'Race-Black') & (~(pat_cyto_data[list(tcols)]>0).any(axis = 1))
    mdu_mask = (pat_cyto_data['Race'] == 'Race-Black') & (pat_cyto_data[list(tcols)]>0).any(axis = 1)
    violin_plots(name, col, check_col, mask)
    oname = vname % {'drug':drug, 'sample':name, 'VL':vlabel}
    plt.savefig(oname, dpi = 200)
    plt.close()
    
    tup = make_imagesc_figure(name, col, check_col, mask, nreps = 10000)
    tup += (name, drug, vlabel)
    trend_results.append(tup)
    oname = hname % {'drug':drug, 'sample':name, 'VL':vlabel}
    plt.savefig(oname, dpi = 200)
    plt.close()
    
pval_res = DataFrame(trend_results, columns = ['m', 'b', 'r^2', 'p-val', 'm-pval', 'Npats', 'Col', 'Drug', 'VL'])


    

# <codecell>


def do_bh(ser):
    
    _, bh_pvals, _, _ = multipletests(ser, alpha = 0.05, method = 'fdr_bh')
    return Series(bh_pvals, index = ser.index)

def do_bon(ser):
    _, bh_pvals, _, _ = multipletests(ser, alpha = 0.05, method = 'bonferroni')
    return Series(bh_pvals, index = ser.index)


_, bon_pvals, _, _ = multipletests(pval_res['p-val'], alpha = 0.05, method = 'bonferroni')
pval_res['BH-qvals'] = pval_res.groupby(['Drug', 'VL'])['p-val'].transform(do_bh)
pval_res['Bon-qvals'] = pval_res.groupby(['Drug', 'VL'])['p-val'].transform(do_bon)

wcols = ['Drug', 'Col', 'Npats', 'm', 'b', 'r^2', 'p-val', 'BH-qvals', 'Bon-qvals']
pval_res[wcols].to_excel('tables/TrendingPvals.xlsx', index=False)

# <codecell>


sens_cohorts = []
for num in range(1,4):
    sens_cohorts.append((str(num), pat_cyto_data['NumTotalVisits'] >= num))

drug_groups = list(pat_cyto_data['Grouping'].dropna().unique())
eqns = [
('HCV_AllSample', (('ALL.Benzodiazepines', 'ALL.Benzodiazepines'),
                  ('ALL.Cannabinoid', 'ALL.Cannabinoid'),
                  ('ALL.Cocaine', 'ALL.Cocaine'),
                  ('ALL.Opiates', 'ALL.Opiates'),
                ('Gender', 'as.factor(Gender)'), 
                ('Race', 'as.factor(Race)'), 
                ('HAART', 'as.factor(HAART)'), 
                ('HCV', 'as.factor(HCV)'))),
]

age_factors = [('AgePast50', 'AgePast50')]
    
group_references = ['PN']
    
    

def linker_fun(tup):
    
    (eqn_name, items), (cohort_name, cohort_mask), cyto, age_factor, group_ref = tup
    data = pat_cyto_data.ix[cohort_mask]
    
    #print eqn_name
    if ('Sample' in eqn_name) and (group_ref != 'PN'):
        return None
    #print 'after', eqn_name
    factor_refs = {
                    'HCV':'False',
                    'Grouping':group_ref,
                    'Race':'Race-Black',
                    'HIVDI':'False',
                  }
    #print 'testing'
    output = RuiRegressTest(data, items+(age_factor,), factor_refs, cyto, debug = False)
    if output is None:
        #print 'wrong'
        return None
    #print 'yes'
    ttables, nitems = output
    ttables['Cytokine'] = cyto
    ttables['Correction'] = eqn_name
    ttables['CohortName'] = cohort_name
    ttables['DrugReference'] = group_ref
    ttables['AgeFactor'] = age_factor[0]
    ttables['NumSamples'] = len(nitems)
    return ttables
    
test_items = product(eqns, sens_cohorts, cytos, age_factors, group_references)
num_to_do = len(eqns)*len(subcohorts)*len(cytos)*len(age_factors)*len(group_references)

sense_test_res = []

#with ProcessPoolExecutor(max_workers = 20) as E:
results = imap(linker_fun, test_items)
for num, res in enumerate(results):
    if res is None:
        #print 'wrongness'
        continue
    #print res.T.to_string()
    #raise KeyError
    sense_test_res.append(res)
    if num % 1000 == 0:
        print num, num_to_do

sens_data = concat([df.reset_index() for df in sense_test_res], axis = 0, ignore_index=True)

# <codecell>

#['Correction', 'CohortName', 'AgeFactor', 'DrugReference', 'Cytokine']
for drug in drugs + ['HCV']:
    res = pivot_table(sens_data[sens_data['index'] == 'p-value'], rows = 'CohortName', cols = 'Cytokine', values = drug)<0.1
    res = res[cytos]
    plt.figure(figsize = (20, 20))
    plt.imshow(res.values, interpolation = 'nearest', cmap = get_cmap('copper_r'))
    _ = plt.xticks(range(len(res.columns)), res.columns, rotation = 90, fontsize = 16);
    _ = plt.yticks(range(len(res.index)), res.index, fontsize = 16);
    xpos = np.arange(len(res.columns))+0.5
    ypos = np.arange(len(res.index))+0.5
    PlottingTools.add_grid(plt.gca(), xpos, ypos, color = 'w', lw = 2)
    plt.ylabel('# Required Visits')
    plt.savefig('figures/HeatMaps/%s_sensitivity_analysis.png' % drug)

# <codecell>

res = pat_cyto_data['NumTotalVisits'].groupby(level = ['Patient ID', 'VisitNum']).max()
for num in range(10):
    print num, (res >= num).sum()

# <codecell>

from scipy.stats import ttest_ind

def group_fun(series):
    if not series.any():
        return 'PN'
    elif series.sum()>1:
        return 'MDU'
    elif series['ALL.Benzodiazepines']:
        return 'PBe'
    elif series['ALL.Cannabinoid']:
        return 'PCa'
    elif series['ALL.Cocaine']:
        return 'PC'
    elif series['ALL.Opiates']:
        return 'PO'
    
    
wanted_cols = ['Latest CD4 count', 'Log10 Viral Load', 'CD4:CD8', 'Latest CD8 count']
res = pat_cyto_data[wanted_cols].groupby(level = ['Patient ID', 'VisitNum']).max()
groupings = (pat_cyto_data[['ALL.Benzodiazepines', 
                            'ALL.Cannabinoid', 
                            'ALL.Cocaine', 
                            'ALL.Opiates']]>0).groupby(level = ['Patient ID', 'VisitNum']).any()
grouping_res = groupings.apply(group_fun, axis = 1)
check_data = res.copy()
check_data['Grouping'] = grouping_res

all_groups = grouping_res.unique()
out_checks = []
for col in wanted_cols:
    for A, B in combinations(all_groups, 2):
        maskA = check_data['Grouping'] == A
        maskB = check_data['Grouping'] == B
        dataA = check_data[col].ix[maskA].dropna().values
        dataB = check_data[col].ix[maskB].dropna().values
        _, pval = ttest_ind(dataA, dataB)
        out_checks.append((col, A, B, pval))

out_df = DataFrame(out_checks, columns = ['Data', 'GroupA', 'GroupB', 'pval'])
out_df.to_excel('tables/group_pval_table.xlsx')

# <codecell>

tA = big_table['Pval'][['Cocaine', 'Cannabinoid']].ix['AA'].ix['HCV_AllSample'].ix[cytos]
fig = make_heatmap_df(-tA.applymap(np.log10), figsize = (2,10), colormap = 'copper_r', grid_kwargs=grid_kwargs)
cbar = plt.colorbar(ticks = -np.log10([1, 0.1, 0.05, 0.01, 0.001, 0.0001]))
plt.clim([0, 4])
cbar.ax.set_yticklabels(['1', '0.1', '0.05', '0.01', '0.001', '0.0001']);
plt.savefig('figures/HeatMaps/Coc_canab_cyto.png')

# <codecell>

tB = big_table['Pval'][['Cocaine', 'Cannabinoid', 'HCV']].ix['AA'].ix['HCV_AllSample'].ix[cytos]
fig = make_heatmap_df(-tB.applymap(np.log10), figsize = (2,10), colormap = 'copper_r', grid_kwargs=grid_kwargs)
cbar = plt.colorbar(ticks = -np.log10([1, 0.1, 0.05, 0.01, 0.001, 0.0001]))
plt.clim([0, 4])
cbar.ax.set_yticklabels(['1', '0.1', '0.05', '0.01', '0.001', '0.0001']);
plt.savefig('figures/HeatMaps/HCV_Coc_canab_cyto.png')

# <codecell>

wanted = [('GroupingCorrection', 'AA', 'CCM'),
          ('HCV_AtSample', 'AA-Restricted', 'Restricted-wLCM'),
          ('HCV_AtSample', 'AA', 'At-Sample'),
          ('HCV_ToSample', 'AA', 'To-Sample'),
          ('HCV_AllSample', 'AA', 'All-Sample')]
cols = ['Cocaine', 'Cannabinoid']
order = ['CCM', 'Restricted-wLCM', 'At-Sample',
         'To-Sample', 'All-Sample']

for typ, colormap in [('Effect', 'RdYlGn'), ('Pval', 'copper_r')]:
    for drug in cols:
        tdict = {}
        for check, cohort, name in wanted:
            tdict[name] =  big_table[typ][drug].ix[cohort].ix[check].ix[cytos]
        tdf = DataFrame(tdict)
        if typ == 'Pval':
            tdf = tdf.applymap(np.log10)
        else:
            tdf = tdf.applymap(np.sign)
        fig = make_heatmap_df(-tdf[order], 
                              figsize = (2,10), 
                              colormap = colormap, 
                              grid_kwargs=grid_kwargs)
        if typ == 'Pval':
            cbar = plt.colorbar(ticks = -np.log10([1, 0.1, 0.05, 0.01, 0.001, 0.0001]))
        else:
            cbar = plt.colorbar()
            
        plt.savefig('figures/HeatMaps/drexel_cyto_heatmap_%s_%s.png' % (drug, typ))
        plt.close()
    
    

# <codecell>

from subprocess import call
import shlex

cmd = 'rsync -rvh /home/will/HIVSystemsBio/NewCytokineAnalysis/ /home/will/Dropbox/Wigdahl\ HIV\ Lab/NewCytokineAnalysis/'
call(shlex.split(cmd))

# <codecell>

tres = pval_res[pval_res['p-val']<0.05]
outs = []
for key, row in tres.iterrows():
    outs.append({
                 'Drug':row['Drug'],
                 'Effect':abs(row['m']),
                 'Direction': 'Up' if row['m'] > 0 else 'Down',
                 'Analysis':'Trending',
                 'Param': row['Col']
                 })
    
checks = ['Cocaine', 'Cannabinoid', 'LVL', 'Age']
tmp = final_table.swaplevel(0, 1, axis = 1)
for c in checks:
    t = tmp[c].ix['AA'].reset_index()
    mask = (t['Pval'] < 0.05) & t['Correction'].map(lambda x: not x.endswith('_NVL'))
    for key, row in t[mask].iterrows():
        outs.append({
                     'Drug':c,
                     'Effect': abs(row['Effect']),
                     'Direction': 'Up' if row['Effect'] > 0 else 'Down',
                     'Analysis':row['Correction'],
                     'Param': row['Cytokine']
                     })
        
DataFrame(outs).to_csv('tables/unpivoted_data.csv')

# <codecell>


# <codecell>


