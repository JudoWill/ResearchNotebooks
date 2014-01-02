# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=1>

# Cannabinoid Effect on HIVD score

# <markdowncell>

# Here I am looking to find the effect of Cannabinoid Use on HIVD score. HIVD score is a value that ranges from 0-12 with anything below 10 marking 'impaired' and >=10 implying not impaired. This test is done at every patient visit along with blood-tests confirming drug use.

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
os.chdir('/home/will/HIVSystemsBio/NewCytokineAnalysis/')

import Rtools
import PlottingTools

# <headingcell level=2>

# Data Description

# <markdowncell>

# I took data from the 1/16/2013 (most recent) cohort data dump. I extracted drug test resuts and HIVD scores along with Date of Visit and Visit Number.

# <codecell>

store = HDFStore('/home/will/HIVReportGen/Data/SplitRedcap/2013-01-16/EntireCohort.hdf')
redcap_data = store['redcap']

t = redcap_data['Event Name'].dropna().apply(lambda x: x.split(' - ')[0])
redcap_data['Patient visit number'] = redcap_data['Patient visit number'].combine_first(t)
redcap_data["Drugs used (choice='Other')"] = redcap_data["Drugs used (choice='Other')"].dropna() == 'Checked'

drug_cols = ['Amphetamines',
             'Barbiturates',
            'Benzodiazepines',
            'Cannabinoid',
            'Cocaine + metabolite',
            'Opiates',
            'Phencyclidine'
            ]
wanted_drug_names = ['Benzodiazepines', 'Cannabinoid', 'Cocaine', 'Opiates']
all_drug_names = ['Benzodiazepines', 'Cannabinoid', 'Cocaine', 'Opiates','Amphetamines','Barbiturates','Phencyclidine']
race_cols = ["Race (choice='Asian')",
 "Race (choice='American Indian/Alaska Native')",
 "Race (choice='Black or African American')",
 "Race (choice='Native Hawaiian or other Pacific Islander')",
 "Race (choice='White')",
 "Race (choice='More than one race')",
 "Race (choice='Unknown')"]

admit_cols = ["Drugs used (choice='Marijuana')",
 "Drugs used (choice='Cocaine (crack, nasal, smoke, inject)')",
 "Drugs used (choice='Heroin (nasal, inject)')",
 "Drugs used (choice='Methamphetamine (smoke, nasal, inject)')",
 "Drugs used (choice='Benzodiazapine (i.e. valium, ativan, xanax, klonipin, etc)')",
 "Drugs used (choice='Narcotics')",
 "Drugs used (choice='Ecstasy')",
 "Drugs used (choice='PCP')",
 "Drugs used (choice='Ritalin')",
 "Drugs used (choice='Other')"]


wanted_cols = ['Patient ID', 'Patient visit number', 'Year of Birth','HIV seropositive date','Current ART status',
                'Date of visit', 'Total Modified Hopkins Dementia Score']+drug_cols+race_cols+admit_cols
wanted_redcap = redcap_data[wanted_cols].dropna(how = 'all')
#wanted_redcap["Drugs used (choice='Other')"] = wanted_redcap["Drugs used (choice='Other')"] == 'Checked'
data = wanted_redcap.rename(columns= {
                                'Patient visit number':'VisitNum',
                                'Date of visit':'Date',
                                'Total Modified Hopkins Dementia Score':'nHIVD',
                                'Cocaine + metabolite': 'Cocaine',
                                "Drug Use and HIV Status":'DrugStatus',
                                'Current ART status':'ART'
                            })
data.sort(['Patient ID', 'VisitNum'], inplace=True)
data = data.set_index(['Patient ID', 'VisitNum'])

get_years = lambda x: x/np.timedelta64(1,'D')/365
drug_names = ['Cocaine', 'Cannabinoid', 'Opiates', 'Benzodiazepines']
data['DrugFree'] = (data[drug_names].dropna() == 0).all(axis = 1)
data['bDate'] = data['Year of Birth'].dropna().map(lambda x: datetime(int(x), 1,1))
data['Age'] = data[['Date', 'bDate']].dropna().apply(lambda x: x['Date'] - x['bDate'], axis = 1).map(get_years)
data['TimeSeropositive'] = data[['Date', 'HIV seropositive date']].dropna().apply(lambda x: x['Date'] - x['HIV seropositive date'], axis = 1).map(get_years)

# <codecell>

race_dict = {"Race (choice='Asian')":'Asian',
             "Race (choice='American Indian/Alaska Native')":'Indian',
             "Race (choice='Black or African American')":'Black',
             "Race (choice='Native Hawaiian or other Pacific Islander')":'Hawaiian',
             "Race (choice='White')":'White',
             "Race (choice='More than one race')":'Multi-Race',
             "Race (choice='Unknown')":'Unknown'}

def fix_race(indf):
    tmp = indf.mean().idxmax()
    return race_dict[tmp]

race_data = data[race_cols].groupby(level = 'Patient ID').apply(fix_race)
data['Race'] = [race_data[pat] for pat, _ in data.index]

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

name_mappings = {'Benzodiazepines':'SBe', 
                'Cannabinoid':'SCa', 
                'Cocaine':'PC',
                'Opiates':'PO',
                'Amphetamines':np.nan,
                'Barbiturates':np.nan,
                'Phencyclidine':np.nan
                }

def niz_groupings(indf):
    
    if len(indf.index) < 3:
        return np.nan
    indf = indf.dropna()
    ever_used = indf.any(axis = 0)
    if not ever_used.any():
        return 'PN'
    all_used = indf.all(axis = 0)
    if all_used.sum() == 1:
        return name_mappings[all_used.idxmax()]
    elif all_used.sum() > 1:
        return 'MDU'
    
    pure_cols = []
    non_pure_cols = []
    for col in indf.columns:
        if count_with_skips(indf[col], 1):
            pure_cols.append(col)
        else:
            non_pure_cols.append(col)
    if ever_used[non_pure_cols].any():
        return np.nan
            
    if len(pure_cols) == 1:
        return name_mappings[pure_cols[0]]
    else:
        return 'MDU'
    
admit_data = data[admit_cols].any(axis = 1).groupby(level = 'Patient ID').any()
drug_names = ['Benzodiazepines', 'Cannabinoid', 'Cocaine', 'Opiates']
niz_groupings = data[all_drug_names].groupby(level = 'Patient ID').apply(niz_groupings)
niz_groupings[(niz_groupings == 'PN') & (admit_data)] = np.nan
data['NizGroupings'] = [niz_groupings[pat] for pat, _ in data.index]

# <codecell>

sample_counts = data['NizGroupings'].value_counts()
patient_counts = data.groupby(level = 'Patient ID')['NizGroupings'].first().value_counts()
print DataFrame({'SampleCounts':sample_counts, 'PatientCounts':patient_counts})

# <codecell>

data = data.drop(['Year of Birth', 'bDate', 'HIV seropositive date']+race_cols+admit_cols, axis = 1)

print data

# <headingcell level=2>

# HIVD Scores

# <codecell>

fig, axes = plt.subplots(2,2, sharey = True, sharex = True, figsize = (10,10))
plt.sca(axes.flatten()[0])
data.groupby(level = 'Patient ID')['nHIVD'].first().hist(bins = range(13))
plt.title('Entry HIVD scores');

plt.sca(axes.flatten()[1])
data.groupby(level = 'Patient ID')['nHIVD'].mean().hist(bins = range(13))
plt.title('Mean HIVD scores');

plt.sca(axes.flatten()[2])
data.groupby(level = 'Patient ID')['nHIVD'].max().hist(bins = range(13))
plt.title('Max HIVD scores');

plt.sca(axes.flatten()[3])
data.groupby(level = 'Patient ID')['nHIVD'].min().hist(bins = range(13))
plt.title('Min HIVD scores');

# <markdowncell>

# We can see from these that all of the different ways of summarizing HIVD scores gives drastically different answers. If we look at the Max score then almost everyone is 'Healthy'. Looking at the Min score we can see that about half of the patients were impaired at some point. I wonder what would happen if we 'smoothed' them.

# <codecell>

def cum_dev(inser):
    #print inser
    #print np.sum(np.abs(np.diff(inser.dropna().values)))
    #raise KeyError
    return np.sum(np.abs(np.diff(inser.dropna().values)))
    
abs_dev = data.groupby(level = 'Patient ID').agg({'nHIVD':cum_dev})
abs_dev = abs_dev.sort('nHIVD')
high_dev_pats = abs_dev.tail(n=4).index

fig, axs = plt.subplots(4,1, sharey=True, figsize = (10,15))
for pat, ax in zip(high_dev_pats, axs.flatten()):
    
    vals = data.ix[pat][['Date', 'nHIVD']].set_index('Date')
    sHIVD = rolling_mean(vals['nHIVD'], 3, min_periods=1, freq = '12M', axis = 0)
    start_date = max(vals.index[0], sHIVD.index[0])
    end_date = min(vals.index[-1], sHIVD.index[-1])
    sHIVD = sHIVD.truncate(before=start_date, after=end_date)
    plt.sca(ax)
    plt.hold(True)
    plt.plot(vals.index, vals['nHIVD'].values, lw = 4, alpha = 0.6, marker = '*')

    plt.plot(sHIVD.index, sHIVD.values, lw = 3, color = 'r', marker = 'o')
    plt.hold(False)
    plt.ylabel('HIVD')
    plt.title(pat)
plt.ylim([0, 12])
plt.savefig('figures/NeuroSpecificFigures/smoothed_HIVD.png', dpi = 200)

# <codecell>

def normalize_by_date(indf):
    tdf = indf.reset_index().set_index('Date').drop(['Patient ID', 'VisitNum'], axis = 1)
    tdf[drug_names + ['DrugFree']] = tdf[drug_names + ['DrugFree']].applymap(float)
    tdata = rolling_mean(tdf[['nHIVD','DrugFree']+drug_names], 3, min_periods=1, freq = '12M', axis = 0)#.reset_index()
    odf = tdf.resample('12M', how = 'last')
    odf[['DrugFree']+drug_names] = tdata[['DrugFree']+drug_names]
    odf['sHIVD'] = tdata['nHIVD']
    try:
        odf['iHIVD'] = tdf['nHIVD'].dropna().values[0]
    except IndexError:
        odf['iHIVD'] = np.nan
    odf['LsHIVD'] = odf['sHIVD'].shift(1)
    odf['dHIVD'] = rolling_apply(odf['sHIVD'], 2, np.diff)
    #odf['ART'] = tdf['ART'].resample('12M', how = 'last')
    odf = odf.reset_index()
    odf.index.name = 'Year'
    return odf

date_data = data.groupby(level = 'Patient ID').apply(normalize_by_date).reset_index()
print date_data

# <headingcell level=2>

# Cytokine Related Dataset

# <markdowncell>

# First, I'll look at the subset of the cohort that was in the previous cytokine analysis. For consistency I'll use the Drug-Status groupings from the previous analysis. We'll branch out into a more expansive dataset and a more 

# <codecell>

norm_cyto_data = read_csv('CytoRawDataNorm.csv', 
                            sep = '\t')
known_pat_data = read_csv('CytoPatData.csv', 
                            sep = '\t')
#print known_pat_data

wanted_pats = set(known_pat_data['Patient ID'])
groupings = dict(zip(known_pat_data['Patient ID'].values, known_pat_data['Grouping'].values))
cyto_dataset = date_data.ix[date_data['Patient ID'].map(lambda x: x in wanted_pats)]
cyto_dataset['DrugStatus'] = cyto_dataset['Patient ID'].map(lambda x: groupings[x])

print cyto_dataset

# <codecell>

cyto_dataset['DrugStatus'].value_counts()

# <codecell>

drug_groups = ['PN', 'SCa', 'PC', 'SBe','MDU']
tdata = []
for group in drug_groups:
    mask = cyto_dataset['DrugStatus'] == group
    neuro_data = cyto_dataset['sHIVD'][mask]
    tdata.append(neuro_data.dropna().values)

plt.boxplot(tdata);
plt.xticks(range(1,len(drug_groups)+1), drug_groups)
plt.ylim([0,12])
plt.ylabel('sHIVD');
plt.savefig('figures/NeuroSpecificFigures/cyto-sHIVD.png')

# <markdowncell>

# Here we can see that there may be a difference between the different drug groupings. It would seem that Cannabinoids have less variation then Non-Users. However, this figure has all data points lumped into one boxplot. Would it look different if we grouped by year?

# <codecell>

from itertools import product

fig, axs = plt.subplots(2,2, sharex=True, sharey=True, figsize=(10,10))

ax = axs.flatten()[0]
plt.sca(ax)
tdata = []
for group in drug_groups:
    mask = cyto_dataset['DrugStatus'] == group
    neuro_data = cyto_dataset['iHIVD'][mask]
    tdata.append(neuro_data.dropna().values)
plt.boxplot(tdata);
plt.xticks(range(1,len(drug_groups)+1), drug_groups)
plt.ylim([0,12])
plt.ylabel('sHIVD');
plt.title('Intake')


for year, ax in enumerate(axs.flatten()[1:], 1):
    plt.sca(ax)
    year_data = date_data['Year'] == year
    tdata = []
    for group in drug_groups:
        mask = cyto_dataset['DrugStatus'] == group
        neuro_data = cyto_dataset['sHIVD'][mask]
        tdata.append(neuro_data.dropna().values)

    plt.boxplot(tdata);
    plt.xticks(range(1,len(drug_groups)+1), drug_groups)
    plt.ylim([0,12])
    plt.ylabel('sHIVD');
    plt.title('HIVD - year %i' % (year))
plt.savefig('figures/NeuroSpecificFigures/cyto-YearHIVD.png')

# <markdowncell>

# Looking at this small dataset I cannot see a difference between Cannabinoid and Non-Users. If anything it looks like Benzos have a protective effect.

# <markdowncell>

# This is another way of looking at the same data. Here I've plotted the HIVD score afer X-years vs the Intake Score. You can see that there is a strong linear relationship between the intake HIVD score and the Year-1 and Year-2 score, by year 4 there just isn't enough data. It does look like there is a difference in slopes by Year 2 with the Pure-Cocaine (red) and Pure-Cannabinoid (green) having a smaller slope then the Non-users (grey) and MDU (blue). It seems so odd that Non-Users are on-par with MDU when you would expect the opposite! I wonder if I group people by intake scores and then by drug-use if this effect will go away.

# <codecell>

vals = []
for grouping in cyto_dataset['DrugStatus'].unique():
    good_mask = cyto_dataset['DrugStatus'] == grouping
    tdata = cyto_dataset[good_mask]
    for pat, rows in tdata.groupby('Patient ID'):
        cum_diff = rolling_apply(rows['sHIVD'], 2, np.diff)
        tmp = DataFrame({'sHIVD': rows['sHIVD'], 'dHIVD':cum_diff})
        tmp['DrugStatus'] = grouping
        tmp['Age'] = rows['Age']
        tmp['TimeSeropositive'] = rows['TimeSeropositive']
        tmp['ART'] = rows['ART']
        tmp['iHIVD'] = rows['iHIVD']
        tmp['Patient ID'] = pat
        tmp['Year'] = rows['Year']
        tmp['Race'] = rows['Race']
        vals.append(tmp.reset_index())
deltaData = concat(vals, axis = 0, ignore_index=True)
print deltaData.dropna().head(n=10).to_string()

# <codecell>

data['Race'].unique()

# <codecell>

import statsmodels.api as sm
from patsy import dmatrices, dmatrix
from copy import deepcopy
result_objs = []

group_data = deltaData.dropna()
y, X = dmatrices('sHIVD ~ Year+iHIVD+TimeSeropositive+Age+C(DrugStatus, Treatment("PN"))+C(ART, Treatment("naive"))+C(Race, Treatment("White"))', group_data, return_type = 'dataframe')
res = sm.GLM(y,X).fit()
print res.summary()
print res.pvalues
print 
print res.params[res.pvalues < 0.05]
group_data['cHIVD'] = res.fittedvalues
result_objs.append(('Cytokine', 'Grouping', deepcopy(res)))

# <markdowncell>

# We can see that there is a HUGE correlation with the initial HIVD score and has a positive effect on the nHIVD.

# <codecell>

fig, axs = plt.subplots(2,2, sharex=True, sharey=True, figsize = (10, 10))
for year, ax in enumerate(axs.flatten(), 1):
    if year == 4:
        ymask = True
        title = 'All Years'
    else:
        title = 'Score at year %i' % year
        ymask = group_data['Year'] == year
    tdata = []
    for group in drug_groups:
        gmask = group_data['DrugStatus'] == group
        tdata.append(group_data['cHIVD'][ymask & gmask].dropna().values)
    plt.sca(ax)
    plt.boxplot(tdata)
    plt.xticks(range(1,len(drug_groups)+1), drug_groups, rotation=90)
    plt.title(title)
    plt.ylabel('Adjusted HIVD')
plt.savefig('figures/NeuroSpecificFigures/ctyo-Adj-HIVD-Group.png', dpi = 200)

# <markdowncell>

# Well, it certainly looks like Can users have an increased HIVD score. I wonder if this holds true in the 'linear contribution' model.

# <codecell>

def tmp_rolling(inser):
    return -rolling_apply(inser, 2, np.diff, min_periods=2)
    

#Since we don't want the Year-0 and Year > 5 patients skewing the data
linear_data = cyto_dataset[(cyto_dataset['Year']>0) & (cyto_dataset['Year']<5)].dropna()
print linear_data.drop(['Date', 'Opiates'], axis  =1).head().to_string()

# <codecell>

y, X = dmatrices('sHIVD ~ iHIVD+Cannabinoid+Cocaine+Year+TimeSeropositive+Age+C(ART, Treatment("naive"))+C(Race, Treatment("White"))', 
                    linear_data, return_type = 'dataframe')

res = sm.GLM(y,X).fit()
print res.summary()
print res.pvalues, '\n'
print res.params[res.pvalues < 0.05]
linear_data['cHIVD'] = res.fittedvalues
result_objs.append(('Cytokine', 'LCM', deepcopy(res)))

# <markdowncell>

# So we see a huge effect of Cannabinoid use! Using Cannabinoids causes a 0.77 increase in HIVD score.

# <codecell>

from scipy.stats import linregress
import re
from types import StringType
robj = re.compile('[\d.]+')
def get_mid(obj):
    nums = robj.findall(str(obj))
    if len(nums) != 2:
        return np.nan
    return float(nums[1])

fig, axs = plt.subplots(2,3, sharey=True, figsize = (15, 10))
groups = [('DrugFree', np.arange(0,1.1,0.1), 0.1),
            ('Cannabinoid', np.arange(0,1.1,0.1), 0.1), 
            ('Cocaine', np.arange(0,1.1,0.1), 0.1),
            ('TimeSeropositive', np.arange(0,30,5), 5),
            ('ART', ['on', 'off', 'naive', 'non-adherent'], None),
            ('Age', np.arange(20,70, 10), 10),]
for ax, (group, bins, width) in zip(axs.flatten(), groups):
    
    tdata = linear_data[[group, 'cHIVD']].dropna()
    
   
    if type(bins[0]) != StringType:
        pos = Series(map(get_mid, cut(tdata[group], bins)), 
                index = tdata.index)
        tdata['Pos'] = pos
        tdata = tdata.dropna().sort('Pos')
        m, b, r, p, _ = linregress(tdata[group].values, y = tdata['cHIVD'].values)
        x = np.linspace(bins[0]-width, bins[-1]+width)
        y = m*x + b
    else:
        tdata['Pos'] = tdata[group]
    #tdata.boxplot(column = 'cHIVD', by = 'Pos', ax = ax)
    tmp = []
    for b in bins:
        tmp.append(tdata['cHIVD'][tdata['Pos'] == b].values)
    plt.sca(ax)
    plt.hold(True)
    if type(bins[0]) != StringType:
        plt.plot(x, y, lw=10, alpha=0.2, color = 'r')
        plt.boxplot(tmp, positions = bins, widths=width*0.7)
        plt.xlim([bins[0]-width, bins[-1]+width])
    else:
        plt.boxplot(tmp)
        plt.xticks(range(1, len(bins)+1), bins)
    plt.title(group)
    plt.ylabel('Adj-HIVD')
    plt.hold(False)
       
plt.savefig('figures/NeuroSpecificFigures/cyto-Adj-HIVD-Trends.png', dpi = 200)
    
    
    

# <codecell>

vals = []
for grouping in date_data['NizGroupings'].unique():
    good_mask = date_data['NizGroupings'] == grouping
    tdata = date_data[good_mask]
    for pat, rows in tdata.groupby('Patient ID'):
        cum_diff = rolling_apply(rows['sHIVD'], 2, np.diff)
        tmp = DataFrame({'sHIVD': rows['sHIVD'], 'dHIVD':cum_diff})
        tmp['DrugStatus'] = grouping
        tmp['Age'] = rows['Age']
        tmp['TimeSeropositive'] = rows['TimeSeropositive']
        tmp['ART'] = rows['ART']
        tmp['iHIVD'] = rows['iHIVD']
        tmp['Patient ID'] = pat
        tmp['Year'] = rows['Year']
        tmp['Race'] = rows['Race']
        vals.append(tmp.reset_index())
large_deltaData = concat(vals, axis = 0, ignore_index=True)
print large_deltaData.dropna().head(n=10).to_string()

# <codecell>

import statsmodels.api as sm
from patsy import dmatrices, dmatrix
group_data = large_deltaData.dropna()
y, X = dmatrices('sHIVD ~ Year+iHIVD+TimeSeropositive+Age+C(DrugStatus, Treatment("PN"))+C(ART, Treatment("naive"))+C(Race, Treatment("White"))', group_data, return_type = 'dataframe')
res = sm.GLM(y,X).fit()
print res.summary()
print res.pvalues
print 
print res.params[res.pvalues < 0.05]
group_data['cHIVD'] = res.fittedvalues
result_objs.append(('Entire', 'Grouping', deepcopy(res)))

# <codecell>

fig, axs = plt.subplots(2,2, sharex=True, sharey=True, figsize = (10, 10))
for year, ax in enumerate(axs.flatten(), 1):
    if year == 4:
        ymask = True
        title = 'All Years'
    else:
        title = 'Score at year %i' % year
        ymask = group_data['Year'] == year
    tdata = []
    for group in drug_groups:
        gmask = group_data['DrugStatus'] == group
        tdata.append(group_data['cHIVD'][ymask & gmask].dropna().values)
    plt.sca(ax)
    plt.boxplot(tdata)
    plt.xticks(range(1,len(drug_groups)+1), drug_groups, rotation=90)
    plt.title(title)
    plt.ylabel('Adjusted HIVD')
plt.savefig('figures/NeuroSpecificFigures/large-Adj-HIVD-Group.png', dpi = 200)

# <codecell>

large_linear_data = date_data.dropna()
y, X = dmatrices('sHIVD ~ iHIVD+Cannabinoid+Cocaine+Year+TimeSeropositive+Age+C(ART, Treatment("naive"))+C(Race, Treatment("White"))', 
                    large_linear_data, return_type = 'dataframe')

res = sm.GLM(y,X).fit()
print res.summary()
print res.pvalues, '\n'
print res.params[res.pvalues < 0.05]
large_linear_data['cHIVD'] = res.fittedvalues
result_objs.append(('Entire', 'LCM', deepcopy(res)))

# <codecell>

fig, axs = plt.subplots(2,3, sharey=True, figsize = (15, 10))
groups = [('DrugFree', np.arange(0,1.1,0.1), 0.1),
            ('Cannabinoid', np.arange(0,1.1,0.1), 0.1), 
            ('Cocaine', np.arange(0,1.1,0.1), 0.1),
            ('TimeSeropositive', np.arange(0,30,5), 5),
            ('ART', ['on', 'off', 'naive', 'non-adherent'], None),
            ('Age', np.arange(20,70, 10), 10),]
for ax, (group, bins, width) in zip(axs.flatten(), groups):
    
    tdata = large_linear_data[[group, 'cHIVD']].dropna()
    
   
    if type(bins[0]) != StringType:
        pos = Series(map(get_mid, cut(tdata[group], bins)), 
                index = tdata.index)
        tdata['Pos'] = pos
        tdata = tdata.dropna().sort('Pos')
        m, b, r, p, _ = linregress(tdata[group].values, y = tdata['cHIVD'].values)
        x = np.linspace(bins[0]-width, bins[-1]+width)
        y = m*x + b
    else:
        tdata['Pos'] = tdata[group]
    #tdata.boxplot(column = 'cHIVD', by = 'Pos', ax = ax)
    tmp = []
    for b in bins:
        tmp.append(tdata['cHIVD'][tdata['Pos'] == b].values)
    plt.sca(ax)
    plt.hold(True)
    if type(bins[0]) != StringType:
        plt.plot(x, y, lw=10, alpha=0.2, color = 'r')
        plt.boxplot(tmp, positions = bins, widths=width*0.7)
        plt.xlim([bins[0]-width, bins[-1]+width])
    else:
        plt.boxplot(tmp)
        plt.xticks(range(1, len(bins)+1), bins)
    plt.title(group)
    plt.ylabel('Adj-HIVD')
    plt.hold(False)
       
plt.savefig('figures/NeuroSpecificFigures/large-Adj-HIVD-Trends.png', dpi = 200)

# <codecell>

nindex = []
ndata = []
for cohort, anal, res in result_objs:
    nindex.append((cohort, anal, 'pvalues'))
    nindex.append((cohort, anal, 'effects'))
    ndata.append(res.pvalues.copy())
    ndata.append(res.params.copy())
    
mi = MultiIndex.from_tuples(nindex, names = ['Cohort', 'Analysis', 'Type'])
df = DataFrame(ndata, index = mi).reset_index()

# <codecell>

df.reset_index().T.to_excel('figures/NeuroSpecificFigures/regression_results.xlsx')

# <codecell>

def normalize_by_date_nosample(indf):
    tdf = indf.reset_index().set_index('Date').drop(['Patient ID', 'VisitNum'], axis = 1)
    tdf[drug_names + ['DrugFree']] = tdf[drug_names + ['DrugFree']].applymap(float)
    tdata = rolling_mean(tdf[['nHIVD','DrugFree']+drug_names], 2, min_periods=1, freq = '6M', axis = 0)#.reset_index()
    odf = tdf.resample('6M', how = 'last')
    odf[['DrugFree']+drug_names] = tdata[['DrugFree']+drug_names]
    odf['sHIVD'] = tdata['nHIVD']
    try:
        odf['iHIVD'] = tdf['nHIVD'].dropna().values[0]
    except IndexError:
        odf['iHIVD'] = np.nan
    odf['LsHIVD'] = odf['sHIVD'].shift(1)
    odf['dsHIVD'] = rolling_apply(odf['sHIVD'], 2, np.diff)
    odf['LnHIVD'] = odf['nHIVD'].shift(1)
    #odf[]
    #odf['ART'] = tdf['ART'].resample('12M', how = 'last')
    odf = odf.reset_index()
    odf.index.name = 'Year'
    return odf

new_date_data = data.groupby(level = 'Patient ID').apply(normalize_by_date_nosample).reset_index()
print new_date_data

# <codecell>

tmp_data = new_date_data.dropna(subset = ['nHIVD','LnHIVD','dsHIVD','sHIVD', 'Cannabinoid', 'Cocaine', 'TimeSeropositive', 'Age', 'ART'])
tmp_data = tmp_data[tmp_data['NizGroupings'] != 'PO']
tmp_data = tmp_data[tmp_data['Year']<6]

y, X = dmatrices('sHIVD ~ iHIVD+Year+Cocaine+Cannabinoid+TimeSeropositive+Age+ART', 
                    tmp_data, return_type = 'dataframe')
res = sm.GLM(y, X).fit()
#print tmp_data
print res.summary()

# <codecell>

def resolve_date_last_used(inval):
    
    if '/' in inval:
        parts = [x for x in inval.split('/')]
        if len(parts) == 2:
            month = int(parts[0])
            year = parts[1]
            day = None
        elif len(parts) == 3:
            day = int(parts[1])
            month = int(parts[0])
            year = parts[2]
        else:
            print parts, inval
            raise ValueError
    elif '-' in inval:
        parts = [x for x in inval.split('-')]
        if len(parts) == 2:
            month = int(parts[1])
            year = parts[0]
            day = None
        elif len(parts) == 3:
            day = int(parts[2])
            month = int(parts[1])
            year = parts[0]
        else:
            print parts, inval
            raise ValueError
    else:
        day = None
        month = None
        year = inval
        
    if len(year) == 2:
        if int(year) > 20:
            year = '19'+year
        else:
            year = '20'+year
    freq = 'D'
    if day is None:
        day = 1
        freq = 'M'
    if month is None:
        month = 1
        freq = 'Y'
    
    #try:
    #print year, month, day, freq
    guess_date = Period(year = int(year), month = month, day = day, freq = freq)
    #except ValueError:
    #    print inval
        #rise ValueError
    #    return np.nan
    return guess_date
get_days = lambda x: x/np.timedelta64(1,'D')
drug_started_data = read_csv('DrugStartedData.tsv', sep = '\t', parse_dates=['Date of visit', 'HIV seropositive date']).set_index(['Patient ID', 'VisitNum'])
drug_started_data['DateLastUsedDrugs'] = drug_started_data['Date last used drugs'].dropna().map(resolve_date_last_used)
drug_started_data['DateLastUsedDrugs'] = drug_started_data['DateLastUsedDrugs'].dropna().map(lambda x: x.to_timestamp('D', how='s'))
drug_started_data['DaysSinceLastDrugUse'] = (drug_started_data['Date of visit'] - drug_started_data['DateLastUsedDrugs']).dropna().map(get_days)

# <codecell>

tested_pos = data[all_drug_names].dropna().any(axis = 1)
pos_tests, last_drug = tested_pos.align(drug_started_data['DaysSinceLastDrugUse'].dropna(), join = 'inner')
last_drug_year = last_drug/365
fig, axs = plt.subplots(2,1, figsize=(10,10))

plt.sca(axs.flatten()[0])
plt.hold(True)
last_drug.hist(bins = np.linspace(0, 30, 30), ax = plt.gca())
last_drug[pos_tests].hist(bins = np.linspace(0, 30, 30), ax = plt.gca())
plt.legend(['Neg', 'Pos'], 'upper right')
plt.xlabel('Days since last drug use - Admit')
plt.ylabel('#Visits')
plt.hold(False)

plt.sca(axs.flatten()[1])
plt.hold(True)
last_drug_year.hist(bins = np.linspace(0,10, 20), ax = plt.gca())
last_drug_year[pos_tests].hist(bins = np.linspace(0,10, 20), ax = plt.gca())
plt.legend(['Neg', 'Pos'], 'upper right')
plt.xlabel('Years since last drug use - Admit')
plt.ylabel('#Visits')
plt.hold(False)

plt.savefig('figures/NeuroSpecificFigures/lying_admitters.png')

# <codecell>

cols = [('Initial CD4 count (cells/uL)','CD4','Date of initial CD4 count'),
        ('Nadir CD4 count (cells/uL)', 'CD4', 'Date of nadir CD4 count'),
        ('Latest CD4 count (cells/uL)', 'CD4', 'Date of latest CD4 count'),
        ('Initial CD8 count (cells/uL)','CD8', 'Date of initial CD8 count'),
        ('Nadir CD8 count (cells/uL)','CD8', 'Date of nadir CD8 count'),
        ('Latest CD8 count (cells/uL)','CD8', 'Date of latest CD8 count'),
        ('Initial viral load (copies/mL)','VL', 'Date of initial viral load'),
        ('Peak viral load (copies/mL)','VL', 'Date of peak viral load'),
        ('Latest viral load','VL', 'Date of latest viral load'),
        ('Amphetamines', 'Test-Amphetamines'),
        ('Barbiturates', 'Test-Barbiturates'),
        ('Benzodiazepines', 'Test-Benzodiazepines'),
        ('Cannabinoid', 'Test-Cannabinoid'),
        ('Cocaine + metabolite', 'Test-Cocaine'),
        ('Opiates','Test-Opiates'),
        ('Phencyclidine', 'Test-Phencyclidine'),
        ("Race (choice='Asian')",'Asian'),
        ("Race (choice='American Indian/Alaska Native')",'Indian'),
        ("Race (choice='Black or African American')", 'Black'),
        ("Race (choice='Native Hawaiian or other Pacific Islander')", 'Hawaiian'),
        ("Race (choice='White')", 'White'),
        ("Race (choice='More than one race')",'Multi-Race'),
        ("Race (choice='Unknown')", 'Unknown'),
        ("Drugs used (choice='Marijuana')",'Admit-Cannabinoid'),
        ("Drugs used (choice='Cocaine (crack, nasal, smoke, inject)')", 'Admit-Cocaine'),
        ("Drugs used (choice='Heroin (nasal, inject)')",'Admit-Opiates'),
        ("Drugs used (choice='Methamphetamine (smoke, nasal, inject)')",'Amphetamines'),
        ("Drugs used (choice='Benzodiazapine (i.e. valium, ativan, xanax, klonipin, etc)')",'Admit-Benzodiazepines'),
        ("Drugs used (choice='Narcotics')",'Admit-Opiates'),
        ("Drugs used (choice='Ecstasy')",'Admit-Phencyclidine'),
        ("Drugs used (choice='PCP')",'Admit-Phencyclidine'),
        ("Drugs used (choice='Ritalin')",'Admit-Ritalin'),
        ("Drugs used (choice='Other')", 'Admit-Other'),
        ('Year of Birth', 'Year of Birth'),
        ('HIV seropositive date', 'HIV seropositive date'),
        ('Current ART status', 'ART'),
        ('Total Modified Hopkins Dementia Score', 'HIVD'),
        ('Drug Use and HIV Status', 'Drug Use and HIV Status'),
        ('Age first used drug', 'Age first used drug'),
        ('VisitNum', 'VisitNum')]
tmp_list = []
for tup in cols:
    if len(tup) == 2:
        col, typ = tup
        date_col = 'Date of visit'
    else:
        col, typ, date_col = tup
    tmp = redcap_data[['Patient ID', col, date_col]].dropna()
    tmp['Type'] = typ
    tmp = tmp.rename(columns={col:'Measure', date_col:'Date'})
    tmp_list.append(tmp)

drug_started_data = read_csv('DrugStartedData.tsv', sep = '\t', parse_dates=['Date of visit', 'HIV seropositive date'])
drug_started_data['DateLastUsedDrugs'] = drug_started_data['Date last used drugs'].dropna().map(resolve_date_last_used)
drug_started_data['DateLastUsedDrugs'] = drug_started_data['DateLastUsedDrugs'].dropna().map(lambda x: x.to_timestamp('D', how='s'))
drug_started_data['DaysSinceLastDrugUse'] = (drug_started_data['Date of visit'] - drug_started_data['DateLastUsedDrugs']).dropna().map(get_days)

tmp = drug_started_data[['Patient ID', 'Date of visit', 'DaysSinceLastDrugUse']]
tmp['Type'] = 'DaysSinceLastDrugUse'
tmp_list.append(tmp.rename(columns = {'Date of visit':'Date', 'DaysSinceLastDrugUse':'Measure'}))
    
clinical_params = concat(tmp_list, axis = 0, ignore_index=True)
clinical_params = pivot_table(clinical_params, rows = ['Patient ID', 'Date'], cols = 'Type', values='Measure', aggfunc='first')
clinical_params.head(n=10)

# <codecell>

from datetime import timedelta
def safe_mean(inser):
    try:
        return inser.mean()
    except TypeError:
        return np.nan
        
def safe_max(inser):
    try:
        return inser.max()
    except TypeError:
        return np.nan

def safe_min(inser):
    try:
        return inser.min()
    except TypeError:
        return np.nan
        
def most_common(inser):
    nser = inser.dropna().value_counts()
    try:
        return nser.index[0]
    except IndexError:
        return np.nan

first_visit = redcap_data[['Patient ID', 'Date of visit']].groupby('Patient ID').first()['Date of visit']
last_visit = redcap_data[['Patient ID', 'Date of visit']].groupby('Patient ID').last()['Date of visit']
date_seropositive = redcap_data[['Patient ID', 'HIV seropositive date']].groupby('Patient ID').agg(most_common)

test_cols = ['Test-Amphetamines',
             'Test-Barbiturates',
             'Test-Benzodiazepines',
             'Test-Cannabinoid',
             'Test-Cocaine',
             'Test-Opiates',
             'Test-Phencyclidine']

admit_cols = ['Admit-Benzodiazepines',
             'Admit-Cannabinoid',
             'Admit-Cocaine',
             'Admit-Opiates',
             'Admit-Other',
             'Admit-Phencyclidine',
             'Admit-Ritalin']
other_admit_cols = ['DaysSinceLastDrugUse']
clinical_cols = ['VL', 'CD4', 'CD8', 'HIVD']

param_groups = [('mean-', safe_mean, test_cols+admit_cols+clinical_cols+other_admit_cols),
                ('max-', safe_max, clinical_cols+other_admit_cols),
                ('min-', safe_min, clinical_cols),
                ('', most_common, ['ART'])]


def calc_params(indf):
    
    pat_id = indf.index[0][0]
    indf = indf.ix[pat_id]
    fday = min(first_visit.ix[pat_id], last_visit.ix[pat_id])
    lday = max(first_visit.ix[pat_id], last_visit.ix[pat_id])
    
    pers = date_range(fday, lday, freq = '12M')
    indf = indf.truncate(fday-timedelta(days = 1), lday)
    
    out_data = []
    for prepend, func, cols in param_groups:
        ndf = indf[cols].resample('12M', closed = 'right', label = 'right', how = func)
        ncols = dict([(name, prepend+name) for name in cols])
        ndf.rename(columns = ncols, inplace = True)
        out_data.append(ndf.copy())
    odf = concat(out_data, axis = 1)
    #print 'date', date_seropositive.ix[pat_id]
    #print 'series', Series(odf.index)
    #print date_seropositive.ix[pat_id].values - np.array(list(odf.index), np.dtype('<M8[ns]'))
    if date_seropositive.ix[pat_id].dropna():
        odf['DaysSinceSeropositive'] = np.array(list(odf.index), np.dtype('<M8[ns]')) - date_seropositive.ix[pat_id].values
    else:
        odf['DaysSinceSeropositive'] = np.nan
        odf['DaysSinceSeropositive'] = odf['DaysSinceSeropositive'].astype('datetime64[ns]')
    return odf

nres = []
for pat, group in clinical_params.groupby(level = 'Patient ID'):
    tres = calc_params(group).reset_index()
    tres['Patient ID'] = pat
    tres['VisitYear'] = range(len(tres))
    nres.append(tres.copy())

res = concat(nres, axis = 0, ignore_index = True)
res = res.set_index(['Patient ID', 'VisitYear'])
#pat_df['Cannabinoid'].resample('12M', closed = 'right', label = 'right', how = safe_mean)#.truncate(fday, lday)

# <codecell>

ban_cols = set(['mean-Test-Amphetamines',
            'mean-Test-Barbiturates',
            'mean-Test-Benzodiazepines',
            'mean-Test-Cocaine',
            'mean-Test-Opiates',
            'mean-Test-Phencyclidine',
            'mean-Test-Cannabinoid',
            'mean-Admit-Benzodiazepines',
            'mean-Admit-Cocaine',
            'mean-Admit-Cannabinoid',
            'mean-Admit-Opiates',
            'mean-Admit-Other',
            'mean-Admit-Phencyclidine',
            'mean-Admit-Ritalin'])
require_cols = [(set(), 'PN'),
                (set(['mean-Test-Cocaine', 'mean-Admit-Cocaine']), 'Cocaine'),
                (set(['mean-Test-Cannabinoid', 'mean-Admit-Cannabinoid']), 'Cannabinoid'),
                (set(['mean-Test-Cocaine', 'mean-Admit-Cocaine', 'mean-Test-Cannabinoid', 'mean-Admit-Cannabinoid']), 'Co-Ca')]

tmp_drug_cols = []
for allowed, name in require_cols:
    zero_cols = list(ban_cols - allowed)
    one_cols = list(allowed)
    take_nothing_but_alowed = (res[zero_cols].fillna(0)==0).groupby(level = 'Patient ID').all().all(axis = 1)
    take_allowed = (res[one_cols].fillna(1)==1).groupby(level = 'Patient ID').all().all(axis = 1)
    
    wanted_pats = (take_nothing_but_alowed & take_allowed)
    wanted_pats.name = name
    tmp_drug_cols.append(wanted_pats.copy())

pure_pats = concat(tmp_drug_cols, axis = 1)
print pure_pats.sum()

# <codecell>

years = range(4)
cols = [('PN', 'k'), 
        ('Cocaine', 'b'), 
        ('Cannabinoid' ,'g'), 
        ('Co-Ca', 'r')]
plt.figure(figsize = (10,10))
plt.hold(True)
tdata = merge(res[['min-HIVD']].reset_index(), pure_pats,
                left_on = 'Patient ID', right_index = True,
                how = 'inner')
for col, color in cols:
    wdata = tdata[tdata[col]]
    for pat, group in wdata.groupby('Patient ID'):
        plt.plot(group['VisitYear'], group['min-HIVD'], color = color, alpha = 0.5, lw = 10)
    

# <codecell>


