# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=1>

# New Cannabinoid Analysis

# <markdowncell>

# Here I am looking for the effect of Cannabinoid use on the progression of HIV related dementia. This is currently being measured by the Total Modified Hopkins Dementia Score (HIVD-score). This is a non-trivial analysis due to a few confouding issues:
# 
#  - HIVD score is a catagorical variable between 1 and 12.
#  - The precision and accuracy of this score is difficult to assess.
#  - The HIVD score was not gathered at every examination.
#  - Examinations are not equally spaced.
#  - Pure cannabinoid use is rare in our cohort.
#  - The difference between _confirmed_ drug-use (by blood test) and _admitted_ drug-use.
#  - The varying sevarity of HIV disease among these individuals.

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

# Data Extraction

# <markdowncell>

# I'm currently using the 1/16/2013 version of the Redcap database.

# <codecell>

store_loc = os.path.join('home', 'will', 'HIVReportGen', 'Data',
                        'SplitRedcap', '2013-01-16',
                        'EntireCohort.hdf')
store = HDFStore('/'+store_loc)
redcap_data = store['redcap']

#fix the Visit Number issue
t = redcap_data['Event Name'].dropna().apply(lambda x: x.split(' - ')[0])
vnum_field = 'Patient visit number'
redcap_data[vnum_field] = redcap_data[vnum_field].combine_first(t)

#Fix the Other-drug issue
ofield = "Drugs used (choice='Other')"
redcap_data[ofield] = redcap_data[ofield].dropna() == 'Checked'

# <markdowncell>

# Here are the columns I'm extracting from the database. I'm treating the CD4, CD8, and VL counts in a different way from our 'traditional' method. Since the database includes dates associated with each  datapoint I'm keeping them. This allows me to do a better job interpretting the clinical paramaters that should be associated with each visit. It also lets me treat things in a 'visit-independent' method. It also lets me use time-based sampling methods which I'll mention later.

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
        ("Race (choice='Asian')",'Race-Asian'),
        ("Race (choice='American Indian/Alaska Native')",'Race-Indian'),
        ("Race (choice='Black or African American')", 'Race-Black'),
        ("Race (choice='Native Hawaiian or other Pacific Islander')", 'Race-Hawaiian'),
        ("Race (choice='White')", 'Race-White'),
        ("Race (choice='More than one race')",'Race-Multi-Race'),
        ("Race (choice='Unknown')", 'Race-Unknown'),
        ("Drugs used (choice='Marijuana')",'Admit-Cannabinoid'),
        ("Drugs used (choice='Cocaine (crack, nasal, smoke, inject)')", 'Admit-Cocaine'),
        ("Drugs used (choice='Heroin (nasal, inject)')",'Admit-Opiates'),
        ("Drugs used (choice='Methamphetamine (smoke, nasal, inject)')",'Admit-Amphetamines'),
        ("Drugs used (choice='Benzodiazapine (i.e. valium, ativan, xanax, klonipin, etc)')",'Admit-Benzodiazepines'),
        ("Drugs used (choice='Narcotics')",'Admit-Opiates'),
        ("Drugs used (choice='Ecstasy')",'Admit-Phencyclidine'),
        ("Drugs used (choice='PCP')",'Admit-Phencyclidine'),
        ("Drugs used (choice='Ritalin')",'Admit-Ritalin'),
        ("Drugs used (choice='Other')", 'Admit-Other'),
        ('Current Alcohol use', 'Admit-Alcohol'),
        ('Year of Birth', 'Year of Birth'),
        ('HIV seropositive date', 'HIV seropositive date'),
        ('Current ART status', 'ART'),
        ('Total Modified Hopkins Dementia Score', 'HIVD'),
        ('Drug Use and HIV Status', 'Drug Use and HIV Status'),
        ('Age first used drug', 'Age first used drug'),
        ('Patient visit number', 'VisitNum'),
        ('Gender', 'Gender'),
        ]
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

get_days = lambda x: x/np.timedelta64(1,'D')
    
clinical_params = concat(tmp_list, axis = 0, ignore_index=True)
clinical_params = pivot_table(clinical_params, 
                                rows = ['Patient ID', 'Date'], 
                                cols = 'Type', values='Measure', 
                                aggfunc='first')

# <codecell>

clinical_params['ART-on'] = clinical_params['ART'] == 'on'
clinical_params['ART-off'] = clinical_params['ART'] == 'off'
clinical_params['ART-nonadherent'] = clinical_params['ART'] == 'non-adherent'
clinical_params['ART-naive'] = clinical_params['ART'] == 'naive'
age_func = lambda x: x['Date'].year - x['Year of Birth']
ages = clinical_params['Year of Birth'].reset_index().apply(age_func, axis = 1)
ages.index = clinical_params.index
clinical_params['Age'] = ages

has_date_mask = clinical_params['HIV seropositive date'].notnull()
clinical_params['YearsSeropositive'] = np.nan
years_sero_func = lambda x: x['Date'].year - x['HIV seropositive date'].year
years_sero_data = clinical_params['HIV seropositive date'].reset_index().dropna().apply(years_sero_func, axis = 1)
years_sero_data.index = has_date_mask[has_date_mask].index
clinical_params['YearsSeropositive'][has_date_mask] = years_sero_data

clinical_params['Admit-Alcohol'] = clinical_params['Admit-Alcohol'] == 'Yes'

clinical_params['AoT-Alcohol'] = clinical_params['Admit-Alcohol']
clinical_params['AoT-Benzodiazepines'] = clinical_params[['Test-Benzodiazepines', 'Admit-Benzodiazepines']].mean(axis = 1)>0
clinical_params['AoT-Cannabinoid'] = clinical_params[['Test-Cannabinoid', 'Admit-Cannabinoid']].mean(axis = 1)>0
clinical_params['AoT-Cocaine'] = clinical_params[['Test-Cocaine', 'Admit-Cocaine']].mean(axis = 1)>0
clinical_params['AoT-Opiates'] = clinical_params[['Test-Opiates', 'Admit-Opiates']].mean(axis = 1)>0
clinical_params['AoT-Other'] = clinical_params['Admit-Other']
clinical_params['AoT-Phencyclidine'] = clinical_params[['Test-Phencyclidine', 'Admit-Phencyclidine']].mean(axis = 1)>0
clinical_params['AoT-Ritalin'] = clinical_params['Admit-Ritalin']
clinical_params['AoT-Amphetamines'] = clinical_params[['Test-Amphetamines', 'Admit-Amphetamines']].mean(axis = 1)>0


drop_cols = ['Drug Use and HIV Status', 'Year of Birth', 
            'HIV seropositive date', 'Age first used drug']
clinical_params = clinical_params.drop(drop_cols, axis = 1)

# <codecell>

admit_cols =['Admit-Alcohol', 'Admit-Benzodiazepines', 
            'Admit-Cannabinoid', 'Admit-Cocaine', 
            'Admit-Opiates', 'Admit-Other', 
            'Admit-Phencyclidine', 'Admit-Ritalin', 
            'Admit-Amphetamines'] 
test_cols = ['Test-Amphetamines', 'Test-Barbiturates', 
            'Test-Benzodiazepines', 'Test-Cannabinoid', 
            'Test-Cocaine', 'Test-Opiates', 'Test-Phencyclidine']
therapy_cols = ['ART-on', 'ART-off', 'ART-nonadherent', 'ART-naive']
clinical_cols = ['CD4', 'CD8', 'VL', 'HIVD']
patient_cols = ['Age', 'YearsSeropositive']
race_cols = ['Race-Asian', 'Race-Black', 'Race-Hawaiian', 
            'Race-Indian', 'Race-Multi-Race', 'Race-Unknown', 'Race-White']

fill_cols = admit_cols+therapy_cols+clinical_cols+patient_cols+race_cols
filled_data = clinical_params[fill_cols].groupby(level = 'Patient ID').fillna(method = 'pad')

date_data = clinical_params.combine_first(filled_data)
visit_data = date_data.dropna(subset = ['VisitNum']).reset_index().set_index(['Patient ID', 'VisitNum'])

# <headingcell level=2>

# Checking Patient Definitions

# <markdowncell>

# The hardest thing to do here is to determine which patients should be in the analysis. Since we're looking _logitudinally_ its important that the patients have all relevant datapoints at _EVERY_ visit we analyze. I also believe that the visits should be consecutive. I don't think every patient needs to start with R00 (in-fact this would remove a significant portion of the patients) but it should be N consecutive visits.
# 
# We also have to consider the Admit and Test drug data. We've gotten any number of comments about whether to include admit data. Or whether to exclude patients based on admit data. This can work to both _increase_ or _decrease_ our dataset sizes. There are numerous patients who always admit to Cocaine (for example) but occasonally test negative. Should they be excluded from the 'Pure Cocaine' group because they happened to skip thier rock the day before the exam? Some of this is resolved with our wLCM discussed later. However, since the wLCM is still a novel idea, its helpful to compare it to a more understood method.
# 
# For consistency with our previous analysis I'm limiting the data to people who checked African-American at ALL visits under consideration.

# <codecell>

def check_pat(indf, empty_cols, req_cols, num_visits, allow_skips):
    
    checked_data = []
    
    valid_rows = range(len(indf.index)-num_visits)
    for n in valid_rows:
        tdata = indf.iloc[n:n+num_visits]
        neg_mask = (tdata[list(empty_cols)] == 0).all(axis = 1)
        pos_mask = (tdata[list(req_cols)] == 1).all(axis = 1)
        missing_hivd = tdata['HIVD'].isnull().sum()
        is_aa = tdata['Race-Black'].all()
        wanted = neg_mask & pos_mask & is_aa
        
        if wanted.all():
            checked_data.append((missing_hivd, tdata.copy()))
    if checked_data:
        miss, data = min(checked_data, key = lambda x: x[0])
        if miss == 0:
            return data
        elif allow_skips and miss == 1:
            return data
        
            
    

def get_drug_data(num_visits, req_type, allow_skips):
    patient_defs = [('NonUser', []),
                    ('Pure-Canab', ['%s-Cannabinoid' % req_type]),
                    ('Pure-Coc', ['%s-Cocaine'% req_type]),
                    ('Co-Ca', ['%s-Cannabinoid'% req_type, '%s-Cocaine'% req_type])]

    drug_cols = set([col for col in visit_data.columns if col.startswith(req_type+'-')])
    wanted_pats = []

    for name, req_cols in patient_defs:
        empty_cols = drug_cols - set(req_cols)
        for pat, group in visit_data.groupby(level = 'Patient ID'):
            if len(group) >= num_visits:
                nres = check_pat(group, list(empty_cols), 
                                req_cols, num_visits, 
                                allow_skips)
                if nres is not None:
                    nres['Grouping'] = name
                    nres['ConVisit'] = range(num_visits)
                    wanted_pats.append(nres.reset_index())
    drug_data = concat(wanted_pats)
    return drug_data

# <codecell>

check_visits = [2, 3, 4, 5]
groups = ['NonUser', 'Pure-Canab', 'Pure-Coc', 'Co-Ca']
test_types = ['Test', 'Admit', 'AoT']
skip_defs = [True, False]
tdata = []

for num_visits, test_type, skip in product(check_visits, test_types, skip_defs):
    drug_data = get_drug_data(num_visits, test_type, skip)
    drug_data['AllowedSkip'] = skip
    drug_data['TestType'] = test_type
    drug_data['RequiredVisits'] = num_visits
    tdata.append(drug_data.copy())
    
group_data = concat(tdata, axis = 0, ignore_index = True)

# <codecell>

group_counts = pivot_table(group_data, cols = ['AllowedSkip', 'RequiredVisits'], 
                            rows = ['TestType', 'Grouping'], values='Patient ID', 
                            aggfunc = lambda x: len(x.unique())).fillna(0)
order = list(product(['Admit', 'Test', 'AoT'], 
                        ['NonUser', 'Pure-Canab', 
                            'Pure-Coc', 'Co-Ca']))
print group_counts.ix[order]

# <markdowncell>

# So we can see that our dataset size is drastically reduced by requiring the HIVD score at each visit. Even allowing a single skip doesn't help much and it will make the understanding much more complicated. You can also see that there is a HUGE difference whether we use Admit, Test or Admit-or-Test. I'm going to use the 3-visit data and show that the answer can drasctically change depending on our definitions.

# <codecell>

wanted_groups = group_data[group_data['RequiredVisits']==3]
mean_hivds = pivot_table(wanted_groups, rows = ['TestType', 'Grouping'], 
                            cols = ['AllowedSkip', 'ConVisit'], values = 'HIVD', 
                            aggfunc = lambda x: np.mean(x))

std_hivds = pivot_table(wanted_groups, rows = ['TestType', 'Grouping'], 
                            cols = ['AllowedSkip', 'ConVisit'], values = 'HIVD', 
                            aggfunc = lambda x: np.std(x))
print mean_hivds.ix[order][False].to_string(float_format=lambda x: '%10.2f' % x)

# <codecell>

colors = {
'NonUser':'k',
'Pure-Coc':'r',
'Pure-Canab':'g',
'Co-Ca':'b'
}

fig, axs = plt.subplots(2,3, figsize = (10,10), sharex = True, sharey = True)
order = product(enumerate(skip_defs), enumerate(test_types))
for (snum, skip), (tnum, test_type) in order:
    plt.sca(axs[snum, tnum])
    for offset, group in zip([-0.1, -0.05, 0, 0.05], groups):
        mean_data = mean_hivds.ix[test_type].ix[group][skip]
        err_data = std_hivds.ix[test_type].ix[group][skip]
        xvals = np.arange(3)+offset
        plt.errorbar(xvals, mean_data, yerr = err_data, 
                        color = colors[group], lw = 2, alpha = 0.7)
    plt.xlim([-0.5, 2.5])
    plt.xticks([0,1,2])
    plt.xlabel('Visit')
    plt.ylabel('HIVD')
    title = 'With-Skip' if skip else 'No-Skip'
    title += ' ' + test_type
    plt.title(title)

# <markdowncell>

# Above is the progression of HIVD scores under varying definitions of:
# 
# - Non-Users : Black
# - Pure-Canab : Green
# - Pure-Coc : Red
# - Coc-Ca : Blue
# 
# You can see that the story changes dramatically based on the different definitions of the groups. In the Admit and the AoT data it would appear that Canabinoids has a protective effect. However, you have to remember that this is based off of _one patient_. If you look at the Test data there are at least 3 patient in the Cannabinoid group but there doesn't seem to be an effect. It also appears in all of these cases that the MDU group has a higher HIVD score then Non-users. Also counter-intuitive. Also, patient sizes this small it would seem that its difficult to get a consistent story.
# 
# Since we need a figure for Niz's thesis I'm picking the 'Test with skip' since this is the only one that actually has more then one Pure-Canab patients.

# <codecell>

thesis_mask = (group_data['RequiredVisits']==3) & \
                (group_data['AllowedSkip'] == True) & \
                (group_data['TestType'] == 'Test')
nwanted = group_data[thesis_mask]
plot_data = pivot_table(wanted_groups, rows = ['Grouping', 'Patient ID'], 
                            cols = 'ConVisit', values = 'HIVD', 
                            aggfunc = lambda x: np.mean(x))

# <codecell>

fig, axs = plt.subplots(4,1, figsize = (4,10), sharex = True, sharey = True)
for ax, group in zip(axs.flatten(), groups):
    medians = plot_data.ix[group].median()
    ax.boxplot(plot_data.ix[group].fillna(medians).values, bootstrap=1000)
    
    ax.set_title(group + ' N:%i' % len(plot_data.ix[group]))
    ax.set_ylabel('HIVD')
    
ax.set_xlabel('Visit');
plt.savefig('figures/NeuroSpecificFigures/rawHIVD_grouped.png')

# <markdowncell>

# When we correct for confounders:
# 
# - Age
# - Gender
# - ART
# - Drug-Use Stattus
# 

# <codecell>

tmp_data = plot_data.reset_index()
tmp = []
for idx, row in tmp_data.iterrows():
    for col in [0,1,2]:
        tmp.append((row['Patient ID'], 
                    row['Grouping'],
                    col,
                    row[col]))
out_data = DataFrame(tmp, columns = ['Patient ID', 'Grouping', 'Visit', 'HIVD'])
pat_info = visit_data[['ART', 'Age', 'YearsSeropositive', 'Gender']].groupby(level = 'Patient ID').first()
fit_data = merge(out_data, pat_info, 
                 left_on = 'Patient ID',
                 right_index = True, 
                 how = 'left').dropna()
fit_data

# <codecell>

import statsmodels.api as sm
from patsy import dmatrices, dmatrix

eqn = 'HIVD ~ C(Visit, Treatment(0))+Age+Gender+YearsSeropositive+C(Grouping, Treatment("NonUser"))+C(ART, Treatment("naive"))'
y, X = dmatrices(eqn, fit_data, return_type = 'dataframe')
res = sm.GLM(y,X).fit()
fit_data['iHIVD'] = res.fittedvalues
res.summary()


# <codecell>

fig, axs = plt.subplots(4,1, figsize = (4,10), sharex = True, sharey = True)
nplot_data = pivot_table(fit_data, 
                         rows = ['Grouping', 'Patient ID'],
                         cols = 'Visit',
                         values = 'iHIVD')
for ax, group in zip(axs.flatten(), groups):
    medians = nplot_data.ix[group].median()
    ax.boxplot(nplot_data.ix[group].fillna(medians).values)
    
    ax.set_title(group + ' N:%i' % len(plot_data.ix[group]))
    ax.set_ylabel('aHIVD')
    
ax.set_xlabel('Visit');
plt.savefig('figures/NeuroSpecificFigures/adjHIVD_grouped.png')

# <markdowncell>

# The main difficulty is our requirement of having an HIVD score at every visit, even allowing a missing visit doesn't improve our numbers. A hidden complication in the above analysis is that by using a 'visit centric' analysis we hide the complication that our visits are unequally spaced. To deal with both of these complications I'm proposing a 'mean sampling method'.

# <headingcell level=2>

# Mean Sampling Method

# <markdowncell>

# This sampling method uses the idea that our samples are unequally spaced but can be _regularized_ by taking an Exponentially Weighted Moving Average (EWMA). As demonstrated below on a few particularly variable patients.

# <codecell>

def do_ewma(inser):
    tser = inser.reset_index(level = 'Patient ID', drop = True).map(float)
    res = ewma(tser,freq = '6M', span = 3)
    return res

def do_mean(inser):
    tser = inser.reset_index(level = 'Patient ID', drop = True).map(float)
    res = rolling_mean(tser,2, freq = '12M')
    return res

hivd_data = date_data.dropna(subset = ['VisitNum'])['HIVD']

ewma_smoothed_HIVD = hivd_data.groupby(level = 'Patient ID').apply(do_ewma)
ewma_smoothed_HIVD.name = 'sHIVD'
rolling_smoothed_HIVD = hivd_data.groupby(level = 'Patient ID').apply(do_mean)
rolling_smoothed_HIVD.name = 'sHIVD'

example_pats = ['A0014', 'A0132', 'A0110', 'A0062']
fig, axs = plt.subplots(4,1, sharey=True, figsize = (10,8))
for ax, pat in zip(axs.flatten(), example_pats):
    plt.sca(ax)
    pat_ewma_smooth = ewma_smoothed_HIVD.ix[pat].reset_index()
    pat_rolling_smooth = rolling_smoothed_HIVD.ix[pat].reset_index()
    pat_raw = date_data.dropna(subset = ['VisitNum']).ix[pat]['HIVD'].reset_index()
    
    plt.plot_date(pat_ewma_smooth['Date'], pat_ewma_smooth['sHIVD'], 'r.-', 
                    label = 'EWMA', lw = 5, alpha = 0.7)
    plt.plot_date(pat_rolling_smooth['Date'], pat_rolling_smooth['sHIVD'], 'g.-', 
                    label = 'Rolling', lw = 5, alpha = 0.7)
    plt.plot_date(pat_raw['Date'], pat_raw['HIVD'], 'b.-', 
                    label = 'Raw', lw = 5, alpha = 0.7)
    plt.title(pat)
    plt.ylim([0, 13])
    plt.ylabel('HIVD')
plt.legend(loc = 'best');
plt.tight_layout()

# <markdowncell>

# Here the idea is to use a sliding window average in which nearby samples are weighted more then distant samples. This has particular attractiveness in our situtation due to the unequal spacing of our sampling. A simple rolling mean considers all points within the window around the sampling point to be equal, even if one is 2 weeks away and the other two are more then 9 months in either direction. The EWMA takes this distiction into account.
# 
# In this case I'm sampling every 6 months using an "18-month Exponetially Weighted Sampling". This means the measurements up to 18-months from the sample contribute 50% to the average and it decreases exponentially after that.
# 
# Now I don't have to worry about the unequal sampling or missing values. I'm going to continue with the wLCM.
# 
# For Niz's thesis/paper I mae the same figure but only included the Raw and EWMA lines. It seems less confusing.

# <codecell>

example_pats = ['A0014', 'A0132', 'A0110', 'A0062']
fig, axs = plt.subplots(4,1, sharey=True, figsize = (10,8))
for ax, pat in zip(axs.flatten(), example_pats):
    plt.sca(ax)
    pat_ewma_smooth = ewma_smoothed_HIVD.ix[pat].reset_index()
    pat_rolling_smooth = rolling_smoothed_HIVD.ix[pat].reset_index()
    pat_raw = date_data.dropna(subset = ['VisitNum']).ix[pat]['HIVD'].reset_index()
    
    plt.plot_date(pat_ewma_smooth['Date'], pat_ewma_smooth['sHIVD'], 'r.-', 
                    label = 'EWMA', lw = 5, alpha = 0.7)
    plt.plot_date(pat_raw['Date'], pat_raw['HIVD'], 'b.-', 
                    label = 'Raw', lw = 5, alpha = 0.7)
    plt.title(pat)
    plt.ylim([0, 13])
    plt.ylabel('HIVD')
plt.tight_layout()
plt.savefig('figures/NeuroSpecificFigures/HIVDsampling.png', dpi = 300)
plt.close()

# <headingcell level=2>

# wLCM

# <markdowncell>

# Here I'm going to _summarize_ patients based on thier Admit, Test and Admit-or-Test data. So I will ultimately get one drug value for each patient. Its important to note that I'm *not* modeling the dynamics of drug use. I'm assuming that the drug-testing gives an idea of the 'proportion' of time the patient has a drug in their system. I'm specifically excluding patients that Admit-or-Test to any drug other than Cannabinoids or Cocaine. I'm also limiting the analysis to Black/AA patients who are on ART for the entire timeperiod.

# <codecell>

def safe_mean(inser):
    return inser.mean()
    
exclude_drugs = [ 'AoT-Other', 'AoT-Phencyclidine', 
                 'AoT-Amphetamines', 'AoT-Benzodiazepines', 'AoT-Opiates']
check_drugs = ['Test-Cocaine', 'Test-Cannabinoid']

aa_pats = visit_data['Race-Black'].groupby(level = 'Patient ID').all()
non_trans_pats = (visit_data['Gender'] != 'Transgender').groupby(level = 'Patient ID').all()
art_patients = (visit_data['ART']=='on').groupby(level = 'Patient ID').all()
odrug_valid = (visit_data[exclude_drugs]==0).all(axis = 1).groupby(level = 'Patient ID').all()
drug_fracs = visit_data[check_drugs].groupby(level = 'Patient ID').agg(safe_mean)
visit_counts = visit_data['HIVD'].notnull().groupby(level = 'Patient ID').sum()

wanted_mask = aa_pats & non_trans_pats & odrug_valid & (visit_counts >= 3)
wanted_fracs = drug_fracs[wanted_mask]

hivd_visits = merge(visit_data[['ART', 'Age', 'YearsSeropositive', 'Gender', 'Date', 'HIVD']].reset_index(),
                    wanted_fracs,
                    left_on = 'Patient ID',
                    right_index = True,
                    how = 'right').dropna()
hivd_visits['HIVD'] = hivd_visits['HIVD'].map(float)
hivd_visits['NumVisit'] = hivd_visits['VisitNum'].map(lambda x: float(x[1:]))
#hivd_visits['ART'] = hivd_visits['ART'].replace('non-adherent', 'off').replace('naive', 'off')
hivd_visits.rename(columns = {
                              'Test-Cocaine':'Cocaine', 
                              'Test-Cannabinoid':'Cannabinoid'
                              }, inplace = True)

hivd_visits

# <codecell>

eqn = 'HIVD ~ NumVisit+Age+Gender+YearsSeropositive+ART+Cocaine+Cannabinoid'
y, X = dmatrices(eqn, hivd_visits, return_type = 'dataframe')
#print y.head()
#print X.head()
res = sm.GLM(y,X).fit()
hivd_visits['aHIVD'] = res.fittedvalues
res.summary()

# <codecell>

DataFrame({'pvalues':res.pvalues, 'effectsize':res.params}).to_excel('tables/HIVD_effects.xlsx')

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

fig, axs = plt.subplots(1,3, sharey=True, figsize = (15, 5))
groups = [('YearsSeropositive', np.arange(0,30,5), 5),
          ('ART', ['naive', 'on', 'non-adherent', 'off'], None),
          ('Age', np.arange(20,70, 10), 10),]
for ax, (group, bins, width) in zip(axs.flatten(), groups):
    
    tdata = hivd_visits[[group, 'aHIVD']].dropna()
    
   
    if type(bins[0]) != StringType:
        pos = Series(map(get_mid, cut(tdata[group], bins)), 
                index = tdata.index)
        tdata['Pos'] = pos
        tdata = tdata.dropna().sort('Pos')
        m, b, r, p, _ = linregress(tdata[group].values, y = tdata['aHIVD'].values)
        x = np.linspace(bins[0]-width, bins[-1]+width)
        y = m*x + b
    else:
        tdata['Pos'] = tdata[group]
    #tdata.boxplot(column = 'cHIVD', by = 'Pos', ax = ax)
    tmp = []
    for b in bins:
        tmp.append(tdata['aHIVD'][tdata['Pos'] == b].values)
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
    plt.ylim([0, 13])
    plt.ylabel('Adj-HIVD')
    plt.hold(False)
    
plt.savefig('figures/NeuroSpecificFigures/HIVD_trends.png')

# <codecell>

def new_ewma(indf):
    t_df = indf[['Date', 'HIVD', 'aHIVD', 'Cocaine', 'Cannabinoid']].sort('Date')
    t_df.set_index('Date', inplace=True)
    res = ewma(t_df,freq = '6M', span = 3).reset_index(drop = True)
    res.index.name = 'StdVisit'
    res.rename(columns = {'Cocaine':'I-Cocaine', 'Cannabinoid':'I-Cannabinoid'}, 
               inplace = True)
    res['Cocaine'] = indf['Cocaine'].mean()
    res['Cannabinoid'] = indf['Cannabinoid'].mean()
    
    return res

smoothed_data = hivd_visits.groupby('Patient ID').apply(new_ewma)
print smoothed_data.head().to_string()
sum_data = smoothed_data[['Cannabinoid', 'Cocaine']].groupby(level = 'Patient ID').mean()
print 'Only-Cocaine', ((sum_data['Cannabinoid'] == 0) & (sum_data['Cocaine'] > 0)).sum()
print 'Only-Cannabinoid', ((sum_data['Cannabinoid'] > 0) & (sum_data['Cocaine'] == 0)).sum()
print 'Non-User', ((sum_data['Cannabinoid'] == 0) & (sum_data['Cocaine'] == 0)).sum()
print 'Multi-User', ((sum_data['Cannabinoid'] > 0) & (sum_data['Cocaine'] > 0)).sum()

# <codecell>

from mpl_toolkits.mplot3d import Axes3D
from itertools import islice

fig, axs = plt.subplots(1, 3, figsize = (15,5), subplot_kw={'projection':'3d'})

check_cols = [['Cannabinoid'], ['Cocaine'], ['Cannabinoid', 'Cocaine']]

for check_col, ax in zip(check_cols, axs.flatten()):
    
    ax.hold(True)
    ax.set_title('+'.join(check_col))
    ax.set_zlabel('HIVD')
    ax.set_xlabel('Drug-Use')
    ax.set_ylabel('StdVisit')
    ax.set_xlim([0,1])
    
    for pat, row in smoothed_data.groupby(level = 'Patient ID'):
        nrow = row.reset_index()
        drug_per = nrow[check_col].mean(axis = 1)
        df = DataFrame({
                        'sHIVD':nrow['aHIVD'], 
                        'StdVisit':nrow['StdVisit'],
                        'DrugPerc':drug_per
                    }).dropna()
        if len(df)>2:
            
            ax.plot(df['DrugPerc'], df['StdVisit'], df['sHIVD'], 
                    lw = 5, color = [0, 0, drug_per.mean()], alpha = 0.7)
            
fig.tight_layout()
plt.savefig('figures/NeuroSpecificFigures/threeD_cocaine_figure.png', dpi = 300)

# <codecell>

from pylab import get_cmap

def make_crazy_plots(prop, ntitle, cuts, limit):

    tmp = pivot_table(smoothed_data.reset_index(),
                      rows = 'Patient ID',
                      cols = 'StdVisit',
                      values = 'HIVD')
    
    counts = tmp.applymap(np.isnan).sum(axis = 1)
    order = DataFrame({'prop':prop, 'counts':counts})
    order.sort(['prop', 'counts'], inplace = True)
    
    diffs = tmp.T.pct_change(periods=1).T
        
    fig, axs = plt.subplots(2,3, figsize = (10,6))
    for (title, (min_v, max_v)), taxs in zip(cuts, axs.transpose()):
        wpats = (order['prop'] >= min_v) & (order['prop'] <= max_v)
        avg_ax = taxs[1]
        img_ax = taxs[0]
        avg_ax.hold(True)
        img_ax.hold(True)
        pdata, _ = tmp.align(wpats[wpats], axis = 0, join='right')
        img_ax.imshow(pdata, 
                      interpolation='nearest', 
                      aspect = 'auto', 
                      cmap = get_cmap('Reds'),
                      vmin = 0, vmax = 12)
        img_ax.set_xlim([-0.1, limit+0.1])
        img_ax.set_xticklabels([])
        img_ax.set_yticklabels([])
    
        img_ax.set_title(title + ' %s N=%i' % (ntitle, wpats.sum()))
        plt.sca(img_ax)
        
        avg_score = pdata.mean()
        std_score = pdata.std()
        avg_ax.errorbar(np.arange(12), avg_score.values, yerr=std_score.values, lw = 2)
        
        if title == '1-99%':
            
            img_ax.yaxis.tick_right()
            img_ax.set_yticks(np.arange(len(pdata)+0.5))
            img_ax.set_yticklabels(order['prop'][wpats].values)
            img_ax.set_ylim([0.5, wpats.sum()-0.5])
            
            drg_ax = avg_ax.twinx()
            ntmp = pivot_table(smoothed_data.reset_index(),
                      rows = 'Patient ID',
                      cols = 'StdVisit',
                      values = 'I-'+ntitle)
            
            dpdata, _ = ntmp.align(wpats[wpats], axis = 0, join='right')
            mean_use = dpdata.mean()
            #print dpdata.head()
            #raise KeyboardInterrupt
            std_use = dpdata.std()
            errorbar(np.arange(12), 
                     mean_use.values,
                     lw = 2, color = 'g', alpha = 0.7)
            drg_ax.set_ylim([.25,.50])
            #rint ntmp
            #aise KeyError
        avg_ax.set_xlim([-0.1, limit+0.1])
        avg_ax.set_ylim([0, 13])
        avg_ax.set_xticks(range(0,limit+1))
        avg_ax.set_xticklabels([6*x for x in range(0,limit+1)])
        avg_ax.set_xlabel('Months')
        
            
            
        
    axs[0,0].set_ylabel('Patient')
    axs[1,0].set_ylabel('HIVD')

cuts = [('0%', (0,0)), 
        ('100%', (1, 1)),
        ('1-99%', (0.01, 0.99))]
mud_cuts = [('0%', (0,0)), 
        ('1-75%', (0.01, 0.75)), 
        ('76-100%', (0.76, 1))]
    
plots = [('Cannabinoid', smoothed_data[['Cannabinoid']].mean(axis = 1).groupby(level = 'Patient ID').max(), cuts),
         ('Cocaine', smoothed_data[['Cocaine']].mean(axis = 1).groupby(level = 'Patient ID').max(), cuts),
         ('MDU', smoothed_data[['Cannabinoid', 'Cocaine']].mean(axis = 1).groupby(level = 'Patient ID').max(), mud_cuts),
         ]
for (name, mask, ncuts), limit in product(plots, [6, 10]):
    make_crazy_plots(mask, name, ncuts, limit)
    plt.tight_layout()
    plt.savefig('figures/NeuroSpecificFigures/%s_%imonths_neuro_time_course.png' % (name, limit), dpi = 300)
    plt.close()

plt.figure(figsize = (10, 10))
plt.imshow(np.tile(np.arange(0, 13), (5, 1)), 
           interpolation='nearest', 
           aspect = 'auto', 
           cmap = get_cmap('Reds'),
           vmin = 0, vmax = 12)
plt.colorbar()
plt.savefig('figures/NeuroSpecificFigures/neuro_time_course_colorbar.png', dpi = 300)
plt.close()

# <codecell>

cyto_data = read_csv('CytoRawData.csv', sep='\t')
mcyto = cyto_data.groupby(['Patient ID']).mean().drop(['SampleNum'], axis=1)
cytos = mcyto.columns

# <codecell>

def safe_mean(ser):
    return np.mean(ser)

pat_cols = [('max-VL', lambda x:x['VL'].groupby(level='Patient ID').max()),
            ('min-VL', lambda x:x['VL'].groupby(level='Patient ID').min()),
            ('mean-VL', lambda x:x['VL'].groupby(level='Patient ID').agg(safe_mean)),
            ('max-CD4', lambda x:x['CD4'].groupby(level='Patient ID').max()),
            ('min-CD4', lambda x:x['CD4'].groupby(level='Patient ID').min()),
            ('mean-CD4', lambda x:x['CD4'].groupby(level='Patient ID').agg(safe_mean))]

extra_data = {}
for name, func in pat_cols:
    extra_data[name] = func(visit_data)
extra_pat_data = DataFrame(extra_data)

npat_data = concat(extra_pat_data.align(mcyto, axis=0, level='Patient ID', join='right'), axis=1)
check_cols = npat_data.columns

# <codecell>

from scipy.stats import ttest_ind

def nonnan_count(ser):
    return ser.notnull().sum()

cyto_hivd_data = concat(smoothed_data.align(npat_data, join='inner', axis=0, level='Patient ID'), axis=1)

extra_data = []
for cyto in check_cols:
    plt.figure(figsize = (5,5))
    plt.hold(True)
    col = 'aHIVD'
    med_val = cyto_hivd_data[cyto].groupby(level='Patient ID').max().median()
    no_coke_mask = cyto_hivd_data['Cocaine'] == 0
    high_pats = (cyto_hivd_data[cyto] >= med_val) & no_coke_mask
    low_pats = (cyto_hivd_data[cyto] < med_val) & no_coke_mask
    wpats = high_pats | low_pats

    diff_scoring = cyto_hivd_data[[col]][wpats].reset_index()
    plot_data = pivot_table(diff_scoring, 
                            rows='Patient ID', 
                            cols='StdVisit', 
                            values=col)
    plot_data['ishigh'] = high_pats.groupby(level='Patient ID').any()
    
    low_data = cyto_hivd_data[[col]][low_pats].dropna()
    high_data = cyto_hivd_data[[col]][high_pats].dropna()
    _, pval = ttest_ind(low_data.values, high_data.values)
    extra_data.append((cyto, 'pval', pval[0]))
    extra_data.append((cyto, 'low_count', low_pats.groupby(level='Patient ID').all().sum()))
    extra_data.append((cyto, 'high_count', high_pats.groupby(level='Patient ID').all().sum()))
    extra_data.append((cyto, 'cutoff', med_val))
    
    
    plots = plot_data.groupby('ishigh').agg(['mean', 'std', nonnan_count]).swaplevel(1,0, axis=1)
    order = [(False, 'r'), 
             (True, 'b')]
    for num, (row, color) in enumerate(order):
        mask = plots['nonnan_count'].ix[row] >= 5
        pdata = plots['mean'].ix[row][mask]
        edata = plots['std'].ix[row][mask]
        plt.errorbar(np.array(pdata.index) + 0.2*num, 
                    pdata.values, 
                    yerr = edata.values, 
                    color = color)
        
    plt.title(cyto + ' p: %0.2e' % pval)
    plt.xlabel('StdVisits')
    plt.ylabel('aHIVD')
    plt.savefig('figures/NeuroSpecificFigures/cytosplits/%s.png' % cyto.replace('.', ''), 
                dpi = 300)
    plt.close()
    

# <codecell>

tdf = DataFrame(extra_data, columns = ['Cyto', 'Type', 'Val'])
out_pvals = pivot_table(tdf, rows = 'Cyto', cols = 'Type', values='Val')
out_pvals.to_excel('tables/neuro_cytokine_split_pvals.xlsx')

# <headingcell level=2>

# Predicting NeuroCog in Time-Chunks

# <codecell>

def roll_slopes(indf):
    
    slopes = []
    for a, b in zip(range(1, len(indf)), range(len(indf))):
        slope = get_slopes(indf.iloc[[b,a]])
        slopes.append(slope)
    slopes.append(np.nan)
    out = Series(slopes, 
                 index=indf.reset_index()['VisitNum'], 
                 name='HIVD-m')
    
    return out
    

def get_slopes(indf):
    
    ts = indf['Date'].iloc[1] - indf['Date'].iloc[0]
    ns = indf['HIVD'].iloc[1] - indf['HIVD'].iloc[0]
    return ns/ts.days
    

hivd_slopes = visit_data[['Date', 'HIVD']].dropna().groupby(level='Patient ID').apply(roll_slopes)
pred_data = merge(DataFrame(hivd_slopes), visit_data,
                  left_index = True,
                  right_index = True)
pred_data['Gender-M'] = pred_data['Gender'] == 'M'

# <codecell>

binary_cols = ['ART-naive', 'ART-nonadherent', 'ART-off', 'ART-on', 
               'Admit-Alcohol', 'Admit-Amphetamines', 'Admit-Benzodiazepines', 
               'Admit-Cannabinoid', 'Admit-Cocaine', 'Admit-Opiates', 'Admit-Other', 
               'Admit-Phencyclidine', 'Admit-Ritalin', 
               'AoT-Alcohol', 'AoT-Amphetamines', 'AoT-Benzodiazepines', 
               'AoT-Cannabinoid', 'AoT-Cocaine', 'AoT-Opiates', 'AoT-Other', 
               'AoT-Phencyclidine', 'AoT-Ritalin', 'Gender-M', 
               'Race-Asian', 'Race-Black', 'Race-Hawaiian', 'Race-Indian', 
               'Race-Multi-Race', 'Race-Unknown', 'Race-White', 
               'Test-Amphetamines', 'Test-Barbiturates', 'Test-Benzodiazepines', 
               'Test-Cannabinoid', 'Test-Cocaine', 'Test-Opiates', 'Test-Phencyclidine']
               
cont_cols = ['Age', 'CD4', 'CD8', 'HIVD', 'VL', 'YearsSeropositive']
bset = set(binary_cols)
ccols = set(cont_cols)
check_cols = binary_cols+cont_cols

pvals = []
for col in check_cols:
    #print col
    tmp = pred_data[['HIVD-m', col]].dropna()
    if col in ccols:
        tmp[col] = tmp[col] > tmp[col].median()
    try:
        _, pval = ttest_ind(tmp['HIVD-m'][tmp[col]==True].values, tmp['HIVD-m'][tmp[col]==False].values)
        pvals.append((np.log10(pval), col))
    except ZeroDivisionError:
        pass
    
    
pvals.sort()
print pvals[:10]
order = [col for _, col in pvals]

# <codecell>

def my_score(xt, xp):
    score = mean_absolute_error(xt, xp)
    return -score

# <codecell>

from sklearn.metrics import mean_absolute_error, r2_score, mean_squared_error, zero_one_loss
from sklearn.svm import SVC
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier, ExtraTreesClassifier
from sklearn.dummy import DummyRegressor, DummyClassifier
from sklearn.neighbors import KNeighborsRegressor, KNeighborsClassifier
from sklearn.feature_selection import SelectKBest, f_regression
from sklearn.pipeline import Pipeline
from sklearn.grid_search import GridSearchCV, IterGrid
from sklearn.cross_validation import LeavePLabelOut, KFold, LeaveOneOut, cross_val_score, StratifiedKFold
from itertools import islice

models = [('Naive', DummyClassifier(), {}),
          ('GradBoost', GradientBoostingClassifier(), {}),
          ('DecisionTree', DecisionTreeClassifier(), {'model__criterion':['gini', 'entropy']}),

          ('RandomForest', RandomForestClassifier(), {
                                                      'model__criterion':['gini', 'entropy'],
                                                      'model__n_estimators':[10, 20, 50]
                                                      }),
          ('ExtraTrees', ExtraTreesClassifier(), {
                                                  'model__criterion':['gini', 'entropy'],
                                                  'model__n_estimators':[10, 20, 50]
                                                  }),
          ('KNN', KNeighborsClassifier(warn_on_equidistant=False,
                                       weights='distance'), {
                                                             'model__n_neighbors':[5,10,20]
                                                             }),
          ('SVR', SVC(), {'model__kernel':['linear', 'rbf']})]

num = 50
ncols = order[:num]
tmp = pred_data[['HIVD-m']+ncols].dropna()

X = tmp[ncols].values
y = tmp['HIVD-m'].values < 0
num_select = 20
default_params = {'select__k':range(5, 40, 10)}
for name, model, extra_params in models:
    
    pipe = Pipeline([('select', SelectKBest(k = num_select, 
                                            score_func = f_regression)),
                     ('model', model)])
    extra_params.update(default_params)
    print name, len(list(IterGrid(extra_params)))
    grid = GridSearchCV(pipe, extra_params,
                        cv = StratifiedKFold(y, n_folds=2),
                        n_jobs = 50)
    grid.fit(X,y)
    scores = cross_val_score(grid.best_estimator_, X, y, 
                             score_func=zero_one_loss, 
                             cv = LeaveOneOut(len(y)),
                             n_jobs = 50)
    
    print name, np.mean(scores), grid.best_params_

# <codecell>


# <codecell>


# <codecell>

from subprocess import call
import shlex

cmd = 'rsync -rvh /home/will/HIVSystemsBio/NewCytokineAnalysis/ /home/will/Dropbox/Wigdahl\ HIV\ Lab/NewCytokineAnalysis/'
call(shlex.split(cmd))

# <codecell>


