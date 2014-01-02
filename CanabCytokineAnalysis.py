# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=1>

# Canabinoid Cytokine Analysis

# <markdowncell>

# This analysis will look at the effect of Cannabinoid use on Cytokine Profiles. This effect is complicated due to the poly-abuse background.

# <headingcell level=2>

# Data Extraction

# <codecell>

from __future__ import division
import os, os.path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('/home/will/PySeqUtils/')
sys.path.append('/home/will/PatientPicker/')
os.chdir('/home/will/HIVSystemsBio/MoreCytokineAnalysis/')

# <codecell>

raw_cyto_data = pd.read_csv('../NewCytokineAnalysis/CytoRawData.csv', sep = '\t')

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
        return pd.Series(np.nan, index=inds)
    indf = indf.dropna()
    ever_used = indf.any(axis = 0)
    if not ever_used.any():
        return pd.Series('PN', index=inds)
    all_used = indf.all(axis = 0)
    if all_used.sum() == 1:
        return pd.Series(name_mappings[all_used.idxmax()], index=inds)
    elif all_used.sum() > 1:
        return pd.Series('MDU', index=inds)
    
    pure_cols = []
    non_pure_cols = []
    for col in indf.columns:
        if count_with_skips(indf[col], 1):
            pure_cols.append(col)
        else:
            non_pure_cols.append(col)
    if ever_used[non_pure_cols].any():
        return pd.Series(np.nan, index=inds)
            
    if len(pure_cols) == 1:
        return pd.Series(name_mappings[pure_cols[0]], index=inds)
    else:
        return pd.Series('MDU', index=inds)

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


def get_days_since_bl(ser):
    
    fdate = ser.dropna().min()
    diff_dates = (ser-fdate).apply(safe_days)
    return diff_dates


def guess_race(indf):
    race_cols = [col for col in indf.columns if col.startswith('Race-')]
    race = indf[race_cols].sum().idxmax()
    vnums = [v for _, v in indf.index]
    rser = pd.Series([race]*len(indf), index = vnums)
    return rser

def convert_haart(inval):
    tdict = {'on': 'cH',
             'non-adherent': 'dH',
             'off': 'dH',
             'naive': 'nH'}
    return tdict.get(inval, np.nan)


cols = {'Age':redcap_data.groupby(level = 0)['Age'].transform('min'),
        'NumTotalVisits': redcap_data.groupby(level = 0)['Age'].transform(len),
        'Latest CD4 count': redcap_data['Latest CD4 count (cells/uL)'],
        'Current Alcohol use': redcap_data['Current Alcohol Use'],
        'Current Tobacco use': redcap_data['Current Tobacco Use'],
        'Days since baseline': redcap_data.groupby(level = 0)['Date Of Visit'].transform(get_days_since_bl),
        'Gender': redcap_data['Gender'],
        'Grouping': niz_groupings,
        'Race': redcap_data.groupby(level=0).apply(guess_race),
        'HAART': redcap_data['Current ART status'].map(convert_haart),
        'Hepatitis C status (HCV)': redcap_data.groupby(level=0)['Hepatitis C status (HCV)'].transform(pd.expanding_max),
        'HIVD score': redcap_data['TMHDS'],
        'HIVD.I': redcap_data['TMHDS']<10,
        'Hepatitis B status (HBV)': redcap_data.groupby(level=0)['Hepatitis B status (HBV)'].transform(pd.expanding_max),
        'Latest CD8 count': redcap_data['Latest CD8 count (cells/uL)'],
        'Nadir CD4 count': redcap_data.groupby(level=0)['Nadir CD4 count (cells/uL)'].transform(pd.expanding_min),
        'Nadir CD8 count': redcap_data.groupby(level=0)['Nadir CD8 count (cells/uL)'].transform(pd.expanding_min),
        'Peak viral load': redcap_data.groupby(level=0)['Peak viral load (copies/mL)'].transform(pd.expanding_max),
        'Latest Viral load': redcap_data['Latest viral load'],
        'Years Seropositive': redcap_data['Years Seropositive'],
        'TOSample.Benzodiazepines': redcap_data.groupby(level=0)['Test-Benzodiazepine'].transform(pd.expanding_mean),
        'TOSample.Cannabinoid': redcap_data.groupby(level=0)['Test-Cannabinoid'].transform(pd.expanding_mean),
        'TOSample.Cocaine': redcap_data.groupby(level=0)['Test-Cocaine'].transform(pd.expanding_mean),
        'TOSample.Opiates': redcap_data.groupby(level=0)['Test-Opiates'].transform(pd.expanding_mean),
        'ALL.Benzodiazepines': redcap_data.groupby(level=0)['Test-Benzodiazepine'].transform('mean'),
        'ALL.Cannabinoid': redcap_data.groupby(level=0)['Test-Cannabinoid'].transform('mean'),
        'ALL.Cocaine': redcap_data.groupby(level=0)['Test-Cocaine'].transform('mean'),
        'ALL.Opiates': redcap_data.groupby(level=0)['Test-Opiates'].transform('mean'),
        'ATSample.Benzodiazepines': redcap_data['Test-Benzodiazepine'],
        'ATSample.Cannabinoid': redcap_data['Test-Cannabinoid'],
        'ATSample.Cocaine': redcap_data['Test-Cocaine'],
        'ATSample.Opiates': redcap_data['Test-Opiates']}

known_pat_data = pd.DataFrame(cols)

# <codecell>

known_pat_data

# <codecell>

list(redcap_data.columns)

# <codecell>


