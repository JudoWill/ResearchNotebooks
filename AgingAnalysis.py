# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=1>

# Aging Analysis of Cytokine and SNP data

# <markdowncell>

# This analysis will look at how aging effects HIV-1 disease progression. This will include looking at things like clinical parameters, LTR SNPs, Cytokine Profiling, and NeuroCog Impairment.

# <headingcell level=2>

# Data Extraction

# <codecell>

from __future__ import division
import os, os.path
import numpy as np
import pandas as pd
from patsy import dmatrices
from patsy.contrasts import Treatment
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

sys.path.append('/home/will/PySeqUtils/')
sys.path.append('/home/will/PatientPicker/')
os.chdir('/home/will/AgingAnalysis/')

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
    rser = pd.Series([race]*len(indf), index = vnums)
    return rser

def convert_haart(inval):
    tdict = {'on': 'cH',
             'non-adherent': 'dH',
             'off': 'dH',
             'naive': 'nH'}
    return tdict.get(inval, np.nan)

pat_groups = redcap_data.groupby(level = 0)



cols = {'Age':redcap_data.groupby(level = 0)['Age'].transform('min'),
        'NumTotalVisits': redcap_data.groupby(level = 0)['Age'].transform(len).map(safe_float),
        'CD4': redcap_data['Latest CD4 count (cells/uL)'].map(safe_float),
        'Alcohol': redcap_data['Current Alcohol Use'],
        'Tobacco': redcap_data['Current Tobacco Use'],
        'DaysSinceBaseline': pat_groups['Date Of Visit'].transform(get_days_since_bl),
        'Gender': redcap_data['Gender'],
        'Grouping': niz_groupings,
        'Race': pat_groups.apply(guess_race).map(lambda x: x.replace('-', '')),
        'HAART': redcap_data['Current ART status'].map(convert_haart),
        'HCV': pat_groups['Hepatitis C status (HCV)'].transform(pd.expanding_max),
        'HIVD': redcap_data['TMHDS'].map(safe_float),
        'HIVDI': redcap_data['TMHDS']<10,
        'HBV': pat_groups['Hepatitis B status (HBV)'].transform(pd.expanding_max).map(safe_float),
        'CD8': redcap_data['Latest CD8 count (cells/uL)'].map(safe_float),
        'NadirCD4': pat_groups['Nadir CD4 count (cells/uL)'].transform(pd.expanding_min).map(safe_float),
        'NadirCD8': pat_groups['Nadir CD8 count (cells/uL)'].transform(pd.expanding_min).map(safe_float),
        'PeakLVL': pat_groups['Peak viral load (copies/mL)'].transform(pd.expanding_max).map(safe_float).map(np.log10),
        'LVL': redcap_data['Latest viral load'].map(safe_float).map(np.log10),
        'YearsSeropositive': redcap_data['Years Seropositive'].map(safe_float),
        'TOSampleBenzodiazepines': pat_groups['Test-Benzodiazepine'].transform(pd.expanding_mean).map(safe_float),
        'TOSampleCannabinoid': pat_groups['Test-Cannabinoid'].transform(pd.expanding_mean).map(safe_float),
        'TOSampleCocaine': pat_groups['Test-Cocaine'].transform(pd.expanding_mean).map(safe_float),
        'TOSampleOpiates': pat_groups['Test-Opiates'].transform(pd.expanding_mean).map(safe_float),
        'ALLBenzodiazepines': pat_groups['Test-Benzodiazepine'].transform('mean').map(safe_float),
        'ALLCannabinoid': pat_groups['Test-Cannabinoid'].transform('mean').map(safe_float),
        'ALLCocaine': pat_groups['Test-Cocaine'].transform('mean').map(safe_float),
        'ALLOpiates': pat_groups['Test-Opiates'].transform('mean').map(safe_float),
        'ATSampleBenzodiazepines': redcap_data['Test-Benzodiazepine'].map(safe_float),
        'ATSampleCannabinoid': redcap_data['Test-Cannabinoid'].map(safe_float),
        'ATSampleCocaine': redcap_data['Test-Cocaine'].map(safe_float),
        'ATSampleOpiates': redcap_data['Test-Opiates'].map(safe_float)}

known_pat_data = pd.DataFrame(cols)

# <codecell>

import Rtools
cytos = sorted(['IL.8','VEGF','IL.1beta',
        'G.CSF','EGF','IL.10','HGF',
        'FGF.basic','IFN.alpha','IL.6',
        'IL.12','Rantes','Eotaxin',
        'GM.CSF','MIP.1beta',
        'MCP.1','IL.5','IL.13', 'IFN.gamma','TNF.alpha',
        'IL.RA','IL.2','IL.7','IP.10',
        'IL.2R','MIG','IL.4','IL.15',
        'IL.17','MIP.1alpha']) + ['Th1', 'Th2']

raw_cyto_data = pd.read_csv('/home/will/HIVSystemsBio/NewCytokineAnalysis/CytoRawData.csv', 
                            sep = '\t').groupby(['Patient ID', 'VisitNum', 'SampleNum']).first().applymap(safe_float)

raw_cyto_data['Th1'] = raw_cyto_data[['IFN.gamma', 'IL.2', 'TNF.alpha']].sum(axis=1)
raw_cyto_data['Th2'] = raw_cyto_data[['IL.4', 'IL.5', 'IL.10']].sum(axis=1)
#norm_cyto_data = Rtools.quantile_norm_with_R(raw_cyto_data[cytos].applymap(safe_float))
raw_cyto_data['Th1Th2'] = raw_cyto_data['Th1']/raw_cyto_data['Th2']
raw_cyto_data = raw_cyto_data.groupby(level=[0,1]).median()

# <codecell>

from functools import partial

known_pat_data['CD4CD8'] = known_pat_data['CD4'] / known_pat_data['CD8']

known_pat_data['HIVHealthy'] = (known_pat_data['LVL'] <= 2) & \
                                (known_pat_data['CD4'] >= 250)
known_pat_data['OlderThan50'] = known_pat_data['Age']>50
known_pat_data['OlderThan60'] = known_pat_data['Age']>60
known_pat_data['AgePast50'] = (known_pat_data['Age'] - 50).map(partial(max, 0))
known_pat_data['AgePast60'] = (known_pat_data['Age'] - 60).map(partial(max, 0))

# <codecell>

anal_data = pd.merge(raw_cyto_data.reset_index(), known_pat_data.reset_index(),
                     on = ['Patient ID', 'VisitNum']).set_index(['Patient ID', 'VisitNum'])

# <codecell>

anal_data['YearsSeropositive'].describe()

# <codecell>

pd.pivot_table(anal_data, rows='OlderThan50', values=cytos+['HIVD']).T
               

# <headingcell level=2>

# Cytokine Stats

# <codecell>

from statsmodels.graphics.boxplots import beanplot
from itertools import chain, combinations
from statsmodels.graphics.regressionplots import plot_fit

def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))

def unique_powerset(iterable):
    seen = set()
    for tup in powerset(iterable):
        fz = frozenset(tup)
        if fz not in seen:
            yield tup
            seen.add(fz)

def most_common(inser):
    try:
        return inser.value_counts().index[0]
    except IndexError:
        return np.nan

def make_agg(indata, confounders):
    agg_dict = {}
    for col in confounders:
        if indata[col].dtype == 'O':
            agg_dict[col] = most_common
        else:
            agg_dict[col] = 'mean'
    return agg_dict

def extract_cols(indata, cyto, cols):
    
    extract_cols = [cyto]+cols
    tmp_data = indata[extract_cols]
    agg_dict = make_agg(indata, extract_cols)
    tmp_agg = tmp_data.groupby(level=[0,1]).agg(agg_dict).dropna()
    tmp_agg = tmp_agg.rename(columns={cyto:'y'})
    eqn = 'y ~ ' + ' + '.join(cols)
    return eqn, tmp_agg
    


def make_50_binary(indata, cyto, ax, confounders):
    eqn, tmp_cyto = extract_cols(indata, cyto, ['OlderThan50']+confounders)
    y, X = dmatrices(eqn, tmp_cyto, return_type='dataframe')
    try:
        model = sm.OLS(y, X).fit()
    except:
        print indata
        print confounders
        print eqn
        print tmp_cyto
        print X
        print y
        raise TypeError
    
    if ax is not None:
        boxes = [tmp_cyto['y'][tmp_cyto['OlderThan50']], tmp_cyto['y'][~tmp_cyto['OlderThan50']]]
        try:
            beanplot(boxes, labels=['>50', '<50'], ax=ax)
        except:
            ax.boxplot(boxes)
    
    return model, model.f_pvalue, model.pvalues['OlderThan50[T.True]'], model.params['OlderThan50[T.True]']


def make_35_binary(indata, cyto, ax, confounders):
    
    valid_mask = (indata['Age'] >= 50) | (indata['Age'] <= 35)
    
    eqn, tmp_cyto = extract_cols(indata[valid_mask], cyto, 
                                     ['OlderThan50']+confounders)
    y, X = dmatrices(eqn, tmp_cyto, return_type='dataframe')
    model = sm.OLS(y, X).fit()
    boxes = [tmp_cyto['y'][tmp_cyto['OlderThan50']], tmp_cyto['y'][~tmp_cyto['OlderThan50']]]
    
    if ax is not None:
        try:
            beanplot(boxes, labels=['>50', '<35'], ax=ax)
        except:
            ax.boxplot(boxes)
    
    return model, model.f_pvalue, model.pvalues['OlderThan50[T.True]'], model.params['OlderThan50[T.True]']


def make_age_linear(indata, cyto, ax, confounders):
    eqn, tmp_cyto = extract_cols(indata, cyto, 
                                     ['Age']+confounders)
    y, X = dmatrices(eqn, tmp_cyto, return_type='dataframe')
    model = sm.OLS(y, X).fit()
    num = (num for num, col in enumerate(X.columns) if col=='Age').next()
    
    if ax is not None:
        plot_fit(model, num, ax=ax)
    
    return model, model.f_pvalue, model.pvalues['Age'], model.params['Age']


def check_more_data(indata, check_func, cyto, confounders, num_extras = [0.5, 1.0, 2.0, 5.0]):
    
    for num in num_extras:
        num_choose = int(num*len(indata))
        pats = [p for p, v in indata.index]
        tdata = indata.copy()
        for e in range(num_choose):
            rep_pat = np.random.choice(pats)
            pat = indata.ix[rep_pat]
            pat.index = pd.MultiIndex.from_tuples([(rep_pat+str(e), v) for v in pat.index])
            tdata = pd.concat([tdata, pat], axis=0)
        _, model_p, age_p, age_e = check_func(tdata, cyto, None, confounders)
        yield num_choose, model_p, age_p, age_e
        

# <headingcell level=3>

# Patient Groups

# <codecell>

all_selectors = set(['HCV', 'LVL', 'HAART', 'Race', 'Grouping'])
restrict_dict = {
                 'HCV':anal_data['HCV']==False,
                 'LVL':anal_data['LVL']<=2,
                 'HAART':anal_data['HAART']=='cH',
                 'Race':anal_data['Race']=='RaceBlack',
                 'Grouping':anal_data['Grouping']=='PN'
                 }
anal_methods = [('50|50 split', make_50_binary),
                ('35|50 split', make_35_binary),
                ('Linear', make_age_linear)]
checks = cytos+['CD4', 'LVL', 'HIVD']

num_subsets = len(list(unique_powerset(all_selectors)))
print num_subsets

# <codecell>

anal_data['Race']

# <codecell>

base_confounders = ['YearsSeropositive']


results = []
for cyto in checks:
    print cyto
    max_val = anal_data[cyto].max()*1.1
    for row_num, selectors in enumerate(unique_powerset(all_selectors)):
        print selectors
        fig, axs = plt.subplots(1, 3, figsize=(15,5), sharey=True)
        set_sels = set(selectors)
        confounders = base_confounders + list(all_selectors-set_sels)
        wanted_mask = anal_data[cyto].notnull()
        for sel in selectors:
            wanted_mask &= restrict_dict[sel]
        wanted_pats = anal_data[wanted_mask]
        sel_bool = [s in set_sels for s in sorted(all_selectors)]
        for ax, (pname, pfunc) in zip(axs.flatten(), anal_methods):
            model, model_p, age_p, age_e = pfunc(wanted_pats, cyto, ax, confounders)
            results.append(sel_bool+[cyto, pname, model_p, age_p, age_e, 0])
            ax.set_title(pname + ' p=%f' % age_p)
            if ax.is_first_col():
                ax.set_ylabel(cyto)
            if ax.is_last_col():
                ax.yaxis.set_label_position("right")
                ax.set_ylabel(', '.join(confounders))
            for num_extra, model_p, age_p, age_e in check_more_data(wanted_pats, pfunc, cyto, confounders):
                results.append(sel_bool+[cyto, pname, model_p, age_p, age_e, num_extra])
    
        if len(selectors):
            sels = '_'.join(sorted(selectors))
        else:
            sels = 'None'
        fname = cyto + '_' + sels + '.png'
        fig.savefig('draft_figures/aging_cyto_dump/'+fname)
        plt.close(fig)
        
    


# <codecell>

print model.params

# <codecell>

res = pd.DataFrame(results, columns = sorted(all_selectors)+['Cyto', 'AnalName', 'ModelP', 'AgeP', 'NumExtra'])
mask = (res['AgeP']<0.05) & (res['ModelP']<0.05) & (res['NumExtra'] == 0)
res[mask].groupby(['Cyto', 'AnalName', 'NumExtra'])['AgeP'].min()

# <codecell>


more_counts = pd.pivot_table(res[mask], rows='Cyto', cols='AnalName', 
                             values='NumExtra', aggfunc='min')
more_counts.to_excel('needed_pats.xlsx')

# <codecell>

import statsmodels.formula.api as smi
from scipy.stats import ttest_ind, ks_2samp
from patsy import dmatrices
from statsmodels.graphics.boxplots import beanplot, violinplot

fig = make_50_binary(wanted_cyto)
fig.savefig('draft_figures/Binary_50_cytos.png', dpi=500)

# <codecell>


fig = make_35_binary(wanted_cyto)
fig.savefig('draft_figures/Binary_35_cytos.png', dpi=500)
    

# <codecell>

fig = make_age_linear(wanted_cyto)
fig.savefig('draft_figures/Age_cytos.png', dpi=500)

# <codecell>

larger_wanted_pats = pat_cyto_data['HCV']==False
#larger_wanted_pats &= pat_cyto_data['VL']<=100
larger_wanted_pats &= pat_cyto_data['HAART']=='cH'
larger_wanted_pats &= pat_cyto_data['Race']=='Black/AA'
larger_wanted_pats &= (pat_cyto_data['Grouping']=='PN')|(pat_cyto_data['Grouping']=='PC')
print larger_wanted_pats.sum()/4
larger_wanted_cyto = pat_cyto_data[larger_wanted_pats]

# <codecell>


fig = make_50_binary(larger_wanted_cyto, confounders=['YearsSeropositive', 'Grouping'])
fig.savefig('draft_figures/Binary_50_cytos_withPC.png', dpi=500)

# <codecell>

fig = make_35_binary(larger_wanted_cyto, confounders=['YearsSeropositive', 'Grouping'])
fig.savefig('draft_figures/Binary_35_cytos_withPC.png', dpi=500)

# <codecell>

fig = make_age_linear(larger_wanted_cyto, confounders=['YearsSeropositive', 'Grouping'])
#fig.savefig('draft_figures/Age_cytos_withPC.png', dpi=500)

# <codecell>


