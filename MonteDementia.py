# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

from __future__ import division
import numpy as np
import pandas as pd
import sys
import os
sys.path.append('/home/will/PatientPicker/')
import LoadingTools
from itertools import chain, islice
os.chdir('/home/will/Dropbox/MonteDementia/')

# <codecell>

redcap_data = LoadingTools.load_redcap_data()

# <codecell>

wcols = ['Age',
         'IsMale',
         'Race-Asian',
         'Race-Indian',
         'Race-Black',
         'Race-Hawaiian',
         'Race-White',
         'Race-Multiple',
         'Race-Unknown',       
         'Admit-Cannabinoid',
         'Admit-Cocaine',
         'Admit-Heroin',
         'Admit-Amphetamines',
         'Admit-Benzodiazapine',
         'Admit-Narcotics',
         'Admit-Ecstasy',
         'Admit-PCP',
         'Admit-Ritalin',
         'Admit-Other',
         'Admit-None',
         'Test-Amphetamines',
         'Test-Barbiturates',
         'Test-Benzodiazepine',
         'Test-Cannabinoid',
         'Test-Cocaine',
         'Test-Opiates',
         'Test-Phencyclidine',
         'Hepatitis C status (HCV)',
         'Diabetes',
         'Hypertension',
         'TMHDS',
         'HAART-Naive',
         'HAART-Non-Adherent',
         'HAART-Off',
         'HAART-On',
         'HAART-Missing']
date_col = 'Date Of Visit'
redcap_data['IsMale'] = redcap_data['Gender'] == 'Male'
other = [('CD4', 'Initial CD4 count (cells/uL)', 'Date of initial CD4 count'),
         ('CD4', 'Nadir CD4 count (cells/uL)', 'Date of nadir CD4 count'),
         ('CD4', 'Latest CD4 count (cells/uL)', 'Date of latest CD4 count'),
         ('CD8', 'Initial CD8 count (cells/uL)', 'Date of initial CD8 count'),
         ('CD8', 'Nadir CD8 count (cells/uL)', 'Date of nadir CD8 count'),
         ('CD8', 'Latest CD8 count (cells/uL)', 'Date of latest CD8 count'),
         ('VL', 'Initial viral load (copies/mL)', 'Date of initial viral load'),
         ('VL', 'Peak viral load (copies/mL)', 'Date of peak viral load'),
         ('VL', 'Latest viral load', 'Date of latest viral load')]

baseline_date = redcap_data.groupby('Patient ID')['Date Of Visit'].agg(lambda x: x.dropna().min()).to_dict()
reshaped_data = []
checks = zip(wcols, wcols,[date_col]*len(wcols))+other
for _, pat_row in redcap_data.iterrows():
    for out_col, check_col, date in checks:
        reshaped_data.append({
                              'Patient ID': pat_row['Patient ID'],
                              'Date': pat_row[date],
                              'Column': out_col,
                              'Value': float(pat_row[check_col])
                              })

# <codecell>

out_data = pd.pivot_table(pd.DataFrame(reshaped_data), 
                          rows = ['Patient ID', 'Date'], 
                          cols = 'Column',
                          values = 'Value',
                          aggfunc = 'mean')
out_data['VL'] = np.log10(out_data['VL'])
out_data['CD4'] = np.log10(out_data['CD4'])

# <codecell>

test_cols = ['Test-Amphetamines',
         'Test-Barbiturates',
         'Test-Benzodiazepine',
         'Test-Cannabinoid',
         'Test-Cocaine',
         'Test-Opiates',
         'Test-Phencyclidine']

def get_pure(indf):
    
    is_pure = []
    for col in test_cols:
        if indf[col].all():
            is_pure.append(col)
    if len(is_pure) == 1:
        return is_pure[0]
    elif len(is_pure) == 2:
        return 'MDU'
    else:
        vals = indf.values.flatten()
        good_vals = vals[~np.isnan(vals)]
        if np.all(good_vals == 0):
            return 'PN'
        else:
            return ''



groups = out_data[test_cols].groupby(level = 0).apply(get_pure)

# <codecell>

boxes

# <codecell>

from statsmodels.graphics.api import violinplot
gdict = groups.to_dict()

out_data['Grouping'] = [gdict[pid] for pid, _ in out_data.index]
check_col = 'VL'
wanted = ['PN', 'Test-Cocaine', 'Test-Cannabinoid', 'MDU']
boxes = []
labels = []
for group in wanted:
    df = out_data[out_data['Grouping'] == group]
    if group in wanted:
        labels.append(group)
        vals = df[check_col].dropna().values
        boxes.append(vals[np.isfinite(vals)])

haart_mask = out_data['HAART-On'] >= 0.5
df = out_data[haart_mask]
labels.append('cH')
vals = df[check_col].dropna().values
boxes.append(vals[np.isfinite(vals)])

df = out_data[~haart_mask]
labels.append('dH')
vals = df[check_col].dropna().values
boxes.append(vals[np.isfinite(vals)])




plt.figure(figsize = (10, 10))
ax = plt.subplot(111)
_=violinplot(boxes, labels = labels, ax = ax)
ax.set_ylim([0, 8])
plt.savefig('VLfigure.png', dpi = 500)

# <codecell>

from functools import partial

ewma_func = partial(pd.ewma, span=3, freq = '6M')


out = []
for pid, indf in out_data.reset_index().groupby('Patient ID'):
    ndf = indf.set_index('Date')
    tmp = ewma_func(ndf)
    tmp['Patient ID'] = pid
    tmp['pTMHDS'] = tmp['TMHDS'].shift(1)
    out.append(tmp.reset_index().copy())
    
tout = pd.concat(out, axis = 0, ignore_index=True)
tout['IntakeDate'] = tout['Patient ID'].map(lambda x: baseline_date[x])
tout['DaysSinceBL'] = (tout['Date']-tout['IntakeDate']).map(lambda x: x.astype('timedelta64[D]').astype(int))/365
trans_data = tout[tout['DaysSinceBL'] >= 0]

# <codecell>

plot_cols = ['VL', 'CD4', 'TMHDS']
check_pats = ['A0073', 'A0107', 'A0041']
trans_pat_data = trans_data.groupby(['Patient ID', 'Date']).first()
xlims = {}
xticks = {}
for pat in check_pats:
    fdate = baseline_date[pat]
    sampled = out_data.ix[pat].truncate(before=fdate)
    estimated = trans_pat_data.ix[pat].truncate(before=fdate)
    if len(sampled) < 5:
        continue
    fig, axs = plt.subplots(3,1, sharex=True, figsize=(10,10))
    for col, ax in zip(plot_cols, axs.flatten()):
        ts = sampled[col].dropna()
        es = estimated[col].dropna()
        ax.plot_date(ts.index, ts, color = 'g', marker = 'd', linestyle = '-', markersize = 20, linewidth = 5, alpha = 0.8)
        ax.plot_date(es.index, es, color = 'b', marker = '.', linestyle = '-', markersize = 25, linewidth = 5, alpha = 0.8)
        if col == 'TMHDS':
            ax.set_ylim([0, 12.5])
            ax.set_ylabel(col)
        elif col == 'CD4':
            ax.set_ylim([2, 3])
            ax.set_ylabel('log(' + col + ')')
        elif col == 'VL':
            ax.set_ylabel('log(' + col + ')')
            ax.set_ylim([1, 4])
        
        if ax.is_first_row():
            ax.set_title(pat)
    plt.tight_layout()
    xticks[pat] = ax.get_xticks()
    xlims[pat] = ax.get_xlim()
    plt.savefig('sample_fig_%s.png' % pat, dpi = 500)
    plt.close()
    

# <codecell>

import statsmodels.api as sm

admit_cols = ['Admit-Cannabinoid',
         'Admit-Cocaine',
         'Admit-Heroin',
         'Admit-Amphetamines',
         'Admit-Benzodiazapine',
         'Admit-Narcotics',
         'Admit-Ecstasy',
         'Admit-PCP',
         'Admit-Ritalin',
         'Admit-Other',
         'Admit-None']
test_cols = [#'Test-Barbiturates',
         #'Test-Benzodiazepine',
         'Test-Cannabinoid',
         'Test-Cocaine',
         'Test-Opiates']
pat_cols = ['Age',
         'IsMale']
race_cols = ['Race-Asian',
         'Race-Indian',
         'Race-Black',
         'Race-Hawaiian',
         'Race-White',
         'Race-Multiple',
         'Race-Unknown']

haart_cols = ['HAART-Naive',
         'HAART-Non-Adherent',
         'HAART-Off',
         'HAART-On',
         'HAART-Missing']
clinical_cols = ['VL',
         'CD4', ]#'pTMHDS'
         
mask = trans_data['HAART-On']==1
mask &= trans_data['Race-Black']==1
X = trans_data[pat_cols+clinical_cols+test_cols+['DaysSinceBL', 'Admit-None']][mask]
y = trans_data['TMHDS'][mask]

tmask = X.applymap(np.isnan).any(axis=1) | y.isnull()
nX = X[~tmask]
ny = y[~tmask]

res = sm.GLM(ny,nX).fit()
res.summary()

# <codecell>

from sklearn.preprocessing import Normalizer, Imputer
from sklearn.pipeline import Pipeline
from sklearn.linear_model import LinearRegression, BayesianRidge, LarsCV
from sklearn.neighbors import KNeighborsRegressor
from sklearn.cross_validation import cross_val_score
from sklearn.cross_validation import Bootstrap
from sklearn.metrics import mean_absolute_error, make_scorer
from sklearn.ensemble import AdaBoostRegressor, RandomForestRegressor
from sklearn.dummy import DummyRegressor
from sklearn.svm import SVR

regressors = [('Dummy', DummyRegressor()),
              ('Adaboost', AdaBoostRegressor()),
              ('linear', LinearRegression()),
              ('RF', RandomForestRegressor()),
              ('ridge', BayesianRidge()),
              ('KNN', KNeighborsRegressor()),
              ('SVR', SVR())]


for name, reg in regressors:
    pipe = Pipeline(steps=[('Norm', Normalizer()),
                           ('Regress', reg)])

    scores = cross_val_score(pipe, nX.values, ny.values, scoring=make_scorer(mean_absolute_error),
                             cv=Bootstrap(len(ny), n_iter=10, train_size=0.6))
    print name, scores.mean()

# <codecell>

from scipy.stats import norm
cdf = norm(x=0, scale=0.5)
reg = AdaBoostRegressor()
reg.fit(nX.values, ny.values)
pred_cols = pat_cols+clinical_cols+test_cols+['DaysSinceBL', 'Admit-None']
for pat in check_pats:
    fdate = baseline_date[pat]
    
    sampled = out_data.ix[pat].truncate(before=fdate)
    estimated = trans_pat_data.ix[pat].truncate(before=fdate)
    est_data = estimated[pred_cols].dropna()
    predicted = pd.Series(reg.predict(est_data.values), est_data.index)
    #predicted = predicted.combine_first(estimated['TMHDS']+cdf.rvs(estimated['TMHDS'].values.shape))
    plt.figure(figsize = (10, 5))
    ax = plt.subplot(111)
    ts = sampled['TMHDS'].dropna()
    es = estimated['TMHDS'].dropna()
    
    ax.plot_date(ts.index, ts, color = 'g', marker = 'd', linestyle = '-', markersize = 20, linewidth = 5, alpha = 0.8)
    ax.plot_date(es.index, es, color = 'b', marker = '.', linestyle = '-', markersize = 25, linewidth = 5, alpha = 0.8)
    ax.plot_date(est_data.index, predicted, color = 'r', marker = '*', linestyle = '-', markersize = 25, linewidth = 5, alpha = 0.8)
    ax.set_title(pat)
    ax.set_ylim([0, 12.5])
    ax.set_xlim(xlims[pat])
    ax.set_xticks(xticks[pat])
    plt.tight_layout()
    plt.savefig('pred_fig_%s.png' % pat, dpi = 500)

# <codecell>

res.pvalues

# <codecell>


