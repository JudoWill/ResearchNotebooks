# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

from __future__ import division
from pandas import *
import os, os.path
import sys
import numpy as np

sys.path.append('/home/will/HIVReportGen/AnalysisCode/')
sys.path.append('/home/will/PySeqUtils/')

# <codecell>

store = HDFStore('/home/will/HIVReportGen/Data/SplitRedcap/2013-01-16/EntireCohort.hdf')

# <codecell>

redcap_data = store['redcap']
seq_data = store['seq_data']

t = redcap_data['Event Name'].dropna().apply(lambda x: x.split(' - ')[0])
t.unique()
redcap_data['Patient visit number'] = redcap_data['Patient visit number'].combine_first(t)

# <codecell>

wanted_cols = ['Patient ID', 'Patient visit number', 'Date of visit', 'Latest CD4 count (cells/uL)', 'Latest viral load', 'Current ART status']
wanted_redcap = redcap_data[wanted_cols]

data = merge(wanted_redcap, seq_data[['LTR-bin-align']],
            left_on = ['Patient ID', 'Patient visit number'],
            right_index = True, how = 'inner').dropna()
data = data.rename(columns= {
                                'Patient visit number':'VisitNum',
                                'Date of visit':'Date',
                                'Latest CD4 count (cells/uL)':'CD4',
                                'Latest viral load':'VL',
                                'Current ART status':'ART'
                            })
print data

# <codecell>

data['ART'].unique()

# <codecell>

def QuantitateChange(seq_list):
    
    delta = 0
    for a, b in zip(seq_list, seq_list[1:]):
        delta =+ (np.array(a[200:-100]) == np.array(b[200:-100])).sum()
    return delta/len(seq_list)

# <codecell>

def ProcessPatient(df):
    
    if len(df) < 2:
        return Series({
                    'NiaveVisits':np.nan,
                    'ARTVisits':np.nan,
                    'NumWellControlled':np.nan,
                    'LongestWellControlled':np.nan,
                    'NumUnControlled':np.nan,
                    'ControlledVariation':np.nan,
                    'UnControlledVariation':np.nan})
    
    num_niave = (df['ART'] == 'naive').sum()
    num_on = (df['ART'] == 'on').sum()
    
    well_controlled = (df['VL'] <= 100) & (df['CD4'] >= 250)
    num_well = well_controlled.sum()
    num_uncontrolled = (~well_controlled).sum()
    
    run = 0
    longest_run = 0
    for well in well_controlled.values:
        if well:
            run += 1
        else:
            longest_run = max(longest_run, run)
            run = 0
        
    longest_well_controlled = longest_run
    
    well_controlled_variation = 0
    if num_well > 1:
        well_controlled_variation = QuantitateChange(df['LTR-bin-align'][well_controlled])
    
    uncontrolled_variation = 0
    if num_uncontrolled > 1:
        uncontrolled_variation = QuantitateChange(df['LTR-bin-align'][~well_controlled])
        
        
        
    return Series({
                    'NiaveVisits':num_niave,
                    'ARTVisits':num_on,
                    'NumWellControlled':num_well,
                    'LongestWellControlled':longest_well_controlled,
                    'NumUnControlled':num_uncontrolled,
                    'ControlledVariation':well_controlled_variation,
                    'UnControlledVariation':uncontrolled_variation})
    

ndata = data.groupby('Patient ID').apply(ProcessPatient)
print ndata

# <codecell>

niave_to_art = (ndata['NiaveVisits']>0) & (ndata['ARTVisits']>1)

print ndata.ix[niave_to_art].to_string()

uncontrolled_to_well_controlled = (ndata['LongestWellControlled'] >= 2) & (ndata['NumUnControlled'] >= 1)
print ndata.ix[uncontrolled_to_well_controlled].to_string()

# <codecell>

max_day = ndata['DateFromImpaired'].max()
min_day = ndata['DateFromImpaired'].min()

time_axis = np.arange(min_day, max_day)

good_max = ndata['DateFromImpaired'][ndata['DateFromImpaired']>0].median()
good_min = ndata['DateFromImpaired'][ndata['DateFromImpaired']<0].median()
print good_max, good_min

# <codecell>

npats = len(ndata['Patient ID'].unique())
ltr_len = len(data['LTR-bin-align'].values[0])
img_mat = np.ones((len(time_axis), ltr_len, npats))*np.nan


for depth, (pat, df) in enumerate(ndata.groupby('Patient ID')):
    
    ndf = df.copy()
    ndf.sort('DateFromImpaired')
    #ltr_mat = np.vstack(ndf['LTR-bin-align'].values)
    seqs = list(ndf['LTR-bin-align'].values)
    nseqs = [seqs[0]] + seqs + [seqs[-1]]
    times = list(ndf['DateFromImpaired'].values)
    ntimes = [min_day] + times + [max_day]
    
    for start_time, stop_time, ltr in zip(ntimes, ntimes[1:], nseqs):
        start_ind = int(start_time - min_day)
        stop_ind = int(stop_time - min_day)
        #print ntimes, start_time, start_ind, stop_time, stop_ind, ltr
        #raise KeyError
        img_mat[start_ind:stop_ind,:,depth] = ltr

    #print ltr_mat
    #raise KeyError
        
    
    

# <codecell>

from scipy.stats import ttest_ind
t_img = nanmean(img_mat, axis = 2)
pre_mask = time_axis>0

pvals = []
for n in range(ltr_len):
    pre_data = t_img[pre_mask, n]
    post_data = t_img[~pre_mask, n]
    if len(post_data)>5 and len(pre_data)>5:
        _, pval = ttest_ind(pre_data[pre_data == pre_data], 
                            post_data[post_data == post_data])
    else:
        pval = np.nan
    pvals.append(pval)
pvals = np.array(pvals)

# <codecell>

from pylab import get_cmap
plt.figure(figsize = (10,10))
good_min = -24
good_max = 24
plot_mask = (time_axis > good_min) & (time_axis < good_max)

tmp = pvals < 1e-30
n_img = t_img.copy()
n_img = n_img[plot_mask, :]
n_img = n_img[:, tmp]
ticks = np.where(tmp)

plt.imshow(np.flipud(n_img), 
                interpolation = 'nearest', 
                cmap = get_cmap('gray'), 
                aspect = 'auto', vmin = 0, vmax = 1.0,
                extent = (0, len(ticks[0]), good_max, good_min))
plt.xticks(np.arange(0, len(ticks[0]), 20), ticks[0][::20], rotation=90)
plt.colorbar()
plt.ylabel('Months from Last Impaired visit');
plt.xlabel('LTR pos');
plt.title('Significant LTR')

# <codecell>

tmp = np.abs(np.diff(t_img, axis = 0)).sum(axis = 0)
#print tmp
plt.hist(tmp, bins = 50)

# <codecell>

t_img.shape

# <codecell>


