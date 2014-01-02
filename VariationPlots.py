# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=1>

# Examination of HIV Variation

# <markdowncell>

# In an effort to examine the amount of genetic variation from longitudinal visits of well-controlled patients. We are determining the number of mutations that occur in the LTR over a set of consecutive visits in which the patient has maintained a Viral-Load <100 copies/mL and varying levels of CD4 counts.

# <codecell>

from __future__ import division
from pandas import *
import os, os.path
import sys
import numpy as np

sys.path.append('/home/will/HIVReportGen/AnalysisCode/')
sys.path.append('/home/will/PySeqUtils/')
os.chdir('/home/will/HIVVariation/')

# <codecell>

from GeneralSeqTools import call_muscle

# <headingcell level=2>

# Data Extraction

# <markdowncell>

# Using the Redcap and sequence data up until 1/16/2013.

# <codecell>

store = HDFStore('/home/will/HIVReportGen/Data/SplitRedcap/2013-01-16/EntireCohort.hdf')
redcap_data = store['redcap']
seq_data = store['seq_data']

t = redcap_data['Event Name'].dropna().apply(lambda x: x.split(' - ')[0])
t.unique()
redcap_data['Patient visit number'] = redcap_data['Patient visit number'].combine_first(t)


# <codecell>

def VisitType(row):
    
    if row['ART'] == 'naive':
        return 'naive'
    elif (row['CD4'] >= 250) & (row['VL']<=100):
        return 'controlled'
    else:
        return 'wild'
    

# <codecell>

wanted_cols = ['Patient ID', 'Patient visit number', 
                'Date of visit', 
                'Latest CD4 count (cells/uL)', 'Date of latest CD4 count',
                'Latest viral load', 'Date of latest viral load',
                'Current ART status']
wanted_redcap = redcap_data[wanted_cols]
data = merge(wanted_redcap, seq_data[['LTR']],
            left_on = ['Patient ID', 'Patient visit number'],
            right_index = True, how = 'inner')
data = data.rename(columns= {
                                'Patient visit number':'VisitNum',
                                'Date of visit':'Date',
                                'Latest CD4 count (cells/uL)':'CD4',
                                'Date of latest CD4 count':'CD4-Date',
                                'Latest viral load':'VL',
                                'Date of latest viral load':'VL-Date',
                                'Current ART status':'ART'
                            }).dropna()
mask = data['Patient ID'] == 'A0139'
data = data.drop(mask[mask].index, axis = 0)
data.sort(['Patient ID', 'VisitNum'], inplace=True)
print 'Valid samples from Redcap/Sequencing'
print data

# <codecell>

hxb2_ltr = """TGGAAGGGCTAATTTACTCCCAAAAAAGACAAGATATCCTTGATCTGTGGGTC
TACCACACACAAGGCTACTTCCCTGATTGGCAGAACTACACACCAGGGCCAGG
GATCAGATATCCACTGACCTTTGGATGGTGCTTCAAGCTAGTACCAGTTGAGC
CAGAGAAGGTAGAAGAGGCCAATGAAGGAGAGAACAACAGCTTGTTACACCCT
ATGAGCCTGCATGGGATGGAGGACCCGGAGAAAGAAGTGTTAGTGTGGAAGTT
TGACAGCCGCCTAGCATTTCATCACATGGCCCGAGAGCTGCATCCGGAGTACT
ACAAGGACTGCTGACATCGAGCTTTCTACAAGGGACTTTCCGCTGGGGACTTT
CCAGGGAGGCGTGGCCTGGGCGGGACTGGGGAGTGGCGAGCCCTCAGATGCTG
CATATAAGCAGCTGCTTTTTGCCTGTACTGGGTCTCTCTGGTTAGACCAGATC
TGAGCCTGGGAGCTCTCTGGCTAACTAGGGAACCCACTGCTTAAGCCTCAATA
AAGCTTGCCTTGAGTGCTTCAAGTAGTGTGTGCCCGTCTGTTGTGTGACTCTG
GTAACTAGAGATCCCTCAGACCCTTTTAGTCAGTGTGGAAAATCTCTAGCA""".replace('\n', '')


def align_seq_ser(seq_series):
    
    nseqs = [(i, s) for i, s in zip(seq_series.index, seq_series.values)]
    nseqs += [('hxb2', hxb2_ltr)]
    daln = dict(call_muscle(nseqs))
    aln = [daln[str(s)] for s, _ in nseqs]
    aln_ser = Series(aln[:-1], seq_series.index)
    return aln_ser

data['LTR-align'] = data.groupby('Patient ID')['LTR'].apply(align_seq_ser)

# <codecell>

wanted_seq_cols = [108, 115, 120, 160, 168, 181, 251, 293] #1-based!!


def get_wanted_seq_cols(seq_series):
    
    nseqs = [(i, s) for i, s in zip(seq_series.index, seq_series.values)]
    nseqs += [('hxb2', hxb2_ltr)]
    daln = dict(call_muscle(nseqs))
    aln = [daln[str(s)] for s, _ in nseqs]
    outs = [[] for _ in range(len(aln)-1)]
    hxb2pos = 0
    for tup in zip(*aln):
        if tup[-1] != '-':
            hxb2pos += 1 #1-based!
        if hxb2pos in wanted_seq_cols:
            for out, let in zip(outs, tup):
                out.append(let)
    
    out_ser = Series(outs, seq_series.index)
    return out_ser
data['SNPCols'] = data.groupby('Patient ID')['LTR'].apply(get_wanted_seq_cols)

# <codecell>

data['LTR'].map(len).hist(bins = 50);
plt.xlabel('LTR Sequence Lengths');
plt.title('Our Data');

# <markdowncell>

# Some of our sequences are much smaller then the whole LTR. This may drastically effect the Mut/LTR calculation. For example if the sequence is only ~50 NT and it has 2 mutations then when I scale this up it may not represent the true number of mutations 'per LTR'. So I'm excluding any sequences smaller then 100NT.

# <codecell>

def calc_variation(seq1, seq2, scale):
    
    nmuts = 0
    nlets = 0
    for a, b in zip(seq1, seq2):
        if (a != '-') and (b != '-'):
            nlets += 1
            nmuts += a != b
    if nlets < 100:
        return np.nan
    frac = nmuts/nlets
    
    if scale == 'LTR':
        return frac * 600
    else:
        return scale * frac
    
    
def calc_cum_variation(aln_series, scale = 'LTR'):
    
    var = []
    sA = aln_series.values[0]
    for sB in aln_series.values:
        var.append(calc_variation(sA, sB, scale))
    return Series(var, index = aln_series.index)
    
    
def calc_series_variation(aln_series, scale = 'LTR'):
    
    var = [np.nan]
    for sA, sB in zip(aln_series.values, aln_series.values[1:]):
        var.append(calc_variation(sA, sB, scale))
    return Series(var, index = aln_series.index)

def calc_time_series_variation(var_time_series, scale = 'LTR', time_scale = 365):
    
    items = list(var_time_series.iterrows())
    if len(items) == 1:
        #rint 'short'
        return Series(np.nan, index = var_time_series.index)
    times = [np.nan]
    for (_, rA), (t_, rB) in zip(items, items[1:]):
        delta = (rA['Date'] - rB['Date']).days
        times.append(abs(delta)/time_scale)
    res = var_time_series['MutsPerLTR']/Series(times, index = var_time_series.index)
    #print 'new',res
    return res
    
    

data['MutsPerLTR'] = data.groupby('Patient ID')['LTR-align'].apply(calc_series_variation)
data['CumMutsPerLTR'] = data.groupby('Patient ID')['LTR-align'].apply(calc_cum_variation)
ndata = data.groupby('Patient ID', as_index = False).apply(calc_time_series_variation)
data['MutsPerLTRPerYear'] = ndata
    

# <codecell>

from collections import deque


    
def pat_yield_runs(df, min_run, check_funs, 
                    controlled, max_run = None, 
                    art_status = 'on'):
    
    crun = deque(maxlen = max_run)
    for _, row in df.iterrows():
        goods = [fun(row) for fun in check_funs]
        valid_pat = all(goods) == controlled
        valid_art = row['ART'] == art_status
        valid_mut = ~np.isnan(row['MutsPerLTRPerYear'])
        if valid_pat and valid_art and valid_mut:
            crun.append(row.copy())
        else:
            if len(crun) >= min_run:
                yield concat([DataFrame(it) for it in crun], 
                              axis = 1, ignore_index=True).T
            crun = []
        if max_run and len(crun) == max_run:
            yield concat([DataFrame(it) for it in crun], 
                            axis = 1, ignore_index=True).T


# <headingcell level=2>

# Figures

# <codecell>

def linker(col, val, direc, row):
    if val is None:
        return True
    if direc == 'gt':
        return row[col] >= val
    else:
        return row[col] < val

def make_boxplots(cd4_cuts, vl_cuts, wanted_col, ylabel):
    
    fig, axes = plt.subplots(len(cd4_cuts),len(vl_cuts), 
                            figsize = (10,15), 
                            sharey = True, 
                            sharex = True)
    
    for ax, (cd4_cut, vl_cut) in zip(axes.flatten(), product(cd4_cuts, vl_cuts)):
    
        cd_control = partial(linker, 'CD4', cd4_cut, 'gt')
        vl_control = partial(linker, 'VL', vl_cut, 'lt')
    
        check_funs = [cd_control, vl_control]
        checks = [('Controlled', True, 'on'),
                  ('Wild', False, 'on'),
                  ('Naive-C', True, 'naive'),
                  ('Naive-W', False, 'naive')]
        outs = []
        for name, control_val, art_status in checks:
            odata = []
            for _, df in data.groupby('Patient ID'):
                for run in pat_yield_runs(df, 3, check_funs, 
                                            control_val, 
                                            art_status=art_status):
                    odata += list(run[wanted_col].dropna())
            outs.append(deepcopy(odata))
        plt.sca(ax)
        plt.boxplot(outs, bootstrap = 1000, sym = '')
        plt.title('VL: %s, CD4: %s' % (vl_cut, cd4_cut))
        
        plt.xticks([1,2,3,4], ['Controlled', 'Wild', 'Naive-C', 'Naive-W'], rotation=90)
        plt.ylabel(ylabel)

# <headingcell level=3>

# Patient Information

# <codecell>

from itertools import product
from functools import partial
from copy import deepcopy

cd4_cuts = [0, 200, 250, 500]
vl_cuts = [50, 100, 200]

make_boxplots(cd4_cuts, vl_cuts, 'CD4', 'CD4 Count')
plt.savefig('boxplot_cd4.png')

# <codecell>

data['Log-VL'] = np.log10(data['VL'])
make_boxplots(cd4_cuts, vl_cuts, 'Log-VL', 'log10(VL)')
plt.savefig('boxplot_viral.png')

# <markdowncell>

# The viral loads of the 'Wild' Patients vary quite a bit. From the Cutoff all the way to 10000.

# <codecell>

make_boxplots(cd4_cuts, vl_cuts, 'MutsPerLTRPerYear', 'Muts/LTR/Year')
plt.savefig('boxplot_muts.png')

# <headingcell level=3>

# Longitudinal Variation

# <codecell>

def get_months(dt):
    return dt.astype('timedelta64[D]')/np.timedelta64(24*30, 'h')


def make_long_plots(cd4_cuts, vl_cuts, plot_col = 'MutsPerLTRPerYear'):
    
    fig, axes = plt.subplots(len(cd4_cuts),len(vl_cuts), 
                            figsize = (10,10), 
                            sharey = True, 
                            sharex = True)
    axes[0,0].set_ylabel('Muts/LTR/Year')
    for ax, (cd4_cut, vl_cut) in zip(axes.flatten(), product(cd4_cuts, vl_cuts)):
    
        cd_control = partial(linker, 'CD4', cd4_cut, 'gt')
        vl_control = partial(linker, 'VL', vl_cut, 'lt')
    
        check_funs = [cd_control, vl_control]
        checks = [('Controlled', True),
                  ('Wild', False)]
        outs = []
        plt.sca(ax)
        plt.title('VL: %s, CD4: %s' % (vl_cut, cd4_cut))
        plt.hold(True)
        for name, control_val in checks:
            odata = []
            for _, df in data.groupby('Patient ID'):
                for run in pat_yield_runs(df, 3, check_funs, control_val):
                    tdelta = (run['Date'] - run['Date'].min()).map(get_months)
                    tmuts = run[plot_col]
                    odata.append((tdelta.copy(), tmuts.copy()))
            colors = np.linspace(0.1,1,len(odata))
            
            for (tims, muts), color in zip(odata, colors):
                if name == 'Controlled':
                    ocolor = (color, 0, 0)
                else:
                    ocolor = (0, 0, color)
                try:
                    plt.plot(tims, muts, color = ocolor)
                except TypeError:
                    print tims
                    print muts
                    raise TypeError
        plt.xlim([0,36])
        
        plt.hold(False)
                
make_long_plots(cd4_cuts, vl_cuts)       
plt.savefig('long_var_plots.png')

# <markdowncell>

# This figure shows the lognitudinal change of the LTR in each patient. Red patients are 'Wild' and Blue paitents are 'Controlled'. We can see from this figure that there is no difference between the two groups. Everyone has change.

# <codecell>

make_long_plots(cd4_cuts, vl_cuts, plot_col = 'CumMutsPerLTR')
plt.savefig('cum_long_var_plots.png')

# <markdowncell>

# The figure above shows the Cumulative variation of the sequences. At each point I've plotted the variation from the _initial_ visit. I still don't see a difference between controlled (blue) and wild (red).

# <headingcell level=3>

# Mutation Correlation with VL

# <codecell>

wanted_cols = ['Patient ID', 'Patient visit number', 
                'Date of visit', 
                'Latest CD4 count (cells/uL)', 'Date of latest CD4 count',
                'Latest viral load', 'Date of latest viral load',
                'Current ART status']
wanted_redcap = redcap_data[wanted_cols]
new_data = merge(wanted_redcap, seq_data[['LTR']],
            left_on = ['Patient ID', 'Patient visit number'],
            right_index = True, how = 'inner')
new_data = new_data.rename(columns= {
                                'Patient visit number':'VisitNum',
                                'Date of visit':'Date',
                                'Latest CD4 count (cells/uL)':'CD4',
                                'Date of latest CD4 count':'CD4-Date',
                                'Latest viral load':'VL',
                                'Date of latest viral load':'VL-Date',
                                'Current ART status':'ART'
                            })
mask = new_data['Patient ID'] == 'A0139'
new_data = new_data.drop(mask[mask].index, axis = 0)
#drop out patients which have short LTRS
new_data['LTR'][new_data['LTR'].fillna('').map(len)<400] = np.nan

aln_data = new_data[['Patient ID', 'VisitNum', 'LTR', 'Date']].dropna()
aln_data['LTR-align'] = aln_data.groupby('Patient ID')['LTR'].apply(align_seq_ser)

aln_data['MutsPerLTR'] = aln_data.groupby('Patient ID')['LTR-align'].apply(calc_series_variation)
aln_data['MutsPerLTR'][aln_data['MutsPerLTR']> 200] = np.nan
aln_data['ACumMutsPerLTR'] = aln_data.groupby('Patient ID')['LTR-align'].apply(calc_cum_variation)
aln_data['CumMutsPerLTR'] = aln_data.groupby('Patient ID')['MutsPerLTR'].cumsum()
ndata = aln_data.groupby('Patient ID', as_index = False).apply(calc_time_series_variation)
aln_data['MutsPerLTRPerYear'] = ndata

new_seq_data = merge(new_data, aln_data.drop(['Date', 'LTR'], axis = 1), 
                     left_on = ['Patient ID', 'VisitNum'],
                     right_on = ['Patient ID', 'VisitNum'],
                     how = 'left')


# <codecell>

vl_data = new_seq_data[['Patient ID', 'VL', 'VL-Date']].rename(columns = {'VL-Date':'Date'})
cd4_data = new_seq_data[['Patient ID', 'CD4', 'CD4-Date']].rename(columns = {'CD4-Date':'Date'})
var_data = new_seq_data[['Patient ID', 'ART', 'Date', 'CumMutsPerLTR', 'MutsPerLTRPerYear']]

long_data = concat([vl_data, cd4_data, var_data], axis = 0, ignore_index = True)
long_data.sort(['Patient ID', 'Date'], inplace=True)
long_data = long_data.groupby(['Patient ID', 'Date']).first()

# <codecell>

import matplotlib.dates as mdates
from matplotlib.dates import num2date, date2num

def date_fmt(inum):
    return num2date(inum).strftime("%m-%Y")

def plot_pat_data(pat_data, vl_ax):
    
    colors = {
              'on':'k',
              'off':'r',
              'naive':'b'
              }
    tdata = pat_data[['CumMutsPerLTR', 'VL', 'ART']].dropna(how = 'all').reset_index()
    tdata['DateNums'] = tdata['Date'].map(date2num)
    vl_data = tdata[['DateNums', 'VL']].dropna()
    vl_ax.plot(vl_data['DateNums'].values, 
               vl_data['VL'].values, 'r-',
               lw = 4, alpha = 0.9)
    vl_ax.set_yscale('log')
    for tl in vl_ax.get_yticklabels():
        tl.set_color('r')

    cum_ax = vl_ax.twinx()
    cum_data = tdata[['DateNums', 'CumMutsPerLTR']].dropna()
    cum_ax.plot(cum_data['DateNums'].values, 
                cum_data['CumMutsPerLTR'].values, 'g.-', 
                lw = 4, alpha = 0.9)
    for tl in cum_ax.get_yticklabels():
        tl.set_color('g')

    vl_ax.set_xticklabels(map(date_fmt, cum_ax.get_xticks()))
    vl_ax.set_ylabel('Viral Load')
    cum_ax.set_ylabel('Cumulative Muts')
    
    art_data = tdata[['DateNums', 'ART']].dropna()
    
    height = cum_ax.get_ylim()[1]*0.05
    for r1, r2 in zip(art_data.index, art_data.index[1:]):
        xy = (art_data.ix[r1]['DateNums'], 0)
        width = art_data.ix[r2]['DateNums'] - art_data.ix[r1]['DateNums']
        color = colors[art_data.ix[r1]['ART']]
        cum_ax.add_patch(Rectangle(xy,width,height, facecolor = color, alpha = 0.5))
    
    plt.gcf().autofmt_xdate()

# <codecell>

highest_vl = long_data['VL'].groupby(level = 'Patient ID').max()
frac_controlled_vl = (long_data['VL']<200).groupby(level = 'Patient ID').mean()
has_naive = (long_data['ART'] == 'naive').groupby(level = 'Patient ID').any()
has_art = (long_data['ART'] == 'on').groupby(level = 'Patient ID').any()
num_seqs = long_data['CumMutsPerLTR'].notnull().groupby(level = 'Patient ID').sum()

pat_props = DataFrame({
                       'Max-VL':highest_vl,
                       'FracCont':frac_controlled_vl,
                       'has_naive':has_naive,
                       'has_art':has_art,
                       'NumSeqs':-num_seqs
                       })

# <codecell>

always_controlled = pat_props['has_art'] & (pat_props['Max-VL'] <= 200)
controlled_pats = pat_props[always_controlled].sort('NumSeqs')
agree_cont_pat = controlled_pats.index[9]
dis_cont_pat = controlled_pats.index[2]

# <codecell>

wanted_pats = (pat_props['has_naive']==0) & (pat_props['NumSeqs']<-3)
controlled_pats = pat_props[wanted_pats].sort('Max-VL')
agree_wild_pat = controlled_pats.index[-2]
dis_wild_pat = controlled_pats.index[-4]

# <codecell>

wanted_pats = (pat_props['has_naive']==1) & (pat_props['NumSeqs']<-3)
controlled_pats = pat_props[wanted_pats].sort('Max-VL')
agree_naive_pat = controlled_pats.index[-2]
dis_naive_pat = controlled_pats.index[-1]

# <codecell>

wanted_pats = [agree_cont_pat, dis_cont_pat, 
               agree_wild_pat, dis_wild_pat, 
               agree_naive_pat, dis_naive_pat]

fig, axs = plt.subplots(3, 2, figsize = (10, 10))

for ax, pat in zip(axs.flatten(), wanted_pats):
    
    plot_pat_data(long_data.ix[pat], ax)
    ax.set_title(pat)
    
plt.tight_layout()
#plt.savefig('VL_var_plots.png')

# <codecell>

first_visits = new_data[['Patient ID', 'Date']][new_data['VisitNum'] == 'R00'].set_index('Patient ID')['Date']
def time_deltas(indf):
    print indf
    pid = indf.index[0][0]
    try:
        fdate = first_visits.ix[pid]
        deltas = (indf-fdate).map(get_months)
    except KeyError:
        return Series([np.nan]*len(indf), index = indf.index)
    return deltas
    
a, b = first_visits.align(long_data['CumMutsPerLTR'], level = 'Patient ID', join = 'right')
long_data['FirstVisit'] = a
tmp = long_data[['FirstVisit']].dropna()
tmp_i = tmp.index
tmp = tmp.reset_index()
tmp.index = tmp_i
deltaT = (tmp['Date']-tmp['FirstVisit'])
_, deltaT = long_data.align(deltaT, axis = 0, join = 'left')
long_data['DeltaT'] = deltaT

# <codecell>

from datetime import date
from scipy.stats import linregress

t = long_data[['CumMutsPerLTR', 'DeltaT']].dropna().reset_index(drop = True).sort('DeltaT')

t['StdDate'] = t['DeltaT'] + date.today()
t = t[t['DeltaT'].notnull()]
date_muts = t.groupby('StdDate')['CumMutsPerLTR'].median()
res = ewma(date_muts, span = 50, freq = '30d')
res.name = 'Mean-Mutations'
mean_muts = res.reset_index()
mean_muts['DeltaT-m'] = (mean_muts['StdDate'] - date.today()).map(get_months)

m, b, r, _, _ = linregress(t['DeltaT'].map(get_months).values,
                           t['CumMutsPerLTR'].values)
mean_muts['Regress'] = m*mean_muts['DeltaT-m']+b
print r**2

# <codecell>

long_data['DeltaT-m'] = long_data['DeltaT'].map(get_months)
low_vl = long_data['VL'].fillna(method = 'pad')<=100
high_cd = long_data['CD4'].fillna(method = 'pad')>=200
long_data['gART'] = long_data['ART'].groupby(level = 'Patient ID').fillna(method = 'pad')

scatter_data = long_data[(low_vl & high_cd) | (long_data['gART'] == 'naive')]
scatter_data = scatter_data[['CumMutsPerLTR', 'gART', 'DeltaT-m']].dropna()
colors = {
          'on':'g',
          'off':'r',
          'naive':'b',
          'non-adherent':'r'
          }

scatter_data['Color'] = scatter_data['gART'].map(lambda x: colors[x])
plt.figure(figsize = (8,5))
plt.scatter(scatter_data['DeltaT-m'], scatter_data['CumMutsPerLTR'], 
            c = list(scatter_data['Color'].values),
            alpha = 0.8, edgecolor = list(scatter_data['Color'].values))
plt.ylim([0, 180])
plt.xlim([0, 75])
plt.hold(True)
plt.plot(mean_muts['DeltaT-m'], mean_muts['Mean-Mutations'], 'k', lw = 2)
plt.ylabel('Cumulative Mutations')
plt.xlabel('Months in the study')
plt.hold(False)
plt.savefig('scatter_variation.png', dpi = 200)

# <codecell>


# <codecell>


# <markdowncell>

# From these figures it looks like there is a roughly equal amount of variation when you look at well controlled and uncontrolled patients. We can also guess that in general there are bursts of genetic variation which wanes over time. To examine this I'm going to look at consecutive pairs of visits (instead of requiring 3+ visits) and then compare the results of consecutive well-controlled visits to consecutive un-controlled visits.

# <headingcell level=3>

# Paired Visits Table

# <codecell>

from collections import defaultdict

dfdict = defaultdict(list)
for (cd4_cut, vl_cut) in product(cd4_cuts, vl_cuts):
    
    cd_control = partial(linker, 'CD4', cd4_cut, 'gt')
    vl_control = partial(linker, 'VL', vl_cut, 'lt')
    
    check_funs = [cd_control, vl_control]
    checks = [('Controlled', True, 'on'),
              ('Wild', False, 'on'),
              ('Naive-C', True, 'naive'),
              ('Naive-W', False, 'naive')]
    
    for name, control_val, art_status in checks:
        for pat, df in data.groupby('Patient ID'):
            for run in pat_yield_runs(df, 2, check_funs, 
                                        control_val, 
                                        art_status=art_status, 
                                        max_run=2):
                items = run['MutsPerLTRPerYear'].dropna()
                dfdict['MutsPerLTRPerYear'] += list(items)
                dfdict['Patient ID'] += [pat]*len(items)
                dfdict['Name'] += [name]*len(items)
                dfdict['VLCut'] += [vl_cut]*len(items)
                dfdict['CD4Cut'] += [cd4_cut]*len(items)
                dfdict['ArtStatus'] += [art_status]*len(items)
                
                
paired_data = DataFrame(dfdict)
print paired_data.head()

# <codecell>

agg_pairs = pivot_table(paired_data, rows = ['CD4Cut', 'VLCut', 'Patient ID'], 
                        cols= 'Name', values = 'MutsPerLTRPerYear', 
                        aggfunc=np.mean)
agg_pairs.describe()

# <markdowncell>

# Even under the 'pair-wise' model there are very few Naive-Controlled patients. So I'm going to exclude them from the rest of the analysis.

# <headingcell level=3>

# Variation Boxplot

# <codecell>

fig, axes = plt.subplots(len(cd4_cuts),len(vl_cuts), 
                            figsize = (10,15), 
                            sharey = True, 
                            sharex = True)
    
order = ['Controlled', 'Wild', 'Naive-W']
for ax, (cd4_cut, vl_cut) in zip(axes.flatten(), product(cd4_cuts, vl_cuts)):
    if cd4_cut is None:
        cd4_mask = paired_data['CD4Cut'].isnull()
    else:
        cd4_mask = paired_data['CD4Cut'] == cd4_cut
    vl_mask = paired_data['VLCut'] == vl_cut
    tdata = paired_data[vl_mask & cd4_mask]
    masks = [tdata['Name'] == o for o in order]
    
    
    plt.sca(ax)
    outs = [tdata[mask]['MutsPerLTRPerYear'] for mask in masks]
    plt.boxplot(outs, bootstrap = 1000, sym='')
    plt.title('VL: %s, CD4: %s' % (vl_cut, cd4_cut))
    plt.ylabel('Muts/LTR/Year')
    plt.xticks([1,2,3,4], order)
plt.tight_layout()
plt.savefig('pairwise_varaiation.png')

# <markdowncell>

# In order to examine whether the consecutive visits have a different number of mutations per 100bp I made the above boxplots. There are 76 'wild' patients, 114 controlled and 12 naive patients. There clearly no differences in the number of mututations on different therapy regmines. The bottom plot normalizes the changes by the sampling time, just in case there was an issue with one group going longer between visits. Clearly not an issue.

# <codecell>

wanted_seq_cols = [108, 115, 120, 160, 168, 181, 251, 293]
ref_val = ['A', 'A', 'C', 'C', 'G', 'A', 'G', 'G']

names = ['%i-%s' % (c, v) for c,v in zip(wanted_seq_cols, ref_val)]
def check_data(series):
    
    out = []
    for let, wlet in zip(series['SNPCols'], ref_val):
        if let == '-':
            out.append(np.nan)
        else:
            out.append(float(let == wlet))
    #print out
    return Series(out, index=names)

snp_data = data.apply(check_data, axis = 1)
ndata = concat([data, snp_data], axis = 1)

# <codecell>

dfdict = defaultdict(list)
for (cd4_cut, vl_cut) in product(cd4_cuts, vl_cuts):
    
    cd_control = partial(linker, 'CD4', cd4_cut, 'gt')
    vl_control = partial(linker, 'VL', vl_cut, 'lt')
    
    check_funs = [cd_control, vl_control]
    checks = [('Controlled', True, 'on'),
              ('Wild', False, 'on'),
              ('Naive-C', True, 'naive'),
              ('Naive-W', False, 'naive')]
    
    for name, control_val, art_status in checks:
        for pat, df in ndata.groupby('Patient ID'):
            for run in pat_yield_runs(df, 1, check_funs, control_val, art_status=art_status, max_run=1):
                items = run.dropna(subset = ['MutsPerLTRPerYear'])
                #dfdict['MutsPerLTRPerYear'] += list(items[['MutsPerLTRPerYear']])
                dfdict['Name'] += [name]*len(items)
                dfdict['VLCut'] += [vl_cut]*len(items)
                dfdict['CD4Cut'] += [cd4_cut]*len(items)
                dfdict['ArtStatus'] += [art_status]*len(items)
                for col in ndata.columns:
                    dfdict[col] += list(items[col])
single_data = DataFrame(dfdict)

# <codecell>

def min_num_mean(vals):
    mval = 10
    if len(vals) < mval:
        return np.nan
    else:
        return np.mean(vals)

snp_names = ['%i-%s' % (c, v) for c,v in zip(wanted_seq_cols, ref_val)]
snp_frac = pivot_table(single_data, rows = ['CD4Cut', 'VLCut', 'Name'], 
                         values = snp_names, 
                        aggfunc=min_num_mean).dropna()
snp_frac.to_csv('snp_grouped_fractions.csv')
snp_frac

# <codecell>

def check_single_transition(dates, snp_col):

    
    deltas = dates - dates.min()
    begin_state = snp_col.values[0]
    out_state = 'CBloss' if begin_state else 'CBgain'
    for delta, snp_val in zip(deltas.values, snp_col.values):
        if snp_val != begin_state:
            return delta.days, out_state
    return np.nan, np.nan

def check_run_transition(run):
    
    times = []
    changes = []
    wnames = []
    for snp_name in snp_names:
        trun = run[['Date', snp_name]].dropna()
        if len(trun.index) > 1:
            d, c = check_single_transition(trun['Date'], trun[snp_name])
            times.append(d)
            changes.append(c)
            wnames.append(snp_name)
        
    out_df = DataFrame({
                        'TransTimes':times,
                        'StateChanges':changes,
                        'SNPNames': wnames
                        })
    return out_df

# <codecell>

trans_list = []
for (cd4_cut, vl_cut) in product(cd4_cuts, vl_cuts):
    
    cd_control = partial(linker, 'CD4', cd4_cut, 'gt')
    vl_control = partial(linker, 'VL', vl_cut, 'lt')
    
    check_funs = [cd_control, vl_control]
    checks = [('Controlled', True, 'on'),
              ('Wild', False, 'on'),
              ('Naive-C', True, 'naive'),
              ('Naive-W', False, 'naive')]
    
    for name, control_val, art_status in checks:
        for pat, df in ndata.groupby('Patient ID'):
            for run in pat_yield_runs(df, 2, check_funs, 
                                        control_val, art_status=art_status):
                trans_data = check_run_transition(run)
                if trans_data['StateChanges'].notnull().any():
                    trans_data['Patient ID'] = pat
                    trans_data['Name'] = name
                    trans_data['VLCut'] = vl_cut
                    trans_data['CD4Cut'] = cd4_cut
                    trans_list.append(trans_data.dropna())
all_trans = concat(trans_list, axis = 0, ignore_index=True)                

# <codecell>

from functools import partial
def min_num_mean(mval, vals):
    if len(vals) < mval:
        return np.nan
    else:
        return np.mean(vals)

single_support = pivot_table(all_trans, rows = ['CD4Cut', 'VLCut', 'Name', 'StateChanges'], 
                                     cols = 'SNPNames', values = 'TransTimes', 
                                    aggfunc=partial(min_num_mean,1)).dropna(how = 'all')
single_support.to_csv('trans_times_single.csv')
single_support

# <markdowncell>

# Looking at the correlation between CD4 and number of mutations. I do not really see any relationship between CD4 count and number of mutations. The outliers up there worry me a little, but there are roughly equal numbers of each type.

# <codecell>

multi_support = pivot_table(all_trans, rows = ['CD4Cut', 'VLCut', 'Name', 'StateChanges'], 
                                     cols = 'SNPNames', values = 'TransTimes', 
                                    aggfunc=partial(min_num_mean,2)).dropna(how = 'all')
multi_support.to_csv('trans_times_multi.csv')
multi_support

# <codecell>


