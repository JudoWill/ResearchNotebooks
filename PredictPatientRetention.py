# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=1>

# Predicting Patient Retention Rates

# <markdowncell>

# Here I am looking for a simple method to predict which patients are likely to return. My idea is to look at the average time between visits across all patients and then across this specific patient.

# <codecell>

from __future__ import division
import os, os.path
import sys
import pandas as pd
import numpy as np

sys.path.append('/home/will/HIVReportGen/AnalysisCode/')
sys.path.append('/home/will/PySeqUtils/')
sys.path.append('/home/will/PatientPicker/')
import LoadingTools
#os.chdir('/home/will/HIVVariation/')

# <codecell>

store = pd.HDFStore('/home/will/HIVReportGen/Data/SplitRedcap/2013-01-16/EntireCohort.hdf')
redcap_data = store['redcap']
seq_data = store['seq_data']

t = redcap_data['Event Name'].dropna().apply(lambda x: x.split(' - ')[0])
t.unique()
redcap_data['VisitNum'] = redcap_data['Patient visit number'].combine_first(t)
redcap_data['VisitNum'][redcap_data['VisitNum']=='A03'] = 'R03'

# <codecell>

fields = ['Patient ID', 'VisitNum', 'Date of visit']
data = redcap_data[fields].rename(columns = {'Date of visit':'Date'})

# <headingcell level=2>

# Data Description

# <markdowncell>

# Here I define my return or not-return patients. In this case I'm defining every patient that _actually_ returned as True and every patient that has been more the 365\*2 days since thier last visit (using 1/16/2013 as the 'current date') as False. If its been less then two years the I exclude that visit from the analysis.

# <codecell>

from datetime import datetime
def get_diff_days(inser):
    return np.diff(inser)/(1e9*60*60*24)

def get_visit_diffs(inpat):
    
    nvisit = pd.DataFrame({
                        'Date':[datetime(2013,1,16)],
                        'VisitNum':['RN']
                        })
    
    ndata = pd.concat([inpat, nvisit], axis = 0, ignore_index=True)
    ndata.sort('Date', inplace=True)
    ndata['DiffDate'] = pd.rolling_apply(ndata['Date'], 2, get_diff_days)
    return ndata.set_index('VisitNum').drop('Patient ID', axis = 1)

odata = data.groupby('Patient ID').apply(get_visit_diffs).dropna()
print odata.head(n=30)

# <codecell>

from scipy.stats import norm
cohort_level_data = odata.groupby(level=1)['DiffDate'].agg({'std':'std', 
                                                            'mean':'mean', 
                                                            'count':'count'})
print cohort_level_data

# <markdowncell>

# Above is a table of the average times between visits. This only includes patients that _actually returned_ for the R0X visit. You can see that for the first few visits the average is well above the 6-month goal but it levels out around R05.
# 

# <headingcell level=2>

# Prediction

# <markdowncell>

# Here I'm builing a cohort-level 'Surivial Function'. In this case I'm using mean and std from between-visit times for all patients at each visit. I'm assuming a Gaussian Distribution. Then I build a Patient-Level 'Survival Function' based on thier between-visit times. For patients with less than 3 visits I built a SF from all patients.

# <codecell>

cohort_norm = {}
for key, row in cohort_level_data.iterrows():
    cohort_norm[key] = norm(loc = row['mean'], scale = row['std'])

# <codecell>

pat_mu = odata['DiffDate'].drop('RN', axis=0, level=1).mean()
pat_std = odata['DiffDate'].drop('RN', axis=0, level=1).std()

def add_sf_data(inpat):
    if len(inpat) > 3:
        mu = inpat['DiffDate'].mean()
        st = inpat['DiffDate'].std()
        obj = norm(loc=mu, scale=st)
    else:
        obj = norm(loc=pat_mu, scale=pat_std)
    
    inpat['CohortSF'] = np.nan
    inpat['PatSF'] = np.nan
    inpat['Returned'] = True
    inpat['Returned'].iloc[-1] = np.nan if inpat['DiffDate'].iloc[-1] < 365*2 else False

    for key, row in inpat.iterrows():
        inpat['CohortSF'].ix[key] = cohort_norm[key[1]].sf(row['DiffDate'])
        inpat['PatSF'].ix[key] = obj.sf(row['DiffDate'])
    return inpat
        
ndata = odata.groupby(level=0).apply(add_sf_data)
print ndata[['DiffDate', 'CohortSF', 'PatSF', 'Returned']].head(n = 30)

# <markdowncell>

# Guessing at how to combine these two Survival Functions is difficult. As such I'm using the SKLearn package to build a DecisionTree and a Naive-Bayes predictor using ONLY THESE PARAMETERS. This has the advantage of not directly biasing any future selection from these results. I'm also comparing this to a simple DummyClassifier which will guess that all patients return.

# <codecell>

from sklearn.metrics import auc_score, zero_one_loss, roc_curve
from sklearn.tree import DecisionTreeClassifier, export_graphviz
from sklearn.cross_validation import cross_val_score, Bootstrap
from sklearn.naive_bayes import GaussianNB
from sklearn.svm import SVC
from sklearn.ensemble import AdaBoostClassifier
from sklearn.dummy import DummyClassifier
from sklearn.linear_model import LogisticRegression

X = ndata.dropna()[['CohortSF', 'PatSF']].values
y = ndata.dropna()['Returned'].values

# <codecell>

import matplotlib.pyplot as plt
from collections import defaultdict

classifiers = [(DecisionTreeClassifier(), 'DecisionTree', 'r'),
               (GaussianNB(), 'NaiveBayes', 'g'),
               (AdaBoostClassifier(), 'Adaboost', 'c'),
               (LogisticRegression(), 'Logistic', 'k'),
               (DummyClassifier(), 'Dummy', 'b')]
plt.figure(figsize = (10,10))

losses = defaultdict(float)
nreps = 5
for train, test in Bootstrap(len(y), nreps):
    
    for pred, name, color in classifiers:
        fitted = pred.fit(X[train, :], y[train])
        pred_y = fitted.predict_proba(X[test, :])
        fpr, tpr, _ = roc_curve(y[test], pred_y[:,1])
        plt.plot(fpr, tpr, color, label=name)
        losses[name] += zero_one_loss(y[test], fitted.predict(X[test, :]))/nreps
    
plt.xlabel('False-Positive-Rate')
plt.ylabel('True-Positve-Rate')
plt.legend(['DecisionTree', 'NaiveBayes', 'Adaboost', 'Logistic', 'Dummy'], 'lower right');

# <markdowncell>

# The ROC curve is a commno method to look at the effectiveness of a classifier. It measures the trade-off between True-Positves and False-Negatives. The larger the Area Under the Curve the better the score. A random coin flip would have an area of 0.5 (the blue line).

# <codecell>

scores = numpy.array([losses[name] for name in ['DecisionTree', 'NaiveBayes', 'Adaboost', 'Logistic', 'Dummy']])
plt.bar([1, 2, 3,4,5], scores, width = 0.5)
plt.ylabel('%Miss-classified')
plt.xticks([1.25,2.25,3.25,4.25,5.25], ['DecisionTree', 'NaiveBayes', 'Adaboost', 'Logistic', 'Dummy']);

# <markdowncell>

# This graph shows the effectiveness of each of the three methods. The y-axis represents the fraction of miss-classified samples (averaged over 5 trials). We can see that the DecisionTree has only a 3% likelihood of mis-classifying a patient as return or not return. We can use this classifier to prioritize which patients we call for return visits.

# <codecell>

def expand_sf_data(inpat):
    dates = pd.date_range('1/1/2013', periods = 30, freq = 'M')
    if len(inpat) > 3:
        mu = inpat['DiffDate'].mean()
        st = inpat['DiffDate'].std()
        obj = norm(loc=mu, scale=st)
    else:
        obj = norm(loc=pat_mu, scale=pat_std)
    
    outdata = pd.DataFrame(columns = ['CohortSF', 'PatSF', 'LastVisit', 'DiffDays'],  
                           index = pd.Index(dates, name = 'Date'))
    
    try:
        ldate = inpat.iloc[-2]['Date']
        lvisit = inpat.index[-2][1]
    except IndexError:
        lvisit = 'R01'
        ldate = inpat.iloc[0]['Date']
    
    outdata['LastVisit'] = lvisit
    
    for date in dates:
       
        diff_date = (date - ldate).days

        outdata.ix[date]['CohortSF'] = cohort_norm[lvisit].sf(diff_date)
        outdata.ix[date]['PatSF'] = obj.sf(diff_date)
        outdata.ix[date]['DiffDays'] = diff_date
        
    return outdata

edata = odata.groupby(level=0).apply(expand_sf_data)

# <codecell>

X = ndata.dropna()[['CohortSF', 'PatSF']].values
y = ndata.dropna()['Returned'].values
predictor = AdaBoostClassifier().fit(X,y)

edata['LikelyReturn'] = predictor.predict(edata[['CohortSF', 'PatSF']].values)

# <codecell>

date_count = edata.groupby(level = 'Date')['LikelyReturn'].sum()
date_count.plot(figsize=(15,10))
plt.title('Returnable Cohort Size')
plt.xlabel('Starting Date')
plt.ylabel('Patients Likely To Return')

# <codecell>

print date_count.diff().mean(), 'Patients lost per month wasted!'

# <markdowncell>

# The above figure shows the number of patients that are predicted to return if called given a particular starting date. We can see that the longer we wait the less patients we can keep for 'longitudinal' natures. From this graph we can estimate that we lose ~1.5 patients per week that we don't start!

# <headingcell level=2>

# Make Call-Back Sheets

# <markdowncell>

# Here I want to make a set of sheets for the clinic to use to re-call patients. For each month I'll make a list (sorted by likelihood) of patients who are likely to return. I'll also mark which patients have 3+ which should be seen by the neurologist.

# <codecell>


pred_pat_data = edata.swaplevel(0, 1)
pred_pat_data['ProbReturn'] = predictor.predict_proba(pred_pat_data[['CohortSF', 'PatSF']].values)[:,1]

out_file = '/home/will/RedcapQC/CallbackPats/HIV_DrexelMed_GeneticAnalysisCohort_recall_list.xlsx'
writer = pd.ExcelWriter(out_file)
sheet_name = '%i-%i'
for tdate, rows in pred_pat_data.groupby(level=0):
    
    if tdate > datetime.now():
        srows = rows.sort(['ProbReturn', 'LastVisit'], ascending=False)
        srows['Neuro'] = ''
        srows['Neuro'][srows['LastVisit']>='R03'] = 'Prefered Neuropsych visit'
        rem_days = srows['DiffDays'] % 180
        month_mask = (rem_days < 31)
        
        tmp = srows[['Neuro']][month_mask].reset_index()
        print tdate, month_mask.sum()
        ndate = tdate.to_datetime()
        tmp[['Patient ID', 'Neuro']].to_excel(writer, 
                                              sheet_name % (ndate.year, ndate.month), 
                                              index=False)
writer.save()

# <codecell>

pd.DataFrame().to_excel?

# <codecell>

pd.DataFrame().to_excel

# <codecell>

tmp[['Patient ID', 'Neuro']].to_excel

# <codecell>

data['Date'].map(lambda x: x.month).value_counts()

# <codecell>


