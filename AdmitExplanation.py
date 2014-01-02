# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import pandas as pd
import sys
import os, os.path

sys.path.append('/home/will/PatientPicker/')
import LoadingTools

# <codecell>

redcap= LoadingTools.load_redcap_data().set_index(['Patient ID', 'VisitNum'])

# <codecell>

admit_no_drugs = [('Current-Drug-Use-NO', redcap['Current Drug use']=='No'),
                  ('Current-Drug-Use-NEVER', redcap['Current Drug use']=='Never'),
                  ('Date-Stopped-Drug-Use', redcap['Date Stopped Drug Use']<redcap['Date Of Visit']),
                  ('Drug-Use-And-HIV-Status-BEFORE', redcap['Drug Use And HIV Status']=='Used after HIV+')]
test_cols = [col for col in redcap.columns if col.startswith('Test-')]
admit_cols = [col for col in redcap.columns if (col.startswith('Admit-') and ('None' not in col))]
ever_test = redcap[test_cols].groupby(level='Patient ID').agg('any')

admit_no_drug_df = pd.concat([redcap[admit_cols], pd.DataFrame(dict(admit_no_drugs))], axis = 1)

# <codecell>

tmp = pd.concat(admit_no_drug_df.align(ever_test, axis=0, level='Patient ID'), 1)
tmp.iloc[100:105].T

# <codecell>

def fix_num(num):
    if num < 1:
        return -1/num
    else:
        return num

res = []
tmp['Test-Anything'] = tmp[test_cols].any(axis=1)
ntest_cols = [col for col in tmp.columns if col.startswith('Test-')]
nadmit_cols = [col for col in tmp.columns if not col.startswith('Test-')]
for col in nadmit_cols:
    say_yes = tmp[col]
        
    yes_frac = tmp[ntest_cols][say_yes].mean()
    no_frac = tmp[ntest_cols][~say_yes].mean()
    odds_r = yes_frac/no_frac
    
    res.append(odds_r.to_dict())
    res[-1]['Col'] = col
    
nres = pd.DataFrame(res).set_index('Col').applymap(fix_num)

# <codecell>

print nres

# <codecell>

#nres.to_excel('/home/will/DrugStuff/admit_explanations.xlsx')

# <codecell>

from sklearn.cross_validation import cross_val_score, StratifiedShuffleSplit
from sklearn.naive_bayes import GaussianNB
from sklearn.linear_model import LogisticRegression
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import AdaBoostClassifier


X = tmp[nadmit_cols]
y_mat = tmp[ntest_cols]

classifiers = [('Naive-Bayes', GaussianNB()),
               ('Logistic-Regression', LogisticRegression()),
               ('Decision-Tree', DecisionTreeClassifier()),
               ('Adaboost', AdaBoostClassifier()),
               ]
tres = []
for col in ntest_cols:
    if len(y_mat[col].dropna().unique()) < 2:
        continue
    for name, cl in classifiers:
        
        vals= cross_val_score(cl, X.values, y=y_mat[col].values, 
                        cv = StratifiedShuffleSplit(y_mat[col].values),
                        scoring = 'roc_auc')
        tres.append({'Col':col, 
                     'Predictor': name, 
                     'Accuracy':np.mean(vals)})
tdf = pd.pivot_table(pd.DataFrame(tres),rows = 'Col', cols = 'Predictor', values='Accuracy')

# <codecell>

tdf

# <codecell>


