# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import os, os.path
import sys
sys.path.append('/home/will/PatientPicker/')
import LoadingTools

# <codecell>

import glob
import pandas as pd

files = glob.glob('/home/will/HIVReportGen/Data/PatientFasta/*LTR.fasta')
#redcap_data = LoadingTools.load_redcap_data()
has_ltr = []
for f in files:
    fname = f.rsplit('/',1)[-1]
    pid, vnum, _ = fname.split('-')
    has_ltr.append({
                    'PatientID':pid,
                    'VisitNum':vnum,
                    'HasLTR':'has_ltr'
                    })
has_ltr = pd.DataFrame(has_ltr).groupby(['PatientID', 'VisitNum']).first()

# <codecell>

redcap_data = LoadingTools.load_redcap_data().groupby(['Patient ID', 'VisitNum']).first()

# <codecell>

jean_comments = pd.read_csv('/home/will/HIVVariation/ProblemPCRpatientsamples.csv').groupby(['Patient ID', 'VisitNum']).first()

# <codecell>

left, right = jean_comments.align(has_ltr, join='outer')
ltr_comments = left.copy()
ltr_comments['HasLTR'] = ltr_comments['HasLTR'].combine_first(right['HasLTR'])

# <codecell>

red_data = pd.merge(ltr_comments, redcap_data,
                    left_index=True, 
                    right_index=True,
                    how='outer')
red_data

# <codecell>

group_key = 'HasLTR'
check_cols = ['Latest CD4 count (cells/uL)', 'Latest CD8 count (cells/uL)', 'LVL']

# <codecell>

red_data['LVL'].describe()

# <codecell>

from statsmodels.graphics.boxplots import violinplot
import numpy as np
red_data['LVL'] = red_data['Latest viral load'].map(np.log10)
fig, axs = plt.subplots(3,1, figsize=(10,10))

for col, ax in zip(check_cols, axs.flatten()):
    boxes = []
    labels = []
    for key, group in red_data.groupby(group_key):
        labels.append(key)
        print key, len(group[col].dropna().unique())
        boxes.append(group[col].dropna())
    #iolinplot(boxes, ax=ax, labels=labels)
    ax.boxplot(boxes)
    
    

# <codecell>

list(red_data.columns)

# <codecell>


