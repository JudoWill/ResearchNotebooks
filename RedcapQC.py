# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import pandas as pd
import os, os.path
import numpy as np
os.chdir('/home/will/RedcapQC/')

# <codecell>

raw_data = pd.read_csv('/home/will/HIVReportGen/Data/RedcapDumps/HIVAIDSGeneticAnalys_DATA_LABELS_2013-01-16_1211.csv')
raw_data.rename(columns={raw_data.columns[0]:'Patient ID'}, inplace=True)

# <codecell>

from collections import Counter

def get_common(inser):
    counts = Counter(inser.values)
    return counts.most_common(1)[0][0]

ident_cols = ['Gender', 'Year of Birth', 'Ethnicity']
ident_cols += [col for col in raw_data.columns if col.startswith('Race ')]
ident_cols += [col for col in raw_data.columns if col.startswith('Exposure ')]


#PatID, VisitNum, VisitDate, Field, ActualValue, ExpectedValue
common_fixing = []
all_count = 0
for pat_id, group in raw_data.groupby('Patient ID'):
    
    for field in ident_cols:
        if len(group[field].dropna()) == 0:
            common = np.nan
        else:
            common = get_common(group[field].dropna())
            
        for _, row in group[['Event Name', field, 'Date of visit']].iterrows():
            all_count += 1
            if row[field] != common:
                expected = 'missing' if row[field] != row[field] else common
                if expected == 'missing':
                    common_fixing.append({
                                          'Patient ID':pat_id,
                                          'VisitNum':row['Event Name'].split()[0],
                                          'VisitDate':row['Date of visit'],
                                          'FieldToFix':field,
                                          'ActualValue':'missing',
                                          'ExpectedValue':'a valid value'
                                          })
                else:
                    common_fixing.append({
                                          'Patient ID':pat_id,
                                          'VisitNum':row['Event Name'].split()[0],
                                          'VisitDate':row['Date of visit'],
                                          'FieldToFix':field,
                                          'ActualValue':row[field],
                                          'ExpectedValue':expected
                                          })


# <codecell>

from dateutil import parser
from copy import deepcopy

date_fields = [('CD4', 'Date of latest CD4 count'),
               ('CD8', 'Date of latest CD8 count'),
               ('VL', 'Date of latest viral load')]

#PatID, VisitNum, VisitDate, Field, ActualValue, ExpectedValue
param_fixing = []
for _, row in raw_data.iterrows():
    
    for name, field in date_fields:
        tdict = {
                 'Patient ID':row['Patient ID'],
                 'VisitNum':row['Event Name'].split()[0],
                 'ActualValue':None,
                 'ExpectedValue':None
                 }
        date_dict = {}
        for nfield in [field, 'Date of visit']:
            all_count += 1
            if row[nfield] != row[nfield]:
                date_dict[nfield] = None
                ndict = {
                         'VisitDate':row['Date of visit'],
                         'FieldToFix':nfield,
                         'ActualValue':'Missing!!',
                         'ExpectedValue':'Anything'
                     }
                param_fixing.append(deepcopy(tdict).update(ndict))
            else:
                try:
                    date_dict[nfield] = parser.parse(row[nfield])
                except ValueError:
                    date_dict[nfield] = None
                    ndict = {
                         'VisitDate':row['Date of visit'],
                         'FieldToFix':nfield,
                         'ActualValue':row[nfield],
                         'ExpectedValue':'Invalid Date Format!!'
                         }
                    param_fixing.append(deepcopy(tdict).update(ndict))
            
        vdate, fdate = date_dict['Date of visit'], date_dict[field]
        
            
        if (vdate is not None) and (fdate is not None):
            dif_date = abs((vdate - fdate).days/30)
            if dif_date > 3:
                param_fixing.append({
                                 'Patient ID':row['Patient ID'],
                                 'VisitNum':row['Event Name'].split()[0],
                                 'VisitDate':row['Date of visit'],
                                 'FieldToFix':name  + ' too distantly measured!',
                                 'ActualValue':'%i months' % dif_date,
                                 'ExpectedValue':row[field]
                                 })
        
        

# <codecell>

paired_cols = [('Hepatitis B status (HBV)','Year diagnosed HBV positive'),
               ('Hepatitis C status (HCV)','Year diagnosed HCV positive'),
               ('Cytomegalovirus (CMV)','Year diagnosed CMV positive'),
               ('Human Papillomavirus (HPV)','Year diagnosed HPV positive'),
               ('Herpes Simplex Virus Type 1 (HSV 1)','Year diagnosed HSV 1 positive'),
               ('Herpes Simplex Virus Type 2 (HSV 2)','Year diagnosed HSV 2 positive'),
               ('Tuberculosis','Year diagnosed tuberculosis positive'),
               ('Hypertension','Year diagnosed with hypertension'),
               ('Diabetes','Year diagnosed with diabetes'),
               ('Elevated lipids','Year diagnosed with elevated lipids'),
               ('Asthma','Year diagnosed with asthma'),
               ('Chronic obstructive pulmonary disease (COPD)','Year diagnosed with COPD')]


for pat_id, group in raw_data.groupby('Patient ID'):
    
    for data_col, date_col in paired_cols:
        group[data_col] = group[data_col].replace('ND', np.nan)
        group[date_col] = group[date_col].replace('ND', np.nan)
            
        for field in [data_col, date_col]:
            if len(group[field].dropna()) == 0:
                common = np.nan
            else:
                common = get_common(group[field].dropna())
            if common != common:
                continue
            for _, row in group[['Event Name', field, 'Date of visit']].iterrows():
                all_count += 1
                if row[field] != common:
                    if row[field] != row[field]:
                        actual = 'missing'
                    else:
                        actual = row[field]
                    common_fixing.append({
                                          'Patient ID':pat_id,
                                          'VisitNum':row['Event Name'].split()[0],
                                          'VisitDate':row['Date of visit'],
                                          'FieldToFix':field,
                                          'ActualValue':actual,
                                          'ExpectedValue':common
                                          })


    
    

# <codecell>

from operator import itemgetter

sort_fields = ['Patient ID', 'VisitNum', 'VisitDate', 'FieldToFix']
ordered_fixes = sorted([row for row in (param_fixing + common_fixing) if row is not None], key = itemgetter(*sort_fields))

# <codecell>


fields = ['Patient ID', 'VisitNum', 'VisitDate', 'FieldToFix', 'ActualValue', 'ExpectedValue']

df = pd.DataFrame(ordered_fixes)
#df[fields].to_excel('/home/will/HIVVariation/redcap_QC.xlsx', index = False)

# <codecell>

for pat, rows in df.groupby('Patient ID'):
    rows.to_excel('PatsToFix/%s-%i.xlsx' % (pat, len(rows)), index=False)
    

# <codecell>

num_pats = len(list(df.groupby('Patient ID')))

# <codecell>

num_pats/5

# <codecell>

num_pats

# <codecell>

print all_count, float(len(ordered_fixes))/float(all_count)

# <codecell>


