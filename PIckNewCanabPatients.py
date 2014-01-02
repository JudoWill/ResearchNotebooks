# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import pandas as pd
import sys
import os, os.path

sys.path.append('/home/will/PatientPicker/')

# <codecell>

import LoadingTools

# <codecell>

redcap_data = LoadingTools.load_redcap_data()


# <codecell>

test_cols = [col for col in redcap_data.columns if col.startswith('Test-')]


tmp = redcap_data[['Patient ID', 'Event Name']+test_cols].groupby('Patient ID').agg('count')
mask = (tmp[test_cols] == 0).all(axis = 1)

print mask.sum()
print tmp[mask]['Event Name'].sum()
print tmp[mask]['Event Name'].max()

# <codecell>

old_data = pd.read_csv('/home/will/HIVSystemsBio/NewCytokineAnalysis/CytoRawData.csv', sep='\t')

# <codecell>

old_pats = set(old_data['Patient ID'].unique())
already_done = redcap_data['Patient ID'].map(lambda x: x in old_pats)
nredcap_data = redcap_data[~already_done]
nredcap_data['IsMale'] = nredcap_data['Gender'] == 'Male'
nredcap_data['IsFemale'] = nredcap_data['Gender'] == 'Female'
print nredcap_data

# <codecell>

def safe_sum(inser):
    return inser.sum()

def safe_count(inser):
    return len(inser.dropna())

def safe_first(inser):
    return inser.first()

admit_cols = [col for col in nredcap_data.columns if col.startswith('Admit-')]
test_cols = [col for col in nredcap_data.columns if col.startswith('Test-')]
haart_cols = [col for col in nredcap_data.columns if col.startswith('HAART-')]
race_cols = [col for col in nredcap_data.columns if col.startswith('Race-')]
other_cols = ['Hepatitis C status (HCV)', 'IsMale', 'IsFemale', 'Age']
agg_dict = {'Event Name':safe_count,
            'TMHDS':safe_count,
            }
for col in admit_cols+test_cols+haart_cols+race_cols+other_cols:
    agg_dict[col] = safe_sum
    
agg_dict['Age'] = 'first'

wcols = ['Event Name', 'TMHDS']+haart_cols+test_cols+admit_cols+race_cols+other_cols
rename_dict = {'Event Name':'NumVisits', 'TMHDS':'Num-TMHDS'}
pat_data = nredcap_data.groupby('Patient ID').agg(agg_dict)[wcols].rename(columns=rename_dict)
for col in pat_data.columns:
    pat_data[col] = pat_data[col].map(float)

# <codecell>

num_visit = pat_data['NumVisits'] >= 3
num_tmhds = pat_data['Num-TMHDS'] >= 3
num_race = pat_data['Race-Black'] > 0
canab = pat_data['Test-Cannabinoid'] > 0
odrug_cols = [col for col in test_cols if 'Cannabinoid' not in col]
odrugs = pat_data[odrug_cols].sum(axis=1)==0
wanted_pats = pat_data[num_visit&num_tmhds&canab&odrugs&num_race]

# <codecell>

mask = nredcap_data['Test-Cannabinoid'] == True
first_canab = nredcap_data[mask][['Patient ID', 'VisitNum']]
canab_samples = pd.merge(wanted_pats, first_canab,
                         left_index=True,
                         right_on='Patient ID',
                         how='inner')[['Patient ID', 'VisitNum']]
canab_samples['Experiment'] = 'Canabinoid'
canab_samples['Priority'] = 0

# <codecell>

benj_x4 = [('A0017', 'R02', 'BenjX4'),
           ('A0107', 'R05', 'BenjX4'),
           ('A0208', 'R00', 'BenjX4'),
           ('A0403', 'R01', 'BenjX4')]
benj_samples = pd.DataFrame(benj_x4, columns=['Patient ID', 'VisitNum', 'Experiment']).set_index('Patient ID')
benj_samples['Priority'] = 1

# <codecell>

greg_x4 = [('A0004', 'R07', 'GregX4'),
           ('A0367', 'R05', 'GregX4'),
           ('A0023', 'R01', 'GregX4')]
#skipped A0208-R02 and A0017-R01 because they are long-visit
greg_samples = pd.DataFrame(greg_x4, columns=['Patient ID', 'VisitNum', 'Experiment']).set_index('Patient ID')
greg_samples['Priority'] = 2

# <codecell>

num_visit = pat_data['NumVisits'] >= 3
num_tmhds = pat_data['Num-TMHDS'] >= 3
num_race = pat_data['Race-Black'] > 0
on_haart = pat_data['HAART-On']==pat_data['NumVisits']
old_people = pat_data['Age']>=50
old_pats = pat_data[num_visit&num_tmhds&num_race&old_people&on_haart]

mask = nredcap_data['TMHDS'].notnull()
first_tmhds = nredcap_data[mask].groupby('Patient ID')[['VisitNum']].first()
old_samples = pd.merge(old_pats, first_tmhds,
                       right_index=True,
                       left_index=True,
                       how='inner')[['VisitNum', 'NumVisits']]
old_samples['Experiment'] = 'Aging'
old_samples['Priority'] = '3'
old_samples = old_samples.sort('NumVisits', ascending=False)[['VisitNum', 'Experiment', 'Priority']]
old_samples

# <codecell>

all_samples = pd.concat([canab_samples,
                         benj_samples.reset_index(),
                         greg_samples.reset_index(),
                         old_samples.reset_index(),
                         ], axis=0, ignore_index=True)
out_pats = all_samples.sort(['Priority', 'Patient ID'])
out_pats.drop('Priority', axis=1).to_excel('/home/will/RedcapQC/NewCytokineAnalysis_allvisits.xlsx', index=False)

# <codecell>

out_pats

# <codecell>


