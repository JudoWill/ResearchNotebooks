# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import sys
sys.path.append('/home/will/PatientPicker/')
import pandas as pd

# <codecell>

import LoadingTools

store = pd.HDFStore('/home/will/HIVReportGen/Data/SplitRedcap/2013-01-16/EntireCohort.hdf')
pat_data = store['redcap']

# <codecell>

drug_cols = ['Current Tobacco use',
             'Current Alcohol use',
             'Cocaine + metabolite',
             'Amphetamines',
             'Barbiturates',
             'Benzodiazepines',
             'Cannabinoid',
             'Opiates',
             'Phencyclidine',
             "Drugs used (choice='Ritalin')"]
gender_cols = ['Gender']
age_cols = ['Age', 'Calc-Years-Seropositive']
haart_cols = ['Current ART status']
race_cols = [col for col in pat_data.columns if col.startswith('Race')]
eth_cols = ['Ethnicity']

wanted_cols = drug_cols+gender_cols+age_cols+haart_cols+race_cols+eth_cols+['Patient ID', 'Patient visit number']
spat_data = pat_data[wanted_cols].sort(['Patient ID', 'Patient visit number'])

# <codecell>

intake_mask = spat_data['Patient visit number'] == 'R00'
intake_data = spat_data[intake_mask].set_index('Patient ID').drop(['Patient visit number'], axis = 1)

pc_pats = ['A0022','A0025','A0029','A0039',
            'A0040','A0047','A0056','A0058',
            'A0068','A0083','A0091','A0106',
            'A0124','A0136','A0142','A0151',
            'A0175','A0181','A0191','A0208',
            'A0209','A0262','A0284','A0313',
            'A0379','A0388','A0427']
pn_pats = ['A0010','A0017','A0032','A0078',
           'A0100','A0159','A0195','A0206',
           'A0217','A0220','A0223','A0238',
           'A0239','A0240','A0242','A0255',
           'A0258','A0280','A0294','A0321',
           'A0339','A0356','A0363','A0376',
           'A0380','A0393','A0397','A0405',
           'A0406','A0415','A0440','A0447',
           'A0456']

pc_data = intake_data.ix[pc_pats]
pn_data = intake_data.ix[pn_pats]

# <codecell>

from itertools import product
from copy import deepcopy

def calc_gender(indata):
    tdict = {
             'Male':(indata['Gender']=='Male').sum(),
             'Female':(indata['Gender']=='Female').sum(),
             }
    return tdict

def calc_race(indata):
    wcol = "Race (choice='White')"
    bcol = "Race (choice='Black or African American')"
    ucol = "Race (choice='Unknown')"
    ocols = ["Race (choice='Asian')",
             "Race (choice='American Indian/Alaska Native')",
             "Race (choice='Native Hawaiian or other Pacific Islander')",
             "Race (choice='More than one race')"]
 
    tdict = {
             'White': indata[wcol].sum(),
             'Black/AA': indata[bcol].sum(),
             'Unknown': indata[ucol].sum(),
             'Other': indata[ocols].any(axis=1).sum()
             }
    return tdict

def calc_drug_use(indata):
    
    groups = [('Tobacco', indata['Current Tobacco use']=='Yes'),
              ('Alcohol', indata['Current Alcohol use']=='Yes'),
              ('Cocaine', indata['Cocaine + metabolite']),
              ('Cannabinoids', indata['Cannabinoid']),
              ('Methamphetamines', indata['Amphetamines']),
              ('Benzodiazepines', indata['Benzodiazepines']),
              ('Narcotics', indata['Opiates']),
              ('Ritalin', indata["Drugs used (choice='Ritalin')"])]
    tdict = {}
    for key, col in groups:
        tdict[key] = col.sum()
    return tdict

def calc_eth(indata):
    
    cols = ['Not Hispanic or Latino', 'Hispanic or Latino']
    tdict = {}
    for col in cols:
        tdict[col] = (indata['Ethnicity']==col).sum()
    tdict['Unknown'] = indata['Ethnicity'].isnull().sum()
    return tdict


def calc_haart(indata):
    
    col = 'Current ART status'
    tdict = {
             'cH': (indata[col] == 'on').sum(),
             'dH': ((indata[col] == 'off') | (indata[col] == 'non-adherent')).sum(),
             'nH': (indata[col] == 'naive').sum(),
             }
    return tdict

def calc_age(indata):
    
    return {'':indata['Age'].mean()}

def calc_sero(indata):
    
    return {'':indata['Calc-Years-Seropositive'].mean()}

anal_list = [('Gender', calc_gender),
             ('Race', calc_race),
             ('Ethnicity', calc_eth),
             ('Drug Use', calc_drug_use),
             ('HAART', calc_haart),
             ('Age', calc_age),
             ('Years Seropsitive', calc_sero)]

groups = [('All', intake_data), 
          ('PN', pn_data), 
          ('PC', pc_data)]
out_groups = []
for (gname, group), (varname, func) in product(groups, anal_list):
    
    odict = func(group)
    for key, val in odict.items():
        out_groups.append({
                           'group':gname,
                           'varname':varname,
                           'itemname':key,
                           'itemval':val
                           })
        
    
    

# <codecell>

ndf = pd.DataFrame(out_groups)
table_data = pd.pivot_table(ndf, 
                            rows = ['varname', 'itemname'], 
                            cols = 'group',
                            values = 'itemval')
print table_data
table_data.to_excel('/home/will/quickBrain/pat_groups.xlsx')

# <codecell>

drug_cols = ['Current Tobacco use',
             'Current Alcohol use',
             'Cocaine + metabolite',
             'Amphetamines',
             'Barbiturates',
             'Benzodiazepines',
             'Cannabinoid',
             'Opiates',
             'Phencyclidine']
intake_data['Ethnicity'].unique()

# <codecell>

pat_data['Patient ID'].unique()

# <codecell>


