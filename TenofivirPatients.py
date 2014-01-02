# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

from pandas import *
import os, os.path
import sys

sys.path.append('/home/will/HIVReportGen/AnalysisCode/')
store = HDFStore('/home/will/HIVReportGen/Data/SplitRedcap/2013-01-16/EntireCohort.hdf')

# <codecell>

redcap_data = store['redcap']
seq_data = store['seq_data']
visit_data = store['visit_redcap']
pat_data = store['pat_redcap']

# <codecell>

ofields = ['Latest viral load', 'Latest CD4 count (cells/uL)', 'Total Modified Hopkins Dementia Score']
wanted_fields = ['CalcAge', 'Gender', 'Drug User Classification', 'Hepatitis C status (HCV)', 'Predicted-R5']
seq_fields = ['LTR', 'Vpr', 'Tat', 'V3']

have_seq = seq_data[seq_fields].apply(lambda x: x.notnull()).fillna(False)
pat_fields = visit_data
all_fields = concat([pat_fields, have_seq], axis = 1)
all_fields['Predicted-R5'] = all_fields['Predicted-R5']>=0.8

# <codecell>

def check_for_tenof(df):
    wanted_drugs = ["Current ART (choice='%s')" % d for d in ['TDF', 'Truvada', 'Atripla']] 
    start_niave = df['Current ART status'][0] == 'naive'
    on_therapy = (df['Current ART status'] == 'on').any()
    
    on_wanted = df[wanted_drugs].any().any()
    return start_niave & on_therapy & on_wanted



wanted_drugs = ["Current ART (choice='%s')" % d for d in ['TDF', 'Truvada', 'Atripla']] 
tdata = all_fields[['Current ART status'] + wanted_drugs]
isTenof = tdata.groupby(level = 0).apply(check_for_tenof)

# <codecell>

def tenof_visits(df):
    
    wanted_drugs = ["Current ART (choice='%s')" % d for d in ['TDF', 'Truvada', 'Atripla']] 
    start_niave = df['Current ART status'][0] == 'naive'
    on_therapy = (df['Current ART status'] == 'on').any()
    
    on_wanted = df[wanted_drugs].any().any()
    if start_niave & on_therapy & on_wanted:
        t = df[wanted_drugs].any(axis = 1)
        tmp = iter(zip(t.index, t.values))
        
        before_first_ten_visit = None
        last_ten_visit = None
        for (pat, visit), on_ten in tmp:
            if ~on_ten:
                before_first_ten_visit = visit
            else:
                break
        
        for (pat, visit), on_ten in tmp:
            if on_ten:
                last_ten_visit = visit
            else:
                break
        if last_ten_visit is None:
            last_ten_visit = visit
        return Series([before_first_ten_visit, last_ten_visit], index = ['LastNaive', 'LastTen'])
    else:
        return Series([np.nan, np.nan], index = ['LastNaive', 'LastTen'])
    

wanted_drugs = ["Current ART (choice='%s')" % d for d in ['TDF', 'Truvada', 'Atripla']] 
tdata = all_fields[['Current ART status'] + wanted_drugs]
wanted_samples = tdata.groupby(level = 0).apply(tenof_visits)

# <codecell>

def check_well_controlled(df):
    
    ndf = df.dropna()
    cd4_good = ndf['Latest CD4 count (cells/uL)']>=250
    vl_good = ndf['Latest viral load']<=100
    on_drugs = ndf['Current ART status'] == 'on'
    
    valid_visits = (on_drugs & cd4_good & vl_good)
    consecutive_visits = 0
    best_num = 0
    for n in valid_visits.values:
        if n:
            consecutive_visits += 1
        else:
            best_num = max(best_num, consecutive_visits)
            consecutive_visits = 0
            
    best_num = max(consecutive_visits, best_num)
    
    return best_num

tdata = all_fields[['Current ART status', 'Latest CD4 count (cells/uL)', 'Latest viral load']]
num_controlled = tdata.groupby(level = 0).apply(check_well_controlled)    
    

# <codecell>

morbids = [
'Hepatitis B status (HBV)',
'Hepatitis C status (HCV)',
'Cytomegalovirus (CMV)',
'Human Papillomavirus (HPV)',
'Herpes Simplex Virus Type 1 (HSV 1)',
'Herpes Simplex Virus Type 2 (HSV 2)',
'Tuberculosis',
'Hypertension',
'Diabetes',
'Elevated lipids',
'Asthma',
'Chronic obstructive pulmonary disease (COPD)']

drug_use = [
'Amphetamines',
'Barbiturates',
'Benzodiazepines',
'Cannabinoid',
'Cocaine + metabolite',
'Opiates',
'Phencyclidine']

art = [
"Current ART (choice='AZT')",
"Current ART (choice='ABC')",
"Current ART (choice='DVL')",
"Current ART (choice='ATV')",
"Current ART (choice='T-20')",
"Current ART (choice='3TC')",
"Current ART (choice='TDF')",
"Current ART (choice='SAQ')",
"Current ART (choice='AMP')",
"Current ART (choice='FPV')",
"Current ART (choice='DDI')",
"Current ART (choice='FTC')",
"Current ART (choice='RTV')",
"Current ART (choice='LPV/r')",
"Current ART (choice='DDC')",
"Current ART (choice='EFV')",
"Current ART (choice='NFL')",
"Current ART (choice='TPV')",
"Current ART (choice='D4T')",
"Current ART (choice='NVP')",
"Current ART (choice='IDV')",
"Current ART (choice='DRV')",
"Current ART (choice='Combivir')",
"Current ART (choice='Trizivir')",
"Current ART (choice='Kaletra')",
"Current ART (choice='Epzicom')",
"Current ART (choice='Truvada')",
"Current ART (choice='Atripla')",
"Current ART (choice='Other')",
"Current ART (choice='none')",
"Current ART (choice='ND')",
]

def check_morbid(df):
    ndf = df.dropna()
    return df.mean()

morbid_res = all_fields[morbids].groupby(level = 0).agg(check_morbid)
drug_res = all_fields[drug_use].groupby(level = 0).agg(check_morbid)
art_res = all_fields[art].groupby(level = 0).agg(check_morbid)
other_res = all_fields[seq_fields+ofields].groupby(level = 0).agg(check_morbid)

# <codecell>

outdata = concat([wanted_samples, DataFrame({'NumControlled':num_controlled}), morbid_res, drug_res, other_res, art_res], axis = 1)
outdata[isTenof].to_csv('/home/will/tmpstuf/tenofivir_pats.csv', sep = '\t')

# <codecell>


# <codecell>


