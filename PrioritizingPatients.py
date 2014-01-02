# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

from pandas import *
import os, os.path
import sys

sys.path.append('/home/will/HIVReportGen/AnalysisCode/')

# <codecell>

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

# <codecell>

have_seq = seq_data[seq_fields].apply(lambda x: x.notnull()).fillna(False)
pat_fields = visit_data
all_fields = concat([pat_fields, have_seq], axis = 1)
all_fields['Predicted-R5'] = all_fields['Predicted-R5']>=0.8

# <codecell>

def check_fun(df):
    wanted_drugs = ["Current ART (choice='%s')" % d for d in ['TDF', 'Truvada', 'Atripla']] 
    start_niave = df['Current ART status'][0] == 'naive'
    on_therapy = (df['Current ART status'] == 'on').any()
    
    on_wanted = df[wanted_drugs].any().any()
    return start_niave & on_therapy & on_wanted



wanted_drugs = ["Current ART (choice='%s')" % d for d in ['TDF', 'Truvada', 'Atripla']] 
tdata = all_fields[['Current ART status'] + wanted_drugs]
res = tdata.groupby(level = 0).apply(check_fun)


# <codecell>

all_fields.index.names = ['Patient ID', 'Visit Number']
output = merge(all_fields[[]].reset_index(), DataFrame({'result':res}), left_on = 'Patient ID', right_index = True)
print output[['Patient ID', 'Visit Number', 'result']].head(n = 20).to_string()

# <codecell>

output.to_csv('/home/will/tmpstuf/drugged_data.csv')

# <codecell>


all_fields.fillna(False).to_csv('/home/will/HIVSystemsBio/NewPatientInfo_extreme.csv')

# <codecell>

ols?

# <codecell>

mask = redcap_data['Patient ID'] == 'A0008'
ofields = ['Latest viral load', 'Latest CD4 count (cells/uL)', 'Total Modified Hopkins Dementia Score']
other_fields = ['Gender', 'Current ART status', 'Age', 'Hepatitis C status (HCV)', 'Hepatitis B status (HBV)', 'Years seropositive', 'HIV seropositive date']
race_fields = ["Race (choice='Asian')",
"Race (choice='American Indian/Alaska Native')",
"Race (choice='Black or African American')",
"Race (choice='Native Hawaiian or other Pacific Islander')",
"Race (choice='White')",
"Race (choice='More than one race')",
"Race (choice='Unknown')",
]

drug_fields = [
'Amphetamines',
'Barbiturates',
'Benzodiazepines',
'Cannabinoid',
'Cocaine + metabolite',
'Opiates',
'Phencyclidine']

print redcap_data[['Patient visit number', 'Date of visit']+ other_fields][mask].to_string(), '\n\n\n\n'
print redcap_data[['Patient visit number', 'Date of visit']+ ofields][mask].to_string(), '\n\n\n\n'
print redcap_data[['Patient visit number', 'Date of visit']+ race_fields][mask].T.to_string(), '\n\n\n\n'
print redcap_data[['Patient visit number', 'Date of visit']+ drug_fields][mask].to_string(), '\n\n\n\n'

# <codecell>

t = redcap_data['Event Name'].apply(lambda x: int(x.split(' - ')[0][1:]))
t.unique()
redcap_data['VisitNum'] = redcap_data['VisitNum'].combine_first(t)

# <codecell>


t = all_fields['Event Name'].dropna().apply(lambda x: int(x.split(' - ')[0][1:]))

all_fields['VisitNum'] = all_fields['VisitNum'].combine_first(t)

# <codecell>

all_fields['Drug User Classification'].unique()

# <codecell>

drug_fields = [
'Cocaine + metabolite',
'Amphetamines',
'Barbiturates',
'Benzodiazepines',
'Cannabinoid',
'Opiates',
'Phencyclidine']

drug_fields[1:]

# <codecell>

drug_fields = [
'Cocaine + metabolite',
'Amphetamines',
'Barbiturates',
'Benzodiazepines',
'Cannabinoid',
'Opiates',
'Phencyclidine']

admit_fields = [
"Drugs used (choice='Marijuana')",
"Drugs used (choice='Cocaine (crack, nasal, smoke, inject)')",
"Drugs used (choice='Heroin (nasal, inject)')",
"Drugs used (choice='Methamphetamine (smoke, nasal, inject)')",
"Drugs used (choice='Benzodiazapine (i.e. valium, ativan, xanax, klonipin, etc)')",
"Drugs used (choice='Narcotics')",
"Drugs used (choice='Ecstasy')",
"Drugs used (choice='PCP')",
"Drugs used (choice='Ritalin')",
"Drugs used (choice='Other')"]

tmp = all_fields[drug_fields + admit_fields +['LTR']].reset_index()

def check_PN(df):
    any_pos = df[drug_fields].any().any()
    any_admit = df[admit_fields].any().any()
    return (any_admit | any_pos)

def check_PC(df):
    pos_coc = df[drug_fields[0]].any()
    pos_other = df[drug_fields[1:]].any().any()
    return pos_coc and ~pos_other

def check_mdu(df):
    num_pos = df[drug_fields].any().sum()
    return num_pos > 1
    
def check_ltr(df):
    return df['LTR'].values[-1]

#print tmp
checks = {'LTR': check_ltr,
          'PN': check_PN,  
           'PC': check_PC,         
           'MDU': check_mdu,}
nchecks = list(checks.items())
res = []
valid_visits = tmp['Visit Number']=='A'
for visit in range(10):
    visit_str = 'R%02i' % visit
    visit_mask = tmp['Visit Number'] == visit_str
    valid_visits |= visit_mask
    res.append(('#Patients', visit_str, visit_mask.sum()))
    
    ntmp = tmp.ix[valid_visits]
    pats = ntmp.groupby('Patient ID')
    for pat, ndf in pats:
        for name, func in checks.items():
            nres = func(ndf)
            print nres
            raise KeyError
    
#df = DataFrame(res, columns = ['Header', 'VisitNum', 'Value'])
#res = pivot_table(df, rows = ['VisitNum'], cols='Header', values= 'Value')
#print res

# <codecell>

tmp = read_csv('/home/will/HIVSystemsBio/NewCytokineAnalysis/CytoPatData.csv', sep = '\t')
wanted_pats = tmp['Patient ID']

wanted_data = {}
wanted_visits = dict([(p, v) for p,v in zip(tmp['Patient ID'].values, tmp['VisitNum'].values)])
for key, group in redcap_data.groupby('Patient ID'):
    if key in wanted_visits:
        vname = wanted_visits[key]
        wnum = int(vname[1:])
        wdata = group['VisitNum']<= wnum
        res = group[drug_fields].ix[wdata].mean(axis = 0)
        wanted_data[key] = res

print wanted_data.keys()[:5]
drug_mean = DataFrame(wanted_data).T.rename(columns = dict([(col, 'TOSample-'+col) for col in drug_fields]))

drug_mean.ix[wanted_pats].to_csv('/home/will/HIVSystemsBio/NewCytokineAnalysis/ToSampledrug.csv')

# <codecell>

from itertools import groupby
import csv


def missing_test(visit_nums, visit_dates, check_ser):
    
    for v, date, val in zip(visit_nums, visit_dates, check_ser):
        if val != val:
            yield v, date, 'Missing Value', 1
    
def consistency_test(visit_nums, visit_dates, check_ser):
    
    #print t
    if len(check_ser.dropna().unique())>1:
        for v, date, val in zip(visit_nums, visit_dates, check_ser):
            yield v, date, 'Inconsitent Value', 1

def diagnose_test(visit_nums, visit_dates, check_ser, debug = False):
    
    tmp = DataFrame({'Visit':visit_nums, 'Date':visit_dates, 'Check':check_ser}).dropna()
    #print tmp
    tmp.sort(columns = 'Date')
    
    is_sick = False
    for _, row in tmp.iterrows():
        if (row['Check'] == False) and (is_sick == True):
            yield row['Visit'], row['Date'], 'Inconsistent Diagnosis', 1
        is_sick |= row['Check']==1
    
    

def nearby_date(check_dates, visit_dates):
    
    (check_dates - visit_dates).weeks


with open('/home/will/tmpstuf/test_smells.csv') as handle:
    junk = handle.next()
    check_rules = [row for row in csv.reader(handle, delimiter = '\t') if row[3].strip()]
    
messages = []
for patid, df in redcap_data.groupby('Patient ID'):
    for col, report_col, _, testfun in check_rules:
        
        if (testfun == 'consistency_test') or (testfun == 'date_consistency'):
            msgs = list(consistency_test(df['Patient visit number'], df['Date of visit'], df[col]))
        elif testfun == 'diagnose_test':
            #if col == 'Hepatitis C status (HCV)':
                #print col, df[col]
                #print len(list(diagnose_test(df['Patient visit number'], df['Date of visit'], df[col], debug = True)))
                #raise KeyError
            msgs = list(diagnose_test(df['Patient visit number'], df['Date of visit'], df[col]))
            
        else:
            msgs = list(missing_test(df['Patient visit number'], df['Date of visit'], df[col]))
            
        for v, date, msg, score in msgs:
            messages.append((col, report_col, patid, v, date, msg, score))

# <codecell>

tmp = DataFrame(messages, columns = ['Colname', 'Grouping', 'Patient ID', 'Visit', 'VisitDate', 'Message', 'Wrongness'])

print tmp.head(n= 100).to_string()


# <codecell>

res = pivot_table(tmp, rows = 'VisitDate', cols = 'Message', values = 'Wrongness', aggfunc=np.sum)
#res['Inconsitent Value'].dropna()
plt.figure(figsize = (10,10))
rolling_mean(res, 30, min_periods=2).plot(ax = plt.gca())

# <codecell>

tmp.groupby(['Patient ID']).sum().min()

# <codecell>

redcap_data['Hepatitis C status (HCV)'].dropna()

# <codecell>


