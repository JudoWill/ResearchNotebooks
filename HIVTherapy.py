# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=4>

# IMPORT REDCAP DATA AND THE COLUMNS REGARDING CURRENT AND PAST ART TREATMENT

# <codecell>

from pandas import HDFStore #This is a file storage format for large collections of data
store = HDFStore('/home/will/HIVReportGen/Data/BaseRedcap/HIVAIDSGeneticAnalys_DATA_LABELS_2013-01-16_1211.hdf')
redcap_data = store['redcap']
store.close()

# <codecell>

print redcap_data
print redcap_data.columns

# <codecell>

therapy_data = redcap_data[["Current ART (choice='AZT')","Current ART (choice='ABC')","Current ART (choice='DVL')","Current ART (choice='ATV')",
"Current ART (choice='T-20')","Current ART (choice='3TC')","Current ART (choice='TDF')","Current ART (choice='SAQ')","Current ART (choice='AMP')",
"Current ART (choice='FPV')","Current ART (choice='DDI')","Current ART (choice='FTC')", "Current ART (choice='RTV')","Current ART (choice='LPV/r')", 
"Current ART (choice='DDC')","Current ART (choice='EFV')","Current ART (choice='NFL')","Current ART (choice='TPV')","Current ART (choice='D4T')",
"Current ART (choice='NVP')","Current ART (choice='IDV')","Current ART (choice='DRV')","Current ART (choice='Combivir')","Current ART (choice='Trizivir')", 
"Current ART (choice='Kaletra')","Current ART (choice='Epzicom')","Current ART (choice='Truvada')","Current ART (choice='Atripla')", "Current ART (choice='Other')",
"Current ART (choice='none')","Current ART (choice='ND')","Past ART (choice='AZT')","Past ART (choice='ABC')","Past ART (choice='DVL')","Past ART (choice='ATV')", 
"Past ART (choice='T-20')","Past ART (choice='3TC')","Past ART (choice='TDF')","Past ART (choice='SAQ')","Past ART (choice='AMP')","Past ART (choice='FPV')", 
"Past ART (choice='DDI')","Past ART (choice='FTC')","Past ART (choice='RTV')","Past ART (choice='LPV/r')","Past ART (choice='DDC')","Past ART (choice='EFV')", 
"Past ART (choice='NFL')","Past ART (choice='TPV')","Past ART (choice='D4T')","Past ART (choice='NVP')","Past ART (choice='IDV')","Past ART (choice='DRV')", 
"Past ART (choice='Combivir')","Past ART (choice='Trizivir')","Past ART (choice='Kaletra')","Past ART (choice='Epzicom')","Past ART (choice='Truvada')",
"Past ART (choice='Atripla')","Past ART (choice='Other')","Past ART (choice='none')","Past ART (choice='ND')"]]
print therapy_data

# <headingcell level=4>

# IMPORT INFORMATION THAT CONVERTS DRUGS TO DRUG CLASSES

# <rawcell>

# First import a file (saved on Will's super computer) that can be used to convert the names of the row colums, and make a new DataFrame(?) which groups the names into drug classes

# <codecell>

#therapy_data["Current ART (choice='T-20')"]

#redcap column headers are too long, want to rename them to refer to the drugs more explictly

new_headers = {
              


test_mask = therapy_data["Current ART (choice='ATV')"] == True
print test_mask

# <codecell>


