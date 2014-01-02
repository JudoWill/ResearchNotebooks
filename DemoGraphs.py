# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

from pandas import read_csv
import os, os.path
import csv
import matplotlib.pyplot as plt
import re
import dateutil.parser

os.chdir('/home/will/HIVReportGen/')

# <codecell>

def extract_YOB(inp):
    try:
        return float(inp.split('-')[0])
    except AttributeError:
        return float(inp)
    except ValueError:
        #print('Bad YOB', inp)
        return np.nan
    
def safe_float(m, default = np.nan):
    try:
        return float(m)
    except:
        return default
    
def feet2meters(height):
    
    if (height == 'ND') or (height != height):
        return np.nan
    try:
        res = re.findall('(\d).\s{0,1}(\d{0,2})\D?', height)
    except TypeError:
        #print(height)
        raise TypeError
    try:
        ft = float(res[0][0])
        inc = safe_float(res[0][1], default = 0.0)
    except IndexError:
        #print(height,res)
        raise IndexError
    except ValueError:
        #print(height, res)
        raise ValueError
        
    tot_inches = ft*12+inc
    meters = tot_inches*0.0254
    if meters > 2:
        print(meters, height, res)
    
    return meters

def checkbox_conv(inp):
    
    if inp != inp:
        return np.nan
    valdict = {
                'checked':True,
                'test positive':True,
                'positive':True,
                'yes':True,
                'unchecked':False,
                'test negative':False,
                'negative':False,
                'no':True}
    return valdict.get(inp.lower(), np.nan)
    

def verbose_parser(inp):
    try:
        return dateutil.parser.parse(inp)
    except:
        return np.nan
        
def fix_col_name(name):
    
    if "(choice='" in name:
        return name.split("(choice='",1)[1][:-2]
    else:
        return name

# <codecell>

from datetime import date, datetime
from pandas import merge
from copy import deepcopy

class PatData(object):
    
    def __init__(self, redcap_file, config_file):
        
        if redcap_file is None:
            return           
            
        self.config_data = read_csv(config_file, sep = '\t')
        self._generate_converter_dict()
        self._generate_agg_dict()
        with open(redcap_file) as handle:
            handle.read(1)
            self.redcap = read_csv(handle, converters=self.conv_dict)
        self.clip_dates()
        
        self.visit_redcap = None
        self.pat_redcap = None
    
    
    def CopyFromOtherData(self, OtherData):
        
              
        self.config_data = OtherData.config_data.copy()
        self.redcap = OtherData.redcap.copy()
        self.conv_dict = deepcopy(OtherData.conv_dict)
        self.date_clip_cols = deepcopy(OtherData.date_clip_cols)
        self.visit_agg = deepcopy(OtherData.visit_agg)
        self.pat_agg = deepcopy(OtherData.pat_agg)
        if OtherData.pat_redcap is not None:
            self.pat_redcap = OtherData.pat_redcap.copy()
        if OtherData.visit_redcap is not None:
            self.visit_redcap = OtherData.visit_redcap.copy()
        if OtherData.all_group is not None:
            self.all_group = OtherData.all_group.copy()
    
    def _generate_converter_dict(self):
        cdict = {
            'DateParser':verbose_parser,
            'extract_YOB':extract_YOB,
            'checkbox_conv':checkbox_conv,
            'safe_float':safe_float
            }
        conv_dict = {}
        date_clip_cols = set()
        for colname, convfun in zip(self.config_data['RawName'].values, self.config_data['ConvertFun'].values):
            cfun = cdict.get(convfun, None)
            if cfun:
                conv_dict[colname] = cfun
            if convfun == 'DateParser':
                date_clip_cols.add(colname)
        
        self.conv_dict = conv_dict
        self.date_clip_cols = date_clip_cols
    
    def _generate_agg_dict(self):
        
        self.visit_agg = {}
        self.pat_agg = {}
        for colname, aggfun in zip(self.config_data['RawName'].values, self.config_data['AggFunc'].values):
            if aggfun == aggfun:
                self.visit_agg[colname] = aggfun
                self.pat_agg[colname] = aggfun
                
    def clip_dates(self):
        
        maxday = datetime.today()
        minday = datetime(1900,1,1)
        for col in self.date_clip_cols:
            self.redcap[col] = self.redcap[col].clip(lower = minday, upper = maxday)
    
    def fix_visits(self):
        
        def fix_v(vis):
            if vis == 'first':
                return 0.0
            try:
                return float(vis[1:])
            except:
                return None
        
        self.redcap['VisitNum'] = self.redcap['Patient visit number'].apply(fix_v)
    
    def CalcAge(self):
        visit_years = self.redcap['Date of visit'].dropna().apply(lambda x:x.year)
        birth_years = self.redcap['Year of Birth'].dropna()
        self.redcap['CalcAge'] = visit_years-birth_years
            
    def CalcGender(self):
        
        self.redcap['GenotypicMale'] = (self.redcap['Gender'] == 'Male') | (self.redcap['Transgender designation'] == 'male to female')
        self.redcap['IdentifiesMale'] = (self.redcap['Gender'] == 'Male') & (self.redcap['Transgender designation'] != 'male to female')
    
        self.redcap['GenotypicFemale'] = (self.redcap['Gender'] == 'Female') | (self.redcap['Transgender designation'] == 'female to male')
        self.redcap['IdentifiesFemale'] = (self.redcap['Gender'] == 'Female') & (self.redcap['Transgender designation'] != 'female to male')
    
    def CalcBMI(self):
        
        self.redcap['Weight-kg'] = self.redcap['Weight'].apply(safe_float)/2.2
        self.redcap['Height-m'] = self.redcap['Height'].apply(feet2meters)
        self.redcap['BMI'] = self.redcap['Weight-kg']/(self.redcap['Height-m']*self.redcap['Height-m'])
    
    def CalcYearsSero(self):
        
        visit_years = self.redcap['Date of visit'].dropna().apply(lambda x:x.year)
        seropos_years = self.redcap['HIV seropositive date'].dropna().apply(lambda x:x.year)
        self.redcap['Calc-Years-Seropositive'] =  visit_years - seropos_years
    
    def CalcExposure(self):
        
        merge_cols = {'Exposure-MSM': ("Exposure Category (choice='Men who have sex with men (MSM)')",
                                        "Exposure Category (choice='MSM and IDU')"),
                        'Exposure-IDU': ("Exposure Category (choice='Injection drug use (IDU)')",
                                        "Exposure Category (choice='MSM and IDU')",
                                        "Exposure Category (choice='Heterosexual and IDU')"),
                        'Exposure-Heterosexual': ("Exposure Category (choice='Heterosexual and IDU')",
                                                   "Exposure Category (choice='Heterosexual')"),
                        'Exposure-Hemophilia':("Exposure Category (choice='Hemophilia')",),
                        'Exposure-Transfusion':("Exposure Category (choice='Blood transfusion')",),
                        'Exposure-Perinatal':("Exposure Category (choice='Perinatal')",)
                    }

        for merged_col, check_cols in merge_cols.items():
            self.redcap[merged_col] = False
            for col in check_cols:
                self.redcap[merged_col] |= self.redcap[col]
    
    def AddGroupNames(self):
        
        self.groupnames = dict(zip([True, False], ['PosGroup', 'NegGroup']))
    
    def CalcAll(self):
        
        self.AddGroupNames()
        self.fix_visits()
        self.CalcAge()
        self.CalcYearsSero()
        self.CalcGender()
        self.CalcBMI()
        self.CalcExposure()
    
    def ProcessVisits(self, visit_recap):
        """A method to subclass. Must return a DataFrame of the wanted visits."""
        return visit_recap
    
    def ProcessPatients(self, pat_redcap):
        """A method to subclass. Must return a DataFrame of the wanted patients."""
        return pat_redcap
    
        
    def ProcessRedcap(self):
        
        gkey = ['Patient ID', 'Patient visit number']
        visit_redcap = self.redcap.groupby(gkey).agg(self.visit_agg)
        self.visit_redcap = self.ProcessVisits(visit_redcap)
        
        gkey = 'Patient ID'
        pat_redcap = self.visit_redcap.groupby(level=gkey).agg(self.pat_agg)
        self.pat_redcap = self.ProcessPatients(pat_redcap)
    
            
    def MakePatientGroups(self):
        
        self.SplitGroups(self.pat_redcap, self.pat_agg)
        

    def MakeVisitGroups(self):
        
        aligned_data, _ = self.visit_redcap.align(self.pat_redcap, 
                                                    level = 'Patient ID',
                                                    join = 'inner')
        
        self.SplitGroups(aligned_data, self.visit_agg)
    
    def AssignGroups(self, aligned_data):
        raise NotImplementedError
        
    
    def SplitGroups(self, aligned_data, agg_dict):
        
        cur_levels = aligned_data.index.names
        aligned_data = aligned_data.reset_index()
        aligned_data['Grouping'] = self.AssignGroups(aligned_data)
        self.all_group = aligned_data.groupby(['Grouping'] + cur_levels).agg(agg_dict)
        
        
    def make_demo_figures(self):
        
        self.AgeHist()
        
        for key, group in self.config_data.groupby('PlotName'):
            if ((group['PlotType'] == 'BarChart').all()) & ((group['DemographicFunction'] == 'ChoiceCount').all()):
                try:
                    self.plot_bar_chart(group['RawName'].values, key)
                except:
                    print('bad on ', key)
            elif ((group['PlotType'] == 'BoxPlot').all()) & ((group['DemographicFunction'] == 'MeanFunc').all()):
                self.make_box_plot(group['RawName'].values, key)
            elif ((group['PlotType'] == 'LogBoxPlot').all()) & ((group['DemographicFunction'] == 'MeanFunc').all()):
                self.make_log_box_plot(group['RawName'].values, key)
            elif ((group['PlotType'] == 'BarChart').all()) & ((group['DemographicFunction'] == 'IntegerCount').all()):
                print(key, group['RawName'].values)
                self.make_integer_bar(group['RawName'].values[0], key)
        
        
    def AgeHist(self):
        
        bins = [20,30,40,50,60,70,80]
        fig = plt.figure()

        g1data = Series(np.histogram(self.all_group.ix[True]['CalcAge'].values, bins = bins)[0], index = bins[:-1])
        g2data = Series(np.histogram(self.all_group.ix[False]['CalcAge'].values, bins = bins)[0], index = bins[:-1])

        df = DataFrame({self.groupnames[True]:g1data,
                        self.groupnames[False]:g2data})
        df.plot(kind = 'bar', grid = True)
        plt.xlabel('Age at Visit')
        plt.ylabel('#')
        
        return fig, self.all_group['CalcAge']
    
    
    def plot_bar_chart(self, items, title):


        g1sum = self.all_group.ix[True][items].mean()*100
        g2sum = self.all_group.ix[False][items].mean()*100
        allsum = self.all_group[items].mean()*100
        df = DataFrame({self.groupnames[True]:g1sum, 
                        self.groupnames[False]:g2sum,
                        'All':allsum})
        ncols = dict([(col, fix_col_name(col)) for col in df.index])
        df = df.rename(index=ncols)
        fig = plt.figure()
        df.plot(kind = 'bar', ax = plt.gca(), grid = True)
        plt.title(title)
        plt.ylabel('%')
        
        return fig, self.all_group[items]
    
    def make_box_plot(self, items, title):
        g1items = self.all_group.ix[True][items].reset_index()
        g2items = self.all_group.ix[False][items].reset_index()
        allitems = self.all_group[items].reset_index()
    
        pltdata = [(allitems, 'All'),
                    (g1items, self.groupnames[True]),
                    (g2items, self.groupnames[False])]
        odict = {}
        for item, (data, name) in product(items, pltdata):
            odict[item + '--' + name] = data[item]
        
        fig = plt.figure()
        df = DataFrame(odict)
        df.boxplot(rot = 90, ax = plt.gca())
        plt.title(title)
        plt.ylabel('Value')
        
        return fig, self.all_group[items]
    
    def make_log_box_plot(self, items, title):
        g1items = self.all_group.ix[True][items].reset_index()
        g2items = self.all_group.ix[False][items].reset_index()
        allitems = self.all_group[items].reset_index()
    
        pltdata = [(allitems, 'All'),
                    (g1items, self.groupnames[True]),
                    (g2items, self.groupnames[False])]
        odict = {}
        for item, (data, name) in product(items, pltdata):
            odict[item + '--' + name] = data[item]
        
        fig = plt.figure()
        df = DataFrame(odict)
        df.apply(np.log10).boxplot(rot = 90, ax = plt.gca())
        plt.title(title)
        plt.ylabel('log10(Value)')
        
        return fig, self.all_group[items]
    
    def make_integer_bar(self, col, title):
        
        if len(self.all_group[col].unique()) < 2:
            return None, None
        bins = np.arange(0, self.all_group[col].max()+1)
        g1data = Series(np.histogram(self.all_group.ix[True][col].values, bins = bins)[0], index = bins[:-1])/len(self.all_group.ix[True])
        g2data = Series(np.histogram(self.all_group.ix[False][col].values, bins = bins)[0], index = bins[:-1])/len(self.all_group.ix[False])
        alldata = Series(np.histogram(self.all_group[col].values, bins = bins)[0], index = bins[:-1])//len(self.all_group)

        ndf = DataFrame({'All': alldata*100,
                        self.groupnames[True]: g1data*100,
                        self.groupnames[False]: g2data*100})
        
        fig = plt.figure()
        
        ndf.plot(kind = 'bar', ax = plt.gca(), grid = True)
        plt.title(title)
        
        return fig, self.all_group[col]
        

# <codecell>

class GenderPatData(PatData):
    
    def AddGroupNames(self):
        
        self.groupnames = dict(zip([True, False], ['Male', 'Female']))
    
    def AssignGroups(self, visit_redcap):
        
        return visit_redcap['IdentifiesMale']

# <codecell>

config_file = 'Data/Config/ReportFile.csv'
demo_file = 'Data/RedcapDumps/HIVAIDSGeneticAnalys_DATA_LABELS_2012-12-11_1720.csv'
tmp = GenderPatData(demo_file, config_file)
tmp.CalcAll()

# <codecell>

tmp.ProcessRedcap()
tmp.MakePatientGroups()

# <codecell>

tmp.make_demo_figures()

# <codecell>

class NeuroPatData(PatData):
    
    def AddGroupNames(self):
        
        self.groupnames = dict(zip([True, False], ['No Neuro', 'With Neuro']))
    
    def AssignGroups(self, aligned_data):
        
        return aligned_data["Mental Health Issues (choice='No neurological problems')"]
    
ntmp = NeuroPatData(None, None)
ntmp.CopyFromOtherData(tmp)
ntmp.AddGroupNames()
ntmp.MakePatientGroups()
ntmp.make_demo_figures()

# <codecell>

from itertools import product
items = ['Neurocognitive test',
'MSK Score',
'Psychomotor Speed Score',
'Memory Recall Score',
'Constructional Score',
'Total Modified Hopkins Dementia Score',
]

col = items[-1]


# <codecell>

def safe_float(m):
    try:
        return float(m)
    except:
        return None
demo_data['Weight'].apply(safe_float).hist()

# <codecell>

tmp.redcap

# <codecell>


