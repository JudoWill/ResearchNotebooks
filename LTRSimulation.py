# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

from __future__ import division
from pandas import *
import os, os.path
import numpy as np
from matplotlib import pyplot as plt

os.chdir('/home/will/LTRtfAnalysis/')

# <codecell>

import glob

files = glob.glob('microarray_data/rawdata/*.tsv')
data = DataFrame()
for f in files:
    data = data.combine_first(read_csv(f, sep = '\t', index_col=0))
microdata = data

# <codecell>

wanted = {'1050_at':'cebp',
'4790_at':'nfkb', #also called Rela
'3725_at':'AP-1', #also called FOS
'6667_at':'sp'
}
gene_data = microdata.ix[wanted.keys()].rename(wanted)

# <codecell>



annot_data = read_csv('microarray_data/ExpressionAnnotations.tsv', sep = '\t')
#nnot_data['Sample Name']
annot_data['ColName'] = annot_data['Array Data File'].map(lambda x: x.split('.')[0])
annot_data['CellGroup'] = annot_data['OrganismPart'].combine_first(annot_data['CellType'])
mdata = merge(gene_data.T, annot_data,
                left_index = True, right_on = 'ColName')
#print mdata
#sorted([(key, len(rows)) for key, rows in annot_data.groupby(['CellGroup', 'DiseaseState']).groups.items()], key = lambda x: -x[1])
expression_data = mdata.groupby(['CellGroup', 'DiseaseState']).agg({'AP-1':['mean', 'std'], 'cebp':['mean', 'std'], 'sp':['mean', 'std'], 'nfkb':['mean', 'std']})
wanted_express = expression_data.ix[['B cell', 'CD4+ T cell', 'CD8+ T cell', 'PBMC', 'T cell', 'bone marrow', 'lymph node', 'monocyte', 'bone marrow', 'liver', 'brain']]
print wanted_express.to_string()

# <codecell>

from copy import deepcopy


cebp = u"""A  3   1  4  2  4  2 18 18  0
           C  0   1  1  9  2 15  0  0  6
           G  0   4  6  2 10  0  0  0  2
           T  15 12  7  5  2  1  0  0 10"""


#cebpmot = ModelMotif.load_from_japasr_str(cebp)

sp1 = u"""A  0  0  0  4   2  0  1  0  6  3
          C  32 30 35 27  5 28 31 24 25 26
          G  1  1  0  0  15  1  0  3  0  3
          T  2  4  0  4  13  6  3  8  4  3"""

#sp1mot = ModelMotif.load_from_japasr_str(cebp)
#sp1mot_r = sp1mot.reverse_complement()


nfkb = """A  [ 0  0  0  2 11  5  0  0  0  0  1 ]
          C  [ 0  0  0  0  1  0  5 13 17 18 15 ]
          G  [18 18 18 16  6  2  2  0  0  0  1 ]
          T  [ 0  0  0  0  0 11 11  5  1  0  1 ]"""
#nfkbmot = ModelMotif.load_from_japasr_str(nfkb)


let_probs = {'A':1/4, 'C':1/4, 'G':1/4, 'T':1/4}



# <codecell>

import re
import csv
from collections import defaultdict
from itertools import groupby
from decimal import Decimal

class ModelMotif(object):
    
    def __init__(self):
        self.seq_counts = []
        self.num_obs = None
        self.pwm = []
        self.name = None
    
    def __hash__(self):
        return hash(self.name)
    
    @staticmethod
    def load_from_japasr_str(name, input_str):
        
        num_re = re.compile('\d+')
        let_re = re.compile('\w')
        num_dict = {}
        for line in input_str.split('\n'):
            let = let_re.findall(line)[0][0]
            nums = [int(n) for n in num_re.findall(line)]
            num_dict[let] = nums
        mot = ModelMotif()
        mot.name = name
        
        for let_nums in zip(num_dict['A'], num_dict['C'], num_dict['G'], num_dict['T']):
            mot.seq_counts.append(dict(zip('ACGT', let_nums)))
        mot.num_obs = sum(let_nums)
        
        mot._CalculatePWM()
        return mot
    
    @staticmethod
    def load_from_fasta(fname, mot_name):
        
        seqs = []
        with open(fname) as handle:
            for key, lines in groupby(handle, key = lambda x: x.startswith('>')):
                if key:
                    name = lines.next()[1:].strip()
                else:
                    seq = ''.join(line.strip() for line in lines)
                    seq = seq.replace('-', '')
                    if name == 'K03455':
                        refseq = seq
                    else:
                        seqs.append(seq)
        reflen = len(refseq)
        good_seqs = [seq for seq in seqs if len(seq) == reflen]
        
        mot = ModelMotif()
        mot.name = mot_name
        
        for seqtup in zip(*good_seqs):
            
            mot.seq_counts.append({
                        'A':sum(1 for l in seqtup if l == 'A'),
                        'C':sum(1 for l in seqtup if l == 'C'),
                        'G':sum(1 for l in seqtup if l == 'G'),
                        'T':sum(1 for l in seqtup if l == 'T')})
        mot.num_obs = len(good_seqs)
        mot._CalculatePWM()
        return mot, refseq
    
    
    def _CalculatePWM(self, T = 300):
        
        let_probs = {'A':Decimal('0.25'), 'C':Decimal('0.25'), 'G':Decimal('0.25'), 'T':Decimal('0.25')}
        
        R = Decimal('8.314') #joules per Kelvin
        self.pwm = []
        fallback = Decimal(1/(2*self.num_obs))
        for col in self.seq_counts:
            tmp = {}
            for let in col.keys():
                pf = Decimal(col[let])/self.num_obs
                tmp[let] = R*T*(max(pf/let_probs[let], fallback)).ln()
            self.pwm.append(deepcopy(tmp))
            
    def reverse_complement(self):
        
        rev_mapping = dict(zip('ATCG', 'TAGC'))
        
        revc = ModelMotif()
        revc.num_obs = self.num_obs
        
        revc.seq_counts = []
        for row in self.seq_counts:
            nrow = {}
            for o, c in rev_mapping.items():
                nrow[c] = row[o]
            
            revc.seq_counts.append(deepcopy(nrow))
        revc._CalculatePWM()
        revc.name = self.name + '-R'
        return revc
            
        
    def CalculateKd(self, seq):
        R = Decimal('8.314') #joules per Kelvin
        T = Decimal(300)
        deltaG = Decimal(0)
        for l, col in zip(seq, self.pwm):
            deltaG += col[l]
        #print type(deltaG), type(R), type(T)
        #print (deltaG/R/T).exp()
        return (-deltaG/R/T).exp()
    
    def CalculateKa(self, seq):
        kd = self.CalculateKd(seq)
        #print self.name, seq, kd, 1/kd
        return 1/kd
        
    
    def CalculatePi(self, seq, conc):
        
        ka = self.CalculateKa(seq)
        prob = conc/(ka+conc)
        return prob
        
    
    
    
        
class PromoterModel(object):

    def __init__(self, sites, binding_links):
        
        self.sites = sites
        self.binding_links = binding_links
        self.probs = None
        
    def _calc_site_occup(self, site, site_seqs, prot_concs):
        #print site.name, site_seqs
        
        conc = prot_concs[site.name.split('-')[0]]
        ka = site.CalculateKa(site_seqs[site.name])
        Ci = 0
        coop_sites = []
        for osite_name, effect in self.binding_links[site.name]:
        
            if effect == 1:
                dimer_conc = prot_concs[osite_name.split('-')[0]]*conc
                dimer_ka = self.sites[osite_name].CalculateKa(site_seqs[osite_name])*ka
                coop_sites.append((dimer_conc, dimer_ka))
            else:
                comp_cons = prot_concs[osite_name.split('-')[0]]
                comp_ka = self.sites[osite_name].CalculateKa(site_seqs[osite_name])
                Ci += comp_cons*comp_ka
        #print site.name, ka, conc, coop_sites, Ci
        pi = (1+Ci)/(1+Ci+ka*conc)
        for conc, ka in coop_sites:
            pi *= (1+Ci)/(1+Ci+ka*conc)
        
        return 1-pi
    
    def _calc_all_occupancy(self, site_seqs, prot_concs):
        
        probs = {}
        for site in self.sites.values():
            probs[site.name] = self._calc_site_occup(site, site_seqs, prot_concs)
        return probs
    
    
    def CalcProbs(self, site_seqs, prot_concs, valid_occup, niters = 20):
        
        self.probs = self._calc_all_occupancy(site_seqs, prot_concs)
        #print sorted(self.probs.items())
            
        prob = 0
        for occup in valid_occup:
            tp = 1
            for site, effect in occup.items():
                if effect == 1:
                    tp *= self.probs[site]
                elif effect == -1:
                    tp *= 1-self.probs[site]
            prob += tp
        return prob
            
            
        #print self.probs
    
    @staticmethod
    def load_from_direc(path):
        
        refseqs = {}
        sites = {}
        for f in glob.glob(path + '/*.fasta'):
            site_name = f.rsplit('/',1)[-1].rsplit('.',1)[0]
            mot, refseq = ModelMotif.load_from_fasta(f, site_name)
            refseqs[site_name] = refseq
            sites[site_name] = mot
        #should be site. osite, coop_nature tsv file
        coop_sites = defaultdict(set)
        with open(path + '/coop_binding.tsv') as handle:
            for row in csv.reader(handle, delimiter = '\t'):
                coop_sites[row[0]].add((row[1], int(row[2])))
                coop_sites[row[1]].add((row[0], int(row[2])))
                
        return PromoterModel(sites, coop_sites), refseqs
            
    
    
    
    
        

# <codecell>

model, refseqs = PromoterModel.load_from_direc('motifs/')

rel = 1e-7
#macro_diff_conc = {
#'sp-III':1*rel,
#'sp-II':1*rel,
#'sp-I':1*rel,
#'nfkb-II':2**(8.436198-6.154713)*rel,
#'nfkb-I':2**(8.436198-6.154713)*rel,
#'cebp-I':2**(7.055735-6.154713)*rel,
#}


#tcell_diff_conc = {
#'sp-III':1*rel,
#'sp-II':1*rel,
#'sp-I':1*rel,
#'nfkb-II':2**(9.147623-6.240220)*rel,
#'nfkb-I':2**(9.147623-6.240220)*rel,
#'cebp-I':2**(6.546898-6.240220)*rel,
#}

tmp_dict = defaultdict(lambda : Decimal(rel))


valid_occup = [{
'nfkb-I':1,
'nfkb-II':1,
'sp-I':1,
'sp-II':1,
'sp-III':1
}]
#for nf in ['nfkb-I', 'nfkb-II']:
#    for sp in ['sp-I', 'sp-II', 'sp-III']:
#        valid_occup.append({nf:1, sp:1})



print model.CalcProbs(refseqs, tmp_dict, valid_occup)


# <codecell>

from scipy.stats import norm

obj = norm(loc=5,scale = 2)
obj.rvs()

# <codecell>

from itertools import product, islice, imap, cycle
from concurrent.futures import ProcessPoolExecutor
from scipy.stats import norm

rel = Decimal(3e-9)
data = []
prots = ['cebp', 'nfkb', 'sp']


def concentration_generator(gene_exp, relative_conc, num_items=50000):
    
    pdf_dict = {}
    for key, row in gene_exp.dropna().iterrows():
        pdf_dict[key] = dict([(p, norm(loc=row[p]['mean'], scale = row[p]['std'])) for p in prots])
       
    for key, pdf_dict in islice(cycle(pdf_dict.items()), num_items):
        odict = {}
        for p in prots:
            odict[p] = relative_conc*(Decimal(pdf_dict[p].rvs())**2)
        
        yield key, odict


def linker_fun(tup):
    cell_type, conc_dict = tup
    prob = model.CalcProbs(refseqs, conc_dict, valid_occup)
    
    return [float(conc_dict[p]) for p in prots]+[float(prob), cell_type]
    

with ProcessPoolExecutor(max_workers = 20) as executor:
    conc_items = concentration_generator(wanted_express, rel, num_items=500000)
    data = list(executor.map(linker_fun, conc_items))
    
        
        

# <codecell>

df = DataFrame(data, columns = prots+['Prob', 'CellType'])

# <codecell>

len(df['CellType'].unique())

# <codecell>

from pylab import get_cmap

fig, axes = plt.subplots(5,4, sharex=True, sharey=True, figsize = (20,20))
bins = np.linspace(0,1,100)
for celltype, ax in zip(df['CellType'].unique(), axes.flatten()):
    mask = df['CellType'] == celltype
    ndata = df[mask]['Prob'].values
    H, xedges = np.histogram(ndata, bins = bins, normed = True)
    plt.sca(ax)
    plt.plot(bins[1:], H)
    plt.title(celltype)
    


# <codecell>

mask.sum()

# <codecell>

sp1_conc = 6.143E-1 #ng/mL in PBMCS from https://www.proteinspire.org/MOPED/mopedviews/proteinExpressionDatabase.jsf
sp1_mw = 1.339E-19 #grams (80693 daltons from http://www.uniprot.org/uniprot/P08047)
conc = sp1_conc*(1/1E9)*(1/sp1_mw)*(1/6.02E23)*(1000/1) #moles/L
print conc

# <codecell>

def CalculateKd(pwm, seq, T = 300):
    R = 8.314 #joules per Kelvin
    deltaG = 0
    for l, col in zip(seq, pwm):
        deltaG += col[l]
    return np.exp(-deltaG/R/T)

def CalculatePi(pwm, seq, conc, T = 300):
    kd = CalculateKd(pwm, seq)
    prob = conc/(kd+conc)
    return prob

def CalculatePocc(pwm_dict, coop_binding, base_conc, diff_conc):
    
    ipocc = 1
    for name, seq, coop_factors in make_sites(site_seqs, coop_binding, pwm_dict):
        tconc = base_conc*diff_conc[name]
        #print tconc
        pi = CalculatePi(pwm_dict[name], seq, tconc)
        for cname, cseq in coop_factors:
            oconc = base_conc*diff_conc[cname]
            pi *=  CalculatePi(pwm_dict[cname], cseq, tconc*oconc)
        ipocc *= 1-pi
        
    return 1-ipocc


def make_sites(site_seqs, coop_binding, pwm_dict):
    
    for key, coops in coop_binding.items():
        tsite = (site_seqs[key], pwm_dict[key.split('-')[0]])
        coops = [(site_seqs[site], pwm_dict[site.split('-')[0]]) for site in coops]
        yield tsite, coops
    
    

# <codecell>

site_seqs = {
'cebp-I':'AGCTTTCTACAA',
'nfkb-II':'GGGACTTTCC',
'nfkb-I':'GGGACTTTCC',
'sp-III':'GAGGCGTGGC'
}

sites = [('cebp', site_seqs['cebp-I'], [('nfkb', site_seqs['nfkb-II'])]),
         ('nfkb', site_seqs['nfkb-II'], [('cebp', site_seqs['cebp-I']), ('nfkb', site_seqs['nfkb-I'])]),
         ('nfkb', site_seqs['nfkb-I'], [('nfkb', site_seqs['nfkb-II'])]),
         ('sp1', site_seqs['sp-III'], [])]

cooper_map = {
'cebp-I':set(['nfkb-II']),
'nfkb-II':set(['nfkb-I', 'cebp-I']),
'nfkb-I':set(['nfkb-II']),
'sp-III':set()
}

pwm_dict = {
'cebp':PFMtoPWM(cebp, let_probs, num_items = 18),
'nfkb':PFMtoPWM(nfkb, let_probs, num_items = 18),
'sp':PFMtoPWM(sp1, let_probs, num_items = 18),
}

macro_diff_conc = {
'sp1':1,
'nfkb':2**(8.436198-6.154713),
'cebp':2**(7.055735-6.154713),
}


tcell_diff_conc = {
'sp1':1,
'nfkb':2**(9.147623-6.240220),
'cebp':2**(6.546898-6.240220),
}

base_conc = 0.01

print 'MACRO', CalculatePocc(sites, pwm_dict, base_conc, macro_diff_conc), '\nTCELL', CalculatePocc(sites, pwm_dict, base_conc, tcell_diff_conc)

# <codecell>


# <codecell>


