# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

from __future__ import division
from pandas import DataFrame, Series, merge, read_csv, MultiIndex, Index, concat
from subprocess import check_call
from tempfile import NamedTemporaryFile as NTF
import os, os.path
import numpy as np
from scipy.stats import ttest_ind
from itertools import groupby,combinations, islice
from operator import itemgetter
from Bio import Phylo
import networkx
import sys

from random import shuffle
import csv, shlex, shutil

os.chdir('/home/will/HIVTropism//')
sys.path.append('/home/will/HIVReportGen/AnalysisCode/')
sys.path.append('/home/will/PySeqUtils/')

# <codecell>

from SeqProcessTools import read_pat_seq_data, load_training_seq_data, align_seq_data_frame
from TreeingTools import make_mrbayes_trees, run_bats, get_pairwise_distances, check_distance_pvals
from GeneralSeqTools import fasta_reader, WebPSSM_V3_fasta, yield_chunks
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import glob
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
from itertools import chain

# <codecell>

def simple_translate(inseq):
    seq = Seq(inseq, alphabet=generic_dna)
    return seq.translate().tostring()


seq_files = glob.glob('LANLdata/*.fasta')
seq_data = []
for f in seq_files:
    parts = f.split(os.sep)[-1].split('-')
    prot = parts[1]
    subtype = parts[0]
    with open(f) as handle:
        for name, seq in fasta_reader(handle):
            nseq = ''.join(l for l in seq if l.isalpha())
            if prot != 'LTR':
                nseq = simple_translate(nseq)
            seq_data.append((subtype, name, prot, nseq))
            
pat_ltr_seqs = glob.glob('/home/will/HIVReportGen/Data/PatientFasta/*LTR.fasta')
print len(pat_ltr_seqs)
for f in pat_ltr_seqs:
    with open(f) as handle:
        for name, seq in fasta_reader(handle):
            seq_data.append(('SubB', name, 'LTR', seq))
            
pat_v3_seqs = glob.glob('/home/will/HIVReportGen/Data/PatientFasta/*V3.fasta')
print len(pat_v3_seqs)
for f in pat_v3_seqs:
    with open(f) as handle:
        for name, seq in fasta_reader(handle):
            seq_data.append(('SubB', name, 'V3', seq))

            
seq_df = DataFrame(seq_data, columns=['Subtype', 'Name', 'Prot', 'Seq'])

# <codecell>

v3_seqs = [(name, seq) for subtype, name, prot, seq in seq_data if prot == 'V3']

have_data = []
need_data = set([name for name, _ in v3_seqs])
count = 0
while need_data and count < 5:
    count += 1
    print len(need_data)
    gen = ((name, seq) for name, seq in v3_seqs if name in need_data)
    chunks = yield_chunks(gen, 50)
    with ThreadPoolExecutor(max_workers = 10) as e:
        res = e.map(WebPSSM_V3_fasta, chunks)
        have_data += list(chain.from_iterable(res))
    need_data -= set(row[0] for row in have_data)
            

# <codecell>

import numpy as np
def safe_float(inval):
    try:
        return float(inval)
    except ValueError:
        return np.nan

fields = ['name','score','pred','x4.pct','r5.pct','geno','pos.chg','net.chg','percentile']
pssm_data = DataFrame(have_data, columns=fields)
pssm_data['pred'] = pssm_data['pred'] == '1'

float_cols = [1, 3, 4, 6, 7, 8]
for col in float_cols:
    pssm_data[fields[col]] = pssm_data[fields[col]].map(safe_float)
    
valid_pssm = pssm_data[pssm_data['percentile']<0.95]

# <codecell>

trop_scores = valid_pssm.groupby('name')[['score']].mean()
#trop_scores.to_excel('NewPSSMScores.xls')

# <codecell>

tmp = dict(zip(seq_df['Name'], seq_df['Subtype']))
grouped_seq_df = seq_df.pivot(index = 'Name', columns='Prot', values='Seq').reset_index()
grouped_seq_df = merge(grouped_seq_df, trop_scores, 
                        left_on = 'Name', right_index = True)
grouped_seq_df['Subtype'] = grouped_seq_df['Name'].map(lambda x: tmp[x])
print grouped_seq_df

# <codecell>

from collections import defaultdict
from random import shuffle
def calculate_mutual_info(signal1, signal2, shuf = False, batch=False, **kwargs):
    """Caluculates the Mutual Information shared by two signals.

    Arguements:
    signal1 -- An iterable indicating the first signal
    signal2 -- An iterable indicating the second signal

    Signals MUST be the same length! Items must be hashable!

    Returns:
    Mutual Information -- float"""


    if shuf:
        shuffle(signal1)
        shuffle(signal2)
    if batch:
        res = []
        extra = {}
        for _ in xrange(batch):
            r, extra = calculate_mutual_info(signal1, signal2, shuf=True, want_extra=True, **extra)
            res.append(r)
        return res

    overlap_prob = signal2prob(zip(signal1, signal2))
    signal1_prob = kwargs.get('S1prob',signal2prob(signal1))
    signal2_prob = kwargs.get('S2prob',signal2prob(signal2))

    num_items = len(signal1)
    mut_info = float()


    for (s1, s2), count in overlap_prob.items():
        mut_info += overlap_prob[(s1, s2)]*log(overlap_prob[(s1, s2)]/(signal1_prob[s1]*signal2_prob[s2]))

    if kwargs.get('want_extra', False):
        return mut_info, {'S1prob':signal1_prob, 'S2prob':signal2_prob}
    else:
        return mut_info
    
def signal2prob(signal):
    counts = signal2count(signal)
    num = len(signal)
    for key, val in counts.items():
        counts[key] = val/num
    return counts

def signal2count(signal):
    counts = defaultdict(int)
    for s in signal:
        counts[s] += 1
    return counts

# <codecell>

wanted_seq_data = grouped_seq_df.dropna(subset = ['score'])
align_data = align_seq_data_frame(wanted_seq_data,  '/home/will/HIVReportGen/Data/BlastDB/ConBseqs.txt')

# <codecell>

def decide_tropism(inval):
    if inval < -6.95:
        return 'R5'
    elif inval > -2.88:
        return 'X4'
    return np.nan
cols = ['gp120-seq-align',
        'gp41-seq-align',
        'Nef-seq-align',
        'LTR-seq-align']
align_data['Tropism'] = align_data['score'].map(decide_tropism)
#ask = align_data['Tropism'] == 'X4'
#mp_data = align_data[[col, 'V3-seq-align']][mask].dropna()
#print tmp_data
align_data

# <codecell>

t = align_data.dropna(subset = ['gp120', 'gp41', 'Nef', 'LTR'], how = 'all')
t['Subtype'].value_counts()

# <codecell>


tmp = pivot_table(align_data, rows = ['Subtype', 'Tropism'], values = ['gp120', 'gp41', 'LTR', 'Nef'], aggfunc = 'count')
seq_counts = tmp.drop(['Subtype', 'Tropism'], axis = 1)
seq_counts

#crosstab(align_data['Subtype'][v3mask], )
#align_data.dropna(subset = ['V3'])

# <codecell>

pivot_table?

# <codecell>

def process_tups(tups):
    ocols, tupA, tupB = tups
    batchsize = 100
    posA, colA = tupA
    posB, colB = tupB
    subtype, col, name = ocols
    
    
    tv = calculate_mutual_info(list(colA), list(colB))
    wrongs = []
    numwrongs = 0
    while (numwrongs < 10) and batchsize < 5000:
        wrongs += calculate_mutual_info(list(colA), list(colB), batch = batchsize)
        numwrongs = (np.array(wrongs) > tv).sum()
        batchsize *= 2
    wrongs = np.array(wrongs)
    pval = (wrongs > tv).mean()
    return subtype, col, name, posA, posB, tv, wrongs.mean(), len(wrongs), pval

# <codecell>

from itertools import product, imap
from itertools import dropwhile
from functools import partial

cols = ['gp120-seq-align',
        'gp41-seq-align',
        'Nef-seq-align',
        'LTR-seq-align']
masks = [('All', align_data['V3'].map(len)>0),
          ('R5', align_data['Tropism'] == 'R5'),
          ('X4', align_data['Tropism'] == 'X4')]
subtypes = [('SubBC', (align_data['Subtype'] == 'SubB') | (align_data['Subtype'] == 'SubC')),  
                ('SubC', align_data['Subtype'] == 'SubC'),
                ('SubB', align_data['Subtype'] == 'SubB'), 
            ]

def yield_cols():
    for (subname, submask), col, (trop_name, trop_mask) in product(subtypes, cols, masks):
    
        tmp_data = align_data[[col, 'V3-seq-align']][submask & trop_mask].dropna()
        if len(tmp_data.index) < 20:
            continue
            
        ovals = zip(*tmp_data[col].values)
        v3vals = zip(*tmp_data['V3-seq-align'].values)
        args = product(enumerate(ovals), enumerate(v3vals))
        for tupA, tupB in args:
            yield (subname, col, trop_name), tupA, tupB
    
def check_item(last_res, row):
    
    cols = ['gp120-seq-align',
        'gp41-seq-align',
        'Nef-seq-align',
        'LTR-seq-align']
    names = ['All', 'R5', 'X4']
    
    last_col_num = [num for num, val in enumerate(cols) if val == last_res[1]][0]
    last_mask_num = [num for num, val in enumerate(names) if val == last_res[0]][0]
    last_oval = int(last_res[2])
    last_v3val = int(last_res[3])
    
    (cur_col, cur_name), (posA, _), (posB, _) = row
    #print cur_col, cur_name
    cur_col_num = [num for num, val in enumerate(cols) if val == cur_col][0]
    cur_mask_num = [num for num, val in enumerate(names) if val == cur_name][0]
    
    checks = [(last_col_num, cur_col_num),
                (last_mask_num, cur_mask_num),
                (last_oval, posA),
                (last_v3val, posB)]
    #print checks
    
    for last, cur in checks:
        if last > cur:
            return True

    return False



# <codecell>

results = []#list(csv.reader(open('quick_linkage.csv')))[:-2]
writer = csv.writer(open('subtype_linkage.csv', 'w'))
process_iter = yield_cols()

with ProcessPoolExecutor(max_workers = 30) as e:
    
    res_iter = enumerate(e.map(process_tups, process_iter))
    for num, row in res_iter:
        if (num == 10) or (num == 100) or (num == 1000) or (num % 5000 == 0):
            print num, row
        
        writer.writerow(row)
        results.append(row)

# <codecell>

def fix_fun(row):
    trow = [row[0], row[1], int(row[2]), int(row[3]), 
            float(row[4]), float(row[5]), float(row[6]), float(row[7])]
    return trow

tdf = DataFrame(results, columns = ['Subtype', 'Prot', 'Group', 'TPos', 'Vpos', 'MI', 'nMI', 'Count', 'Pval'])
print tdf.head()

# <codecell>

len(align_data['LTR-seq-align'].dropna().values[0])

# <codecell>

tregions = [('gp120', 462),
            ('gp41', 313),
            ('Nef', 205),
            ('LTR', 630)]
nregions = []
mapper = {}
endpos = 0
for prot, size in tregions:
    mapper[prot] = endpos
    nregions.append((prot, endpos, endpos+size, 0, size, 1, 'product'))
    endpos += size + 50

nregions.append(('V3', 267, 302, 0, 32, 2, 'ROI'))
mapper['V3'] = 267
regions = DataFrame(nregions, columns = ['Region_name', 'Genome_Start', 'Genome_Stop', 'Gene_AA_Start', 'Gene_AA_Stop', 'Frame', 'RegionType'])
regions = regions.set_index('Region_name')
print regions, endpos

# <codecell>

import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.path import Path
from math import pi, sin, cos
from random import randint, shuffle, seed
import pylab
from itertools import cycle
import colorsys

class PyCircos(object):
    def __init__(self, regions, genome_size, figsize = (10,10), fig = None, ax = None):
        #a panda object keyed by the region name
        #it must have the following columns:
        # genome_start, genome_stop, frame
        # gene_name
        # aa_start, aa_stop
        # optional: color
        self.regions = regions 
        if ax:
            self.ax = ax
        if fig:
            self.figure = fig
        if fig is None and ax is None:
            self.figure = figure(figsize = figsize)
            self.ax = plt.gca()
            print 'made defaults'
            

        self.genome_size = genome_size
        self.radius = 10
        self.frame_offset = 0.5
        self._prep_circle()
    
    def _prep_circle(self):
        radius_size = self.radius
        circ=pylab.Circle((0,0),radius=radius_size, fc = 'w')
        lim = radius_size + 4.5*self.frame_offset
        self.ax.set_xlim(-lim,lim)
        self.ax.set_ylim(-lim,lim)
        self.ax.add_patch(circ)
    
    def _get_colors(self, num_colors):
        colors=[]
        seed(50) #so we get a consitent set of random colors!
        for i in np.arange(0., 360., 360. / num_colors):
            hue = i/360.
            lightness = (50 + np.random.rand() * 10)/100.
            saturation = (90 + np.random.rand() * 10)/100.
            colors.append(colorsys.hls_to_rgb(hue, lightness, saturation))
        return colors
    
    def add_labels(self, skip_labels = set(), with_text = True):

        radius_size = self.radius
        frame_offset = self.frame_offset
        #print self.genome_size
        roi_types = ['product', 'ROI']
        colors = self._get_colors(len(self.regions))
        shuffle(colors)
        colors = iter(colors)
        for ty in roi_types:
            tyregions = self.regions[self.regions['RegionType'] == ty]
            for region, row in tyregions.iterrows():
                if region == 'gp160':
                    continue
                if region in skip_labels:
                    _ = colors.next()
                    continue
                nrad = row['Frame']*frame_offset + radius_size
                theta1 = 360*(row['Genome_Start']/self.genome_size)
                theta2 = 360*(row['Genome_Stop']/self.genome_size)
                #print region, theta1, theta2, row
                arc_patch = patches.Arc((0,0), 2*nrad, 2*nrad, 
                                    theta1 = theta1, theta2 = theta2,
                                    edgecolor = colors.next(),
                                    linewidth = 10)
                self.ax.add_patch(arc_patch)
                midangle = (theta1 + abs(theta1-theta2)/2)*(2*pi/360)
                point = [(nrad+0.5)*cos(midangle), (nrad+0.5)*sin(midangle)]
                if with_text:
                    self.ax.text(point[0], point[1], region, fontsize = 18)
                
            
    
    def _make_path(self, source_item, target_item):
        rads_per_item = 2*pi/self.genome_size
        srad = source_item*rads_per_item
        trad = target_item*rads_per_item
        
        svert = [self.radius*cos(srad), self.radius*sin(srad)]
        tvert = [self.radius*cos(trad), self.radius*sin(trad)]
        
        verts = [svert, [0,0], tvert]
        codes = [Path.MOVETO,
                 Path.CURVE3,
                 Path.CURVE3]
        return verts, codes
    
    def add_raw_path(self, genome_start, genome_end, **kwargs):
        
        if 'facecolor' not in kwargs:
            kwargs['facecolor'] = 'none'
        if 'lw' not in kwargs:
            kwargs['lw'] = 2
            
        verts, codes = self._make_path(genome_start, genome_end)
        pt = Path(verts, codes)
        linepatch = patches.PathPatch(pt, **kwargs)
        self.ax.add_patch(linepatch)
    
    def add_raw_with_check(self, genome_start, region_start, genome_end, region_end, **kwargs):
        
        reg = self.regions
        startTF = genome_start >= reg['Genome_Start'][region_start] and genome_start <= reg['Genome_Stop'][region_start]
        endTF = genome_end >= reg['Genome_Start'][region_end] and genome_end <= reg['Genome_Stop'][region_end]
        if startTF and endTF:
            self.add_raw_path(genome_start, genome_end, **kwargs)
        
    
    def add_region_path(self, source_region_name, source_pos, target_region_name, target_pos):
        
        regions = self.regions
        sreg = regions.ix[source_region_name]
        treg = regions.ix[target_region_name]
        source_genome_start = sreg['Genome_Start']#+3*(source_pos-sreg['Gene_AA_Start'])
        target_genome_start = treg['Genome_Start']#+3*(target_pos-treg['Gene_AA_Start'])
        
        self.add_raw_path(source_genome_start, target_genome_start)

        




# <codecell>

#['Subtype', 'Prot', 'Group', 'TPos', 'Vpos', 'MI', 'nMI', 'Count', 'Pval']

subtypes = [ 'SubBC', 'SubB', 'SubC']
prots = ['gp120', 'gp41', 'LTR', 'Nef']
color_map = {'All':'k', 'R5':'r', 'X4':'b'}
cuts = {'All':10, 'R5': 10, 'X4': 5}
ncuts = tdf['Group'].map(lambda x: cuts[x])
groups = ['All', 'R5', 'X4', 'Combined']

all_plots = []
wanted_rows = (tdf['Pval'] == 0) #& (tdf['MI'] > ncuts*tdf['nMI'])
for with_text in [False, True]:
    fig, axs = plt.subplots(3,4, figsize = (15*(5/4),15))
    ax_iter = iter(axs.flatten())
    for subtype in subtypes:
        plotted_args = []
        for group in groups:
            ax = ax_iter.next()
            circ_obj = PyCircos(regions, 1810, ax = ax, fig = fig)
    
            if group != 'Combined':
                gmask = (tdf['Group'] == group) & wanted_rows
                gmask &= tdf['Subtype'] == subtype
                tdata = tdf[gmask].sort('deltaMI').dropna()
                skip_labels = set()
                for prot in prots:
                    pdata = tdata[tdata['Prot'] == prot+'-seq-align']
                    pmask = pdata['deltaMI'].rank() > (len(pdata)-20)
                    if pmask.sum() == 0:
                        skip_labels.add(prot)
                    for _, row in pdata[pmask].iterrows():
                        Cpos = row['TPos'] + mapper[prot]
                        Tpos = row['Vpos'] + mapper['V3']
                        if (prot == 'gp120') and (Cpos > 265) and (Cpos < 305):
                            continue
                        edgecolor = color_map[row['Group']]
                        circ_obj.add_raw_path(Tpos, Cpos, edgecolor = edgecolor, alpha = 0.2)
                        plotted_args.append((Tpos, Cpos, edgecolor, group, subtype, prot))
                        all_plots.append((row['Vpos'], row['TPos'], edgecolor, group, subtype, prot))
            else:
                skip_labels = set()
                for Tpos, Cpos, edgecolor, group, subtype, prot in plotted_args:
                    circ_obj.add_raw_path(Tpos, Cpos, edgecolor = edgecolor, alpha = 0.2)
                #all_plots += plotted_args
                
            circ_obj.add_labels(skip_labels=skip_labels, with_text = with_text)
            ax.set_title(group + ' ' + subtype)
            ax.axis('off')
    fname = 'new_tropism_MI_links'
    if with_text:
        fname += '_with_text'
    plt.savefig(fname + '.png', dpi = 200)
    
plt.close()

# <codecell>

rolling_count?

# <codecell>

from pandas import rolling_sum
def get_group_count(indf):
    tregions = dict([('gp120', 462),
            ('gp41', 313),
            ('Nef', 205),
            ('LTR', 630)])
    nser = Series([0]*tregions[indf['Prot'].values[0]])
    #print nser
    nser[indf['Cpos']] = 1
    res = rolling_sum(nser, 20, min_periods=1)
    out = DataFrame({'RollingCount':res, 'Pos':range(tregions[indf['Prot'].values[0]])})
    out['Prot'] = indf['Prot'].values[0]
    out['Group'] = indf['Group'].values[0]
    out['Subtype'] = indf['Subtype'].values[0]
    return out
    


plotted_df = DataFrame(all_plots, columns = ['V3pos', 'Cpos', 'edgecolor', 'Group', 'Subtype', 'Prot']).sort(['Group', 'Subtype', 'Prot', 'Cpos'])
res = plotted_df.groupby(['Group', 'Subtype', 'Prot'], as_index = False).apply(get_group_count)

out = pivot_table(res, rows = ['Prot', 'Pos'], cols = ['Subtype', 'Group'], values = 'RollingCount').dropna(how = 'all')
wanted = (out>5).any(axis = 1)
out[wanted].to_csv('MIgroupings.csv')
#for key, group in :
#    print key, rolling_count(group['Cpos'], 5)

# <codecell>

tdf['ProtName'] = tdf['Prot'].map(lambda x: x.split('-')[0])
res = pivot_table(tdf[wanted_rows], rows = ['ProtName', 'TPos'], cols = 'Group', values = 'Pval', aggfunc = 'min')
(res == 0).to_excel('mutual_info.xlsx')

# <codecell>

from scipy.stats import fisher_exact
cols = ['gp120-bin-align',
        'gp41-bin-align',
        'Nef-bin-align',
        'LTR-bin-align']
masks = [('R5', align_data['Tropism'] == 'R5'),
          ('X4', align_data['Tropism'] == 'X4')]
subtypes = [('SubBC', (align_data['Subtype'] == 'SubB') | (align_data['Subtype'] == 'SubC')),  
                ('SubC', align_data['Subtype'] == 'SubC'),
                ('SubB', align_data['Subtype'] == 'SubB'), 
            ]

def yield_cols_for_fishers():
    for (subname, submask), col in product(subtypes, cols):
    
        tmp_data = align_data[[col, 'Tropism']][submask].dropna()
        if len(tmp_data.index) < 20:
            continue
        r5_trops = tmp_data['Tropism'] == 'R5'
        ovals = zip(*tmp_data[col].values)
        for pos, tup in enumerate(ovals):
            yield (subname, col.split('-')[0], pos), Series(tup, index = tmp_data.index), r5_trops
            
def process_fishers(intup):
    
    subtype, col, pos = intup[0]
    muts = intup[1]
    trops = intup[2]
    #print muts
    #print trops
    
    table = [[(muts & trops).sum(), (muts & ~trops).sum()],
             [(~muts & trops).sum(), (~muts & ~trops).sum()]]
    odds, pval = fisher_exact(table)
    return subtype, col, pos, odds, pval, trops.sum(), (~trops).sum()

with ProcessPoolExecutor(max_workers = 30) as e:
    fisher_res = list(e.map(process_fishers, yield_cols_for_fishers()))
fisher_df = DataFrame(fisher_res, columns = ['Subtype', 'Prot', 'Pos', 'Odds', 'Pval', 'R5Count', 'X4Count'])
    

# <codecell>

fisher_df = fisher_df.dropna()

# <codecell>

subtypes = ['SubB', 'SubC', 'SubBC', 'Combined']
prots = ['gp120', 'gp41', 'LTR', 'Nef']
color_map = {'SubB':'g', 'SubC':'r', 'SubBC':'b'}
fig, axs = plt.subplots(2,2, figsize = (10,10))
wanted_rows = (fisher_df['Pval'] <= 0.01)
for with_text in [False, True]:
    for subtype, ax in zip(subtypes, axs.flatten()):
        circ_obj = PyCircos(regions, 1810, ax = ax, fig = fig)
        if subtype != 'Combined':
            gmask = (fisher_df['Subtype'] == subtype) & wanted_rows
            skipm = (fisher_df['Subtype'] == subtype) & (fisher_df['X4Count'].fillna(0) > 0)
            skips = set(prots) - set(fisher_df['Prot'][skipm])
        else:
            gmask = wanted_rows.copy()
            skipm = fisher_df['X4Count'].fillna(0) >  0
            skips = set(prots) - set(fisher_df['Prot'][skipm])
        tdata = fisher_df[gmask].dropna()
    
        for _, row in tdata.iterrows():
            prot = row['Prot']
            Cpos = row['Pos']
            if (prot == 'gp120') and (Cpos < 302) and (Cpos > 267):
                continue
            Tpos = 20
            circ_obj.add_raw_path(Tpos+mapper['V3'],Cpos + mapper[prot], edgecolor = color_map[row['Subtype']], alpha = 0.2, lw = 3)
    
        circ_obj.add_labels(skip_labels = skips, with_text=with_text)    
        ax.set_title(subtype)
        ax.axis('off')
    fname = 'new_tropism_fisher_links'
    if with_text:
        fname += '_with_text'
    plt.savefig(fname + '.png', dpi = 200)
    
#plt.savefig('fishers_test_subtype.png')

# <codecell>

from pandas import pivot_table
res = pivot_table(fisher_df, rows = ['Prot', 'Pos'], cols = 'Subtype', values = ['Pval', 'R5Count', 'X4Count'], aggfunc = np.mean)
print res.head(n = 30).to_string()

# <codecell>

res['R5Count']['SubC'] = (res['R5Count']['SubBC']-res['R5Count']['SubB']).combine_first(res['R5Count']['SubC'])
res['X4Count']['SubC'] = (res['X4Count']['SubBC']-res['X4Count']['SubB']).combine_first(res['X4Count']['SubC'])

# <codecell>

mask = (res['Pval'] < 0.01).any(axis = 1)
res[mask].to_csv('fishers_table.csv')

# <codecell>


