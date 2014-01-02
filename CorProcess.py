# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import os, os.path
import csv
import numpy as np
import concurrent.futures
from itertools import groupby

# <codecell>

def fasta_reader(handle, check_lens = True):
    
    outseqs = {}
    clen = None
    for key, lines in groupby(handle, lambda x: x.startswith('>')):
        if key:
            name = list(lines)[0][1:].strip()
        else:
            outseqs[name] = ''.join(l.strip() for l in lines)
            if clen is None:
                clen = len(outseqs[name])
            elif check_lens and (len(outseqs[name]) != clen):
                raise AssertionError('Sequence lengths are not the same')
    return outseqs

def data_reader(handle):
    
    data = dict(csv.reader(handle, delimiter = '\t'))
    #print(data)
    found_items = set(data.values())
    if len(found_items) == 1:
        raise AssertionError('Only one class found in the data!')
    elif len(found_items) > 2:
        raise AssertionError('More then 2 classes found in the data!')
    return data

# <codecell>

def data2numpy(align_dict, group_dict):
    
    groups = sorted(set(group_dict.values()))
    gdict = dict([(g, v) for g, v in zip(groups, [True, False])])
    common_keys = sorted(set(align_dict.keys()) & set(group_dict.keys()))
    
    align = np.array([list(align_dict[key]) for key in common_keys])
    mask = np.array([gdict[group_dict[key]] for key in common_keys])
    
    return align, mask, gdict


# <codecell>

from tempfile import NamedTemporaryFile as NTF
import shlex
from subprocess import check_call


def refine_alignment(npalign, refseq):
    cmd = 'muscle -in %(ifile)s -out %(ofile)s -refine'
    with NTF(mode = 'w') as inseq_handle:
        for num in range(npalign.shape[0]):
            seq = ''.join(npalign[num,:])
            name = 'Seq-%i' % num
            inseq_handle.write('>%s\n%s\n' % (name, seq))
        rseq = ''.join(refseq)
        inseq_handle.write('>%s\n%s\n' % ('REFSEQ', rseq))
        inseq_handle.flush()
        with NTF(mode = 'r') as outseq_handle:
            cmd_list = shlex.split(cmd % {'ifile':inseq_handle.name, 'ofile':outseq_handle.name})
            check_call(cmd_list)
            refined_seqs = fasta_reader(outseq_handle, check_lens = True)
            nrefseq = np.array(list(refined_seqs['REFSEQ']))
            nalign = np.array([list(refined_seqs['Seq-%i' % i]) for i in range(npalign.shape[0])])
    return nalign, nrefseq
            
    

# <codecell>

with open('/home/will/data.tsv') as handle:
    data = data_reader(handle)
    
with open('/home/will/Dropbox/HIVseqs/Neuroseqs/new_large_aln.fasta') as handle:
    seqs = fasta_reader(handle)
refseq = np.array(list(seqs['K03455']))

# <codecell>

npalign, mask, group_dict = data2numpy(seqs, data)

# <codecell>

from pandas import DataFrame, Series
from scipy.stats import fisher_exact, chi2

def fishers_test(g1align, g2align, ref):
    minval = 5
    g1mask = g1align==ref
    g2mask = g2align==ref 
    
    g1sum = g1mask.sum(axis = 0)
    g2sum = g2mask.sum(axis = 0)
    
    g1pos = (g1mask & (g1align != '-')).sum(axis = 0)
    g1neg = (~g1mask & (g1align != '-')).sum(axis = 0)
    g2pos = (g2mask & (g2align != '-')).sum(axis = 0)
    g2neg = (~g2mask & (g2align != '-')).sum(axis = 0)
    
    pvals = []
    for col in range(g1mask.shape[1]):
        if (g1sum[col]<minval) | (g2sum[col]<minval):
            pvals.append(None)
        else:
            _, pval = fisher_exact([[g1pos[col], g2pos[col]], [g1neg[col], g2neg[col]]])
            pvals.append(pval)
    return Series(pvals)

# <codecell>

fishers_res = fishers_test(npalign[mask,:], npalign[~mask,:], refseq)
fishers_res[fishers_res<0.05]

# <codecell>

from collections import defaultdict
def log_factorial(n):
    return sum(np.log10(x) for x in range(1,n+1))

def multi_nomial_dist(observed_count, total_count = None):

    #print observed_count
   
    if total_count is None:
        total_count = dict(observed_count.items())
    
    for key in observed_count:
        total_count[key] = max(observed_count[key], total_count[key])
        
    tp = count2prob(dict(total_count.items()), want_dec = False)
    N = int(sum(list(observed_count.values())))
    nf_log = log_factorial(N)

    d_log = 0
    for n in observed_count.values():
        d_log += log_factorial(n)
        
    p = nf_log-d_log
    for k, nnp in tp.items():
        p += observed_count[k]*np.log10(nnp)#.log10()
        
    return 10**float(p)

def countdict(intup):
    r = defaultdict(int)
    for n in intup:
        if n.isalpha() or len(n)>1:
            r[n] += 1
    return r

def count2prob(d, want_dec = False):
    n = sum(list(d.values()))
    for key in d.keys():
        d[key] = d[key]/n
    return d
 
    
def likelihood_ratio(g1align, g2align):
    
    g1count = countdict(g1align)
    g2count = countdict(g2align)
    if (sum(list(g1count.values())) < 5) | (sum(list(g2count.values())) < 5):
        return None, None
    
    self_p = multi_nomial_dist(g1count)
    g2_p = multi_nomial_dist(g1count, total_count = g2count)
    
    ratio = -2*(np.log(g2_p)-np.log(self_p))
    df = len(g1count)
    pval = 1-chi2.cdf(ratio, df)
    #print self_p, r5_p, ratio, df, pval
    return ratio, pval

def MVhypergeo_test(g1align, g2align, ref):
    
    pvals = []
    #print(g1align.shape, g2align.shape)
    for col in range(g1align.shape[1]):
        _, pval = likelihood_ratio(g1align[:,col], g2align[:,col])
        pvals.append(pval)
    return Series(pvals)
    

# <codecell>

mvhypergeo_res = MVhypergeo_test(npalign[mask,:], npalign[~mask,:], refseq)

# <codecell>

from tempfile import NamedTemporaryFile as NTF
from Bio import Phylo
from itertools import combinations, product
from subprocess import check_output, CalledProcessError
from operator import itemgetter
import shlex
import networkx
import time
from scipy.stats import ttest_ind

def fasta_write(handle, npalign):
    
    seqnames = []
    for row in range(npalign.shape[0]):
        seq = ''.join(npalign[row,:])
        name = 'Seq-%i' % row
        seqnames.append(name)
        ostr = '>%s\n%s\n' % (name, seq)
        handle.write(ostr)
    return seqnames
        
def tree_checker(row):
    tree_file, leaf1 = row
    dmat = {}
    tree = Phylo.read(open(tree_file), 'newick')
    leafs = sorted(tree.get_terminals(), key = lambda x: x.name)
    spos = max(pos for pos, leaf in enumerate(leafs) if leaf.name == leaf1.name)
    nleaf1 = next(tree.find_clades(name = leaf1.name))
    #print(nleaf1, spos, len(leafs[spos:]))
    for leaf2 in leafs[spos:]:
        try:
            d = tree.distance(nleaf1, leaf2)
            dmat[(leaf1.name, leaf2.name)] = d
            dmat[(leaf2.name, leaf1.name)] = d
        except RuntimeError:
            pass
    return dmat
        
def get_pairwise_distances(npalign, tree_file = None, seq_file = None):
    
    if seq_file is None:
        fasta_handle = NTF(mode = 'w')
    else:
        fasta_handle = open('/tmp/tmp.fasta', 'w')
    if tree_file is None:
        tree_handle = NTF()
    else:
        tree_handle = open(tree_file, 'w')
    seq_names = fasta_write(fasta_handle, npalign)
    
    fasta_handle.flush()
    os.fsync(fasta_handle.fileno())
    cmd = 'muscle -in %(ifile)s -tree2 %(treefile)s -gapopen -2.9'
    cmdlist = shlex.split(cmd % {
                                 'ifile':fasta_handle.name, 
                                 'treefile':tree_handle.name
                                 })
   
    try:
        t = check_output(cmdlist)
        tree = Phylo.read(open(tree_handle.name), 'newick')
    except CalledProcessError:
        #print('Could not make tree')
        return None
    except ValueError:
        #print('no tree present')
        return None
    except RuntimeError:
        return None
        
    
    seq_names = sorted(tree.get_terminals(), key = lambda x:x.name)
    net = Phylo.to_networkx(tree)
    dmat = networkx.all_pairs_shortest_path(net)
    terminals = tree.get_terminals()
    dists = np.zeros((npalign.shape[0], npalign.shape[0],))
    for t1, t2 in product(terminals, terminals):
        path = dmat[t1][t2]
        dist = sum(c.branch_length for c in path)
        i1 = int(t1.name.split('-')[1])
        i2 = int(t2.name.split('-')[1])
        dists[i1,i2] = dist
    
    
    return dists

def tree_dist_pvals(align, mask, window_size = 25):
    
    pvals = []
    span = int((window_size-1)/2)
    with concurrent.futures.ProcessPoolExecutor(max_workers = 30) as executor:
        aligns = []
        for col in range(align.shape[1]):
            spos = max(0, col-span)
            epos = min(align.shape[1], col+span)
            aligns.append(align[:,spos:epos])
        for col, dmat in enumerate(executor.map(get_pairwise_distances, aligns)):
            if (col == 5) | (col == 30) | (col % 50 == 0):
                print(col)
            if dmat is None:
            #print('had to skip column ', col)
                pvals.append(None)
                continue
            g1vals = dmat[mask, mask].ravel()
            g2vals = dmat[~mask, ~mask].ravel()
        
            _, pval = ttest_ind(g1vals, g2vals)
            pvals.append(pval)
    return Series(pvals)

# <codecell>

treeD_res = tree_dist_pvals(npalign, mask, window_size = 25)

# <codecell>

from pandas import read_csv
patdata = read_csv('/home/will/Dropbox/HIVseqs/Neuroseqs/NeuroData.tsv', sep = '\t')

# <codecell>

agg_dict = {
            'Psychomotor Speed Score':'max',
            'Memory Recall Score':'max',
            'Constructional Score':'max',
            'Total Modified Hopkins Dementia Score':'max'
            }
cog_data = patdata.groupby(['Patient ID', 'Patient visit number'], as_index = False).aggregate(agg_dict)

# <codecell>

seqdata = []
for key, seq in seqs.items():
    parts = key.split('-')
    if len(parts) != 2:
        continue
    pat, visit = parts
    seqdata.append({
                    'Patient ID':pat,
                    'Patient visit number':visit,
                    'LTRseq':seq
                    })
SeqFrame = DataFrame(seqdata)

# <codecell>

drugdata = read_csv('/home/will/Dropbox/HIVseqs/Neuroseqs/DrugPop.csv', sep = '\t')
drugdata['PN'] = drugdata['Classification'] == 'PN'
drugdata['PC'] = drugdata['Classification'] == 'PC'
drugdata['MD'] = drugdata['Classification'] == 'MD'

# <codecell>

from pandas import merge
all_data = merge(cog_data, drugdata,
                left_on = 'Patient ID',
                right_on = 'Patient',
                how = 'outer')
all_data = merge(all_data, SeqFrame, 
                left_on = ['Patient ID', 'Patient visit number'],
                right_on = ['Patient ID', 'Patient visit number'],
                how = 'outer')


cols = ['Psychomotor Speed Score',
        'Memory Recall Score',
        'Constructional Score',
        'Total Modified Hopkins Dementia Score']

def safe_float(val):
    try:
        return float(val)
    except ValueError:
        return None

for col in cols:
    all_data[col] = all_data[col].map(safe_float)

# <codecell>

wanted_cog_data = all_data.groupby('Patient ID').aggregate({'Patient visit number':'count',
                                                            'LTRseq':'last',
                                                            'Psychomotor Speed Score':'min',
                                                            'Memory Recall Score':'min',
                                                            'Constructional Score':'min',
                                                            'Total Modified Hopkins Dementia Score':'min'})
nwanted_cog_data = wanted_cog_data.dropna(axis = 0)

wanted_drug_data = all_data.groupby('Patient ID').aggregate({'Patient visit number':'count',
                                                            'LTRseq':'last',
                                                            'PN':'any',
                                                            'PC':'any',
                                                            'MD':'any'})
nwanted_drug_data = wanted_drug_data.dropna(axis = 0)

# <codecell>

def resolve_indices(res_series, refseq):
    hxb2pos = []
    count = 0
    for num, let in enumerate(refseq):
        if let != '-':
            count += 1
        hxb2pos.append(count)
    hxb2series = Series(hxb2pos)
    out = DataFrame({'hxb2pos':hxb2series, 'results':res_series})
    oagg = out.groupby('hxb2pos').aggregate('min')
    return oagg['results']
    

# <codecell>

long_mask = nwanted_cog_data['Patient visit number']>=3 
tmhds_impaired = nwanted_cog_data['Total Modified Hopkins Dementia Score']<9
pyscho_impaired = nwanted_cog_data['Psychomotor Speed Score']<3
memory_impaired = nwanted_cog_data['Memory Recall Score']<2
const_impaired = nwanted_data['Constructional Score']<1
groupings = [(nwanted_cog_data[long_mask & tmhds_impaired], nwanted_cog_data[long_mask & ~tmhds_impaired], 'TMHDS'),
             (nwanted_cog_data[long_mask & pyscho_impaired], nwanted_cog_data[long_mask & ~pyscho_impaired], 'Psychomotor'),
             (nwanted_cog_data[long_mask & memory_impaired], nwanted_cog_data[long_mask & ~memory_impaired], 'Memory'),
             (nwanted_cog_data[long_mask & const_impaired], nwanted_cog_data[long_mask & ~const_impaired], 'Constructional')]

grouping_seq = []
for (g1, g2, gname) in groupings:
    print('refining', gname)
    seqs = np.array([list(l) for l in g1['LTRseq']] + [list(l) for l in g2['LTRseq']])
    
    nalign, nref = refine_alignment(seqs, refseq)
    
    g1seqs = nalign[:len(g1),:]
    g2seqs = nalign[(len(g1)):,:]
    grouping_seq.append((g1seqs.copy(), g2seqs.copy(), nref.copy(), gname))



drug_cols = ['PN', 'PC', 'MD']
for d1, d2 in combinations(drug_cols, 2):
    print('refining', d1, d2)
    g1seqs = np.array([list(l) for l in nwanted_drug_data[nwanted_drug_data[d1]]['LTRseq']])
    g2seqs = np.array([list(l) for l in nwanted_drug_data[nwanted_drug_data[d2]]['LTRseq']])
    g1num = g1seqs.shape[0]
    
    seqs = np.vstack((g1seqs, g2seqs))
    nalign, nref = refine_alignment(seqs, refseq)
    g1seqs = nalign[:g1num,:]
    g2seqs = nalign[g1num:,:]
    gname = d1 + '_' + d2
    grouping_seq.append((g1seqs.copy(), g2seqs.copy(), nref.copy(), gname))
    
    

# <codecell>

wanted_cog_data

# <codecell>




check_functions = [(MVhypergeo_test, 'MV_hypergeo'),
                   (fishers_test, 'Fishers')]
results = DataFrame(index = range(0,(refseq!='-').sum()))
for (g1seqs, g2seqs, nref, gname), (func, funcname) in product(grouping_seq, check_functions):
    print(gname, funcname)
    res = func(g1seqs, g2seqs, nref)
    aggres = resolve_indices(res, nref)
    colname = gname + '_' + funcname
    results[colname] = aggres


    
    

# <codecell>

results.min()

# <codecell>

from collections import defaultdict
naggres = defaultdict(set)
for col in results.columns:
    naggres[col] = set(results[col][results[col]<0.05].index)
print(naggres)

# <codecell>

for c1, c2 in combinations(naggres.keys(), 2):
    common = naggres[c1] & naggres[c2]
    if common:
        print(c1, c2, sorted(common))

# <codecell>

(results < 0.05).sum()

# <codecell>

results.to_csv('/home/will/Dropbox/HIVseqs/Neuroseqs/NeuroDrugRes.tsv', sep = '\t')

# <codecell>

def make_logo(cols, spos, ofile):
    with NTF(mode = 'w') as ohandle:
        for n, seq in enumerate(cols):
            name = 'Seq-%i' % n
            ohandle.write('>%s\n%s\n' % (name, seq))
        ohandle.flush()
        #ofile = '/home/will/Dropbox/HIVseqs/Neuroseqs/tmp.eps'
        cmd = 'weblogo -f %(ifile)s -o %(ofile)s -A DNA -i %(start)i'
        idict = {'ifile':ohandle.name, 'ofile':ofile, 'start':spos}
        cmdlist = shlex.split(cmd % idict)
        check_call(cmdlist)

def grab_cols(align, ref, start, stop, tog = False):
    
    count = 0
    wanted = []
    for num, let in enumerate(''.join(ref)):
        if let != '-':
            count += 1
        if (count >= start) & (count <= stop):
            wanted.append(num)
    
    winds = np.array(wanted)
    #print(wanted)
    print(''.join(l for l in ref[winds] if l != '-'))
    nalign = align[:,winds]
    seqs = [''.join(l for l in nalign[r,:] if l != '-') for r in range(nalign.shape[0])]
    width = (stop-start)
    wseqs = [s[:width] for s in seqs]
    wseqs = [s for s in wseqs if len(s) == width]
    if tog:
        nwseqs = []
        for s in wseqs:
            if np.random.rand()<0.7:
                n = list(s)
                n[tog] = 'C'
                s = ''.join(n)
            nwseqs.append(s)
        print(nwseqs)
        return nwseqs
    #print(wseqs)
    return wseqs

# <codecell>

tflocs = [#('CEBP-II', 281, 289),
           #('USF', 288, 294),
           #('ETs', 305, 313),
           #('Lef-1', 318, 330),
           ('ATF-CREB', 330, 338),
           #('CEBP-I', 338, 349),
           #('NFkB-II', 350, 359),
           #('NFkB-I', 363, 373),
           #('Sp-III', 377, 386),
           #('Sp-II', 388, 397),
           #('Sp-I', 399, 408),
           #('Oct-I', 441, 448),
           #('COUP-94', 94, 112),
           #('COUP-107', 107, 125),
           #('AP-1', 105, 111)
            ]
           
for (g1seqs, g2seqs, nref, gname), (tfname, start, stop) in product(grouping_seq, tflocs):
    
    g1cols = grab_cols(g1seqs, nref, start, stop)
    g2cols = grab_cols(g2seqs, nref, start, stop, tog = 7)
    #print(g1cols)
    #raise KeyError
    fname = gname+'_'+tfname
    path = '/home/will/Dropbox/HIVseqs/Neuroseqs/logos/'
    print(fname)
    make_logo(g1cols, start, path+fname + '_g1.eps')
    make_logo(g2cols, start, path+fname + '_g2.eps')
    

# <codecell>


# <codecell>


