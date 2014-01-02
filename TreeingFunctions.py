# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

%load_ext autoreload
%autoreload 2

# <codecell>

from pandas import DataFrame, Series, merge, read_csv, MultiIndex, Index, concat
from subprocess import check_call
from tempfile import NamedTemporaryFile as NTF
import os, os.path
import numpy as np
from scipy.stats import ttest_ind
from itertools import groupby,combinations
from operator import itemgetter
from Bio import Phylo
import networkx

from random import shuffle
import csv, shlex, shutil

os.chdir('/home/will/Dropbox/HIVseqs/')
sys.path.append('/home/will/HIVReportGen/AnalysisCode/')
from SeqProcessTools import read_pat_seq_data, load_training_seq_data, align_seq_data_frame

# <codecell>

import glob

pat_files = glob.glob('/home/will/HIVReportGen/Data/PatientFasta/*.fasta')
pat_seq = read_pat_seq_data(pat_files, '/home/will/HIVReportGen/Data/BlastDB/ConBseqs.txt')

training_files = glob.glob('/home/will/HIVReportGen/Data/TrainingSequences/*.fasta')
training_data = load_training_seq_data(training_files)

align_lanl = align_seq_data_frame(training_data,  '/home/will/HIVReportGen/Data/BlastDB/ConBseqs.txt')
    
    


all_seqs = concat([pat_seq, align_lanl])

# <codecell>

def get_pairwise_distances(seq_series, tree_file = None, seq_file = None):
    
    if seq_file is None:
        fasta_handle = NTF()
    if tree_file is None:
        tree_handle = NTF()
    else:
        tree_handle = open(tree_file, 'w')
    for (pat, visit), seq in zip(seq_series.index, seq_series.values):
        nheader = '%s-%s' % (pat, visit)
        fasta_handle.write('>%s\n%s\n' % (nheader, ''.join(seq)))
    fasta_handle.flush()
    os.fsync(fasta_handle.fileno())
    cmd = 'muscle -in %(ifile)s -tree2 %(treefile)s -gapopen -2.9'
    cmdlist = shlex.split(cmd % {
                                 'ifile':fasta_handle.name, 
                                 'treefile':tree_handle.name
                                 })
    t = check_call(cmdlist)
    tree = Phylo.read(open(tree_handle.name), 'newick')
    seq_names = tree.get_terminals()
    dmat = {}
    for p1, p2 in combinations(seq_names, 2):
        d = tree.distance(p1, p2)
        dmat[(p1.name, p2.name)] = d
        dmat[(p2.name, p1.name)] = d
        
    return dmat

def extract_region(seq_series, start, stop):
    
    nseqs = seq_series.map(lambda x: x[start:stop])
    return nseqs
    

# <codecell>

def check_distance_pvals(mat_data, trop_dict):
    nreps = 500
    frac = 0.5
    g1dist = []
    g2dist = []
    for (key1, key2), dist in mat_data.items():
        if trop_dict[key1] and trop_dict[key2]:
            g1dist.append(dist)
        elif not trop_dict[key1] and not trop_dict[key2]:
            g2dist.append(dist)
    nitems = int(min(frac*len(g1dist), frac*len(g2dist)))
    _, raw_pval = ttest_ind(g1dist, g2dist)
    cor_pvals = []
    for _ in range(nreps):
        shuffle(g1dist)
        shuffle(g2dist)
        _, pval = ttest_ind(g1dist[:nitems], g2dist[:nitems])
        cor_pvals.append(pval)
    return raw_pval, np.mean(cor_pvals), np.mean(g1dist), np.mean(g2dist), np.std(g1dist), np.std(g2dist)

# <codecell>

pssm_data = read_csv('/home/will/HIVReportGen/Data/TrainingSequences/pssm_data.csv', index_col = [0,1])
def decide_tropism(inval):
    if inval < -6.95:
        return True
    elif inval > -2.88:
        return False
    return np.nan
tropism_data = pssm_data['score'].map(decide_tropism).dropna()
trop_dict = {}
for (pat, visit), val in zip(tropism_data.index, tropism_data.values):
    trop_dict[pat+'-'+visit] = val
    
benj_selected = []
with open('BensTropismLabels.csv') as handle:
    reader = csv.DictReader(handle)
    for row in reader:
        trop_dict['%s-%s' % (row['Patient ID'], row['Visit'])] = row['Prediction'] == 'TRUE'
        benj_selected.append((row['Patient ID'], row['Visit']))
benj_selected_index = MultiIndex.from_tuples(benj_selected, names = ['Patient ID', 'Visit number'])

# <codecell>

pure_seqs = MultiIndex.from_tuples([
('A0001','R00'),#R5 XX
('A0107','R05'),#X4 XX
('A0017','R02'),#X4 XX
('AB286955','RN'),#R5 XX
('AB287367','RN'),#R5 XX
('AB480695','RN'),#X4 XX
('AB485642','RN'),#X4 XX
('AB604946','RN'),#X4 XX
('AF042101','RN'),#X4 XX
('AY835766','RN'),#R5 XX
('AY835779','RN'),#X4 XX
('AY352275','RN'),#X4 
('AY970950','RN'),#R5 XX
('DQ358809','RN'),#R5 XX
('EF057102','RN'),#X4 XX
('EF363123','RN'),#R5 XX
('JQ316126','RN'),#X4 XX
('GU647196','RN'),#X4 XX
('DQ990880','RN'),#X4 XX
], names = ['Patient ID', 'Visit number'])

equal_pure_seqs = MultiIndex.from_tuples([
('A0001','R00'),#R5 XX
('A0107','R05'),#X4 XX
('A0017','R02'),#X4 XX
('AB286955','RN'),#R5 XX
('AB287367','RN'),#R5 XX
('AB480695','RN'),#X4 XX
('AB485642','RN'),#X4 XX
('AB604946','RN'),#X4 XX
('AF042101','RN'),#X4 XX
('AY835766','RN'),#R5 XX
('AY835779','RN'),#X4 XX
('AY352275','RN'),#X4 
('AY970950','RN'),#R5 XX
('DQ358809','RN'),#R5 XX
('EF057102','RN'),#X4 XX
('EF363123','RN'),#R5 XX
('JQ316126','RN'),#X4 XX
#('GU647196','RN'),#X4 XX
('DQ990880','RN'),#X4 XX
], names = ['Patient ID', 'Visit number'])




# <codecell>

def make_tree_figure(wanted_seqs, trop_dict, tree_file):
    mat_data = get_pairwise_distances(wanted_seqs, tree_file = tree_file)
    tree = Phylo.read(open(tree_file), 'newick')
    net = Phylo.to_networkx(tree)
    
    node_mapping = {}
    clade = 1
    for node in net.nodes():
        if node.name is None:
            node_mapping[node] = 'Clade-%i' % clade
            clade += 1
        else:
            node_mapping[node] = node.name
    new_net = networkx.relabel_nodes(net, node_mapping)
    
    colors = []
    for node in new_net.nodes():
        if node.startswith('Clade'):
            colors.append('w')
        elif trop_dict[node]:
            colors.append('g')
        elif not trop_dict[node]:
            colors.append('r')
        else:
            print node
    #print colors, len(colors), len(new_net.nodes())
    pos = networkx.graphviz_layout(new_net, 'twopi')
    
    networkx.draw_networkx(new_net, pos, with_labels = False, node_color = colors)
    
    
    
    

# <codecell>

check_regions = [#('Tat-seq-align', 'The Acidic domain', 0, 20), 
                 #('Tat-seq-align', 'Cysteine rich domain', 21, 36),
                 #('Tat-seq-align', 'Core domain', 37, 47),
                 #('Tat-seq-align', 'TAR binding domain', 48, 56),
                 #('Tat-seq-align', 'Domain V-72', 57, 71),
                 #('Tat-seq-align', 'Domain V-86', 57, 85),
                 #('Tat-seq-align', 'Exon II 73', 72, 100),
                 #('Tat-seq-align', 'Exon II 87', 86, 100),
                 #('Tat-seq-align', 'Transactivation', 0, 47),
                 #('Tat-seq-align', 'Co-factor binding', 21, 48),
                 #('Tat-seq-align', 'SP1 binding', 29, 54),
                 #('Tat-seq-align', 'Basic Region', 48, 71),
                 #('Tat-seq-align', 'CEBP binding', 46, 66),
                 #('Tat-seq-align', 'NFAT binding', 0, 25),
                 #('Tat-seq-align', 'DNA-PK binding', 55, 100),
                 #('Vpr-seq-align', 'Nuclear localization', 10, 39),
                 #('Vpr-seq-align', 'Cell Cycle Progression', 14, 34),
                 #('Vpr-seq-align', 'Tansactivation', 13, 21),
                 #('Vpr-seq-align', 'Viron Packaging', 28, 39),
                 #('Vpr-seq-align', 'Nuclear localizations', 53, 74),
                 #('Vpr-seq-align', 'Transactivation', 73, 80),
                 #('Vpr-seq-align', 'G2 arrest', 74, 94),
                 #('Vpr-seq-align', 'DNA binding', 84, 92),
                 ('LTR-seq-align', 'U3', 0, 455),
                 ('LTR-seq-align', 'R', 456, 559),
                 ('LTR-seq-align', 'U5', 560, 612),
                 ('LTR-seq-align', 'TAR', 454, 544),
                 ('LTR-seq-align', 'Integration', 0, 70),
                 ('LTR-seq-align', 'AP1-COUPs', 60, 250),
                 ('LTR-seq-align', 'CEBP-Lef-1s', 280, 330),
                 ('LTR-seq-align', 'SP-sites', 376, 408),
                 ('LTR-seq-align', 'AP1-CREB', 539, 616),
                 ('LTR-seq-align', 'Pre-SP-I', 408, 454),
                 ('LTR-seq-align', 'Pre-SP-I-upstream-half', 408, 431),
                 ('LTR-seq-align', 'Pre-SP-I-downstream-half', 431, 454),
                 ('LTR-seq-align', 'GC-Box', 376, 408),
                 ('LTR-seq-align', 'SP-I', 398, 408),
                 ('LTR-seq-align', 'SP-II', 387, 398),
                 ('LTR-seq-align', 'SP-III', 376, 386),
                 ('LTR-seq-align', 'NfKB-SP-III', 349, 386),
                 ('LTR-seq-align', 'NfKB-II-SP-III', 362, 386),
                 ('LTR-seq-align', 'CEBP-I-NF2', 337, 359),
                 ('LTR-seq-align', 'ATF-CREB-CEBP', 329, 349),
                 ('LTR-seq-align', 'LEF1-CREB', 317, 337),
                 ('LTR-seq-align', 'LEF-1', 317, 330),
                 ('LTR-seq-align', 'ATF-CREB',329, 337),
                 ('LTR-seq-align', 'CEBP-I', 337, 344),
                 ('LTR-seq-align', 'ETS-1', 304, 313),
                 ('LTR-seq-align', 'CEBP-II-USF-1', 280, 294),
                 ('LTR-seq-align', 'AP-I-to-CEBP-II', 221, 280),
                 ('LTR-seq-align', 'AP-I-promixal-half', 221, 251),
                 ('LTR-seq-align', 'CEBP-II-promixal-half', 251, 280),
                 ('LTR-seq-align', 'AP-I', 213, 221),
                 ('LTR-seq-align', 'GRE', 191, 207),
                 ('LTR-seq-align', 'AP-II-to-GRE', 162, 191),
                 ('LTR-seq-align', 'AP-II', 154, 162),
                 ('LTR-seq-align', 'COUP-to-AP-II', 131, 154),
                 ('LTR-seq-align', 'COUP', 93, 131),
                 ('LTR-seq-align', 'Pre-COUP', 0, 93),
                 ('LTR-seq-align', 'Pre-COUP-upstream-half', 0, 45),
                 ('LTR-seq-align', 'Pre-COUP-downstream-half', 45, 93),
                 ('LTR-seq-align', 'NfKB-I', 362, 373),
                 ('LTR-seq-align', 'NfKB-II', 349, 359),
                 ('LTR-seq-align', 'NfKB-I-NfKB-II', 349, 373),
                 ('LTR-seq-align', 'CEBP-I', 337, 349),
                 ('LTR-seq-align', 'CEBP-II', 280, 289),
                 ('LTR-seq-align', 'COUP-I', 116, 130),
                 ('LTR-seq-align', 'COUP-II', 105, 124),
                 ('LTR-seq-align', 'COUP-III', 92, 111),
                 ('LTR-seq-align', 'AP-III', 118, 125),
                 ('LTR-seq-align', 'AP-IV', 103, 110),
                    ]
                 

indexes = [#('All LANL Seqs', tropism_data.index),
           #('BenJ Selected', benj_selected_index),
           #('Benj Pure Seq', pure_seqs),
            ('Benj Pure Equal Seq', equal_pure_seqs), 
]

# <codecell>

from itertools import product


results = []
for (ind_name, inds), (seq_col, name, start, stop) in product(indexes, check_regions):
    wanted = extract_region(all_seqs.ix[inds][seq_col].dropna(), start, stop)
    #print wanted.index
    #print('Treeing')
    prot_name = seq_col.split('-')[0]
    treename = 'fixeddomaintrees/%s-%s-%s-%i.nwk' % (ind_name, prot_name, name, start)
    treename = treename.replace(' ', '-')
    mat_data = get_pairwise_distances(wanted, tree_file=treename)
    #print('Testing')
    raw_p, cor_p, r5mean, x4mean, r5std, x4std = check_distance_pvals(mat_data, trop_dict)
    
    if seq_col.startswith('LTR'):
        start = start-454
        stop = stop-454
    print ind_name, name, start, raw_p, cor_p
    results.append((ind_name, prot_name, name,start+1, stop+1, cor_p, r5mean, x4mean, r5std, x4std))
    #pat_wanted = extract_region(pat_seq.ix[tropism_data.index][seq_col].dropna(), start, stop)
    #fname = '/home/will/BenSeqs/Trees/' + name.replace(' ', '-')
    #plt.figure(figsize = (20,20))
    #make_tree_figure(wanted, trop_dict, fname + '.tree')
    #plt.title(name)
    #plt.savefig( fname + '_pure_seqs.png')
    

# <codecell>

with open('fixed_domain_analysis_with_new_ltr.csv', 'w') as handle:
    writer = csv.writer(handle)
    fields = ['Sequence Set', 'Protein Name', 'Domain Name', 'Region Start', 'Region Stop', 'p-value', 'R5-Mean', 'X4-Mean', 'R5-std', 'X4-std']
    writer.writerow(fields)
    writer.writerows(results)

# <codecell>

#widths = [5,10,15,20,25,30,35,45,50]

#indexes = [('BenJ Selected', benj_selected_index),
#            ('All LANL Seqs', tropism_data.index),
#           ]
            #('Benj Pure Seq', pure_seqs)

#['Sequence Set', 'Protein Name', width, 'midpoint',  'p-value']
#large_results = []
#prots = [('Vpr-seq-align', 'Vpr',range(96))]
    
    
#for (ind_name, inds), width, (seq_col, prot, positions) in product(indexes, widths, prots):
#    
#    for mid in positions:
#        print ind_name, width, prot, mid
#        start = max(int(mid-(width/2)),0)
#        stop = min(int(mid+(width/2)),positions[-1])
#        wanted = extract_region(all_seqs.ix[inds][seq_col].dropna(), start, stop)
#        mat_data = get_pairwise_distances(wanted)
#        raw_p, cor_p, r5mean, x4mean = check_distance_pvals(mat_data, trop_dict)
#        
#        large_results.append((ind_name, prot, width,start+1, stop+1, cor_p, r5mean, x4mean))

# <codecell>

import contextlib
from tempfile import mkdtemp

@contextlib.contextmanager
def tmp_directory(*args, **kwargs):
    """A context manager which changes the working directory to the given
    path, and then changes it back to its previous value on exit.

    """
    path = mkdtemp(*args, **kwargs)
    try:
        yield path + '/'
    finally:
        #shutil.rmtree(path)
        pass

# <codecell>

from StringIO import StringIO
from subprocess import check_output
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import IUPAC

def write_nexus_alignment(seq_series, handle):
    seqs = []
    tmp_handle = StringIO()
    for (pi, vn), seq in zip(seq_series.index, seq_series.values):
        nseq = ''.join(seq).replace('O', '-')
        bseq = SeqRecord(Seq(nseq, alphabet=IUPAC.protein), id = '%s-%s' % (pi, vn))
        seqs.append(bseq)
    SeqIO.write(seqs, tmp_handle, 'nexus')
    tmp_handle.seek(0)
    strdata = tmp_handle.read().replace("'", '')
    handle.write(strdata)
    
def write_mrbayes_commands(handle, alignment_path, output_path):
    
    cmd = """begin mrbayes;
   set autoclose=yes nowarn=yes;
   execute %(align)s;
   prset aamodelpr = mixed;
   mcmc nchains = 3 ngen = 50000 samplefreq=1000 diagnfreq=100000 printfreq=100000 file=%(out)s;
   sump;
   sumt;
end;"""
    
    tdict = {'align':alignment_path, 'out':output_path}
    handle.write(cmd % tdict)
    
def run_mrbayes(cmd_path):
    
    cmd = '/home/will/mb ' + cmd_path 
    check_output(shlex.split(cmd))
    
def reformat_nexus(inhandle, outhandle, trop_dict):
    
    def process_tree_line(line):
        
        parts = line.strip().split()
      
        return 'tree %s [&R] %s\n' % (parts[1], parts[-1])
    
    
    for line in inhandle: #get ri of junk
        if line.strip() == 'begin trees;':
            break
    _ = inhandle.next() #get rid of the 'translate' line
    outhandle.write('#NEXUS\n\n\n')
    
    outhandle.write('begin states;\n')
    for line in inhandle:
        nline = line.strip()
        if nline.startswith('tree'):
            first_tree = process_tree_line(line)
            break
        num, seqname = nline[:-1].split(' ', 1)
        try:
            if trop_dict[seqname.replace('.copy', '')]:
                trop = 'R5'
            else:
                trop = 'X4'
        except KeyError:
            print 'Missing ' + seqname + ' !!'
            trop = 'R5'
            
        outhandle.write('%s %s\n' % (num, trop))
    outhandle.write('End;\n\n')
    outhandle.write('begin trees;\n')
    
    tree_lines = [first_tree] + [process_tree_line(line) for line in inhandle if line.strip() != 'end;']
    for line in tree_lines:
        outhandle.write(line)
    outhandle.write('end;\n')
    
    
def run_bats(formated_nexus_path, nreps = 5000):
    
    cmd = 'java -Xmx3000M -jar /home/will/BaTS_beta_build2.jar single %s %i %i'
    out = check_output(shlex.split(cmd % (formated_nexus_path, nreps, 2)))
    
    handle = StringIO(out)
    for line in handle:
        if line.startswith('Stat'):
            headers = line.strip().split('\t')
            break
    return list(csv.DictReader(handle, fieldnames=headers, delimiter = '\t'))[:-2]

    
    

    
def run_MrBats_analysis(seq_series, trop_dict, tree_file):
    
    
    with tmp_directory(dir = '/home/will/tmpstuf/') as tmpdir:
        align_file = tmpdir + 'seqalign.nxs'
        mrbayes_cmd_file = tmpdir + 'analysis.nxs'
        cons_file = tmpdir + 'seqalign.nxs.con.tre'
        multi_prob = tmpdir + 'seqalign.nxs.trprobs'
        multi_mod = tmpdir + 'seqalign.nxs.trprobs.modified'
        #print align_file
        
        with open(align_file, 'w') as handle:
            #print align_file, len(seq_series)
            write_nexus_alignment(seq_series, handle)
            
        with open(mrbayes_cmd_file, 'w') as handle:
            write_mrbayes_commands(handle, align_file, align_file)
        
        run_mrbayes(mrbayes_cmd_file)
        
        with open(multi_prob) as inhandle:
            with open(multi_mod, 'w') as ohandle:
                reformat_nexus(inhandle, ohandle, trop_dict)
        out = run_bats(multi_mod)
        #out = [{}]
        
        if tree_file:
            shutil.copy(cons_file, tree_file)
        
    return out
    
    
    
    
    

# <codecell>

from copy import deepcopy
from concurrent.futures import ThreadPoolExecutor
from itertools import chain, imap, product

indexes = [#('BenJ Selected', benj_selected_index),
            ('Benj Pure Seq', pure_seqs),
            ('Benj Pure Equal Seq', equal_pure_seqs), 
            #('All LANL Seqs', tropism_data.index),
]

def linker_code(tup):
    
    wanted_seqs, treename, extra_dict = tup
    out = run_MrBats_analysis(wanted_seqs, trop_dict, tree_file = treename)
    final = []
    for row in out:
        row.update(extra_dict)
        final.append(row)
    return final
        
    

#(ind_name, prot_name, name,start+1, stop+1, cor_p, r5mean, x4mean)
new_method_inputs = []
for (ind_name, inds), (seq_col, name, start, stop) in product(indexes, check_regions):
    wanted = extract_region(all_seqs.ix[inds][seq_col].dropna(), start, stop)
    #print wanted.index
    #print('Treeing')
    prot_name = seq_col.split('-')[0]
    treename = 'newdomaintrees/%s-%s-%s-%i.nwk' % (ind_name, prot_name, name, start)
    extra_dict = {
                    'IndName':ind_name,
                    'ProtName':prot_name,
                    'Domain':name,
                    'Start':start,
                    'Stop':stop,
                 }
    treename = treename.replace(' ', '-')
    #print treename, len(wanted)
    if len(wanted)>10:
        new_method_inputs.append((wanted.copy(), treename, deepcopy(extra_dict)))
#raise KeyError
    
#results_so_far = []

with ThreadPoolExecutor(max_workers = 30) as executor:
    res = executor.map(linker_code, new_method_inputs)
    for row in chain.from_iterable(res):
        print row['Domain'], row['IndName'], row['Statistic'] ,row['significance']
        #results_so_far.append(row)

# <codecell>

tmp = DataFrame(results_so_far)
tmp.to_csv('new_method_results.csv')

# <codecell>

from itertools import islice
widths = [5,10,15,20,25,30,35,45,50]
indexes = [('BenJ Selected', benj_selected_index),
                ('All LANL Seqs', tropism_data.index),
            ('Benj Pure Seq', pure_seqs),
           ]
prots = [('Vpr-seq-align', 'Vpr',range(96))]

def take(n, iterable):
    "Return first n items of the iterable as a list"
    return list(islice(iterable, n))

def yield_regions(indexes, widths, prots):
    for (ind_name, inds), width, (seq_col, prot, positions) in product(indexes, widths, prots):
        for mid in positions:
            #print ind_name, width, prot, mid
            start = max(int(mid-(width/2)),0)
            stop = min(int(mid+(width/2)),positions[-1])
            wanted = extract_region(all_seqs.ix[inds][seq_col].dropna(), start, stop)
            extra_dict = {
                    'IndName':ind_name,
                    'ProtName':prot,
                    'Start':start,
                    'Stop':stop,
                    'Mid':mid,
                    'width':width
                 }
    
            if len(wanted)>10:
                yield (wanted.copy(), None, deepcopy(extra_dict))


#['Sequence Set', 'Protein Name', 'width', 'midpoint',  'p-value']
vpr_window_results = []
block_size = 500
with ThreadPoolExecutor(max_workers = 30) as executor:
    
    iterable = yield_regions(indexes, widths, prots)
    block = take(block_size, iterable)
    
    while block:
        res = executor.map(linker_code, block)
        for row in chain.from_iterable(res):
            print row['width'], row['IndName'], row['Statistic'] ,row['significance']
            vpr_window_results.append(row)
        block = take(block_size, iterable)
    
    
    
    
    

    
    

# <codecell>

windowed_df = DataFrame(vpr_window_results)

# <codecell>

windowed_df.head()

# <codecell>

import dendropy

fixtrees = glob.glob('newdomaintrees/*.nwk')
for f in fixtrees:
    if 'Equal' not in f:
        continue
    with open(f) as handle:
        tree = dendropy.Tree.get_from_stream(open(f), 'nexus')
        
    tree.deroot()
    rmnodes = [tree.prune_subtree(t, update_splits = True) for t in tree.leaf_nodes() if t.get_node_str().endswith("copy'")]
    #tree.prune_taxa(rmnodes)
    nf = f.replace('newdomaintrees', 'unrootedtrees-equal')
    with open(nf, 'w') as handle:
        tree.write_to_stream(handle, 'newick')

# <codecell>

tmp = list(tree.leaf_nodes())

# <codecell>

labels = [n.get_node_str() for n in rmnodes]
tree.prune_taxa_with_labels(labels, update_splits = True)

# <codecell>

[tree.prune_subtree(t) for t in tree.leaf_nodes() if t.get_node_str().endswith("copy'")]

# <codecell>

with open('/home/will/tmpstuf/tmpyw6FXN/seqalign.nxs.trprobs') as handle:
    treeL = dendropy.TreeList.get_from_stream(handle, 'nexus')

# <codecell>

print treeL[0].description()

# <codecell>

with open('/home/will/tmpstuf/tmpyw6FXN/tmptree.nxs', 'w') as handle:
    treeL.write_to_stream(handle, 'nexus')

# <codecell>

for tree in treeL:
    for leaf in tree.leaf_iter():
        print str(leaf.taxon)

# <codecell>


