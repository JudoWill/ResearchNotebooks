# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

from __future__ import division
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

from types import ListType
from itertools import combinations, groupby, islice, imap
from collections import Counter
from subprocess import check_output
from operator import itemgetter
from StringIO import StringIO
import csv
import shlex
import sys
import glob
sys.path.append('/home/will/PySeqUtils/')
from GeneralSeqTools import fasta_reader, fasta_writer

# <codecell>



def unique_justseen(iterable, key=None):
    "List unique elements, preserving order. Remember only the element just seen."
    # unique_justseen('AAAABBBCCDAABBB') --> A B C D A B
    # unique_justseen('ABBCcAD', str.lower) --> A B C A D
    return imap(next, imap(itemgetter(1), groupby(iterable, key)))



def check_overlap(tup1, tup2, overlap = 20):
    sa, Ia = tup1
    sb, Ib = tup2
    if ((Ia+len(sa)) - Ib) < overlap:
        return False
    else:
        nS = max(Ia, Ib)
        nE = min(Ia+len(sa), Ib+len(sb))
        A = sa[(nS-Ia):(nE-Ia)]
        B = sb[(nS-Ib):(nE-Ib)]
        #print A, B
        return A == B        

def guess_seq_col(row):
    wlets = 'ATCG-'
    skip_lets = '*HIVREFX:|1234567890'
    
    for num, col in enumerate(row):
        for l in skip_lets:
            if l in col:
                break
        if l in col:
            continue
        if any(l in col for l in wlets):
            print num, col
            return num
    print row
    raise KeyboardInterrupt
    
    
def get_reads(fname):
    if type(fname) is ListType:
        count = 0
        for f in fname:
            for (se, st, ind) in get_reads(f):
                yield se, st, 'join-read-%i' % count
                count += 1
    else:
        if fname.endswith('.bam'):
            cmd = 'samtools view %s' % fname
            out = check_output(shlex.split(cmd))
        else:
            with open(fname) as handle:
                out = ''.join(line for line in handle if not line.startswith('@'))
        reader = csv.reader(StringIO(out), delimiter = '\t')
        #seq_col = guess_seq_col(reader.next())
        getter = itemgetter(9, 3)
        reads = imap(getter, reader)
        for k, (se, st) in enumerate(reads):
            start = int(st)
            if start < 634: #fix the LTR issue
                start += 9086
            yield se, start, 'read-%i' % k
    
    
    
def get_IR_reads(reads):
    n_reads = sorted(reads, key = itemgetter(0,1))
    ir_reads = list(unique_justseen(n_reads, key = itemgetter(0,1)))
    print 'IR reads', len(ir_reads)
    return ir_reads
    

def initialize_graph(ir_reads):
    
    graph = nx.DiGraph()
    for seq, start, ind in ir_reads:
        graph.add_node(ind, seq = seq, start = int(start), end = int(start)+len(seq))
        
    nover = 0
    for (sa, Ia, indA), tup_rows in groupby(combinations(ir_reads,2), lambda tups: tups[0]):
        end_pos = Ia+len(sa)
        for _, (sb, Ib, indB) in tup_rows:
            if end_pos < Ib:
                break
            if check_overlap((sa, Ia), (sb, Ib), overlap=5):
                nover += 1
                graph.add_edge(indA, indB)
    return graph


def add_terminals(graph, ir_reads, start_pos, end_pos):
    
    graph.remove_nodes_from(['source', 'target'])
    
    graph.add_node('source')
    graph.add_node('target')
    
    starts = 0
    targets = 0
    for seq, start, ind in ir_reads:
        end = start+len(seq)
        if (start <= start_pos) & (end > start_pos):
            graph.add_edge('source', ind)
            starts += 1
        if (start < end_pos) and (end >= end_pos):
            graph.add_edge(ind, 'target')
            targets += 1
            
    #print 'Start edges: %i Sink Edges: %i' % (starts, targets)
    return graph


def generate_coverage_map(ir_reads, ax = None, label = None):
    
    pos_count = Counter()
    for seq, start, ind in ir_reads:
        for s in range(start, start+len(seq)):
            pos_count[s] += 1

    x = range(0, max(pos_count.keys()))
    y = [pos_count[p] for p in x]
    
    if ax is None:
        ax = plt.subplot(111)
    ax.plot(x, y, label = label)
    ax.set_ylabel('Coverage')
    
    return ax
    
    

# <codecell>

def simple_paths(graph):
    
    #target_paths = nx.single_source_shortest_path(graph, 'target')
    target_dict = {}
    for node in graph.nodes_iter():
        try:
            target_dict[node] = nx.shortest_path(graph, node, 'target')
        except nx.NetworkXNoPath:
            pass
    #print len(target_dict)
    
    source_dict = nx.single_source_shortest_path(graph, source = 'source')
    #print len(source_dict)
    
    for key in (set(source_dict.keys()) & set(target_dict.keys())):
        #print source_dict[key], target_dict[key]
        yield source_dict[key] + target_dict[key][1:]


    
        
def get_common(col):
    try:
        return Counter(col[col != '-']).most_common(1)[0][0]
    except IndexError:
        return '-'

def adjust_pos(seq, start, hap_begin, hap_end):
    
    
    norm_start = start
    norm_stop = start + len(seq)
    
    ome_start = max(hap_begin, norm_start)
    hap_start = ome_start - hap_begin
    seq_start = max(0, hap_begin-norm_start)
    
    ome_stop = min(norm_stop, hap_end)
    hap_stop = ome_stop - hap_begin
    seq_stop = (hap_stop - hap_start) + seq_start
    
    return seq[seq_start:seq_stop], hap_start, hap_stop

def assemble_haps(paths, graph, hap_begin, hap_end):
    hap_len = hap_end-hap_begin
    haps = set()
    for num, path in enumerate(paths):
        if num % 50 == 0:
            print '%i paths, %i haps' % (num, len(haps))
        seq = np.empty((len(path)-2, hap_len), dtype=str)
        seq[:] = '-'
        for row, node_name in enumerate(path, -1):
            if (node_name != 'source') and (node_name != 'target'):
                node = graph.node[node_name]
                #print node['start'], node['start']+len(node['seq'])
                nseq, hap_start, hap_stop = adjust_pos(node['seq'], node['start'], hap_begin, hap_end)
                try:
                    seq[row, hap_start:hap_stop] = np.array(list(nseq))
                except ValueError:
                    pass
        tmp = ''.join([get_common(seq[:,col]) for col in range(seq.shape[1])])
        #if '-' in tmp:
        #    print 'bad one!'
        #    continue
        #    for row in range(len(path)-2):
        #        print ''.join(seq[row,:])
        #    for row, node_name in enumerate(path, -1):
        #        if (node_name != 'source') and (node_name != 'target'):
        #            node = graph.node[node_name]
        #            #print node['start'], node['start']+len(node['seq'])
        #            nseq, hap_start, hap_stop = adjust_pos(node['seq'], node['start'], hap_begin, hap_end)
        #            print node['seq'], node['start'], nseq, hap_start, hap_stop
        #    
        #    raise KeyboardInterrupt
        haps.add(tmp)
    print len(haps), 'unique haplotypes found!'
    return haps
    

# <codecell>


bases = sorted(['DrexelMed.A0010',  'DrexelMed.A0017', 'DrexelMed.A0019', 'DrexelMed.A0107', 'DrexelMed.A0121',
         'sim_reads', 'DrexelMed.A0017.R02', 'DrexelMed.A0041.R02', 'DrexelMed.A0107.R02', 
         'DrexelMed.A0220'])
methods = ['lastz']#, 'bwa', 'ngm']

fig, axs = plt.subplots(len(bases),2, figsize = (10,20))

for num, bname in enumerate(bases):
    
    for method in methods:
        fname = '/home/will/DeepPipeline/Data/MappingResults/' + bname + '.' + method + '.map'
        all_reads = list(get_reads(fname))
       
        big_ax = axs[num, 0]
        zoom_ax = axs[num, 1]
    
        big_ax = generate_coverage_map(all_reads, ax = big_ax, label = method)
        big_ax.set_xlim([0, 9719])
        big_ax.set_title(bname)
        #big_ax.set_yscale('log')
    
        zoom_ax = generate_coverage_map(all_reads, ax = zoom_ax, label = method)
        zoom_ax.set_xlim([9086, 9719])
        zoom_ax.set_title(bname)
        #zoom_ax.set_yscale('log')
        if zoom_ax.is_first_row():
            zoom_ax.legend()
        
    

fig.tight_layout()
plt.savefig('/home/will/Downloads/HIVDeepSequencingMapping.png')

# <codecell>

sys.path.append('/home/will/DeepPipeline/AnalysisCode/')

# <codecell>

#reload(HapReconTools)

fname = '/home/will/DeepPipeline/Data/ShoRAHruns/DrexelMed.A0010.lastz/tmp.bam'
reads = sorted(HapReconTools.read_from_bam(fname), key = itemgetter(1))
ir_reads = list(HapReconTools.yield_IR_reads(reads))
graph = HapReconTools.generate_hap_graph(ir_reads)

# <codecell>

for tup in islice(graph.in_degree_iter(), 5):
    print tup

# <codecell>

def new_simple_paths(graph):
    
    nodes = graph.in_degree_iter()
    for node, in_degree in nodes:
        if in_degree == 0:
            tmp_tree = nx.dfs_tree(graph, node)
            out_degree = tmp_tree.out_degree()
            for tree_node in tmp_tree.nodes_iter():
                if out_degree[tree_node] == 0:
                    path = nx.shortest_path(graph, node, tree_node)
                    if len(path) > 2:
                        yield path


# <codecell>

def new_assemble_haps(paths, graph):
    haps = set()
    for num, path in enumerate(paths):
        if num % 500 == 0:
            print '%i paths, %i haps' % (num, len(haps))
        hap_begin = graph.node[path[0]]['start']
        hap_end = graph.node[path[-1]]['stop']
        seq = np.empty((len(path), hap_end-hap_begin), dtype=str)
        seq[:] = '-'
        for row, node_name in enumerate(path, 0):
            
            node = graph.node[node_name]
            nseq, hap_start, hap_stop = adjust_pos(node['seq'], node['start'], hap_begin, hap_end)
            
            try:
                seq[row, hap_start:hap_stop] = np.array(list(nseq))    
            except:
                print node['seq'], node['start']
                print nseq, hap_start, hap_stop
                seq[row, hap_start:hap_stop] = np.array(list(nseq))    
                #raise KeyboardInterrupt
            
                
        tmp = ''.join([get_common(seq[:,col]) for col in range(seq.shape[1])])
        haps.add((tmp, hap_begin))
    print len(haps), 'unique haplotypes found!'
    return haps
    
haps = new_assemble_haps(new_simple_paths(graph), graph)
sorted_haps = sorted(haps, key=itemgetter(0))
ir_haps = list(HapReconTools.yield_IR_reads(sorted_haps))
print 'IR haps:', len(ir_haps)

# <codecell>

def new_generate_coverage_map(ir_reads, ax = None, label = None):
    
    pos_count = Counter()
    for seq, start in ir_reads:
        for s in range(start, start+len(seq)):
            pos_count[s] += 1

    x = range(0, max(pos_count.keys()))
    y = [pos_count[p] for p in x]
    
    if ax is None:
        ax = plt.subplot(111)
    ax.plot(x, y, label = label)
    ax.set_ylabel('Coverage')
    
    return ax

new_generate_coverage_map(ir_haps)

# <codecell>

from collections import deque
from copy import deepcopy

def get_read_counts(ir_reads):
    
    read_counts = Counter(seq for seq, _, in ir_reads)
    return read_counts


def log_likelihood(urh, prh, ph):
    
    inner = ph*prh
    inner[inner > 0] = np.log10(inner[inner>0])
    
    return (urh*inner).sum()
    

def initialize(read_counter, haps):
    
    invalid_reads = []
    for read in read_counter.keys():
        keep = False
        for hap in haps:
            if read in hap:
                keep = True
                break
        if not keep:
            invalid_reads.append(read)
    print len(invalid_reads)
    for read in invalid_reads:
        del(read_counter[read])
    print 'valid reads!', sum(read_counter.values())
    sorted_uni_reads = sorted(read_counter.keys())
    sorted_uni_haps = sorted(haps)
    prh = np.zeros((len(sorted_uni_reads), len(sorted_uni_haps)))
    urh = np.zeros((len(sorted_uni_reads), len(sorted_uni_haps)))
    
    for row, read in enumerate(sorted_uni_reads):
        con_haps = [col for col, hap in enumerate(haps) if read in hap]
        if con_haps:
            prh[row, con_haps] = 1/read_counter[read]
            urh[row, con_haps] = read_counter[read]/len(con_haps)
    ph = np.ones((len(sorted_uni_haps),))/len(sorted_uni_haps)
    #ph = Ph(urh)
    
    return sorted_uni_reads, sorted_uni_haps, urh, prh, ph 


def Urh(urh, prh, ph):
    
    ur = urh.sum(axis=1).reshape(-1, 1)
    num = ph*prh
    pr = num.sum(axis=1).reshape(-1, 1)
    nurh = np.zeros_like(urh)
    if np.isnan(nurh).any().any():
        raise KeyboardInterrupt
    return ur*num/pr

def Ph(urh):
    
    ur = urh.sum(axis=1).reshape(-1,1)
    ph = (urh/ur)
    if np.isnan(ph).any().any():
        raise KeyboardInterrupt
    return ph.mean(axis=0)


def estimate_min_bound(ir_reads, hap_begin, hap_end):
    
    p = 0.9 #prob of finding a haplotype
    N = len(ir_reads)
    L = sum(len(seq) for seq, _, in ir_reads)/N #avg read length
    n = hap_end-hap_begin
    
    return -n*np.log(1-p**(1/n))/(N*L)


def do_EM(urh, prh, ph, max_rep = 10000, tol = 0.001):
    
    urh_l = deque([urh], maxlen = 2)
    ph_l = deque([ph], maxlen = 2)
    ll_l = deque([log_likelihood(urh_l[-1], prh, ph_l[-1])], maxlen = 2)
    
    for rep in range(1,max_rep):
        
        urh_l.append(Urh(urh_l[-1], prh, ph_l[-1]))
        ph_l.append(Ph(urh_l[-1]))
        ll_l.append(log_likelihood(urh_l[-1], prh, ph_l[-1]))
        
        if np.log10(rep) == int(np.log10(rep)):
            print ll_l[-1], rep, ph_l[-1].sum(), ph_l[-1].max()
        
        if abs(ll_l[-1] - ll_l[-2]) < tol:
            print 'reached tol!'
            return ph_l[-1]
    print 'exceded limit'
    return ph_l[-1]


def estimate_freqs(read_counts, haps, min_bound):
    
    exclude_haps = set()
    last_len = 0
    for rep in range(len(haps)):
        print rep, 'exluding %i haps but keeping %i' % (len(exclude_haps), len(haps-exclude_haps))
        sorted_uni_reads, sorted_uni_haps, urh, prh, ph = initialize(deepcopy(read_counts), haps-exclude_haps)
        nph = do_EM(urh, prh, ph)
        for p, hap in zip(nph.flatten(), sorted_uni_haps):
            if p < min_bound:
                exclude_haps.add(hap)
        if len(exclude_haps) == last_len:
            print 'done!'
            return sorted(zip(sorted_uni_haps, nph.flatten()), key = lambda x: x[1], reverse=True)
        last_len = len(exclude_haps)
    

# <codecell>

read_counts = get_read_counts(HapReconTools.read_from_bam(fname))
thaps = set(seq for seq, _ in ir_reads)
res = estimate_freqs(read_counts, thaps, estimate_min_bound(ir_reads, 7000, 9000)/100)
print res[0]


    


# <codecell>

from tempfile import NamedTemporaryFile as NTF
from subprocess import check_output, check_call
import shlex
import os
from concurrent.futures import ProcessPoolExecutor
import csv
from StringIO import StringIO
from itertools import islice
from functools import partial

def check_seqs(db_path, seqs):
    cmd = "blastn -db %(db)s -query %(q)s -outfmt '10 qseqid sseqid pident nident' -num_threads 20 -max_target_seqs 1"
    fields = ['SeqA', 'SeqB', 'pident', 'nident']
    dpath =  '/home/will/tmpstuf/haptest/tmpseqs/'
    
    with NTF(suffix='.fa', dir=dpath, delete=False) as check_handle:
        
        fasta_writer(check_handle, seqs)
        check_handle.flush()
        os.fsync(check_handle.fileno())
        
        tdict = {
                 'db':db_path,
                 'q':check_handle.name
                 }
        cmd_list = shlex.split(cmd % tdict)
        out = check_output(cmd_list)
        reader = csv.DictReader(StringIO(out), fieldnames=fields)
        return list(reader)
    

def yield_blocks(iterable, block_size):
    
    block = list(islice(iterable, block_size))
    while block:
        yield block
        block = list(islice(iterable, block_size))
    

def blast_all_v_all(seqsA, seqsB, block_size=20):
        
    dpath = '/home/will/tmpstuf/haptest/tmpseqs/'
    with NTF(suffix='.fa', dir=dpath, delete=False) as db_handle:
        fasta_writer(db_handle, seqsA)
        db_handle.flush()
        os.fsync(db_handle.fileno())
        
        cmd = 'makeblastdb -in %s -dbtype nucl' % db_handle.name
        cmd_list = shlex.split(cmd)
        check_call(cmd_list)
        
        align_func = partial(check_seqs, db_handle.name)
        check_iterable = islice(yield_blocks(iter(seqsB), 200), 20)
        with ProcessPoolExecutor(max_workers=5) as pool:
            res_iter = pool.map(align_func, check_iterable)
            for num, block in enumerate(res_iter):
                print num, len(block)
        
        
    

    
    




# <codecell>

blast_all_v_all(sA, sB)

# <codecell>

with open('/home/will/tmpstuf/haptest/DrexelMed.A0107.R02.fa') as handle:
    sA = list(fasta_reader(handle))
    
with open('/home/will/tmpstuf/haptest/DrexelMed.A0107.fa') as handle:
    sB = list(fasta_reader(handle))

# <codecell>

sA[:5]

# <codecell>


