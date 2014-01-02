# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import sys
sys.path.append('/home/will/PySeqUtils/')
import os, os.path
os.chdir('/home/will/HIVTropism/R5Cluster/')
from itertools import islice
import pandas as pd
from GeneralSeqTools import fasta_reader, call_muscle

# <codecell>

lanl_data = pd.read_csv('results.csv', sep='\t')
lanl_data

# <codecell>

from subprocess import check_call
from GeneralSeqTools import call_muscle
import shutil
import shlex

def split_seqs(inseqs):
    num = int(len(inseqs)/2)
    return inseqs[num:], inseqs[:num]

def take(N, iterable):
    return list(islice(iterable, N))


def recursive_align(seqsA, seqsB, depth=0, max_seqs=1000):
    
    #print depth, len(seqsA), len(seqsB)
    if len(seqsA) > max_seqs:
        print 'recursingA', len(seqsA), depth
        seqsA = recursive_align(*split_seqs(seqsA), 
                                depth=depth+1, max_seqs=max_seqs)
    elif any(len(seqsA[0][1]) != len(s) for _, s in seqsA):
        print 'aligningA', len(seqsA), depth
        seqsA = call_muscle(seqsA)
    else:
        print 'FinishedA', depth
        
    if len(seqsB) > max_seqs:
        print 'recursingB', len(seqsB), depth
        seqsB = recursive_align(*split_seqs(seqsB), 
                                depth=depth+1, max_seqs=max_seqs)
    elif any(len(seqsB[0][1]) != len(s) for _, s in seqsB):
        print 'aligningB', len(seqsB), depth
        seqsB = call_muscle(seqsB)
    else:
        print 'FinishedB', depth
    
    if len(seqsA[0][1]) == len(seqsB[0][1]):
        print 'Easy join', depth
        return seqsA + seqsB
    else:
        print 'Hard join', depth
        return call_muscle(seqsA + seqsB)
        


def large_cluster(iseqs, cluster_size=1000):
    
    running_align = 'running_profile.fasta'
    ifile = 'adding_fasta.fasta'
    
    #do base alignment
    chunk_seqs = take(cluster_size, iseqs)
    muscle_align(chunk_seqs, running_align)
    
    count = len(chunk_seqs)
    print 'Seqs Processed:', count
    chunk_seqs = take(cluster_size, iseqs)
    while chunk_seqs:
        count += len(chunk_seqs)
        muscle_align(chunk_seqs, ifile)
        shutil.move(running_align, running_align + '.tmp')
        muscle_join(running_align + '.tmp', ifile, running_align)
        chunk_seqs = take(cluster_size, iseqs)
        print 'Seqs Processed:', count
        
    shutil.move(running_align, running_align + '.tmp')
    print 'refining!'
    cmd = 'muscle -in %s -out %s -refine' % (running_align + '.tmp', running_align)
    check_call(shlex.split(cmd))
    
    with open(running_align) as handle:
        return list(fasta_reader(handle))
    
    

# <codecell>

V3_aa = 'VQDPTTIQEKESVSREDQGEHLLQEKEIDKHI'
V3_nuc = 'TGTACAAGACCCAACAACAATACAAGAAAAAGAATCCGTATCCAGAGAGGACCAGGGAGAGCATTTGTTACAATAGGAAAAATAGGAAATATGAGACAAGCACATTGT'

# <codecell>

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from collections import Counter, defaultdict
from GeneralSeqTools import fasta_writer, seq_align_to_ref, WebPSSM_V3_fasta
from concurrent.futures import ThreadPoolExecutor
import csv
from itertools import imap


def filter_seq(handle, trans):
    for name, seq in fasta_reader(handle):
        tseq = ''.join(l for l in seq if l.isalpha())
        l = len(tseq)
        if (l > 100) and (l < 120):
            if trans:
                rseq = Seq(tseq, generic_dna).translate()
                yield name, ''.join(l for l in rseq.tostring() if l.isalpha())
            else:
                yield name, tseq
                
def yield_chunks(iterable, chunksize):

    chunk = take(chunksize, iterable)
    while chunk:
        yield chunk
        chunk = take(chunksize, iterable)

fields = ['name','score','pred','x4.pct','r5.pct','geno','pos.chg','net.chg','percentile']
        
with open('hiv-db.fasta') as handle:
    seq_iter = filter_seq(handle, True)
    chunk_iter = yield_chunks(seq_iter, 100)
    with open('LANLTropism.tsv', 'w') as ohandle:
        writer = csv.writer(ohandle, delimiter='\t')
        writer.writerow(fields)
        with ThreadPoolExecutor(max_workers=10) as E:
            for res in E.map(WebPSSM_V3_fasta, chunk_iter):
                writer.writerows(res)
        #with ThreadPoolExecutor(max_workers=10) as E:
        
#    for num, tup in enumerate(seq_iter):
#        pass
#    print 'Total seqs!:', num
                
                
#with open('hiv-db.fasta') as handle:
#    seq_iter = filter_seq(handle, True)
#    with open('V3filter.aa.fasta.raln', 'w') as ohandle:
#        align_iter = seq_align_to_ref(seq_iter, V3_aa, max_workers=20)
#        fasta_writer(ohandle, align_iter)
#print 'finished AA!'
        
#with open('hiv-db.fasta') as handle:
#    seq_iter = filter_seq(handle, False)
#    with open('V3filter.nt.fasta.raln', 'w') as ohandle:
#        align_iter = seq_align_to_ref(seq_iter, V3_nuc, max_workers=20)
#        fasta_writer(ohandle, align_iter)
#print 'finished NUC!'

# <codecell>


