# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import os, os.path
import csv
import sys
from subprocess import check_call, check_output
import shlex
from StringIO import StringIO
from itertools import islice
os.chdir('/home/will/VISPA/exec01/')
datadir = '/home/will/DeepSequencing/testfix/'
sys.path.append('/home/will/PySeqUtils/')

import QuasiTools
import GeneralSeqTools

# <codecell>

ref_seq = datadir + 'RefSequence_single.fasta'
reads = datadir + 'sim_reads.fa'

# <codecell>

def add_seq_to_reads(read_file, out_file):
    with open(read_file) as handle:
        with open(out_file, 'w') as ohandle:
            new_seqs = []
            for name, seq in GeneralSeqTools.fasta_reader(handle):
                new_seqs.append((name + ';' + seq, seq))
            GeneralSeqTools.fasta_writer(ohandle, new_seqs)

# <codecell>

def segemehl_align(ref_seq_file, treated_reads_file, threads = 10):
    cmd = 'segemehl.x -x %(nref_seq)s -d %(ref_seq)s -q %(reads)s -t %(threads)i -K'
    tdict = {
        'ref_seq':ref_seq_file,
        'reads':treated_reads_file,
        'threads':10,
        'nref_seq':ref_seq_file+'.idx'
        }
    ncmd = cmd % tdict
    return check_output(shlex.split(ncmd))

# <codecell>

def extract_align_cols(segemehl_handle, out_handle):
    writer = csv.writer(out_handle, delimiter = '\t')
    for row in csv.reader(segemehl_handle, delimiter = '\t'):
        try:
            writer.writerow((row[1], row[9], row[10], row[11], row[13]))
        except IndexError:
            pass

# <codecell>

def maxi_awk(trimmed_align_file, ref_seq_file, maxi_file):
    cmd = 'java MaxI_awk %(trimmed)s %(maxi)s %(ref_file)s'
    tdict = {
             'trimmed':trimmed_align_file,
             'maxi':maxi_file,
             'ref_file':ref_seq_file
        }
    ncmd = cmd % tdict
    check_call(shlex.split(ncmd))
    
def extend_reference(trimmed_align_file, ref_seq_file, ref_extended_file):
    cmd = 'java CreateReferenceWithI_awk %(trimmed)s %(ref_extended)s %(ref_file)s'
    tdict = {
             'trimmed':trimmed_align_file,
             'ref_extended':ref_extended_file,
             'ref_file':ref_seq_file
        }
    
    ncmd = cmd % tdict
    check_call(shlex.split(ncmd))
    
    
def ConvertSegemehlIraExtended(treated_reads_file, trimmed_align_file, reads_extended, ref_seq_file, ref_maxi_file):
    
    cmd = 'java ConvertSegemehlIraExtended_awk %(treated_reads)s %(trimmed)s %(reads_extended)s %(ref_file)s %(ref_maxi)s'
    
    tdict = {
                'treated_reads':treated_reads_file,
                'trimmed':trimmed_align_file,
                'reads_extended':reads_extended,
                'ref_file':ref_seq_file,
                'ref_maxi':ref_maxi_file
        }
    
    ncmd = cmd % tdict
    check_call(shlex.split(ncmd))
    

# <codecell>

def vispa_assemble(extended_reads_file, assembeled_reads_file, extended_ref_file, nmuts = 2, qsps_muts = 10):
    cmd = 'java -Xms7020m -Xmx7020m -jar vispaAssemb.jar %(ex_reads)s %(ass_reads)s %(extended_ref)s %(muts)i %(qsps_muts)s'
    tdict = {
                'ex_reads':extended_reads_file,
                'ass_reads':assembeled_reads_file,
                'extended_ref':extended_ref_file,
                'muts':nmuts,
                'qsps_muts':qsps_muts,
        }
    ncmd = cmd % tdict
    #print ncmd
    check_call(shlex.split(ncmd))
    out_filename = assembeled_reads_file.rsplit('.',1)[0] + '_I_%i_%i_CNTGS_DIST0.txt' % (nmuts, qsps_muts)
    return out_filename
    
def vispaEM(assembeled_reads_file, cntg_file):
    
    cmd = 'java -Xms7020m -Xmx7020m -jar vispaEM.jar %(ass_reads)s %(cntg_file)s'
    tdict = {
             'ass_reads':assembeled_reads_file,
             'cntg_file':cntg_file
        }
    ncmd = cmd % tdict
    print ncmd
    check_call(shlex.split(ncmd))

# <codecell>

def do_all(ref_seq_file, reads_file):
    
    
    treated_reads_file = reads+'.treated'
    extended_reads_file = treated_reads_file + '.extended'
    align_cols_file = datadir + 'segemehl_trimmed_out.tsv'
    ref_maxi_file = ref_seq_file + '.maxi'
    extended_ref_file = ref_seq_file + '.extended'
    ass_reads = reads + '.ass'
    
    add_seq_to_reads(reads_file, treated_reads_file)
    
    out_align_cols = segemehl_align(ref_seq_file, treated_reads_file)
    with open(align_cols_file, 'w') as handle:
        extract_align_cols(StringIO(out_align_cols), handle)
    
    
    maxi_awk(align_cols_file, ref_seq_file, ref_maxi_file)
    extend_reference(align_cols_file, ref_seq_file, extended_ref_file)
    
    ConvertSegemehlIraExtended(treated_reads_file, align_cols_file, extended_reads_file, ref_seq_file, ref_maxi_file)
    
    cntg_file = vispa_assemble(extended_reads_file, ass_reads, extended_ref_file, nmuts = 2, qsps_muts = 10)
    vispaEM(ass_reads, cntg_file)
    
    
    

# <codecell>

do_all(ref_seq, reads)

# <codecell>

java -Xms7020m -Xmx7020m -jar vispaEM.jar /home/will/DeepSequencing/testfix/sim_reads.fa.ass /home/will/DeepSequencing/testfix/final_output.txt

# <codecell>


