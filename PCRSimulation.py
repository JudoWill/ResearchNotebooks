# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=1>

# PCR Simulation

# <markdowncell>

# This is a simple program to estimate the number of errors that may have been introduced from our PCR procudure. Due to the exponential nature of PCR amplification we cannot use a simple 'item-based' simulation. Instead we need to use a little bit of math and statistics theory. It likely won't be accurate enough to win a Nobel prize but should be enough to give us a rough idea about the potential pitfalls of PCR errors.

# <markdowncell>

# First: Assume there are *n* identical integrated copies of the HIV-1 genome in the PCR reaction.
# 
# Second: Assume that the Polymerase Enzyme introduces errors at a constant rate: $$\mu$$
# 
# Third: Assume that the PCR reaction can be approximated by a doubling of strands at each iteration. Therefore we can represent the number of strands as:
# $$S(i) = n*2^i$$ where _i_ is the current iteration.
# 
# Combining 2 & 3 we can see that the number of mutations introduced at a given site at the *i* th step is:
# $$m(i) = \mu*S(i)$$
# 
# Therefore the number of mutations at a single nucleotide positions after _I_ steps is:
# $$M(I) = m(0) + m(1) + m(2) + ... m(I)$$
# $$M(I) = \sum_{i=0}^{i=I}\mu*S(i)$$
# 
# Finally, the probability of seeing a mutation at a single position after _I_ steps is:
# $$p(I) = M(I)/S(I)$$
# 
# 
# From this we can tell that the initial number of integrated copies *n* does not effect the calculation.

# <codecell>

from __future__ import division
import os, os.path
import numpy as np
from matplotlib import pyplot as plt



# <codecell>

def muts_per_cycle(error_rate, copies_generated):
    
    return error_rate*copies_generated

def muts_at_cycle(error_rate, initial_copies, cycle_num):
    
    muts = 0
    for num in range(1,cycle_num):
        muts += muts_per_cycle(error_rate, initial_copies*2**(num-1))
    return muts

def prob_at_cycle(error_rate, initial_copies, cycle_num):
    
    return muts_at_cycle(error_rate, initial_copies, cycle_num)/(initial_copies*2**(cycle_num))



# <codecell>

cycles = range(10,40,5)
num_initial = range(10,600,50)
error_rate = [10**(-i) for i in range(1,12)]

ER, CN = np.meshgrid(error_rate, cycles)

# <codecell>

print ER

# <codecell>

vectorized_prob = np.vectorize(prob_at_cycle)
cnum = 35
res = vectorized_prob(ER, 10, CN)

# <codecell>

plt.plot(np.log10(res))

# <codecell>


    



def linkfun(st):
    return st.amplify_with_error()

def batch_fun(error_rate, n_items):
    return (np.random.rand(n_items)<error_rate).sum()

class PCRChamber(object):
    
    def __init__(self, num_correct, template_len, mut_strands = []):
        
        self.num_correct_strands = num_correct
        self.template_len = template_len
        self.mut_strands = []
    
    
    def _reaction(self, error_rate, executor = None):
        
        
        batchsize = 10000
        rbatchsize = 500*batchsize
        if self.num_correct_strands < batchsize:
            for _ in range(self.num_correct_strands):
                if uniform(0,1) < error_rate:
                    yield Strand(self.template_len)
                else:
                    self.num_correct_strands += 1
        else:
            need = self.num_correct_strands
            batches = [rbatchsize]*int(need/rbatchsize)+[need%rbatchsize]
            nbatch_fun = partial(batch_fun, error_rate)
            if executor is not None:
                nmistake = sum(executor.map(nbatch_fun, batches))
            else:
                nmistake = sum(imap(nbatch_fun, batches))
            
            self.num_correct_strands += need - nmistake
            for n in range(nmistake):
                yield Strand(self.template_len)
            
            
            
        for strand in self.mut_strands:
            if uniform(0,1) < error_rate:
                yield strand
    
    
    def do_step(self, error_rate, executor = None):
        """Performs a single amplification step."""
        
        strands_with_errors = self._reaction(error_rate, executor=executor)
        
        nstrands = []
        if executor is not None:
            nstrands = executor.map(linkfun, strands_with_errors)
        else:
            nstrands = imap(linkfun, strands_with_errors)
                
        self.mut_strands += list(nstrands)
        
    def make_subset(self, fraction):
        """Creates a new PCR chamber object from a subset of this one."""
        
        num_correct = int(self.num_correct_strands*fraction)
        ostrands = []
        for strand in self.mut_strands:
            if uniform(0,1)<fraction:
                ostrands.append(strand)
        return PCRChamber(num_correct, self.template_len, ostrands)
    
    def find_prevelence(self, seq_len):
        
        counts = defaultdict(int)
        for strand in self.mut_strands:
            for pos in set(strand.muts):
                counts[pos] += 1
        
        fracs = np.array([counts[p] for p in range(seq_len)])/(len(self.mut_strands)+self.num_correct_strands)
        return fracs
            
    
class Strand(object):
    
    def __init__(self, seq_len, muts = []):
        self.muts = muts
        self.seq_len = seq_len
        
    def amplify_with_error(self):

        nmuts = self.muts + [choice(range(self.seq_len))]
        return Strand(self.seq_len, muts = nmuts)
    
    

        

# <codecell>

#orig_strands = [Strand(600) for _ in range(10)]
error_rate = 10e-8
chamber = PCRChamber(10, 600)

with ProcessPoolExecutor(max_workers = 30) as ex:
    for step in range(30):
        print 'doing step', step, chamber.num_correct_strands, len(chamber.mut_strands)
        chamber.do_step(error_rate, executor = ex)
    
    nchamber = chamber.make_subset(0.01)
    for step in range(30):
        print 'doing step', step, nchamber.num_correct_strands, len(nchamber.mut_strands)
        nchamber.do_step(error_rate, executor = ex)
        
    
res = nchamber.find_prevelence(600)
print res

# <codecell>

2**36


# <codecell>


