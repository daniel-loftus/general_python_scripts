#!/usr/bin/env python
# -*- coding: utf-8 -*-

#%%

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import random
import bisect

#%%

#import a bed file containing the 129/cast snps, in mm10
my_file = '/n/whipple_lab/Users/dloftus/sorted.129CastTrack_mm10.bed'
bed = pd.read_csv(my_file, header = None, sep = '\t')

#%%

#get an array of chromosome names 
chroms = np.unique(np.array(bed.iloc[:,0]))

#%%

#construct dictionary of chromosome sizes. Format is {chrom: [min, max]}

sizes = {}

for chrom in chroms:
    
    chrom_snps = bed[bed[0] == chrom]
    
    sizes[chrom] = [(np.min(chrom_snps[1])), (np.max(chrom_snps[1]))]
    

#%%
    

def allelic_read_prob(reps, read_len):
    
    allelic_hits = 0
    
    for i in range(reps):
        
        chrom = random.choice(chroms)
        chrom_min = sizes[chrom][0]
        chrom_max = sizes[chrom][1]
        
        coord = random.randrange(chrom_min, chrom_max)
        
        read_start = coord - read_len/2
        read_end = coord + read_len/2 
        
        chrom_bool = bed[0].isin([chrom])
        chrom_snps = bed[chrom_bool]
        coords = np.array(chrom_snps[1])
        
        lower_bound_i = bisect.bisect_left(coords, read_start)
        upper_bound_i = bisect.bisect_right(coords, read_end, lo = lower_bound_i)
        
        snps = coords[lower_bound_i:upper_bound_i]

        if len(snps) > 0:
            allelic_hits += 1
        
    prob = allelic_hits / reps
    return(prob)
#%%
probs = []
my_range = range(20, 151, 2)
for length in my_range:
    
    temp_prob = allelic_read_prob(1000, length)
    probs.append(temp_prob)
    print(temp_prob)
    
plt.plot(my_range, probs)


#%%

def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'valid') / w    

prob_results = pd.read_csv('W:\\Users\\dloftus\\snps_allelic_prob_vs_read_len_analysis\\probs_range_20_151_reps_1000.txt' , \
                           header = None)
prob_results = np.array(prob_results[0])

#%%
plt.scatter(my_range, prob_results, color = 'black', s = 10)
#plt.plot(moving_average(prob_results, 10))

plt.title('Read Length vs Prob of Crossing SNP')
plt.xlabel('Read Length (bp)')
plt.ylabel('Measured Prob of Crossing SNP')