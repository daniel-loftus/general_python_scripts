# -*- coding: utf-8 -*-
"""
Created on Mon Jan 31 14:26:41 2022

@author: dal1858
"""

#%%

import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
import os 

import math

#%%

my_dir = r'W:\Users\dloftus\capture_hic\expt122_DMR_del_hic\usftp21.novogene.com\analysis\hicup\output\final_hicup_bams\snpsplit_output\merged_bams\final_files_for_analysis'
os.chdir(my_dir)

mat_file = r'DMR_DEL1_42-4_maternal.medium.sorted.cast_del_removed.peg13_locus.subsampled.txt'
pat_file = r'DMR_DEL1_42-4_paternal.medium.sorted.cast_del_removed.peg13_locus.txt'


mat = pd.read_csv(mat_file, sep = '\t', header = None)
pat = pd.read_csv(pat_file, sep = '\t', header = None)

       
#%%
    
#Actual insulation score. Let's work with input_array (10x10)

cell_type = 'DEL'

def insulation_score(input_df, resolution = 10, plot_label = False, plot_color = False, color = 'black'):    
    
    input_array = np.histogram2d(input_df.iloc[:, 3], input_df.iloc[:, 7], bins = 144)[0]
    
    res = resolution
    
    insulation_score = []
    for index, row in enumerate(input_array):
        
        if index <= (len(input_array) - 2 * res):
            
            window = input_array[index:(index + res), (index + res):(index + 2 * res)]
            window_sum = sum(map(sum, window))
            
            #insulation_score.append(math.log2(window_sum))
            insulation_score.append(window_sum)
    
    plt.plot(insulation_score, label = plot_label, color = color)
    plt.legend()
    
    return insulation_score

res_list = [50, 75, 100, 150, 200, 250, 300]

for res in res_list:
    
    my_res = res
    
    mat_ins = insulation_score(mat, plot_label = 'Maternal', color = 'darkred', resolution = int(my_res/10))
    pat_ins = insulation_score(pat, plot_label = 'Paternal', color = 'darkblue', resolution = int(my_res/10))
    plt.title(f'{my_res} kb Sliding Window, {cell_type}')
    plt.savefig(f"insulation_{my_res}_kb_{cell_type}.svg")
    plt.clf()


#%%

log2 = []
for mat_score, pat_score in zip(mat_ins, pat_ins):
    
    log2.append(math.log2(mat_score / pat_score))
    
plt.plot(log2)