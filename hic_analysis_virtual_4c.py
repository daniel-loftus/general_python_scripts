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

#%%

my_dir = r'W:\Users\dloftus\capture_hic\expt122_DMR_del_hic\usftp21.novogene.com\analysis\hicup\output\final_hicup_bams\snpsplit_output\merged_bams\final_files_for_analysis'
os.chdir(my_dir)

mat_file = r'DMR_DEL2_42-10_maternal.medium.sorted.cast_del_removed.peg13_locus.subsampled.txt'
pat_file = r'DMR_DEL2_42-10_paternal.medium.sorted.cast_del_removed.peg13_locus.txt'


mat = pd.read_csv(mat_file, sep = '\t', header = None)
pat = pd.read_csv(pat_file, sep = '\t', header = None)


mat = mat[mat[[3, 7]].notnull().all(1)]
pat = pat[pat[[3, 7]].notnull().all(1)]



#%%

def v4c_avg(hic_df, bin_res = 1, avg_rad = 5, anchor_start = 72547000, anchor_end = 72552000, color = 'black', ylim = False, title = None, xlabel = None, ylabel = None, plot = True, leg_label = None, export_svg = False, svg_name = None):

    '''
    Virtual 4C with a moving average 
    '''
    
    locus = [71850000, 73350000]
    locus_len = locus[1] - locus[0]
    
    anchor1 = hic_df[(hic_df[3] > anchor_start) & (hic_df[3] < anchor_end)]
    anchor2 = hic_df[(hic_df[7] > anchor_start) & (hic_df[7] < anchor_end)]
    
    anchor1_contacts = anchor1[7]
    anchor2_contacts = anchor2[3]
    
    anchor_contacts = anchor1_contacts.append(anchor2_contacts)
    
    hist_coords = np.linspace(locus[0] + avg_rad, locus[1] - avg_rad, round( (locus_len * 10**-3) / bin_res) - avg_rad * 2)
    
    bin_counts_moving_average = []
    for coord in hist_coords:
        
        contacts = [contact for contact in anchor_contacts if (contact > (coord - avg_rad * 10**3)) & (contact < (coord + avg_rad * 10**3))]
        
        bin_counts_moving_average.append(len(contacts))
            
    if plot:
        plt.plot(bin_counts_moving_average, color = color, label = leg_label)
        x = np.linspace(0, len(bin_counts_moving_average), len(bin_counts_moving_average)*10**-2 + 1)
        labels = np.linspace(locus[0] * 10**-6, locus[1] * 10**-6, locus_len/(10**5)).round(2)
        plt.xticks(x, labels, rotation = 45)
        plt.grid(False)
        if title:
            plt.title(title)
        if ylim:
            plt.ylim([0, ylim])
        if xlabel:
            plt.xlabel(xlabel)
        if ylabel:
            plt.ylabel(ylabel)

    results_df = pd.DataFrame()
    results_df['coord'] = hist_coords
    results_df['counts'] = bin_counts_moving_average

    if export_svg:
        plt.savefig(f"{svg_name}.svg")
    
    return results_df





#%%
    
#OVERLAY MAT AND PAT V4C PLOTS. ANCHOR = PEG13
   
plt.figure(figsize=(10,3))

v4c_avg(mat, color = 'red', bin_res = 5, avg_rad = 10, ylim = 600, anchor_start = 72805000, anchor_end = 72815000)
v4c_avg(pat, color = 'blue', bin_res = 5, avg_rad = 10, ylim = 600, anchor_start = 72805000, anchor_end = 72815000)

plt.savefig('v4C_dmr_anchor_del2.svg')

#%%
#OVERLAY MAT AND PAT V4C PLOTS. ANCHOR = Kcnk9
   
plt.figure(figsize=(10,3))

v4c_avg(mat, color = 'red', bin_res = 5, avg_rad = 10, ylim = 600)
v4c_avg(pat, color = 'blue', bin_res = 5, avg_rad = 10, ylim = 600)

plt.savefig('v4C_kcnk9_anchor_del2.svg')
#%%
#OVERLAY MAT AND PAT V4C PLOTS. ANCHOR = Maj enh

plt.figure(figsize=(10,3))
    
v4c_avg(mat, color = 'red', bin_res = 5, avg_rad = 10, anchor_start = 72846000, anchor_end = 72856000, ylim = 400)
v4c_avg(pat, color = 'blue', bin_res = 5, avg_rad = 10, anchor_start = 72846000, anchor_end = 72856000, ylim = 400)

plt.savefig('v4C_maj_enh_anchor_del2.svg')



#%%
#OVERLAY MAT AND PAT V4C PLOTS. ANCHOR = Tc9enh

plt.figure(figsize=(10,3))
    
v4c_avg(mat, color = 'red', bin_res = 5, avg_rad = 10, anchor_start = 73017000, anchor_end = 73029000, ylim = 850)
v4c_avg(pat, color = 'blue', bin_res = 5, avg_rad = 10, anchor_start = 73017000, anchor_end = 73029000, ylim = 850)

plt.savefig('v4C_tc9enh_anchor_del2.svg')