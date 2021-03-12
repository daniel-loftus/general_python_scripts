# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 14:51:07 2020

@author: dal1858
"""


#%%


samples = 5

to_filter = ['N04', 'O22', 'E05', 'C02', 'B22', 'B24', 'E14']

file = "dloftus_2020-01-10 18-07-59_CT010444_baseline_adjust_original.txt"

probe_names = ['Peg13', 'Peg13_2', 'Kcnk9', 'Trappc9_5prime', 'Ago2']

gDNA_count = 2
gDNA_names = ['C6', 'D6']

#%%

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats

#%%

#The name of the tab-delimited text file from the BioRad 384 well qPCR machine. 
#Export 'Well', 'Fluor', and 'Cq'. 

df = pd.read_csv(file, sep="\t", skiprows=19)

#%%

#Filter wells 

for well in to_filter:
    for row in df.index:
        if well == df.loc[row, 'Well']:
            df.loc[row, 'Cq'] = 'N/A'
            
#%%
            
#make fam, vic lists

fam = []
vic = []
for row in df.index:
    if df.loc[row, 'Fluor'] == 'FAM':
        fam.append(df.loc[row, 'Cq'])
    elif df.loc[row, 'Fluor'] == 'VIC':
        vic.append(df.loc[row, 'Cq'])
        
#%%
        
#make dCq list

dCq = []
for row in range(len(fam)): 
    
    try: 
        dCq.append(float(fam[row]) - float(vic[row]))
    except ValueError: 
        dCq.append(None)

#%%
        
#make wells list

wells = []        
for row in df.index: 
    if df.loc[row, 'Well'] not in wells: 
        wells.append(df.loc[row, 'Well'])
        
#%%

#make master dataframe
        
data = pd.DataFrame()
data['Well'] = wells
data['FAM'] = fam
data['VIC'] = vic
data['dCq'] = dCq

#%%

#make a dictionary where each key is a probe and the method is the rows that correspond to that probe

probes = {}
row_groups = [['A', 'B', 'C'], ['D', 'E', 'F'], ['G', 'H', 'I'], ['J', 'K', 'L'], ['M', 'N', 'O']]
for i in range(samples):
    probes[probe_names[i]] = row_groups[i]

#%%

#standard curve. Generates graphs as well as lists containing the slopes and the intercepts 

logRatio = [3,2,1,0,-1,-2,-3]*3
slopes = []
intercepts = []
    
for probe_name, probe_rows in probes.items():
    dCq_temp = []
    for row in probe_rows:
        for k in range(3,16,2):
            if k < 10:
                dCq_temp.append(float(data.loc[data['Well'] == str(row + '0' + str(k)), 'dCq']))
            else:
                dCq_temp.append(float(data.loc[data['Well'] == str(row + str(k)), 'dCq']))               

    df_temp = pd.DataFrame()
    df_temp['dCq_temp'] = dCq_temp
    df_temp['logRatio'] = logRatio
    df_temp = df_temp[~df_temp['dCq_temp'].isnull()]
    slope, intercept, r_value, p_value, std_err = stats.linregress(df_temp.loc[:, 'dCq_temp'], df_temp.loc[:, 'logRatio'])    
    lineX = np.linspace(-6,6,2)
    lineY = slope*lineX + intercept
        
    fig = plt.figure(figsize = (5,20))
    sub = fig.add_subplot(samples, 1, 1)
    
    sub.scatter(dCq_temp, logRatio)
    sub.plot(lineX, lineY)
    
    text = '\n'.join(('R^2 = ' + format(r_value**2, '.3f'), 'Slope = ' + format(slope, '.3f'), 'Intercept = ' + format(intercept, '.3f')))
    textBox = dict(facecolor = 'white')
    sub.text(0.05, 0.95, text, transform=sub.transAxes, fontsize=10,
        verticalalignment='top', bbox = textBox)    
    plt.grid()
    plt.ylim(-4, 4)
    plt.xlim(-6,6)
    plt.title(probe_name)
    plt.xlabel("dCq")
    plt.ylabel("log2(M/P)")

    slopes.append(slope)
    intercepts.append(intercept)
    
#%%

#data analysis 
    
#first makes lists of lists with the probe sets so the indexes match the slope/intercept lists 

probes_list = []
for probe_rows in probes.values():
    probes_list.append(probe_rows)

log2mps = []
mps = []
for well in data.index: 
    
    for probe in range(len(probes_list)): 
        for row in probes_list[probe]: 
            if row in data.loc[well, 'Well']:
                log2mps.append(slopes[probe]*data.loc[well, 'dCq'] + intercepts[probe])

for log2mp in log2mps: 
    mps.append(log2mp**2)

for i in range(24):
    log2mps.append(None)
    mps.append(None)

data['Log2(M/P)'] = log2mps
data['M/P'] = mps

#%%















        
        
        
        
        
        
        