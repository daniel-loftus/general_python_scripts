# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 10:34:25 2020

@author: dal1858
"""

#%%

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import math 
import statistics as stat

#%%
file = "loftus_fam_vic_taqman_redo_2020_10_28_expt40.csv"
df = pd.read_csv(file, sep=',', index_col=False)

#%%

#make dataframes with mat and pat data 

pat = []
mat = []
for i in df.index:
    if df.loc[i,'Target'] == "PAT":
        pat.append(df.loc[i,:])
    elif df.loc[i,'Target'] == "MAT":
        mat.append(df.loc[i,:])
        
mat = pd.DataFrame(mat)
mat = mat.reset_index(drop=True)
pat = pd.DataFrame(pat)
pat = pat.reset_index(drop=True)
#%%
        
#make master dataframe data that contains all required info 

data = pd.DataFrame() 
data['well'] = mat.loc[:, 'Well']
data['droplets'] = mat.loc[:, 'AcceptedDroplets']
for df_temp, string_temp in zip([mat, pat], ['mat', 'pat']):
    data[string_temp + '_conc'] = df_temp.loc[:, 'Concentration']
    data[string_temp + '_max'] = df_temp.loc[:, 'PoissonConfMax']
    data[string_temp + '_min'] = df_temp.loc[:, 'PoissonConfMin']
    data[string_temp + '_pos'] = df_temp.loc[:, 'Positives']
    data[string_temp + '_neg'] = df_temp.loc[:, 'Negatives']


#%%
    
#write a function to plot the mat and pat counts for the given wells 
    
def ddpcrCounts(inputFile, *argv):
    
    import pandas as pd
    import matplotlib.pyplot as plt
    import numpy as np
    
    
    df = pd.read_csv(inputFile, sep=',', index_col=False)
    #make dataframes with mat and pat data 
    
    pat = []
    mat = []
    for i in df.index:
        if df.loc[i,'Target'] == "MAT":
            mat.append(df.loc[i,:])
        elif df.loc[i,'Target'] == "PAT":
            pat.append(df.loc[i,:])
            
    mat = pd.DataFrame(mat)
    mat = mat.reset_index(drop=True)
    pat = pd.DataFrame(pat)
    pat = pat.reset_index(drop=True)

    #make master dataframe data that contains all required info 

    data = pd.DataFrame() 
    data['well'] = mat.loc[:, 'Well']
    data['droplets'] = mat.loc[:, 'AcceptedDroplets']
    for df_temp, string_temp in zip([mat, pat], ['mat', 'pat']):
        data[string_temp + '_conc'] = df_temp.loc[:, 'Concentration']
        data[string_temp + '_max'] = df_temp.loc[:, 'PoissonConfMax']
        data[string_temp + '_min'] = df_temp.loc[:, 'PoissonConfMin']
        data[string_temp + '_pos'] = df_temp.loc[:, 'Positives']
        data[string_temp + '_neg'] = df_temp.loc[:, 'Negatives']    
        
    mat_color = 'darkred'
    pat_color = 'mediumblue'
    
    mat_temp = []
    mat_max_temp = []
    mat_min_temp = []
    mat_nan_index = []

    pat_temp = []
    pat_max_temp = []
    pat_min_temp = []
    pat_nan_index = []
    
    #fill out lists containing the mat and pat counts to be graphed. Also sets up which samples were not called. 
    arg_ind = 0
    for arg in argv:
        
        for i in data.index: 
            if data.loc[i, 'well'] == arg:
                try: 
                    mat_temp.append(float(data.loc[i, 'mat_conc']))
                    mat_max_temp.append(float(data.loc[i, 'mat_max']) - float(data.loc[i, 'mat_conc']))
                    mat_min_temp.append(float(data.loc[i, 'mat_conc']) - float(data.loc[i, 'mat_min']))
                except ValueError:
                    mat_temp.append(0)
                    mat_max_temp.append(0)
                    mat_min_temp.append(0)
                    mat_nan_index.append(arg_ind)
                try: 
                    pat_temp.append(float(data.loc[i, 'pat_conc']))
                    pat_max_temp.append(float(data.loc[i, 'pat_max']) - float(data.loc[i, 'pat_conc']))
                    pat_min_temp.append(float(data.loc[i, 'pat_conc']) - float(data.loc[i, 'pat_min']))
                except ValueError: 
                    pat_temp.append(0)
                    pat_nan_index.append(arg_ind)
                    pat_max_temp.append(0)
                    pat_min_temp.append(0)
        arg_ind += 1
       
    ind = np.arange(len(argv))
    width = 0.35
    
    fig, ax = plt.subplots()
    ax.bar(ind, mat_temp, color = mat_color, width = width, label = 'Mat', yerr = [mat_min_temp, mat_max_temp])
    ax.bar(ind + width, pat_temp, color = pat_color, width = width, label = 'Pat', yerr = [pat_min_temp, pat_max_temp])
    ax.set_xticks(ind + width / 2)
    ax.set_xticklabels(argv)
    leg = ax.legend(loc=(1.02, 0.847))
    leg.get_frame().set_edgecolor('black')
    plt.ylabel("copies/ul")
    plt.title("Allelic counts")
    
    #add "No-call" to missing points
    for mat_nan in mat_nan_index: 
        ax.annotate('No-call', (mat_nan - 0.05, 20), rotation = 90, color = mat_color)
    for pat_nan in pat_nan_index:
        ax.annotate('No-call', (pat_nan + 0.3, 20), rotation = 90, color = pat_color)


#ddpcrCounts('loftus_fam_vic_taqman_2020_10_05_expt36.csv', 'E02', 'F02')
#%%

#outputs bar graphs showing mat and pat proportions 

def ddpcrRatios(*argv):

    mat_color = 'darkred'
    pat_color = 'mediumblue'
    
    mat_prop = []
    pat_prop = []
    nan_index = []
    
    #fill out lists containing the mat and pat counts to be graphed. Also sets up which samples were not called. 
    arg_ind = 0
    for arg in argv:
        
        for i in data.index: 
            try: 
                total = float(data.loc[i, 'mat_conc']) + float(data.loc[i, 'pat_conc'])
            except ValueError:
                total = 0
            if data.loc[i, 'well'] == arg:
                try: 
                    mat_prop.append(float(data.loc[i, 'mat_conc']) / total)
                    pat_prop.append(float(data.loc[i, 'pat_conc']) / total)
                except (ValueError, ZeroDivisionError):
                    mat_prop.append(0)
                    pat_prop.append(0)
                    nan_index.append(arg_ind)
        arg_ind += 1
       
    ind = np.arange(len(argv))
    width = 0.35
    
    fig, ax = plt.subplots()
    ax.bar(ind, mat_prop, color = mat_color, width = width, label = 'Mat')
    ax.bar(ind, pat_prop, color = pat_color, width = width, label = 'Pat', bottom = mat_prop)
    ax.set_xticks(ind)
    ax.set_xticklabels(argv)
    leg = ax.legend(loc=(1.02, 0.847))
    leg.get_frame().set_edgecolor('black')
    plt.ylabel("Mat prop")
    plt.title("Allelic proportions")
    
    #add "No-call" to missing points
    for nan in nan_index: 
        ax.annotate('No-call', (nan - 0.09, 0.07), rotation = 90, color = 'black')

#%%

#same as ddpcrCounts() but takes in kw arguments, where key = well and value = sample name 
    
def ddpcrCountskw(**kwargs):
    
    mat_color = 'darkred'
    pat_color = 'mediumblue'
    
    mat_temp = []
    mat_max_temp = []
    mat_min_temp = []
    mat_nan_index = []

    pat_temp = []
    pat_max_temp = []
    pat_min_temp = []
    pat_nan_index = []
    
    #fill out lists containing the mat and pat counts to be graphed. Also sets up which samples were not called. 
    arg_ind = 0
    for key in kwargs.keys():
        
        for i in data.index: 
            if data.loc[i, 'well'] == key:
                try: 
                    mat_temp.append(float(data.loc[i, 'mat_conc']))
                    mat_max_temp.append(float(data.loc[i, 'mat_max']) - float(data.loc[i, 'mat_conc']))
                    mat_min_temp.append(float(data.loc[i, 'mat_conc']) - float(data.loc[i, 'mat_min']))
                except ValueError:
                    mat_temp.append(0)
                    mat_max_temp.append(0)
                    mat_min_temp.append(0)
                    mat_nan_index.append(arg_ind)
                try: 
                    pat_temp.append(float(data.loc[i, 'pat_conc']))
                    pat_max_temp.append(float(data.loc[i, 'pat_max']) - float(data.loc[i, 'pat_conc']))
                    pat_min_temp.append(float(data.loc[i, 'pat_conc']) - float(data.loc[i, 'pat_min']))
                except ValueError: 
                    pat_temp.append(0)
                    pat_nan_index.append(arg_ind)
                    pat_max_temp.append(0)
                    pat_min_temp.append(0)
        arg_ind += 1
       
    ind = np.arange(len(kwargs))
    width = 0.35
    
    fig, ax = plt.subplots()
    ax.bar(ind, mat_temp, color = mat_color, width = width, label = 'Mat', yerr = [mat_min_temp, mat_max_temp])
    ax.bar(ind + width, pat_temp, color = pat_color, width = width, label = 'Pat', yerr = [pat_min_temp, pat_max_temp])
    ax.set_xticks(ind + width / 2)
    ax.set_xticklabels(kwargs.values(), rotation = 45, ha = 'right', rotation_mode = 'anchor')
    leg = ax.legend(loc=(1.02, 0.847))
    leg.get_frame().set_edgecolor('black')
    plt.ylabel("copies/ul")
    plt.title("Allelic counts")
    
    #add "No-call" to missing points
    for mat_nan in mat_nan_index: 
        ax.annotate('No-call', (mat_nan - 0.11, 20), rotation = 90, color = mat_color)
    for pat_nan in pat_nan_index:
        ax.annotate('No-call', (pat_nan + 0.3, 20), rotation = 90, color = pat_color)
        
#ddpcrRatioskw(G03 = 'Peg13_2 C6 gDNA', H03 = 'Peg13_2 D6 gDNA', G04 = 'Kcnk9 C6 gDNA', H04 = 'Kcnk9 D6 gDNA', G05 = 'Trappc9_5prime C6 gDNA', H05 = 'Trappc9_5prime D6 gDNA', G06 = 'Ago2 C6 gDNA', H06 = 'Ago2 D6 gDNA')


#%%

def ddpcrRatioskw(**kwargs):

    mat_color = 'darkred'
    pat_color = 'mediumblue'
    
    mat_prop = []
    pat_prop = []
    nan_index = []
    
    #fill out lists containing the mat and pat counts to be graphed. Also sets up which samples were not called. 
    arg_ind = 0
    for key in kwargs.keys():
        
        for i in data.index: 
            try: 
                total = float(data.loc[i, 'mat_conc']) + float(data.loc[i, 'pat_conc'])
            except ValueError:
                total = 0
            if data.loc[i, 'well'] == key:
                try: 
                    mat_prop.append(float(data.loc[i, 'mat_conc']) / total)
                    pat_prop.append(float(data.loc[i, 'pat_conc']) / total)
                except (ValueError, ZeroDivisionError):
                    mat_prop.append(0)
                    pat_prop.append(0)
                    nan_index.append(arg_ind)
        arg_ind += 1
       
    ind = np.arange(len(kwargs))
    width = 0.35
    
    fig, ax = plt.subplots()
    ax.bar(ind, mat_prop, color = mat_color, width = width, label = 'Mat')
    ax.bar(ind, pat_prop, color = pat_color, width = width, label = 'Pat', bottom = mat_prop)
    ax.set_xticks(ind)
    ax.set_xticklabels(kwargs.values(), rotation = 45, ha = 'right', rotation_mode = 'anchor')
    leg = ax.legend(loc=(1.02, 0.847))
    leg.get_frame().set_edgecolor('black')
    plt.ylabel("Mat prop")
    plt.title("Allelic proportions")
    
    #add "No-call" to missing points
    for nan in nan_index: 
        ax.annotate('No-call', (nan - 0.09, 0.07), rotation = 90, color = 'black')





#ddpcrCountskw(G03 = 'Peg13_2 C6 gDNA', H03 = 'Peg13_2 D6 gDNA', G04 = 'Kcnk9 C6 gDNA', H04 = 'Kcnk9 D6 gDNA', G05 = 'Trappc9_5prime C6 gDNA', H05 = 'Trappc9_5prime D6 gDNA', G06 = 'Ago2 C6 gDNA', H06 = 'Ago2 D6 gDNA')




#%%


#FUNCTION: PLOTS THE LOG2 RATIO BETWEEN THE TWO ALLELES. ARGUMENTS ARE WELL NAMES IN ALL CAPS (ex: 'A01', 'A02', A03')
        
def allelicLog(*argv):
    
    mat_color = 'darkred'
    pat_color = 'mediumblue'
    equal_color = 'black'
    
    foldChange = []
    nan_index = []
    
    #populate the list 'ratio' with the natural log ratios for each requested well from argv 
    for arg in argv:
        
        for i in data.index:
            if data.loc[i, 'well'] == arg:
                try: 
                    tempValue = math.log(float(data.loc[i, 'mat_conc']) / float(data.loc[i, 'pat_conc']))
                except ValueError:
                    tempValue = 'NaN'
                    nan_index.append(i)
                foldChange.append(tempValue)
    
    
    #remove the no calls 

    ind = np.arange(len(argv))

    width = 0.35
    
    fig, ax = plt.subplots()
    
    
    my_colors = []
    for fc in foldChange:
        if fc > 0:
            my_colors.append(mat_color)
        elif fc < 0:
            my_colors.append(pat_color)
        elif fc == 0:
            my_colors.append(equal_color)
    
    ax.bar(ind, foldChange, color = my_colors, width = width)
    ax.set_xticks(ind)
    ax.set_xticklabels(argv, rotation = 45, ha = 'right', rotation_mode = 'anchor')
    plt.ylabel("Log(Mat/Pat)")
    #plt.title("Allelic proportions")
      
allelicLog('A03', 'B03', 'C03', 'A01', 'B01', 'C01')

#%%

#FUNCTION: TAKES MULTIPLE BIOLOGICAL REPLICATES AND OUTPUTS THE MEAN LOG2(FC) + SEM 
#SYNTAX: log2FC(Peg13 = ['A01', 'B01', 'C01'])

def log2FC(**kwargs):
    
    mat_color = 'darkred'
    pat_color = 'mediumblue'
    equal_color = 'black'
    
    meanFC = []
    semFC = []
    
    nanIndex = []
    
    for wellsList in kwargs.values():
        
        sampleFC = []
        nanCount = 0
        for well in wellsList: 
            for i in data.index:
                if data.loc[i, 'well'] == well:
                    try: 
                        wellFC = math.log2(float(data.loc[i, 'mat_conc']) / float(data.loc[i, 'pat_conc']))
                    except ValueError:
                        wellFC = float('nan')
                        nanIndex.append(i)
                        nanCount += 1
                    sampleFC.append(wellFC)
        
        #calculate mean FC 
        sampleMeanFC = np.nansum(sampleFC)/(len(sampleFC) - nanCount)
        meanFC.append(sampleMeanFC)
        
        #calculate FC SEM 
        try:
            sampleSemFC = stat.stdev(sampleFC)/math.sqrt(len(sampleFC) - nanCount)
        except TypeError:
            sampleSemFC = float('nan')
        semFC.append(sampleSemFC)
        
                
    ind = np.arange(len(kwargs.keys()))
    barWidth = 0.75
    plotWidth = barWidth * 4
    
    fig, ax = plt.subplots(figsize = (plotWidth, 5))
    
    my_colors = []
    for fc in meanFC:
        if fc > 0:
            my_colors.append(mat_color)
        elif fc < 0:
            my_colors.append(pat_color)
        elif fc == 0:
            my_colors.append(equal_color)
    
    ax.bar(ind, meanFC, color = my_colors, width = barWidth, yerr = semFC, capsize = 5)
    ax.set_xticks(ind)
    ax.set_xticklabels(kwargs.keys(), rotation = 45, ha = 'right', rotation_mode = 'anchor')
    ax.set_yticks(np.arange(-10, 10, 1))
    plt.ylabel("Log2(Mat/Pat)")
    plt.title("Log2 Ratio of Allelic Counts")
    plt.ylim(-6, 6)
    plt.style.use('seaborn')
        
log2FC(Peg13_Tc9enh_KD = ['A01', 'B01', 'C01'], Peg13_DMR_KD = ['D01', 'E01', 'F01'], Kcnk9_Tc9enh_KD = ['A03', 'B03', 'C03'], Kcnk9_DMR_KD = ['D03', 'E03', 'F03'], Kcnk9_CTRL1 = ['A04', 'B04', 'C04'])
    
    
    
    
    
    
    