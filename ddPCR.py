# -*- coding: utf-8 -*-
"""
Created on Tue Mar 17 15:38:17 2020

@author: dal1858
"""

#%%

#FUNCTION: TAKES MULTIPLE BIOLOGICAL REPLICATES AND OUTPUTS THE MEAN LOG2(FC) + SEM 
#SYNTAX: log2FC(inputFile, Peg13 = ['A01', 'B01', 'C01'])

def log2FC(inputFile, master_title = None, x_label = None, y_limit = 4, export_svg = False, svg_name = None, **kwargs):
    
    
    import pandas as pd
    import matplotlib.pyplot as plt
    import numpy as np
    import math 
    import statistics as stat
    
    
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
        except ZeroDivisionError:
            sampleSemFC = float('nan')
        except stat.StatisticsError:
            sampleSemFC = float('nan')
        semFC.append(sampleSemFC)
        
                
    ind = np.arange(len(kwargs.keys()))
    barWidth = 0.8
    plotWidth = barWidth * 4
    
    fig, ax = plt.subplots(figsize = (plotWidth, 4))
    
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
    plt.ylim(-y_limit, y_limit)
    
    plt.style.use('seaborn-whitegrid')
    
    if master_title:
        plt.title(master_title)
    else:
        plt.title("Log2 Ratio of Allelic Counts")
    
    if x_label:
        plt.xlabel(x_label)
        
    if export_svg:
        plt.savefig(f"{svg_name}.svg")
        
        
        
        
        
        
        
def counts(inputFile, master_title = None, x_label = None, **kwargs):
    
    import pandas as pd
    import matplotlib.pyplot as plt
    import numpy as np
    import math 
    import statistics as stat
    
    
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
    if master_title:
        plt.title(master_title)
    else:
        plt.title("Allelic counts")
    
    #add "No-call" to missing points
    for mat_nan in mat_nan_index: 
        ax.annotate('No-call', (mat_nan - 0.11, 20), rotation = 90, color = mat_color)
    for pat_nan in pat_nan_index:
        ax.annotate('No-call', (pat_nan + 0.3, 20), rotation = 90, color = pat_color)
        