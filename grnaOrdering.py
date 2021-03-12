# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 14:29:13 2020

@author: dal1858
"""

#%%

#This script outputs a list of oligo sequences from a list of gRNA sequences 

import grnaOligoDesigner as g
import pandas as pd
#%%

grnaFile = pd.read_csv("test.txt", sep="\t")

#%%
df = pd.DataFrame()
seqs = []
names = []
for i in range(len(grnaFile)):
    seqs.append(g.coding(grnaFile.iloc[i,1]))
    seqs.append(g.noncoding(grnaFile.iloc[i,1]))
    
    names.append(grnaFile.iloc[i,0])
    names.append(grnaFile.iloc[i,0] + "_comp")
    
df['Name'] = names
df['Seq'] = seqs