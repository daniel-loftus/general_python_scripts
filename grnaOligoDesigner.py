# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 13:31:22 2020

@author: dal1858
"""

#%%

#The function grnaOligoDesigner() takes a gRNA sequence and outputs the required oligos to be used 
#for the cloning of the sequence into the backbone of pLV hU6-sgRNA hUbC-dCas9-KRAB-T2a-Puro (Addgene #71236) 

from revcomp import revcomp


def coding(grna, construct = 'hU6'):
    if construct == 'mU6':
        return 'TTGTTTG' + grna
    elif construct == 'hU6':
        return 'CACCG' + grna
    elif construct == '7SK':
        return 'CCTCG' + grna
    elif construct == 'H1':
        return 'TCCCA' + grna

def noncoding(grna, construct = 'hU6'):
    if construct == 'mU6':
        return 'AAAC' + revcomp(grna) + 'CAA'
    elif construct == 'hU6':
        return 'AAAC' + revcomp(grna) + 'C'
    elif construct == '7SK':
        return 'AAAC' + revcomp(grna) + 'C'
    elif construct == 'H1':
        return 'AAAC' + revcomp(grna) + 'T'

def grnaDesign(grna, construct = 'hU6'):
    coding_seq = coding(grna, construct)
    noncoding_seq = noncoding(grna, construct)
    print(f'Coding: {coding_seq}')
    print(f'Noncoding: {noncoding_seq}')