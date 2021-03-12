# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 12:43:55 2020

@author: dal1858
"""

#%%

#The function revcomp() takes in a DNA or RNA sequence and outputs its reverse complement DNA strand 

def revcomp(seq):
    revcompseq = ""
    for i in range(len(seq)):
        if (seq[len(seq) - 1 - i] == 'A') or (seq[len(seq) - 1 - i] == 'a'):
            revcompseq += 'T'
        elif (seq[len(seq) - 1 - i] == 'T') or (seq[len(seq) - 1 - i] == 't') or (seq[len(seq) - 1 - i] == 'U') or (seq[len(seq) - 1 - i] == 'u'): 
            revcompseq += 'A'
        elif (seq[len(seq) - 1 - i] == 'G') or (seq[len(seq) - 1 - i] == 'g'):
            revcompseq += 'C'
        elif (seq[len(seq) - 1 - i] == 'C') or (seq[len(seq) - 1 - i] == 'C'):
            revcompseq += 'G'
        else:
            print("Input sequence must only contain A,G,T,C,U")
            break
    
    if len(revcompseq) == len(seq):        
        return(revcompseq)
            
