# -*- coding: utf-8 -*-
"""
Created on Thu Oct 13 15:25:41 2022

@author: dal1858
"""

#%%

'''
This little script will take a DNA sequence and generate the string needed to order
the oligo from IDT as an ASO, with the proper stability modifications.

NOTE: Only works on oligos of 20 NTs
'''

def aso_gen(dna):
        
    if len(dna) == 20:
        aso = '/52MOEr' + dna[0] + '/*/i2MOEr' + dna[1] + '/*/i2MOEr' + dna[2] \
            + '/*/i2MOEr' + dna[3] + '/*/i2MOEr' + dna[4] + '/*' + dna[5] + '*' + dna[6] \
            + '*' + dna[7] + '*' + dna[8] + '*' + dna[9] + '*' + dna[10] + '*' + dna[11] \
            + '*' + dna[12] + '*' + dna[13] + '*' + dna[14] + '*/i2MOEr' + dna[15] \
            + '/*/i2MOEr' + dna[16] + '/*/i2MOEr' + dna[17] + '/*/i2MOEr' + dna[18] \
            + '/*/32MOEr' + dna[19] + '/'
            
        return aso 
    
    else:
        print('The DNA sequence must be 20 NTs!!!')
        

