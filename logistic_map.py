# -*- coding: utf-8 -*-
"""
Created on Fri Mar 19 21:28:58 2021

@author: dal1858
"""

#%%

import numpy as np 
import matplotlib.pyplot as plt 
import random as rand

#%%

def trial(rate, gens_count, pop_start, plot = False):

    '''
    For gens_count generations, simulates the population over time starting 
    with pop_start. Each gen gen_n is a function of gen_n-1 as follows:
        
        gen_n = rate * gen_n-1 * (1 - gen_n-1)
        
    In addition to the (rate * gen_n-1) component, there is also the 
    constraining (1 - gen_n-1) term. Gen size is as a proportion of max 
    possible, and as such has range [0, 1]. 
    
    Rate has effective range [0, 4). 
    
    Output: 
        pop = list of populations sizes by gen
    '''
    gens = np.linspace(0, gens_count, gens_count + 1)
    
    pop = [pop_start]
    for gen in range(gens_count):
        
        pop_new = rate * pop[-1] * (1 - pop[-1])

        pop.append(pop_new)
        
        
    if plot:
        plt.plot(gens, pop)
        plt.xlabel('Generation')
        plt.ylabel('Population Size')
        plt.title('Population Size over Time')
        
    return pop
        
def pop_v_rate(rate_linspace, pop_start = 0.5):
    '''
    Plots the final population vs the growth rate for the given linspace 
    defined range of rates and a starting population. Generates the actual 
    logistic map. 
    '''
    
    final_pop = []
    
    for rate in rate_linspace:
        final = trial(rate = rate, gens_count = rand.randint(100, 1000), pop_start = pop_start, plot = False)[-1]
        if final > 0:
            final_pop.append(final)
        else:
            final_pop.append(0)
    
    plt.scatter(rate_linspace, final_pop, s = 1)
    plt.xlabel('Growth Rate')
    plt.ylabel('Population')
    plt.title('Logistic Map')

#pop_v_rate(rate_linspace = np.linspace(2.5, 3.9, 1000000))        
