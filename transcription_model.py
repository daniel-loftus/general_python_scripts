# -*- coding: utf-8 -*-
"""
Created on Sun Feb 24 21:49:13 2019

@author: dal1858
"""
#%%
import numpy as np
import math as math
import matplotlib.pyplot as plt
import random as rand

t = list(range(100))
R = ['NA']*100
start = 70
ksyn = 10
decay = 0.2

#%%

for i in range(len(R)):
    if i == 0:
        R[i] = start
    else:
        R[i] = R[i - 1] + ksyn - decay*R[i - 1]

#%%%
plt.plot(R)

#%%

Rcal = ['NA']*len(t)
for i in range(len(t)):
    Rcal[i] = start*(ksyn/(decay*start) + math.exp(-decay*t[i]))
plt.plot(Rcal)
plt.plot(R)


#%%
kon = .2
koff = .3
start = 0
ksyn = 10
decay = 0.1

#%%
t = list(range(1000))
state = ['NA']*len(t)
state[0] = 0

R = ['NA']*len(t)

R[0] = start
for i in range(1, len(t)):
    if state[i - 1] == 0:
        R[i] = R[i - 1] - decay*R[i - 1]
        p = rand.randint(1, 100)
        if p > kon*100:
            state[i] = 0
        else:
            state[i] = 1
    if state[i - 1] == 1:
        R[i] = R[i - 1] - decay*R[i - 1] + ksyn
        p = rand.randint(1, 100)
        if p > koff*100:
            state[i] = 1
        else: 
            state[i] = 0
plt.plot(R)

#Clearly this is not all there is to it. It doesn't look like the actual data. 
#There has to be something controlling the burst size besides just the probability 
#of going into the off state. There is also probably a refractory period. 



#%%
kon = .01
koff = .05
start = 0
ksyn = 10
decay = 0.1

#%%
#Diploid
t = list(range(1000))
state = ['NA']*len(t)
state[0] = 0

R1 = ['NA']*len(t)

R1[0] = start
for i in range(1, len(t)):
    if state[i - 1] == 0:
        R1[i] = R1[i - 1] - decay*R1[i - 1]
        p = rand.randint(1, 100)
        if p > kon*100:
            state[i] = 0
        else:
            state[i] = 1
    if state[i - 1] == 1:
        R1[i] = R1[i - 1] - decay*R1[i - 1] + ksyn
        p = rand.randint(1, 100)
        if p > koff*100:
            state[i] = 1
        else: 
            state[i] = 0




state = ['NA']*len(t)
state[0] = 0

R2 = ['NA']*len(t)

R2[0] = start
for i in range(1, len(t)):
    if state[i - 1] == 0:
        R2[i] = R2[i - 1] - decay*R2[i - 1]
        p = rand.randint(1, 100)
        if p > kon*100:
            state[i] = 0
        else:
            state[i] = 1
    if state[i - 1] == 1:
        R2[i] = R2[i - 1] - decay*R2[i - 1] + ksyn
        p = rand.randint(1, 100)
        if p > koff*100:
            state[i] = 1
        else: 
            state[i] = 0

totalR = ['NA']*len(t)
for i in range(len(t)):
    totalR[i] = R1[i] + R2[i]
plt.plot(totalR)