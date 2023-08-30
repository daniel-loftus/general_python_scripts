# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 11:46:12 2020

@author: dal1858
"""

#%%

import matplotlib.pyplot as plt
import scipy.stats as stats
import numpy as np

#%%

psi = [783, 772, 760, 713, 702, 690, 678, 643, 620, 585, 573, 432]
day = [7,8,9,14,15,16,17,20,22,23, 24, 28]

#%%

slope, intercept, r_value, p_value, std_err = stats.linregress(day[:], psi[:])    
lineX = np.linspace(6,30,2)
lineY = slope*lineX + intercept

fig = plt.figure(figsize = (7, 10))
sub = fig.add_subplot(2, 1, 1)
plt.scatter(day, psi)
plt.plot(lineX, lineY)
plt.xlabel("Day (Jan)")
plt.ylabel("PSI")
plt.title("CO2 PSI over time (line omits outlier)")

text = '\n'.join(('R^2 = ' + format(r_value**2, '.3f'), 'Slope = ' + format(slope, '.3f'), 'Intercept = ' + format(intercept, '.3f')))
textBox = dict(facecolor = 'white')
plt.text(0.71, 0.95, text, transform=sub.transAxes, fontsize=10,
verticalalignment='top', bbox = textBox)    
plt.grid()