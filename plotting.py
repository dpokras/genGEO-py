# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 15:48:40 2023

@author: pokra
"""

import CoolProp
from CoolProp.Plots import PropertyPlot

import pandas as pd
import matplotlib.pyplot as plt
HS1 = pd.read_csv('C:/Users/pokra/Documents/ETH/Master Thesis/genGEO_matlab_AGS/dpo_HS1_22_06_23.csv',header=None)
hS1 = pd.read_csv('C:/Users/pokra/Documents/ETH/Master Thesis/genGEO_matlab_AGS/dpo_hS1_22_06_23.csv',header=None)
hS2 = pd.read_csv('C:/Users/pokra/Documents/ETH/Master Thesis/genGEO_matlab_AGS/dpo_hS2_22_06_23.csv',header=None)
PS1 = pd.read_csv('C:/Users/pokra/Documents/ETH/Master Thesis/genGEO_matlab_AGS/dpo_PS1_22_06_23.csv',header=None)
PS2 = pd.read_csv('C:/Users/pokra/Documents/ETH/Master Thesis/genGEO_matlab_AGS/dpo_PS2_22_06_23.csv',header=None)

plot = PropertyPlot('CO2', 'ph')
ax = plot.axis
plt.ylim([2e3,1e5])
plt.xlim([150,550])
plot.calc_isolines()
ax.plot(hS1/1e3, PS1/1e3, 'black',linewidth = 0.4)
# ax.plot(hS2/1e3, PS2/1e3, 'black',linewidth = 0.4)

plot.show()
plot.savefig('C:/Users/pokra/Documents/ETH/Master Thesis/genGEO_matlab_AGS/Figures/Ph_diagram_config1.png', dpi=300)

