# -*- coding: utf-8 -*-
"""
Created on Fri Jul  9 10:05:06 2021

@author: Polina
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import csv

# %% Getting started

"""
improting values .csv from FIJI 
so far need to define file manually

"""

df = pd.read_csv('Test_untouched.csv')
intensities = df["Mean"]

# timepoints = df.values.tolist()


"""
user-defined stuff
will make fancy GUI
or extract from settings
"""

start_point = 600
period_length = 200
# from which timepoint we want to analyze our curves
# and how many ms is one cycle

data = intensities[start_point::]
num_periods = len(data)//period_length


# %% Splitting cycles and normalization
periods = {} # stores data of each period
periods_norm = {}
for ii in range(num_periods):
    periods[ii] = data[ii*num_periods:ii*num_periods+period_length]

""""
convert inside of each period onto dictionary as well
""""
# normalization
per = 0

for jj in periods:
    I_min = np.mean(jj[-15:-5])
    I_max = np.mean(jj[5:15])
    per = per+1
    point = 0
    for kk in jj:
        periods_norm[per][point] = (kk-I_min)/(I_max-I_min)
        point = point + 1
# %% Average curve
# for each timepoint within one period we are finding the average










        