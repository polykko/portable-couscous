# -*- coding: utf-8 -*-
"""
Created on Mon Jul 12 10:36:05 2021

@author: Polina
"""

"""
Pain so far:
    - there's a problem in either normalization or averaging
    (looks like in averaging)
    - determine I_min and I_max
        - can be thrown out of indexation
        - now always between 0 and 1, can be e.g. 0-0.6
        - is the period time correct?
    - determine start point
    - responce time calculation (where goes to plateau?)
    
"""
    
"""
1) import .csv (data type can be adjusted later) as numpy;
2) got period time (ot take it from frequency), start point:
beggining can be ugly due to the calibration, pulse, warming up or whatever;
3) calculated number of periods from data length and user-defined parameters;
4) split each period: dictionary of lists;
5) normalization of intensity values for each period;
7) averaging;
8) determine responce time (still to figure out how);
9) Voila!
"""

import matplotlib.pyplot as plt
import numpy as np


def read_file(input_file):
    dataset = np.loadtxt('input_file.csv', dtype=None, delimiter=',',skiprows=1)
    start_point = 600 # to define
    period_len = 200 # to get from frequency
    dataset_new = dataset[start_point::]
    num_periods = len(dataset_new)//period_len
    
    ## now splitting
    periods = {}
    ## key is the index, value is a list of values
    for ii in range(num_periods):
        periods[ii] = dataset_new[ii*num_periods:ii*num_periods+period_len,1]
    
    return periods

"""
now running over each period and normalizing data,
then take avg of each timepoint of each period
"""

def normalization(periods):
    #periods = read_file('input_file.csv')
    
    periods_norm = {}

    for jj in range(len(periods)): # jj is an index for period
        per = periods[jj] # taking one period number ...

## CAUTION: for some datasets min and max can be first or last values,
## so our I_min & I_max might give and error 
       
    ## "MANUAL" way - take suppesed min/max region
        #I_min = np.mean(per[-15:-5])
        #I_max = np.mean(per[5:15])
    ## "AUTOMATIC" way - taking min/max and around 
        I_min = np.mean(per[np.argmin(per)-1:np.argmin(per)+1])
        I_max = np.mean(per[np.argmax(per)-1:np.argmax(per)+1])
        
        periods_norm[jj] = (per - I_min)/(I_max-I_min)
           
    return periods_norm
    

## need to average each timepoint of all normalized periods
## so far ugly; didn't find a function to do all values ot once
## now summing up all periods' values for each timepoint
## and then dividing on number of periods

def averaging(periods_norm):
    
    average = {}
    
    for ii in range(len(periods_norm[0])):       
        average[ii] = 0
        for jj in range(len(periods_norm)):
            average[ii] = average[ii] + periods_norm[jj][ii]/len(periods_norm)
    # why cannot divide in previous loop - mystery    
    for kk in range(len(average)):
        average[kk] = average[kk]/len(periods_norm)
        
    #print(len(periods_norm)) ## just in case
    
    return average

periods = read_file('input_file.csv')   
periods_norm = normalization(periods)
average = averaging(periods_norm)

for_plotting = average.values()  
plt.plot(for_plotting)
plt.grid()    