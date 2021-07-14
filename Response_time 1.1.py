# -*- coding: utf-8 -*-
"""
Created on Tue Jul 13 14:10:06 2021

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
from scipy.optimize import curve_fit


def read_file(input_file):
    dataset = np.loadtxt('input_file.csv', dtype=None, delimiter=',',skiprows=1)
    start_point = 599 # to define
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
    

## averages each timepoint of all normalized periods

def averaging(periods_norm):
    
    avg = {}
    
    for ii in range(len(periods_norm[0])):      # so over each timepoint 
        for jj in range(len(periods_norm)):
            avg[ii] = np.average(periods_norm[jj][ii])
    return avg


# =============================================================================

#     F(t) = A × (C × exp(–t/t1) + (1 - C) × exp(–t/t2)), where t1 was the time constant

#     of the fast component and t2 was the time constant of the slow component. The

#     percentage of the total magnitude that was associated with the fast component (%t1)

#     was defined as C above.

# =============================================================================

def func(t, t1, t2, a, c):
    return a * (c * np.exp(-t / t1) + (1 - c) * np.exp(-t / t2))

        
def fitting(avg): 
    
    #xdata = (np.array(avg.keys())/1000)[35:135]
    
    ydata = np.array(list(avg.values())[35:135])
    xdata = np.arange(len(ydata))/1000
    
    popt, pcov = curve_fit(func, xdata, ydata, \
                           bounds=([0,0,-np.inf,0],[0.1,0.1,np.inf,1]), maxfev = 500000)
    
    return popt, pcov

periods = read_file('input_file.csv')   
periods_norm = normalization(periods)
avg = averaging(periods_norm)

## popt - array of parameters
## popcov - estimated covariance of popt
popt, pcov = fitting(avg) 

# time = (np.array(avg.keys())/1000)[35:135]
data = np.array(list(avg.values())[35:135])
time = np.arange(len(data))/1000
model = []

for ii in time:
    model.append(func(ii,popt[0],popt[1],popt[2],popt[3]))

plt.figure()
plt.plot(time,data,label='Real data')
plt.plot(time,model,label='Model')
plt.legend()
plt.grid()
plt.show()
#for_plotting = avg.values()  
#plt.plot(for_plotting)
#plt.grid()