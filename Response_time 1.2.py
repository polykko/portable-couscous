# -*- coding: utf-8 -*-
"""
Created on Wed Jul 14 14:48:07 2021

@author: Polina
"""

"""
Pain so far:
    - independently on numer of eriods for averaging,
    average curve is super-noisy
    --> interpolation?? smoothing?
    - automatically get start&end points of averaged curve for fitting
   
"""
## handy to check the averaged curve
#plt.plot(list(avg.keys()), list(avg.values())) 
    
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
import scipy
from scipy.optimize import curve_fit
from scipy import signal
import random
# %matplotlib qt

def smoothTriangle(data, degree):
    triangle=np.concatenate((np.arange(degree + 1), np.arange(degree)[::-1])) # up then down
    smoothed=[]

    for i in range(degree, len(data) - degree * 2):
        point=data[i:i + len(triangle)] * triangle
        smoothed.append(np.sum(point)/np.sum(triangle))
    # Handle boundaries
    smoothed=[smoothed[0]]*int(degree + degree/2) + smoothed
    while len(smoothed) < len(data):
        smoothed.append(smoothed[-1])
    return smoothed

def read_file(input_file):
    PMT_trace = np.load(r"C:\Users\Polina\Desktop\ImPhys\Data + code\PMT_array_2020-08-06_11-34-20.npy",allow_pickle=True)[2:]
    pmt = np.load(r"C:\Users\Polina\Desktop\ImPhys\Data + code\Vp2020-08-06_11-33-01.npy",allow_pickle=True)[7:]
        
    noise = (np.where(PMT_trace < 0.7*np.average(PMT_trace)))[0]
    
    pmt = np.delete(pmt,noise)
    PMT_trace = np.delete(PMT_trace, noise)
    
    for ii in range(len(pmt)-1):
        if (pmt[ii]-pmt[ii+1]) < -0.5:
            start_point = ii
            period1 = ii
            for jj in range(ii+1,len(pmt)-1):
                if (pmt[jj]-pmt[jj+1]) < -0.5:
                    period2 = jj
                    break
            break
    #start_point = 1990 # to define
    #period_len = 4000 # to get from frequency
    period_len = period2-period1
    dataset_new = PMT_trace[start_point::]
    num_periods = len(dataset_new)//period_len
    
    ## now splitting
    periods = {}
    ## key is the index, value is a list of values
    for ii in range(num_periods):
        periods[ii] = dataset_new[ii*num_periods:ii*num_periods+period_len]
    
    return periods, start_point, period_len

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
        avg[ii] = 0
        for jj in range(len(periods_norm)):  # for each period
            #avg[ii] = np.average(periods_norm[jj][ii])
            avg[ii] = avg[ii] + periods_norm[jj][ii]
        avg[ii] = avg[ii] / len(periods_norm)
            
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
    
    ydata = np.array(list(avg.values())[1500:2500])
    xdata = np.arange(len(ydata))/20000
    
    popt, pcov = curve_fit(func, xdata, ydata, \
                           bounds=([0,0,-np.inf,0],[0.1,0.1,np.inf,1]), maxfev = 500000)
    
    return popt, pcov

periods, start_point, period_len = read_file('input_file.csv')   
periods_norm = normalization(periods)
avg = averaging(periods_norm)

## popt - array of parameters
## popcov - estimated covariance of popt
popt, pcov = fitting(avg) 

# time = (np.array(avg.keys())/1000)[35:135]
data = np.array(list(avg.values())[1500:2500])
time = np.arange(len(data))/20000
model = []

for ii in time:
    model.append(func(ii,popt[0],popt[1],popt[2],popt[3]))

plt.figure()
plt.plot(time,data,label='Real data')
plt.plot(time,model,label='Model')
plt.legend()
plt.grid()
plt.title('Fit of data')
plt.xlabel('Time (s)')
plt.ylabel('Intensity')
plt.text(0.1,1, 'Fit formula: F(t) = A × (C × exp(–t/t1) + (1 - C) × exp(–t/t2))\nt1 fast component (ms) = %s\nt2 slow component (ms) = %s\na = %s\nc = %s' % (round(popt[0]*1000,3),round(popt[1]*1000,3),round(popt[2],3),round(popt[3],3)),fontsize='large')
plt.show()

plt.figure()
avg_time = np.array(list(avg.keys()))/20000
plt.plot(avg_time, list(avg.values())) 
plt.title('Averaged normalized curve')
plt.xlabel('Time (s)')
plt.ylabel('Intensity')
plt.grid()
plt.show()

## calculating SNR
avg_val = list(avg.values())

## smoothing (choose method)
avg_val = scipy.signal.savgol_filter(avg_val,51,7) # apply a Savitzky-Golay filter
#avg_val = smoothTriangle(avg_val,5)


avg_mean = np.mean(avg_val)
avg_std = np.std(avg_val)
avg_SNR = avg_mean/avg_std

per = periods_norm[random.randint(0,len(periods_norm))] # to check SNR in a random period 
per_val = list(per)                              # over available ones
per_mean = np.mean(per_val)
per_std = np.std(per_val)
per_SNR = per_mean/per_std

N = avg_SNR/per_SNR
