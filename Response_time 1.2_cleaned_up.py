# -*- coding: utf-8 -*-
"""
Created on Wed Jul 14 14:48:07 2021

@author: Polina
"""

"""
Pain so far:
    - automatically get start & end points of averaged curve for fitting
   
"""
## handy to check the averaged curve
#plt.plot(list(avg.keys()), list(avg.values())) 
    
"""
1) import data as a list;
2) got period time (ot take it from frequency), start point:
3) calculate number of periods from data;
4) split each period and store it as a dictionary of lists;
5) normalize intensity values for each period;
7) averaging;
8) determine response time by curve fitting;
9) Voila!
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import random
# %matplotlib qt

sampling_rate = 20000

# =============================================================================
#   Triangle smoothing function
# =============================================================================
def smoothTriangle(data, degree):
    triangle=np.concatenate((np.arange(degree + 1), np.arange(degree)[::-1]))
    smoothed=[]

    for i in range(degree, len(data) - degree * 2):
        point=data[i:i + len(triangle)] * triangle
        smoothed.append(np.sum(point)/np.sum(triangle))
    # Handle boundaries
    smoothed=[smoothed[0]]*int(degree + degree/2) + smoothed
    while len(smoothed) < len(data):
        smoothed.append(smoothed[-1])
    return smoothed


# =============================================================================
#   Normalization of values in the list
# =============================================================================
def normalize(input_list):
    norma = []
    min_v = min(input_list)
    max_v = max(input_list)
#   ii is an index of value in the list
    for ii in range(len(input_list)):
        norma.append((input_list[ii]-min_v)/(max_v-min_v))
    return norma


# =============================================================================
#   Importing data, excluding outliers (variable "noise")
#   Normalizing voltage signal to determine the starting point - time of start
#   of the rising 
# =============================================================================
def read_file(input_file):
    PMT_trace = np.load(r"C:\Users\Polina\Desktop\ImPhys\Data + code\PMT_array_2020-08-06_11-34-20.npy",allow_pickle=True)[2:]
    voltage_data = np.load(r"C:\Users\Polina\Desktop\ImPhys\Data + code\Vp2020-08-06_11-33-01.npy",allow_pickle=True)[7:]
        
    noise = (np.where(PMT_trace < 0.7*np.average(PMT_trace)))[0]
    
    voltage_data = np.delete(voltage_data,noise)
    PMT_trace = np.delete(PMT_trace, noise)
    
    voltage_norm = normalize(voltage_data)

        
#   determine starting point: period1 and period2 are starting timepoints of
#   two consecutive periods

#   main principle: if there is a big difference between voltage values in two
#   timepoints, it means the pulse signal and start of the period

#   ii is the entire measurement time, but we will break the cycle very soon
#   jj is the timepoint for finding the following period
#   difference between same phases of two consecutive periods is period length
    for ii in range(len(voltage_norm)-1):
        if (voltage_norm[ii]-voltage_norm[ii+1]) <= -0.9:
            start_point = ii
            period1 = ii
            for jj in range(ii+1,len(voltage_norm)-1):
                if (voltage_data[jj]-voltage_norm[jj+1]) <= -0.9:
                    period2 = jj
                    break
            break

    period_len = period2-period1
    dataset_new = PMT_trace[start_point::]
    num_periods = len(dataset_new)//period_len
    
#   now splitting periods and store it in dictionary
#   key(ii) is the index, item is a list of values
    periods = {}
    for ii in range(num_periods):
        periods[ii] = dataset_new[ii*num_periods:ii*num_periods+period_len]
    
#   and then normalize each period, jj is the index
    periods_norm = {}        
    for jj in range(len(periods)):
        periods_norm[jj] = normalize(periods[jj])
    
    return periods_norm, start_point, period_len


# =============================================================================
#   averages each timepoint of all normalized periods
# =============================================================================
def averaging(periods_norm):
    
    avg = {}

#   going over each timepoint within the cycle (ii)
#   and finding an average value for this timepoint over each period (jj)
    for ii in range(len(periods_norm[0])):
        avg[ii] = 0
        for jj in range(len(periods_norm)):
            avg[ii] = avg[ii] + periods_norm[jj][ii]
        avg[ii] = avg[ii] / len(periods_norm)
            
    return avg


# =============================================================================
#     F(t) = A × (C × exp(–t/t1) + (1 - C) × exp(–t/t2)), where t1 was the time
#     constant of the fast component and t2 was the time constant of the slow 
#     component. The percentage of the total magnitude that was associated with
#     the fast component (%t1) was defined as C above.
# =============================================================================

def func(t, t1, t2, a, c):
    return a * (c * np.exp(-t / t1) + (1 - c) * np.exp(-t / t2))


# =============================================================================
#   CAUTION: pain
#   still to determine time which we need to take for fitting
#   so far we suppose that the data is too noisy for automatic determination
# =============================================================================        
def fitting(avg): 
    
    ydata = np.array(list(avg.values())[1500:2500])
    xdata = np.arange(len(ydata))/sampling_rate
    
#   popt - array of parameters
#   popcov - estimated covariance of popt    
    popt, pcov = curve_fit(func, xdata, ydata, \
                           bounds=([0,0,-np.inf,0],[0.1,0.1,np.inf,1]), maxfev = 500000)
    
    return popt, pcov


# =============================================================================
#   determining all the parameters using previous functions
# =============================================================================
periods_norm, start_point, period_len = read_file('input_file.csv')   
avg = averaging(periods_norm)
popt, pcov = fitting(avg) 
data = np.array(list(avg.values())[1500:2500])
time = np.arange(len(data))/sampling_rate
#   modelling a curve using parameters from fitting results
model = []
for ii in time:
    model.append(func(ii,popt[0],popt[1],popt[2],popt[3]))

# =============================================================================
#   plotting data and model
# =============================================================================
plt.figure()
plt.plot(time,data,label='Real data')
plt.plot(time,model,label='Model')
plt.legend()
plt.grid()
plt.title('Fit of data')
plt.xlabel('Time (s)')
plt.ylabel('Intensity')
plt.text(0.0025,0.62, 'Fit formula: F(t) = A × (C × exp(–t/t1) + (1 - C) × exp(–t/t2))\nt1 fast component (ms) = %s\nt2 slow component (ms) = %s\na = %s\nc = %s' % (round(popt[0]*1000,3),round(popt[1]*1000,3),round(popt[2],3),round(popt[3],3)),fontsize='large')
plt.show()

# =============================================================================
#   plotting average curve for comparison
# =============================================================================
plt.figure()
avg_time = np.array(list(avg.keys()))/sampling_rate
avg_val = list(avg.values())

#   choose smoothing method
#avg_val = scipy.signal.savgol_filter(avg_val,51,7) # apply a Savitzky-Golay filter
avg_val = smoothTriangle(avg_val,5)


plt.plot(avg_time, list(avg.values())) 
plt.title('Averaged normalized curve')
plt.xlabel('Time (s)')
plt.ylabel('Intensity')
plt.grid()
plt.show()


# =============================================================================
#   checking whether there is an improvement in SNR after curve averaging
# =============================================================================

def SNR_improvement(avg_val,periods_norm):
#   SNR for averaged curve
    avg_mean = np.mean(avg_val)
    avg_std = np.std(avg_val)
    avg_SNR = avg_mean/avg_std
    
# to check SNR in a random period over available ones
    per = periods_norm[random.randint(0,len(periods_norm))]  
    per_val = list(per)                              
    per_mean = np.mean(per_val)
    per_std = np.std(per_val)
    per_SNR = per_mean/per_std
    
    N = avg_SNR/per_SNR
    
    return N

N = SNR_improvement(avg_val,periods_norm)