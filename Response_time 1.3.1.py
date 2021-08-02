# -*- coding: utf-8 -*-
"""
Created on Tue Jul 27 10:51:57 2021

@author: Polina
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import random
# %matplotlib qt


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
#   There is the possibility to exclude the noise, but it will produce the
#   slight phase shift, thus it is commented out
# =============================================================================
def read_file(input_file):
    PMT_trace = np.load(r"C:\Users\Polina\Desktop\ImPhys\Data + code\PMT_array_2020-08-06_11-34-20.npy",allow_pickle=True)[2:]
    #voltage_data = np.load(r"C:\Users\Polina\Desktop\ImPhys\Data + code\Vp2020-08-06_11-33-01.npy",allow_pickle=True)[7:]
        
    #noise = (np.where(PMT_trace < 0.9*np.average(PMT_trace)))[0]
    
    #voltage_data = np.delete(voltage_data,noise)
    #PMT_trace = np.delete(PMT_trace, noise)       

#   also importing from settings
    sampling_rate = 2e4
    sampling_freq = 5
    period_len = int(sampling_rate//sampling_freq)  # points
    sampling_time = len(PMT_trace)  # milliseconds
    num_periods = sampling_time//period_len 
    
#   splitting periods and store it in 2-dimensional array
    periods = PMT_trace.reshape(num_periods,period_len)
    periods_norm = [normalize(periods[x,:]) for x in range(num_periods)]
    
    return periods_norm, period_len, sampling_rate

# =============================================================================
#     F(t) = A × (C × exp(–t/t1) + (1 - C) × exp(–t/t2)) + D, where t1 is the time
#     constant of the fast component and t2 was the time constant of the slow 
#     component. The percentage of the total magnitude that was associated with
#     the fast component (%t1) was defined as C above. D is an offset.
# =============================================================================

def func(t, t1, t2, a, c, d):
    return a * (c * np.exp(-t / t1) + (1 - c) * np.exp(-t / t2)) + d

        
def fitting(avg):     
    ydata = np.array(avg[0:period_len//2])            #rising phase
    #ydata = np.array(avg[period_len//2]:period_len)  #falling phase
    xdata = np.arange(len(ydata))/sampling_rate
#   popt - array of parameters
#   popcov - estimated covariance of popt    
    popt, pcov = curve_fit(func, xdata, ydata, \
                           bounds=([0,0,-np.inf,0,-np.inf],[0.1,0.1,np.inf,1,np.inf]), \
                               maxfev = 500000)
    return popt, pcov


# =============================================================================
#   check
# =============================================================================
periods_norm, period_len, sampling_rate = read_file('input_file.csv')
y=np.array([np.array(xi) for xi in periods_norm])
avg = [np.mean(y[:,x]) for x in range(period_len)]
#avg = smoothTriangle(avg,5)
#avg = normalize(avg)

popt, pcov = fitting(avg) 
data = np.array(avg[0:period_len//2])               #rising phase
#data = np.array(avg[period_len//2]:period_len)     #falling phase
time = np.arange(len(data))/sampling_rate
#   modelling a curve using parameters from fitting results
model = []
for ii in time:
    model.append(func(ii,popt[0],popt[1],popt[2],popt[3],popt[4]))

# =============================================================================
#   plotting data and model
# =============================================================================
plt.figure()
plt.plot(time,data,label='Real data')
plt.plot(time,model,label='Model')
plt.legend(loc='upper right')
plt.grid()
#plt.title('Fit of data')
plt.xlabel('Time (s)')
plt.ylabel('Intensity')
plt.text(1e-4,np.max(data)+0.005, 'Fit formula: F(t) = A × (C × exp(–t/t1) + (1 - C) × exp(–t/t2))\nt1 fast component (ms) = %s\nt2 slow component (ms) = %s\na = %s\nc = %s'\
             % (round(popt[0]*1000,3),round(popt[1]*1000,3),round(popt[2],3), \
                round(popt[3],3)),fontsize='large')
plt.show()

# =============================================================================
#   plotting average curve for comparison
# =============================================================================
plt.figure()
avg_time = np.array(avg)/sampling_rate
avg_val = avg

#   choose smoothing method
#avg_val = scipy.signal.savgol_filter(avg_val,51,7) # apply a Savitzky-Golay filter
#avg_val = smoothTriangle(avg_val,7)

plt.plot(avg) 
plt.title('Averaged curve of normalized data')
plt.xlabel('Timepoints')
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
    per_val = per                              
    per_mean = np.mean(per_val)
    per_std = np.std(per_val)
    per_SNR = per_mean/per_std
    
    N = avg_SNR/per_SNR
    
    return N

N = SNR_improvement(avg_val,periods_norm)