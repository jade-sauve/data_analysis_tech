# -*- coding: utf-8 -*-
# @Author: jadesauve
# @Date:   2020-04-02 15:46:24
# @Last Modified by:   jadesauve
# @Last Modified time: 2020-04-27 11:42:56
"""

file to read and plot oxygen tranform
also computes and plot the autospectrum
adapted from steves code for float data

"""

import numpy as np
from scipy import stats
import xarray as xr
import os
import sys
import time
import datetime
import matplotlib.pyplot as plt
from netCDF4 import Dataset
pht = os.path.abspath('/Users/jadesauve/Documents/Python/scripts/python_tools/')
if pht not in sys.path:
    sys.path.append(pht)
from toolbox_float import *
from toolbox_math import *
from toolbox_plot import *


################# Parameters ###############
directory_in = '/Users/jadesauve/Coding/data/SOCCOM_HRQC_LIAR_netcdf_20200223_AddLRQC/'  # SOCCOM_testing_data_2, SOCCOM_HRQC_LIAR_netcdf_20200112_AddLRQC
filename = '9101SOOCN_HRQC.nc'

# list of variables to be dropped
drop_var = ['Type', 'Parameters', 'REFERENCE_DATE_TIME', 'Chl_a', 'Chl_a_QFA', 'Chl_a_corr', 'Chl_a_corr_QFA',
            'b_bp700', 'b_bp700_QFA', 'b_bp_corr', 'b_bp_corr_QFA', 'b_bp532', 'b_bp532_QFA', 'POC', 'POC_QFA', 'pH25C',
            'pH25C_QFA', 'JULD']

depref = 10  # reference pressure to be taken as closest available to this number in dbar

# choose how much data wrangling to do : remove mean or remove mean and trend
dataw = 'notrend' #'nomean', 'notrend' 
# choose the sampling rate: raw/unchanged or regular 10/9 days
srate = 'reg' # 'reg', 'raw'
# Choose to plot annual/semi-annul f or peaks
plotfeat = 'an' #'peak','an'
# Choose to plot the fourier amplitdue or the energy density
plotw = 'edens' #'famp', 'edens'

###########################################

# define function

def ds_to_1Darr(varname,ds,srate='reg',dataw='notrend'):
    """
    var is a string, varibale name from the ds 
    # choose how much data wrangling to do : remove mean or remove mean and trend
    dataw = 'nomean' or 'notrend' 
    # choose the sampling rate: raw/unchanged or regular 10/9 days
    srate = reg' or 'raw'

    """
    # select var data
    var = ds[varname]
    # select top 150m
    var = var.isel(N_LEVELS=slice(-100,-1))
    # avg over depth
    var_avg = var.mean(dim='N_LEVELS',skipna=True).values

    # create lists of date in readable format
    date_list = []
    date_obj = []
    for profile in ds.N_PROF:
        date_str = ds.mon_day_yr.values[profile].decode("utf-8")
        time_str = ds.hh_mm.values[profile].decode("utf-8")
        date_list.append(date_str + '/' + time_str)
        date_obj.append(datetime.datetime.strptime(date_str + '/' + time_str, '%m/%d/%Y/%H:%M'))

    # make arrays for the time in days and hours
    time_days=np.zeros(len(date_obj))
    time_days[0] = 0
    for i in range(len(date_obj)-1):
      time_days[i+1] = (date_obj[i+1]-date_obj[0]).days
    time_hours = time_days*24.

    # compute the interval in time between each data point
    inter = np.diff(time_days)
    # inter = np.zeros(len(date_obj)-1)
    # for i in range(len(date_obj)-1):
    #   inter[i] = (date_obj[i+1]-date_obj[i]).days

    # find the index where the first value close to 10 is (there all only 5s near the start)
    ind = np.where(inter >5)[0][0] # MODIFY IF CHANGE FILE

    # compute a new list of dates to have a regular interval
    date_obj_10 = date_obj.copy()
    var_avg_10 = var_avg.copy()
    for t in np.arange(ind-2,0,-2):
        del date_obj_10[t]

    # same for o2
    var_avg_10=np.delete(var_avg_10, np.arange(ind-2,0,-2))

    time_days_10=np.zeros(len(date_obj_10))
    time_days_10[0] = 0
    for i in range(len(date_obj_10)-1):
      time_days_10[i+1] = (date_obj_10[i+1]-date_obj_10[0]).days
    time_hours_10 = time_days_10*24.

    # apply parameter
    if srate == 'raw':
        pass
    elif srate == 'reg':
        var_avg = var_avg_10
        time_days = time_days_10
    else:
        print('Choice of srate is not valid.')

    #remove the nan values from dic and o2 and time
    #inds = np.where(~np.isnan(dic_avg))[0]
    if varname == 'DIC_LIAR':
        ind_end = 149 #MODIFY IF CHANGE FILE
        var_avg = var_avg[:ind_end]
        time_days = time_days[:ind_end]

    # remove the mean from the data
    var_mean = np.nanmean(var_avg)
    var_anom = var_avg - var_mean

    # find the least squares regression line
    slope, intercept, r_value, p_value, std_err = stats.linregress(time_days,var_avg)
    y = intercept + slope*time_days

    if dataw == 'notrend':
        # remove the trend from the data
        var_anom = var_avg - y


    return inter, var_anom, time_days, var_avg


# main program

# read the input file (netcdf)
ds,floatnum, date_list, date_obj = load_nc_to_ds(directory_in, filename, drop_var, depref)

#avg lat and lon
lat_avg = ds.Lat.mean().values
lon_avg = ds.Lon.mean().values

# wrangle data into shape
inter, o2_anom, time_days_o2, o2_avg = ds_to_1Darr('Oxygen',ds)
inter_dic, dic_anom, time_days_dic, dic_avg = ds_to_1Darr('DIC_LIAR',ds)
inter_T, T_anom, time_days_T, T_avg = ds_to_1Darr('Temperature',ds)
inter_s, s_anom, time_days_s, s_avg = ds_to_1Darr('Salinity',ds)
#inter_ml, mlp_anom, time_days_ml, mlp_avg = ds_to_1Darr('mlp',ds)

# compute spectrum
freq, fourier_amp, __, spec, __, w, fourier_combined = spectrum(o2_anom)
freq2, fourier_amp2, __, spec2, __, w2, fourier_combined2 = spectrum(o2_anom[:149])
freq_c, fourier_amp_c, __, spec_c,__,__,fourier_combined_c = spectrum(dic_anom)
freq_s, fourier_amp_s, __, spec_s,__, __,fourier_combined_s = spectrum(s_anom)
freq_T, fourier_amp_T, __, spec_T,__, __,fourier_combined_T  = spectrum(T_anom)
#freq_ml, fourier_amp_ml, fourier_phase_ml, spec_ml, spec_amp_ml, omega_max_ml, omega0_ml = spectrum(mlp_anom)

# plot the time series
title = 'SOCCOM FLoat #'+str(ds.Cruise.values)+', Avg Lat = '+str(lat_avg.round(2))+', Avg Lon = '+str(lon_avg.round(2))
plot_2var_from_1Darr(time_days_dic, o2_avg[:149], dic_avg, title, y1label='Oxygen Concentration', y2label='DIC', xlabel='Time (days)')
plot_2var_from_1Darr(time_days_T, T_avg, s_avg, title=title, y1label='Temperature', y2label='Salinity', xlabel='Time (days)')

# plot the spectrum as a function of frequency
# show the maximum and minimum frequencies

title = 'Oxygen and DIC, sampling = '+srate+', '+dataw
plot_spectrum_2var(freq_c, spec2, spec_c, w2[1], w2[0],title=title, maxy=1e-1, miny=5e-7,label1='Oxygen',label2='DIC')

title = 'Temperature and Salinity, sampling = '+srate+', '+dataw
plot_spectrum_2var(freq_T, spec_T, spec_s, w[1], w[0],title=title, maxy=2.5e-5, miny=1e-12,label1='Temperature',label2='Salinity')

title = 'Temperature and Oxygen, sampling = '+srate+', '+dataw
plot_spectrum_2var(freq_T, spec_T, spec, w[1], w[0],title=title, maxy=2.5e-3, miny=2e-11,label1='Temperature',label2='Oxygen')

title = 'Salinity and Oxygen, sampling = '+srate+', '+dataw
plot_spectrum_2var(freq_s, spec_s, spec, w[1], w[0],title=title, maxy=2.5e-3, miny=1e-12,label1='Salinty',label2='Oxygen')


# # #**** test data   ****#

t = np.arange(0,1000,1)
y_sin = np.sin(t*(2*pi/50))

freq_t, fourier_amp_t, fourier_phase_t, spec_t, spec_amp_t, omega_max_t, omega0_t = spectrum(y_sin)

#plot
fig2=plt.figure(figsize=(7,7))
# max_y = 0.
# min_y = 2500.
plt.plot(freq_t,fourier_amp_t,color='darkkhaki')
plt.ylabel('Fourier Coefficient Magnitude for sin function',fontsize=15)
# plt.ylim(min_y,max_y)
plt.xlabel('$\omega$ (radians/day)',fontsize=15,ha='center')
plt.title('Sin Function, Amp=1')
plt.grid()
plt.show()

# #**** ************ ****# 

# Frequencies of interest for different cases:
# srate = 'raw', freq1 =  0.02208
# srate = 'reg', dataw = 'nomean', freq1 = 0.01778, 0.0623, 0.0445, 0.0297, 0.0237
# srate = 'reg', dataw = 'notrend', freq1 = 0.01778, 0.0445, 0.0356

# Period associated with a frequency of interest
# T1 = 1/(freq1/(2*pi))


# fig2=plt.figure(figsize=(10,7))

# max_y = 2e-3
# min_y = 1e-8

# if plotw == 'famp':
#     plt.loglog(freq,fourier_amp,color='magenta')
#     plt.ylabel('Fourier Amplitude ($\mu$mol/kg)',fontsize=15)
# elif plotw == 'edens':
#     plt.loglog(freq,spec,color='orange',label='Oxygen')
#     plt.loglog(freq_c,spec_c,color='brown',label='DIC')
#     plt.ylabel('Energy Density ($\mu$mol/kg)$^2$/(radian/day))',fontsize=15)

# plt.plot([omega_max,omega_max],[min_y,max_y],'--k')
# plt.text(omega_max-0.14,min_y+(min_y*2),'$\omega_{max}$',fontsize=15,color='blue') #modify

# plt.plot([omega0,omega0],[min_y,max_y],'--k',zorder=10)
# plt.text(omega0+0.0001,min_y+(min_y*2),'$\omega_o$',fontsize=15,color='blue')

# # if plotfeat == 'peak':
# #     # freq_285 = 0.02208
# #     plt.text(freq1,min_y+(min_y*2),' '+str(T1)+' days',color='blue')
# #     plt.plot([freq1,freq1],[min_y,max_y],'--k',zorder=10)
# #elif plotfeat == 'an':
# freq_an = 2*pi*(1/365.25)
# plt.text(freq_an,min_y+(min_y*2),' Annual',color='blue')
# plt.plot([freq_an,freq_an],[min_y,max_y],'--k',zorder=10)

# freq_semAn = 2*pi*(2/365.25)
# plt.text(freq_semAn,min_y+(min_y*2),' Semi-An.',color='blue')
# plt.plot([freq_semAn,freq_semAn],[min_y,max_y],'--k',zorder=10)

# freq_seas = 2*pi*(4/365.25)
# plt.text(freq_seas,min_y+(min_y*2),' Seasonal',color='blue')
# plt.plot([freq_seas,freq_seas],[min_y,max_y],'--k',zorder=10)
# # else:
# #     print('Choice for plotfeat is not valid.')

# plt.ylim(min_y,max_y)
# plt.xlabel('$\omega$ (radians/day)',fontsize=15,ha='center')
# plt.title('Oxygen and DIC, sampling = '+srate+', '+dataw)
# plt.grid()
# plt.legend()
# #plt.savefig(path_out2)
# plt.show()
