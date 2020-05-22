"""
code to filter floatt data using some different filter
"""

# Import modules
from scipy import signal
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import sys
import time
import datetime
from netCDF4 import Dataset

toolpath = '/Users/jadesauve/Documents/Python/scripts/python_tools/'

# Import DIY modules
pht = os.path.abspath(toolpath)
if pht not in sys.path:
    sys.path.append(pht)
from toolbox_math import *
from toolbox_float import *

################# Parameters ###############
directory_in = '/Users/jadesauve/Coding/data/SOCCOM_HRQC_LIAR_netcdf_20200223_AddLRQC/'  # SOCCOM_testing_data_2, SOCCOM_HRQC_LIAR_netcdf_20200112_AddLRQC
filename = '9101SOOCN_HRQC.nc' # '9254SOOCN_HRQC.nc'

# list of variables to be dropped
drop_var = ['Type', 'Parameters', 'REFERENCE_DATE_TIME', 'Chl_a', 'Chl_a_QFA', 'Chl_a_corr', 'Chl_a_corr_QFA',
            'b_bp700', 'b_bp700_QFA', 'b_bp_corr', 'b_bp_corr_QFA', 'b_bp532', 'b_bp532_QFA', 'POC', 'POC_QFA', 'pH25C',
            'pH25C_QFA', 'JULD']

depref = 10  # reference pressure to be taken as closest available to this number in dbar

# choose how much data wrangling to do : remove mean or remove mean and trend
# dataw = 'notrend' #'nomean', 'notrend' 
# # choose the sampling rate: raw/unchanged or regular 10/9 days
# srate = 'reg' # 'reg', 'raw'
# # Choose to plot annual/semi-annul f or peaks
# plotfeat = 'an' #'peak','an'
# # Choose to plot the fourier amplitdue or the energy density
# plotw = 'edens' #'famp', 'edens'

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
ds,floatnum, date_list, date_obj = load_nc_to_ds(directory_in,filename,drop_var, depref)

#avg lat and lon
lat_avg = ds.Lat.mean().values
lon_avg = ds.Lon.mean().values

inter_o2, o2_anom, time_days, o2_avg = ds_to_1Darr('Oxygen',ds)
inter_dic, dic_anom, __, dic_avg = ds_to_1Darr('DIC_LIAR',ds)
inter_T, T_anom, __, T_avg = ds_to_1Darr('Temperature',ds)
inter_s, s_anom, __, s_avg = ds_to_1Darr('Salinity',ds)

o2_anom_cut = o2_anom[:149] 

# Filters
o2_combined = pd.Series(index=time_days,data=o2_avg)

# Moving avg - evenly weighted
ma = o2_combined.rolling(center=True,window=11,min_periods=3).mean()

# boxcar
bc = o2_combined.rolling(center=True,window=11,min_periods=3,win_type='boxcar').mean()

# triangular filter
# win_wid = 10
# window = signal.triang(win_wid)
tri = o2_combined.rolling(center=True,window=11,min_periods=3,win_type='triang').mean()

# Gaussian
gau = o2_combined.rolling(center=True,window=11,min_periods=3,win_type='gaussian').mean(std=3)

# Butterworth
fs = 10  # Sampling frequency
fc = 2  # Cut-off frequency of the filter
w = fc / (fs / 2) # Normalize the frequency
b, a = signal.butter(5, w, 'low')
output = signal.filtfilt(b, a, o2_avg)

#plot
#plt.plot(ma,label='Evenly weighted')
plt.plot(time_days,o2_avg,label='raw')
plt.plot(bc,label='Boxcar, win=11')
plt.plot(tri,label='triangular, win=11')
plt.plot(gau,label='Gaussian, win=11')
plt.plot(time_days, output, label='Butterworth, order=5')
plt.grid()
plt.legend()
plt.ylabel('Dissolved Oxygen')
plt.xlabel('Time (days)')
plt.show()




