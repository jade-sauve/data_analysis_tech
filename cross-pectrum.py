# -*- coding: utf-8 -*-
# @Author: jadesauve
# @Date:   2020-04-15 14:48:46
# @Last Modified by:   jadesauve
# @Last Modified time: 2020-04-28 10:44:21
"""
code to compute the cross-spectrum of two variables avg over 21 frequency bands
adapted from steve riser

"""

# file to read and plot WHOI HOT temperature and salinity cross-spectral ampliude
#
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
filename = '9101SOOCN_HRQC.nc' # '9254SOOCN_HRQC.nc'

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
ds,floatnum, date_list, date_obj = load_nc_to_ds(directory_in,filename,drop_var, depref)

#avg lat and lon
lat_avg = ds.Lat.mean().values
lon_avg = ds.Lon.mean().values

inter_o2, o2_anom, time_days, o2_avg = ds_to_1Darr('Oxygen',ds)
inter_dic, dic_anom, __, dic_avg = ds_to_1Darr('DIC_LIAR',ds)
inter_T, T_anom, __, T_avg = ds_to_1Darr('Temperature',ds)
inter_s, s_anom, __, s_avg = ds_to_1Darr('Salinity',ds)

o2_anom_cut = o2_anom[:149] 

freq,__,__,__,__,w,__ = spectrum(o2_anom)
freq_cut,__,__,__,__,w_cut,__ = spectrum(o2_anom_cut)

cospec_amp_oc ,cross_spec_oc, coherence_sq_oc,*rest,cross_phase_combined_oc,coherence_sq_combined_oc = cospec_coher(o2_anom_cut,dic_anom,freq_cut) 
cospec_amp_ot ,cross_spec_ot, coherence_sq_ot,*rest,cross_phase_combined_ot,coherence_sq_combined_ot = cospec_coher(o2_anom,T_anom,freq) 


# plot coherence
fig1=plt.figure(figsize=(9,7))

max_y = 0
min_y = 1

plt.semilogx(freq,coherence_sq_ot,color='purple')

plt.plot([w[1],w[1]],[min_y,max_y],'--k')
plt.text(w[1]-0.1,0.6,'$\omega_{max}$',fontsize=12,color='firebrick')

#plt.text(1.2,0.88,'n_av = 21')

plt.xlabel('$\omega$ (radians/day)',fontsize=15,ha='center')
plt.ylabel('Squared Coherence $\it{O2}$-$\it{DIC}$',fontsize=15)
plt.grid(which='both')
plt.title('Oxygen and DIC Coherence, sampling = '+srate+', '+dataw)
plt.ylim(min_y,max_y)

plt.show()

# plt.close('all')
# N= len(T_anom)
# dt=10
# t_length =  N * dt
# T_fft_raw = fft.rfft(s_anom,n=N)
# S_fft_raw = fft.rfft(o2_anom,n=N) #m
# cross_spec = conj(T_fft_raw) * S_fft_raw / t_length 
# cross_phase_in_rad = arctan2(imag(cross_spec),real(cross_spec))
# cross_phase = 180 * cross_phase_in_rad / pi  
# cospec_amp = real(cross_spec)   
# T_fft_psd = real(conj(T_fft_raw) * T_fft_raw) / t_length          # PSD estimate for T (see note above about t_length)
# S_fft_psd = real(conj(S_fft_raw) * S_fft_raw) / t_length 
# coherence_sq = cospec_amp**2 / (T_fft_psd * S_fft_psd)
# cross_phase_combined = pd.Series(index=freq,data=cross_phase_in_rad)
# coherence_sq_combined = pd.Series(index=freq,data=coherence_sq) 

# plot phase

fig3 = plt.figure(figsize=(14,7))
plt.semilogx(freq,coherence_sq_ot,c='teal',alpha=0.5,lw=0.5,label='Original')
plt.semilogx(coherence_sq_combined.rolling(center=True,window=11,min_periods=3).mean(),c='teal',lw=1.5,label='Band-averaged over 11 bands')
plt.legend()
# plt.semilogx(*signal.coherence(T,S),c='r',lw=1)  # SciPy version
plt.xlabel('$\omega$ (radians/day)',fontsize=15)
plt.ylabel('Squared coherence between $T-O2$',fontsize=15) #m
plt.xlim([min(freq[1:]),max(freq)])
#plt.ylim([0,1])
plt.grid()
plt.title('Squared coherence of Temperature and Dissolved Oxygen',fontsize=15) #m
plt.show()

fig4 = plt.figure(figsize=(14,7))
cross_phase[cross_phase < 0] += 360   # shift from (-180,180) to (0,360)
plt.semilogx(freq,cross_phase,c='orange',alpha=0.5,lw=0.5,label='Original')
# note: band-averaging when result is in degrees from -180 to 180 requires special formula for averaging circular quantities:
#       https://en.wikipedia.org/wiki/Mean_of_circular_quantities#Mean_of_angles
rolling_intermediate_imag = exp(cross_phase_combined * 1j).apply(imag).rolling(center=True,window=11,min_periods=3).sum()
rolling_intermediate_real = exp(cross_phase_combined * 1j).apply(real).rolling(center=True,window=11,min_periods=3).sum()
rolling_cross_phase_averaged_properly = (180/pi)*arctan2(rolling_intermediate_imag,rolling_intermediate_real)
rolling_cross_phase_averaged_properly[rolling_cross_phase_averaged_properly < 0] += 360      # shift from (-180,180) to (0,360)
plt.semilogx(rolling_intermediate_imag.index,rolling_cross_phase_averaged_properly,c='orange',lw=1.5,label='Band-averaged over 11 bands')
plt.legend()
plt.xlabel('$\omega$ (radians/day)',fontsize=15)
plt.ylabel('Phase between $\S-O2$ (degrees)',fontsize=15) #m
plt.xlim([min(freq[1:]),max(freq)])
plt.ylim([0,360])
plt.grid()
plt.title('Phase for Oxygen and Salinity',fontsize=15); #m




#**** test data   ****#
t = np.arange(0,200*np.pi,0.1)
y_sin = np.sin(t)
y_cos = np.cos(t)

pi = np.pi
omega0_t = 2.*pi/(len(t)*0.1)
omega_max_t = pi/0.1 # aslo f_Ny
freq_t = arange(0,omega_max_t+0.000001,omega0_t)

cospec_amp_t ,cross_spec_t, coherence_sq_t,*rest_t,cross_phase_combined_t,coherence_sq_combined_t = cospec_coher(y_sin,y_cos,freq_t,dt=0.1) 


# fig2=plt.figure(figsize=(7,7))

# max_y = 1e8
# min_y = 1e-5

# plt.loglog(freq_av_t,cospec_amp_t,color='darkkhaki')
# plt.ylabel('((U)(U)/RPD)',fontsize=15)

# plt.plot([omega_max_t,omega_max_t],[min_y,max_y],'--k')
# #plt.text(omega_max-0.14,min_y+(min_y*2),'$\omega_{max}$',fontsize=15,color='blue') #modify

# plt.plot([omega0_t,omega0_t],[min_y,max_y],'--k',zorder=10)
# #plt.text(omega0+0.0001,min_y+(min_y*2),'$\omega_o$',fontsize=15,color='blue')

# plt.ylim(min_y,max_y)
# plt.xlabel('$\omega$ (radians/day)',fontsize=15,ha='center')

# plt.title('Sin and Cos cross-spectrum, Amp=1')
# plt.grid()
# plt.show()

# fig3 = plt.figure(figsize=(14,7))
# plt.semilogx(freq_t,coherence_sq_t,c='k',alpha=0.5,lw=0.5,label='Original')
# plt.semilogx(coherence_sq_combined_t.rolling(center=True,window=3,min_periods=3).mean(),c='k',lw=1.5,label='Band-averaged over 3 bands')
# plt.legend()
# # plt.semilogx(*signal.coherence(T,S),c='r',lw=1)  # SciPy version
# plt.xlabel('$\omega$ (radians/day)',fontsize=15)
# plt.ylabel('Squared coherence between Sin and Cos',fontsize=15)
# plt.xlim([min(freq[1:]),max(freq)])
# #plt.ylim([0,1])
# plt.grid()
# plt.title('Squared coherence is 1.0 because Sin and Cos are antiphase.',fontsize=15)
# plt.show()

#**** ************ ****# 

# begin plotting the spectrum

freq1 = 0.01685
T1 = 1/(freq1/(2*pi))

fig1=plt.figure(figsize=(7,7))

max_y = 20
min_y = 2e-3 

plt.loglog(freq,cospec_amp,color='purple')

plt.plot([omega_max,omega_max],[min_y,max_y],'--k')
plt.text(omega_max-0.1,min_y+(min_y),'$\omega_{max}$',fontsize=12,color='firebrick')

plt.plot([omega0,omega0],[min_y,max_y],'--k',zorder=10)
plt.text(omega0+0.0001,min_y+(min_y*2),'$\omega_o$',fontsize=15,color='blue')

#plt.text(0.01,2*min_y,'n_av = '+str(n_av))
if plotfeat == 'peak':
    plt.text(freq1,min_y+(min_y*2),' '+str(round(T1,2))+' days',color='blue')
    plt.plot([freq1,freq1],[min_y,max_y],'--k',zorder=10)
elif plotfeat == 'an':
    freq_an = 2*pi*(1/365.25)
    plt.text(freq_an,min_y+(min_y*2),' Annual',color='blue')
    plt.plot([freq_an,freq_an],[min_y,max_y],'--k',zorder=10)

    freq_semAn = 2*pi*(2/365.25)
    plt.text(freq_semAn,min_y+(min_y*2),' Semi-Annual',color='blue')
    plt.plot([freq_semAn,freq_semAn],[min_y,max_y],'--k',zorder=10)
else:
    print('Choice for plotfeat is not valid.')

plt.xlabel('$\omega$ (radians/day)',fontsize=15,ha='center')
plt.ylabel(' (($\mu$mol/kg)($\mu$mol/kg)/RPD)',fontsize=15)
plt.ylim(min_y,max_y)
plt.grid(which='both')

plt.title('Oxygen and DIC cross-spectrum, sampling = '+srate+', '+dataw)

plt.show()


# begin plotting the coherence and phase
#
# first the coherence
#
fig1=plt.figure(figsize=(9,7))

max_y = 1.5
min_y = 0.5

plt.semilogx(freq_av,coh_sq,color='purple')

plt.plot([omega_max,omega_max],[min_y,max_y],'--k')
plt.text(omega_max-0.1,0.6,'$\omega_{max}$',fontsize=12,color='firebrick')

#plt.text(1.2,0.88,'n_av = 21')

plt.xlabel('$\omega$ (radians/day)',fontsize=15,ha='center')
plt.ylabel('Squared Coherence $\it{O2}$-$\it{DIC}$',fontsize=15)
plt.grid(which='both')
plt.title('Oxygen and DIC Coherence, sampling = '+srate+', '+dataw)
plt.ylim(min_y,max_y)

plt.show()


#
# now the phase
#
fig2=plt.figure(figsize=(9,7))

max_y = -180.
min_y = 180. 

plt.semilogx(freq_av,cospec_phase,color='orange')

plt.plot([omega_max,omega_max],[min_y,max_y],'--k')
plt.text(omega_max-0.1,min_y+10,'$\omega_{max}$',fontsize=12,color='firebrick')

plt.text(1.2,130.,'n_av = 21')

plt.xlabel('$\omega$ (radians/day)',fontsize=15,ha='center')
plt.ylabel('Phase $\it{O2}$-$\it{DIC}$, degrees',fontsize=15)
plt.ylim(min_y,max_y)
plt.grid(which='both')
plt.title('Oxygen and DIC Phase, sampling = '+srate+', '+dataw)

plt.show()


# spectrum of dic


zz=np.fft.rfft(dic_2,n=nn)/nn
fourier_amp=np.sqrt((np.real(zz)**2+np.imag(zz)**2))
fourier_phase=180.*np.arctan2(np.imag(zz),np.real(zz))/pi
spec=np.real(zz*np.conj(zz))/(2.*pi*nn*dt)
spec_amp=(np.absolute(zz))**2/(2.*pi*nn*dt)


fig2=plt.figure(figsize=(7,7))

max_y = 2e-3
min_y = 2.e-7

plt.loglog(freq,spec,color='darkkhaki')
plt.ylabel('Energy Density ($\mu$mol/kg)$^2$/(radian/day))',fontsize=15)

plt.plot([omega_max,omega_max],[min_y,max_y],'--k')
plt.text(omega_max-0.14,min_y+(min_y*2),'$\omega_{max}$',fontsize=15,color='blue') #modify

plt.plot([omega0,omega0],[min_y,max_y],'--k',zorder=10)
plt.text(omega0+0.0001,min_y+(min_y*2),'$\omega_o$',fontsize=15,color='blue')

freq_an = 2*pi*(1/365.25)
plt.text(freq_an,min_y+(min_y*2),' Annual',color='blue')
plt.plot([freq_an,freq_an],[min_y,max_y],'--k',zorder=10)

freq_semAn = 2*pi*(2/365.25)
plt.text(freq_semAn,min_y+(min_y*2),' Semi-Annual',color='blue')
plt.plot([freq_semAn,freq_semAn],[min_y,max_y],'--k',zorder=10)


plt.ylim(min_y,max_y)
plt.xlabel('$\omega$ (radians/day)',fontsize=15,ha='center')

plt.title('DIC, sampling = reg, notrend')
plt.grid()
#plt.savefig(path_out2)
plt.show()


#
# define a function to band-average the spectrum to eliminate some noise
#
# def band_average(fft_var1,fft_var2,frequency,n_av):

#     # fft_var1 and fft_var2 are the inputs computed via fft
#     # they can be the same variable or different variables
#     # n_av is the number of bands to be used for smoothing (nice if it is an odd number)
#     # this function is limnited to 100,000 points but can easily be modified

#     nmax=100000

#     # define some variables and arrays
#     n_spec=len(fft_var1)
#     n_av2=int(n_av//2+1)
#     spec_amp_av=np.zeros(nmax)
#     spec_phase_av=np.zeros(nmax)
#     freq_av=np.zeros(nmax)

#     # average the lowest frequency bands first (with half as many points in the average)
#     sum_low_amp=0.
#     sum_low_phase=0.
#     count=0
#     spectrum_amp=np.absolute(fft_var1*np.conj(fft_var2))
#     spectrum_phase=np.angle(fft_var1*np.conj(fft_var2),deg=True)

#     for i in range(0,n_av2):
#         sum_low_amp+=spectrum_amp[i]
#         sum_low_phase+=spectrum_phase[i]

#     spec_amp_av[0]=sum_low_amp/n_av2
#     spec_phase_av[0]=sum_low_phase/n_av

#     # compute the rest of the averages

#     for i in range(n_av2,n_spec-n_av,n_av):
#         count+=1
#         spec_amp_est=np.mean(spectrum_amp[i:i+n_av])
#         spec_phase_est=np.mean(spectrum_phase[i:i+n_av])
#         freq_est=frequency[i+n_av//2]
#         spec_amp_av[count]=spec_amp_est
#         spec_phase_av[count]=spec_phase_est
#         freq_av[count]=freq_est

#     # contract the arrays
#     spec_amp_av=spec_amp_av[0:count]
#     spec_phase_av=spec_phase_av[0:count]
#     freq_av=freq_av[0:count]
#     return spec_amp_av,spec_phase_av,freq_av,count  

# # select O2 data
# o2 = ds.Oxygen
# dic = ds.DIC_LIAR

# # select top 150m
# o2 = o2.isel(N_LEVELS=slice(-100,-1))
# dic = dic.isel(N_LEVELS=slice(-100,-1))

# # avg over depth
# o2_avg = o2.mean(dim='N_LEVELS',skipna=True).values
# dic_avg = dic.mean(dim='N_LEVELS',skipna=True).values

# # create lists of date in readable format
# date_list = []
# date_obj = []
# for profile in ds.N_PROF:
#     date_str = ds.mon_day_yr.values[profile].decode("utf-8")
#     time_str = ds.hh_mm.values[profile].decode("utf-8")
#     date_list.append(date_str + '/' + time_str)
#     date_obj.append(datetime.datetime.strptime(date_str + '/' + time_str, '%m/%d/%Y/%H:%M'))

# # make arrays for the time in days and hours
# time_days=np.zeros(len(date_obj))
# time_days[0] = 0
# for i in range(len(date_obj)-1):
#   time_days[i+1] = (date_obj[i+1]-date_obj[0]).days
# time_hours = time_days*24.

# # compue the interval in time between each data point
# inter = np.zeros(len(date_obj)-1)
# for i in range(len(date_obj)-1):
#   inter[i] = (date_obj[i+1]-date_obj[i]).days

# # find the index where the first value close to 10 is (there all only 5s near the start)
# ind = np.where(inter >5)[0][0]

# # compute a new list of dates to have a regular interval
# date_obj_10 = date_obj.copy()
# o2_avg_10 = o2_avg.copy()
# dic_avg_10 = dic_avg.copy()
# for t in np.arange(ind-2,0,-2):
#     del date_obj_10[t]

# # same for o2
# o2_avg_10=np.delete(o2_avg_10, np.arange(ind-2,0,-2))
# dic_avg_10=np.delete(dic_avg_10, np.arange(ind-2,0,-2))

# # check the interval - only 9 and 10 now
# # inter_10 = np.zeros(len(date_obj_10)-1)
# # for i in range(len(date_obj_10)-1):
# #   inter_10[i] = (date_obj_10[i+1]-date_obj_10[i]).days

# time_days_10=np.zeros(len(date_obj_10))
# time_days_10[0] = 0
# for i in range(len(date_obj_10)-1):
#   time_days_10[i+1] = (date_obj_10[i+1]-date_obj_10[0]).days
# time_hours_10 = time_days_10*24.

# # apply parameter
# if srate == 'raw':
#     pass
# elif srate == 'reg':
#     o2_avg = o2_avg_10
#     dic_avg = dic_avg_10
#     time_days = time_days_10
# else:
#     print('Choice of srate is not valid.')

# #remove the nan values from dic and o2 and time
# #inds = np.where(~np.isnan(dic_avg))[0]
# ind_end = 149
# o2_avg = o2_avg[:ind_end]
# dic_avg = dic_avg[:ind_end]
# time_days = time_days[:ind_end]

# # remove the mean from the data
# o2_mean = np.nanmean(o2_avg)
# o2_2 = o2_avg - o2_mean

# dic_mean = np.nanmean(dic_avg)
# dic_2 = dic_avg - dic_mean

# # find the least squares regression line
# slope, intercept, r_value, p_value, std_err = stats.linregress(time_days,o2_avg)
# y_o2 = intercept + slope*time_days

# slope, intercept, r_value, p_value, std_err = stats.linregress(time_days,dic_avg)
# y_dic = intercept + slope*time_days

# if dataw == 'notrend':
#     # remove the trend from the data
#     diff = o2_avg - y_o2
#     o2_2 = diff

#     diff = dic_avg - y_dic
#     dic_2 = diff

# determine the length of the data and the resulting number of Fourier estimates
# nn = len(o2_avg) # number of sample n
# mm = int(nn/2+1) # number of Fourier estimates N/2


# # determine the frequencies for the spectrum
# if srate == 'raw':
#     dt = 8 # time interval of sampling [days]
# elif srate == 'reg':
#     dt = 10 

# pi=np.pi
# omega0=2.*pi/(nn*dt)
# omega_max = pi/dt # aslo f_Ny

# freq=np.zeros(mm)
# for i in range(0,mm):
#     freq[i]=i*omega0

# compute the fft of the input data 
# the fft will yield a set of complex numbers
# also compute their complex conjugates
# multiply these together to get the spectrum
#
# n_av=1
# o2_fft=np.fft.rfft(o2_2,n=nn)/nn
# dic_fft=np.fft.rfft(dic_2,n=nn)/nn
# dic_star=np.conj(dic_fft)
# #cospec_amp,cospec_phase,freq_av,count=band_average(o2_fft,dic_fft,freq,n_av) # should it be dic_star?

# o2_spec,o2_phase,freq_av,count=band_average(o2_fft,o2_fft,freq,n_av)
# dic_spec,dic_phase,freq_av,count=band_average(dic_fft,dic_fft,freq,n_av)
# cospec_amp,cospec_phase,freq_av,count=band_average(o2_fft,dic_star,freq,n_av)

# coh_sq=cospec_amp**2/(o2_spec*dic_spec)

# n_av_t = 1
# n_t = len(t) 
# m_t = int(n_t/2+1)
# dt_t = 0.1
# omega0_t = 2.*pi/(n_t*dt_t)
# omega_max_t = pi/dt_t # aslo f_Ny

# freq_t=np.zeros(m_t)
# for i in range(0,m_t):
#     freq_t[i]=i*omega0_t #/2*pi

# sin_fft=np.fft.rfft(y_sin,n=n_t)#/n_t
# cos_fft=np.fft.rfft(y_cos,n=n_t)#/n_t
# cos_star = np.conj(cos_fft)

# sin_spec,sin_phase,freq_av_t,count_t=band_average(sin_fft,sin_fft,freq_t,n_av_t)
# cos_spec,cos_phase,freq_av_t,count_t=band_average(cos_fft,cos_fft,freq_t,n_av_t)
# cospec_amp_t,cospec_phase_t,freq_av_t,count_t=band_average(sin_fft,cos_star,freq_t,n_av_t)

# coh_sq_t=cospec_amp_t**2/(sin_spec*cos_spec)