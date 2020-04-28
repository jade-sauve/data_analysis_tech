# -*- coding: utf-8 -*-
# @Author: jadesauve
# @Date:   2020-04-15 14:50:08
# @Last Modified by:   jadesauve
# @Last Modified time: 2020-04-15 14:51:12
"""
code to compute the squared coherence between two variables avg over 21 frequency bands

adapted from steve riser

"""
# file to read HOT T and S data from WHOTS mooring and to
# estimate coherence and phase between T and S
#
import numpy as np
import os
import time
import matplotlib.pyplot as plt
from netCDF4 import Dataset
#
# read and plot WHOI HOT mooring data
#
# define the netcdf reading function
#
def import_data(file_name):   
    '''This function works regardless of the file called as long as a list of 
    variables can be identified easily in the file.'''
    
    data_netcdf = Dataset(file_name, mode = 'r')
    data = {}
#    
    for vname in list(data_netcdf.variables):
        data[str(vname)] = data_netcdf.variables[vname][:] 
#            
    data_netcdf.close()
#            
    return data
#
# define a function to band-average the spectrum to eliminate some noise
#
def band_average(fft_var1,fft_var2,frequency,n_av):
#
# this is a function to estimate the autospectrum or co-spectrum for 2 variables
# fft_var1 and fft_var2 are the inputs computed via fft
# they can be the same variable or different variables
# n_av is the number of bands to be used for smoothing (nice if it is an odd number)
# this function is limnited to 100,000 points but can easily be modified
#
    nmax=100000
#
# define some variables and arrays
#
    n_spec=len(fft_var1)
    n_av2=int(n_av//2+1)
    spec_amp_av=np.zeros(nmax)
    spec_phase_av=np.zeros(nmax)
    freq_av=np.zeros(nmax)
#
# average the lowest frequency bands first (with half as many points in the average)
# 
    sum_low_amp=0.
    sum_low_phase=0.
    count=0
    spectrum_amp=np.absolute(fft_var1*np.conj(fft_var2))
    spectrum_phase=np.angle(fft_var1*np.conj(fft_var2),deg=True)
#
    for i in range(0,n_av2):
        sum_low_amp+=spectrum_amp[i]
        sum_low_phase+=spectrum_phase[i]
#
    spec_amp_av[0]=sum_low_amp/n_av2
    spec_phase_av[0]=sum_low_phase/n_av
#
# compute the rest of the averages
#
    for i in range(n_av2,n_spec-n_av,n_av):
        count+=1
        spec_amp_est=np.mean(spectrum_amp[i:i+n_av])
        spec_phase_est=np.mean(spectrum_phase[i:i+n_av])
        freq_est=frequency[i+n_av//2]
        spec_amp_av[count]=spec_amp_est
        spec_phase_av[count]=spec_phase_est
        freq_av[count]=freq_est
#
# contract the arrays
#
    spec_amp_av=spec_amp_av[0:count]
    spec_phase_av=spec_phase_av[0:count]
    freq_av=freq_av[0:count]
    return spec_amp_av,spec_phase_av,freq_av,count    
#    
# main program
#
# define the input and output files
#
path_in='/Users/riser/Desktop/ocean.569A/datasets/OS_WHOTS_201606_D_MICROCAT-025m.nc'
path_out1='/Users/riser/Desktop/ocean.569A/HOT.coherence.jpg'
path_out2='/Users/riser/Desktop/ocean.569A/HOT.phase.jpg'
#
# read the input file (netcdf)
#
HOT_data=import_data(path_in)
#
# determine the length of the data and the resulting number of spectral estimates
#
nn=len(HOT_data['TIME'])
print ('nn=',nn)
mm=int(nn/2+1)
#
# define the data arrays
#
time_meas=np.zeros(nn)
time_meas_days=np.zeros(nn)
temp_meas=np.zeros(nn)
sal_meas=np.zeros(nn)
tt=np.zeros(nn)
ss=np.zeros(nn)
freq=np.zeros(mm)
#
# parse the time, temperature, and salinity data from the input file

for i in range(0,nn):
   time_meas[i]=24.*(float(HOT_data['TIME'][i])-float(HOT_data['TIME'][0]))
   time_meas_days[i]=(float(HOT_data['TIME'][i])-float(HOT_data['TIME'][0]))
   temp_meas[i]=float(HOT_data['TEMP'][i])
   sal_meas[i]=float(HOT_data['PSAL'][i])
#
# remove the temperature and salinity means from the data
#
tt=temp_meas-np.mean(temp_meas)
ss=sal_meas-np.mean(sal_meas)
# 
# determine the frequencies for the spectrum
#
delt=0.00208333
T_length=nn
pi=np.pi
omega0=2.*pi/(T_length*delt)
#
for i in range(0,mm):
    freq[i]=i*omega0
#
# compute the fft of the input data (temperature and salinity here)
# the fft will yield a set of complex numbers
# also compute their complex conjugates
# use the function band_average to get the spectral estimates
# estimate the coherence and phase using the results of the function
#
n_av=21
zz1=np.fft.rfft(tt,n=nn)
zz2=np.fft.rfft(ss,n=nn)
zz2_star=np.conj(zz2)
temp_spec,temp_phase,freq_av,count=band_average(zz1,zz1,freq,n_av)
salt_spec,salt_phase,freq_av,count=band_average(zz2,zz2,freq,n_av)
cospec_amp,cospec_phase,freq_av,count=band_average(zz1,zz2_star,freq,n_av)
coh_sq=cospec_amp**2/(temp_spec*salt_spec)
#
# begin plotting the coherence and phase
#
# first the coherence
#
fig1=plt.figure(figsize=(9,7))
plt.ylim(0.,1.)
plt.semilogx(freq_av,coh_sq,color='purple')
plt.xlabel('$\omega$ (radians/day)',fontsize=15,ha='center')
plt.ylabel('Squared Coherence $\it{T}$-$\it{S}$',fontsize=15)
freq_nyquist=pi/delt
freq_T=2.*pi/(nn*delt)
plt.plot([freq_nyquist,freq_nyquist],[0.,1.],'--k')
plt.text(8.e2,0.2,'$\omega_{max}$',fontsize=12,color='firebrick')
plt.text(1.2,0.88,'n_av = 21')
plt.grid(which='both')
plt.title('MBARI M1 Mooring')
plt.show()
plt.savefig(path_out1)
plt.show()
#
# now the phase
#
fig2=plt.figure(figsize=(9,7))
plt.ylim(-180.,180.)
plt.semilogx(freq_av,cospec_phase,color='orange')
plt.xlabel('$\omega$ (radians/day)',fontsize=15,ha='center')
plt.ylabel('Phase $\it{T}$-$\it{S}$, degrees',fontsize=15)
freq_nyquist=pi/delt
freq_T=2.*pi/(nn*delt)
plt.plot([freq_nyquist,freq_nyquist],[-180.,180.],'--k')
plt.text(8.e2,-110.,'$\omega_{max}$',fontsize=12,color='firebrick')
plt.text(1.2,130.,'n_av = 21')
plt.grid(which='both')
plt.title('MBARI M1 Mooring')
plt.savefig(path_out2)
plt.show()