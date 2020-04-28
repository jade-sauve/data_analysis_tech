# -*- coding: utf-8 -*-
# @Author: jadesauve
# @Date:   2020-04-15 14:34:59
# @Last Modified by:   jadesauve
# @Last Modified time: 2020-04-16 18:30:57
"""
Steve's code

"""

# code to read the MBARI 19 year T/S dataset and plot the T and S spectra
#
import numpy as np
import matplotlib.pyplot as plt
import random
#
nn=156885
mm=int(nn/2+1)
nsk=6
eps=0.002
#
# define the input paths
#
path_in='/Users/riser/Desktop/ocean.569A/datasets/MBARI.M1.txt'
#
# begin reading the data using readline
#
data_in=open(path_in,'r')
#
# define some arrays
#
time=np.zeros(nn)
month=np.zeros(nn)
tt=np.zeros(nn)
tt_nomean=np.zeros(nn)
ss=np.zeros(nn)
ss_nomean=np.zeros(nn)
freq=np.zeros(mm)
#
# skip the header stuff
#
for i in range(0,nsk):
    data_read=data_in.readline()
#
# initialize counters and begin reading one line at a time
# strip off the end-of-line character, then split the line into
# the things we want
# edit for bad or missing data; here we'll just do a quick fill
# that is good enough for demonstration purposes
#
count_tt=0
count_ss=0
#
for i in range(0,nn):
   data_read=data_in.readline()
   raw_data0=data_read.strip('\n')
   raw_data1=raw_data0.replace(',',' ')
   rr=raw_data1.split()
   month[i]=float(rr[0][0:2])
   time[i]=i/(24.*365)
   tt[i]=rr[6]
   ss[i]=rr[5]
   if tt[i] > 100:
       tt[i]=tt[77832]
   if ss[i] > 36. or ss[i] < 31.5:
       count_ss+=1
       ss[i]=ss[i-1]*(1.+eps*(random.random()-0.5))
#
data_in.close()
#
# remove the mean from temperature and salinity prior to estimating the spectra
#
tt_nomean=tt-np.mean(tt)
ss_nomean=ss-np.mean(ss)
#
fs=14
#
# define the spectral parameters
#
delt=1./24.
T_length=nn
pi=np.pi
omega0=2.*pi/(T_length*delt)
#
# define the frequency array
#
for i in range(0,mm):
    freq[i]=i*omega0
#
# compute the fft of the input data (temperature in this case)
# the fft will yield a set of complex numbers
# also compute their complex conjugates
# multiply these together to get the spectrum
#
zz_tt=np.fft.rfft(tt_nomean,n=nn)
fourier_amp=np.sqrt((np.real(zz_tt)**2+np.imag(zz_tt)**2))
fourier_phase=180.*np.arctan2(np.imag(zz_tt),np.real(zz_tt))/pi
spec_tt=np.real(zz_tt*np.conj(zz_tt))/(2.*pi*T_length*delt)
spec_tt_amp=(np.absolute(zz_tt))**2/(2.*pi*T_length*delt)
#
# plot the temperature spectrum as a function of frequency
# show the maximum and minimum frequencies
#
fig3=plt.figure(figsize=(7,7))
plt.ylim(0.01,6.e5)
plt.loglog(freq,spec_tt,color='darkkhaki')
plt.xlabel('$\omega$ (radians/day)',fontsize=15,ha='center')
plt.ylabel('Energy Density ($^o$C$^2$)/(radian/day))',fontsize=15)
freq_nyquist=pi/delt
freq_T=2.*pi/(nn*delt)
plt.plot([freq_nyquist,freq_nyquist],[0.01,6.e5],'--k')
plt.text(28.,5.e4,'$\omega_{max}$',fontsize=12,color='blue')
plt.plot([freq_T,freq_T],[0.01,6.e5],'--k',zorder=10)
plt.text(1.2e-3,0.3,'$\omega_o$',fontsize=12,color='blue')
#
# show the the M2 and K1 tide and add a grid
#
plt.text(5.,1.5e2,'K1',color='blue')
plt.text(10.,2.5e1,'M2',color='blue')
plt.text(5./365.,1.7e5,'1 yr',color='blue')
plt.grid()
plt.show()
#
# now do the same for the salinity spectrum
#
zz_ss=np.fft.rfft(ss_nomean,n=nn)
fourier_amp=np.sqrt((np.real(zz_ss)**2+np.imag(zz_ss)**2))
fourier_phase=180.*np.arctan2(np.imag(zz_ss),np.real(zz_ss))/pi
spec_ss=np.real(zz_ss*np.conj(zz_ss))/(2.*pi*T_length*delt)
spec_ss_amp=(np.absolute(zz_ss))**2/(2.*pi*T_length*delt)
#
# plot the salinity spectrum as a function of frequency
# show the maximum and minimum frequencies
#
fig4=plt.figure(figsize=(7,7))
plt.ylim(1.e-4,1.e4)
plt.loglog(freq,spec_ss,color='dodgerblue')
plt.xlabel('$\omega$ (radians/day)',fontsize=15,ha='center')
plt.ylabel('Energy Density (PSU$^2$)/(radian/day))',fontsize=15)
freq_nyquist=pi/delt
freq_T=2.*pi/(nn*delt)
plt.plot([freq_nyquist,freq_nyquist],[1.e-4,1.e4],'--k')
plt.text(28.,5.e3,'$\omega_{max}$',fontsize=12,color='purple')
plt.plot([freq_T,freq_T],[1.e-4,1.e4],'--k',zorder=10)
plt.text(1.2e-3,0.03,'$\omega_o$',fontsize=12,color='purple')
#
# show the the annual peak and add a grid
#
plt.text(5./365.,3.e3,'1 yr',color='purple')
plt.grid()
#plt.savefig(path_out4)
plt.show()
#