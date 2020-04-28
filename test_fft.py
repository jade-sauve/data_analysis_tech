import numpy as np
import matplotlib.pyplot as plt
#
# test the fourier transform and inverse transform in the numpy.fft library
# define a simple function, compute and plot its fft,
# then compute and plot the inverse fft of the fft
# this is a good way to test the fft normalization
#
pathout1='/Users/riser/Desktop/ocean.569A/fft.test.out.jpg'
pathout2='/Users/riser/Desktop/ocean.569A/inverse.fft.test.out.jpg'
#
# set the number of points and number of frequencies
#
nn=1000
mm=int(nn/2+1)
#
# define the arrays
#
tt0=np.zeros(nn)
tt1=np.zeros(nn)
yy0=np.zeros(nn)
yy1=np.zeros(nn)
zz=np.zeros(mm)
freq_rad=np.zeros(mm)
freq=np.zeros(mm)
#
# define some constants
#
pi=np.pi
amp=1.0
n_fac=50.
omega0=2.*pi/nn
#
# set the time and amplitude coordinates
# 
for i in range(0,nn):
    tt0[i]=i # t
    yy0[i]=amp*np.sin(2.*pi*tt0[i]/n_fac)    #y_sin
#
# set the frequency values
#
for i in range(0,mm):
    freq_rad[i]=i*omega0
    freq[i]=freq_rad[i]/(2.*pi)
#
# carry out the fourier transform
# convert from real-imaginary parts to amplitude and phase
#
zz=np.fft.rfft(yy0,n=nn)/nn
zz_real=zz.real
zz_imag=zz.imag
zz_mag=np.sqrt(zz_real**2+zz_imag**2)
zz_phase=np.arctan(zz_imag/zz_real)
#
# plot the magnitude of the transform as a function of frequency
#
sz=15
fig=plt.figure(figsize=(10,5))
plt.xlim([0,0.5])
plt.ylim([0,1])
plt.plot(freq_rad,zz_mag,color='firebrick')
plt.plot([0,pi],[0.5,0.5],'--k')
plt.plot([0.12566,0.12566],[0,1],'--k')
plt.xlabel('$\omega$ (cycles)  (continues to $\pi$)',fontsize=sz)
plt.ylabel('Fourier coefficient magnitude (FFT/$\it{n}$)',fontsize=0.85*sz)
plt.text(0.132,0.55,'$\omega$=0.12566, $\it{f}$=0.02')
plt.grid(color='blue')
plt.title('FFT of sin(2$\pi$$\it{t}$/50)')
#plt.savefig(pathout1)
plt.show()
#
# now test the inverse transform
#
yy1=nn*np.fft.irfft(zz,n=nn)
#
fig=plt.figure(figsize=(10,5))
plt.xlim([0,1000])
plt.ylim([-1,1])
plt.plot(yy1)
plt.plot([0,1000],[0,0],'--m')
per=50.
plt.plot([per,per],[-1,1],'--k')
plt.xlabel('Time',fontsize=sz)
plt.ylabel('Amplitude',fontsize=sz)
plt.text(35,-1.12,'50',fontsize=15)
plt.title('Inverse FFT of [FFT of sin(2$\pi$$\it{t}$/50)]')
#plt.savefig(pathout2)
plt.show()