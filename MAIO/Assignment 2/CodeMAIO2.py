#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 12 13:01:42 2017

@author: victoronink
code for the second MAIO assignment
"""

import numpy as np  
import matplotlib.pyplot as plt  
import cmath as cm
from scipy import signal as sig

def DiscreteFourier(series):
    # Function to carry out a discrete Fourier transformation of a timeseries
    N = len(series)
    Fourier = []
    FourierAbs = []
    PreFac=1/np.sqrt(2*np.pi)
    for m in range(N):
        Fou = 0
        for n in range(N):
            Fou += series[n] * cm.exp(- 2j * np.pi * m * n / N)
        Fourier.append(PreFac*Fou/N)
        FourierAbs.append(PreFac*abs(Fou)/N)
    Spectrum=np.abs(Fourier)**2
    Spectrum=Spectrum[0:len(Spectrum)/2]
    return Fourier, Spectrum, FourierAbs

def SeriesInt(series,step):
    # Function to carry out an integration of a time series
    [N,tot] = [len(series),0]
    for i in range(N):
        tot+=series[i]*step
    return tot

def DeltaTimeSeries(length,step,deltapoint):
    # Function to generate a timeseries with a delta function peak
    # all of length, step and deltapoint have to be in terms of kyr
    # length = length of the time series
    # step = time step
    # deltapoint = time point of the non-zero delta function value
    N=length/step
    series=np.zeros((int(N),2))
    for i in range(int(N)):
        series[i,0]=i*step
        if i == int(deltapoint/step):
            series[deltapoint/step,1]=N
    return series

def SinTimeSeries(length,step,amplitude,wavelength):
    # Function to generate a sinusoidal time series
    # everything in terms of kyr 
    # length = length of time series
    # step = time step
    # amplitude = amplitude of sinusoidal oscillation
    # wavelength = wavelength of sinusoidal oscillation
    N=length/step
    series=np.zeros((int(N),2))
    for i in range(int(N)):
        series[i,0]=i*step
        series[i,1]=amplitude*np.sin((2.*np.pi/wavelength)*(i*step))
    return series

def SawToothSeries(length,step,amplitude,wavelength,wobble,ran_amp):
    # Function to generate a sawtooth time series
    # everything in terms of kyr
    # length = length of time series
    # step = time step
    # amplitude = amplitude of the sawtooth function
    # wavelength = wavelength of sawtooth function
    # wobble= variation in the wavelength
    # ran_amp = standard deviation of random perturbation on amplitude
    N=length/step
    series=np.zeros((int(N),2))
    for i in range(int(N)):
        series[i,0]=i*step
        w=np.random.uniform(-wobble,wobble)
        if ran_amp==0:
            r=0
        else:
            r=np.random.normal(0,ran_amp)
        series[i,1]=(amplitude+r)* \
            sig.sawtooth(2.*np.pi/(wavelength+w)*(i*step))
    return series

#Computation
# part=0 => delta series generation
# part=1 => sin series generation
# part=2 => sawtooth series generation
part=input('Que? 0 is delta series, 1 is sin series and 2 is sawtooth series ')
plt.figure(1)
plt.xlabel('Time (kyr)')
plt.ylabel('Output')
plt.figure(2)
plt.xlabel('Frequency ($kyr^{-1}$)')
plt.ylabel('Amplitude')
plt.figure(3)
plt.xlabel('Wavenumber ($kyr^{-1}$)')
plt.ylabel('Amplitude')
plt.figure(4)
plt.xlabel('Frequency ($kyr^{-1}$)')
plt.ylabel('Spectral Density')
if part==0: 
    #Delta Series Computation
    DeltaSeries=DeltaTimeSeries(600,1,300)
    FourierSeries=DiscreteFourier(DeltaSeries[:,1])
    plt.figure(1)
    plt.plot(DeltaSeries[:,0],DeltaSeries[:,1])
    plt.title('Dirac-Delta Time Series')
    #plt.savefig("DiracDeltaFunction.pdf")
    plt.figure(2)
    plt.plot(DeltaSeries[:,0],FourierSeries[0])
    plt.title('Dirac-Delta Series Fourier Transform')
    #plt.savefig("DiracDeltaFourier.pdf")
    plt.figure(3)
    plt.plot(DeltaSeries[:,0],FourierSeries[2])
    plt.title('Dirac-Delta Series Fourier Transform Absolute')
    #plt.savefig("DiracDeltaFourierAbs.pdf")
    plt.figure(4)
    plt.plot(FourierSeries[1])
    plt.title('Dirac-Delta Time Series Power Spectrum')
    plt.xlim([0,0.05])
    #plt.savefig("DiracDeltaPowerSpec.pdf")
elif part==1:
    #Sin series computation
    SinSeries=SinTimeSeries(1000,1,1,100)
    FourierSeries=DiscreteFourier(SinSeries[:,1])
    SinInt=SeriesInt(FourierSeries[1],1)
    plt.figure(1)
    plt.plot(SinSeries[:,0],SinSeries[:,1])
    plt.title('Sinusoidal Time Series')
    #plt.savefig('SinSeriesFunction.pdf')
    plt.figure(2)
    plt.plot(SinSeries[:,0],FourierSeries[0])
    plt.title('Sinusoidal Time Series Fourier Transform')
    #plt.savefig('SinSeriesFourier.pdf')
    plt.figure(3)
    plt.plot(SinSeries[:,0],FourierSeries[2])
    plt.title('Sinusoidal Time Series Fourier Transform Absolute')
    #plt.savefig('SinSeriesFourierAbs.pdf')
    plt.figure(4)
    plt.plot(FourierSeries[1])
    plt.title('Sinusoidal Time Series Power Spectrum')
    #plt.savefig('SinSeriesPowerSpec.pdf')
elif part==2:
    #Sawtooth series computation
    SawtoothSeries=SawToothSeries(1000,1,1,100,20,1)
    FourierSeries=DiscreteFourier(SawtoothSeries[:,1])
    SawInt=SeriesInt(FourierSeries[1],1)
    plt.figure(1)
    plt.plot(SawtoothSeries[:,0],SawtoothSeries[:,1])
    plt.title('Sawtooth Time Series')
    #plt.savefig('SawtoothFunctionWobbleNoise.pdf')
    plt.figure(2)
    plt.plot(SawtoothSeries[:,0],FourierSeries[0])
    plt.title('Sawtooth Time Series Fourier Transform')
    #plt.savefig('SawtoothFourierWobbleNoise.pdf')
    plt.figure(3)
    plt.plot(SawtoothSeries[:,0],FourierSeries[2])
    plt.title('Sawtooth Time Series Fourier Transform Absolute')
    #plt.savefig('SawtoothFourierAbsWobbleNoise.pdf')
    plt.figure(4)
    plt.plot(FourierSeries[1])
    plt.title('Sawtooth Time Series Power Spectrum')
    #plt.savefig('SawtoothPowerSpecWobbleNoise.pdf')
else:
    print('Wrong answer, run again.')