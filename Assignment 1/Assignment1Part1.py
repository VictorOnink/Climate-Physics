#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  6 11:06:01 2017

@author: victoronink
First MAIO assignment for the analysis of AR and MA processes
"""

import numpy as np
import matplotlib.pyplot as plt

# Define the constants of the AR(2) process
alpha_1=[0]
alpha_2=[0.75]
points=100

#Function for determining the autocorrelation for a given lag tau
def autocorrelation(tau,data):
    mu=np.mean(data)
    var=np.var(data)
    acor=0
    A=1./(len(data)*var)
    for i in range(len(data)-tau):
        acor=acor+(data[i]-mu)*(data[i+tau]-mu)
    return A*acor

#set the length 100 vector for the generated data for each of the three runs
data=np.ones((points,len(alpha_1)))
#set vector for the autocovariance
autocor=np.zeros((len(data[:,0])-1,len(alpha_1)))
#white noise generator
noise_gen=1 #if zero, then noise will not be generated, if 1, then it will be
mean_wn=0
std_wn=1
if noise_gen==0:
    wnoise=np.zeros((points,len(alpha_1)))
elif noise_gen==1:
    wnoise=np.random.normal(mean_wn, std_wn,points)

#loops to generate the data based on the given values of alpha according to
#the AR(2) scheme
for i in range(len(alpha_1)):
    for j in range(len(data[:,i])):
        if j==0:
            data[0,i]=data[0,i]+wnoise[0]
        elif j==1:
            data[1,i]=alpha_1[i]*data[0,i]+wnoise[1]
        elif j>1:
            data[j,i]=alpha_1[i]*data[j-1,i]+alpha_2[i]*data[j-2,i]+wnoise[j]
    
    for tau in range(len(autocor[:,i])):
        autocor[tau,i]=autocorrelation(tau,data[:,i])

#time for plotting
plt.figure(1)
plt.plot(data)
plt.title('AR(2) Function Output')

plt.figure(2)
plt.plot(autocor)
plt.title('AR(2) Autocorrelation')