#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 25 14:10:13 2017

@author: victoronink

Numerical Simulation of the Molenkamp Test with:
    a) Forward Euler Scheme
    b) Two Step Matsuno Scheme
    c) Lax-Wendroff Scheme
    d) Spectral Method
"""

import numpy as np
import matplotlib.pyplot as plt

#Set Parameters
time_max=2.5e4 # seconds
dt= 2e2# seconds
T= int(time_max/dt) # number of time steps
length=2.5e3 # kilometers
dx=2.5e1 # kilometers
N=int(length/dx) # number of spatial grids
u=1e-2 # km s^-1

#Initial starting arrangement
t=0
x=np.arange(0,length,dx)
conc=np.zeros((T,N))
for i in range(N):
    if 1125<=i*dx: 
        if i*dx<=1375:
            conc[0,i]=1
        
#loop for computing the Forward Euler approach
for t in range(0,T-1):
    for i in range(N):
        conc[t+1,i]=conc[t,i]-u*dt*(conc[t,i]-conc[t,i-1])/dx
    for i in range(N):
        conc[t+1,i]=conc[t,i]-u*dt*(conc[t+1,i]-conc[t+1,i-1])/dx

#plot Forward-Euler approach
plt.figure(1)
plt.plot(x,conc[0,:],x,conc[-1,:])
plt.title('Matsuno Scheme')

    