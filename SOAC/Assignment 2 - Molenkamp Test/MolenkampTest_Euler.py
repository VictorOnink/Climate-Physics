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
import cmath as cm

#Set Parameters
length=2500 # kilometers
dx=25 # kilometers
N=int(length/dx) # number of spatial grids
u=1e-2 # km s^-1
time_max=int(length/u)#2.5e4 # seconds
dt= 2e2# seconds
T= int(time_max/dt) # number of time steps
pi=np.pi

Cour=u*dt/dx

#Initial starting arrangement
t=0
x=np.arange(0,length,dx)
C_eul=np.zeros((T,N))
C_mats=np.zeros((T,N))
C_MW=np.zeros((T,N))
C_S=np.zeros((T,N))
C_RK=np.zeros((T,N))
RK1=np.zeros((T,N))
RK2=np.zeros((T,N))
RK3=np.zeros((T,N))
RK4=np.zeros((T,N))

for i in range(N):
    if 1125<=i*dx: 
        if i*dx<=1375:
            C_eul[0,i]=1
            C_mats[0,i]=1
            C_MW[0,i]=1
            C_S[0,i]=1
            

       
#loop for computing the Forward Euler approach
for t in range(0,T-1):
    for i in range(N):
        C_eul[t+1,i]=C_eul[t,i]-Cour*(C_eul[t,i]-C_eul[t,i-1])

#loop for Matsuno approach
for t in range(0,T-1):
    for i in range(N):
        C_mats[t+1,i]=C_mats[t,i]-Cour*(C_mats[t,i]-C_mats[t,i-1])
    for j in range(N):
        C_mats[t+1,j]=C_mats[t,j]-Cour*(C_mats[t+1,j]-C_mats[t+1,j-1])

#loop for the Lax-Wendroff Scheme
for t in range(0,T-1):
    for i in range(N):
        C_MW[t+1,i]=(C_MW[t,i]-Cour*(C_MW[t,(i+1)%N]-C_MW[t,i-1])/2+(Cour**2)/2
        *(C_MW[t,(i+1)%N]-2*C_MW[t,i]+C_MW[t,i-1]))

#loop for the spectral method
a_k=np.ones((T,N/2))*(0+0j)
b_k=np.ones((T,N/2))*(0+0j)
c_k=np.ones((T,N/2))*(0+0j)
for t in range(0,T-1):
    for k in range(N/2):
        if t==0:
            if k==0:
                a_k[t,k]=1./N*np.sum(C_S[t,:])
            else:
                for j in range(N):
                    a_k[t,k]+=2./N*C_S[t,j]*np.cos(-2*pi*k*j/N)
                    b_k[t,k]+=2./N*C_S[t,j]*np.sin(-2*pi*k*j/N)
            c_k[t,k]=a_k[t,k]+1j*b_k[t,k]
            a_k[t+1,k]=a_k[t,k]+dt*(2*pi*u*k/length)*b_k[t,k]
            b_k[t+1,k]=b_k[t,k]-dt*(2*pi*u*k/length)*a_k[t,k]
            a_k[t+1,k]=a_k[t,k]+dt*(2*pi*u*k/length)*b_k[t+1,k]
            b_k[t+1,k]=b_k[t,k]-dt*(2*pi*u*k/length)*a_k[t+1,k]
            c_k[t+1,k]=a_k[t+1,k]+1j*b_k[t+1,k]            
        else:
            a_k[t+1,k]=a_k[t,k]+dt*(2*pi*u*k/length)*b_k[t,k]
            b_k[t+1,k]=b_k[t,k]-dt*(2*pi*u*k/length)*a_k[t,k]
            a_k[t+1,k]=a_k[t,k]+dt*(2*pi*u*k/length)*b_k[t+1,k]
            b_k[t+1,k]=b_k[t,k]-dt*(2*pi*u*k/length)*a_k[t+1,k]
            c_k[t+1,k]=a_k[t+1,k]+1j*b_k[t+1,k]
for t in range(0,T):
    for xx in range(N):
        for kk in range(N/2):
            C_S[t,xx]+=np.real(c_k[t,kk]*cm.exp(2j*pi*kk*dx*xx/length))


#%%
#plot approaches
plt.figure(1)
plt.plot(x,C_eul[0,:],'b',label='Initial')
plt.plot(x,C_eul[-1,:],'r',label='Forward-Euler')
plt.plot(x,C_mats[-1,:],'k',label='Matsuno')
plt.plot(x,C_MW[-1,:],'g',label='Lax-Wendroff')
plt.plot(x,C_S[-1,:],color='orange',label='Spectral')
plt.legend()
plt.title('One-Dimensional Molenkamp Test')
plt.xlabel('Distance (km)')
plt.ylabel('Concentration')
#plt.savefig('MolenkampTestdx25dt10.png')