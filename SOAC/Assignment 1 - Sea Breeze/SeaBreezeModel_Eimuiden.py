"""
Code for the modelling of the sea-breeze in Eimuiden 
on the 7th and 8th of May, 1976

"""
import numpy as np  # http://www.numpy.org
import matplotlib.pyplot as plt   # http://matplotlib.org
import math as M  # https://docs.python.org/2/library/math.html
import pandas as pd

# PARAMETER VALUES
dt = 360.0 # time step in s
tmax = 480 # max. number of time steps in an integration 1440/60; 360/240
omega = 0.000072792  # angular velocity Earth [s^-1]
lat = 52 # latitude of IJmuiden in degrees
pi = np.pi
phi = 0.0 # phase of the pressure gradient 
ro = 1.25 # density
A = 0.00065 # Pa/m
phase = 3*pi/2# phase of surface pressure gradient in time
f = 2 * omega *  M.sin(lat * pi/180)  # Coriolis parameter
C1= - (A/(f*ro)) * ( (M.pow(omega,2) /(M.pow(f,2) - M.pow(omega,2))) + 1)
C3 = A * omega /(ro*(M.pow(f,2) - M.pow(omega,2)))
lam = 0.00001 # Rayleigh damping coefficient

#  DEFINE time, u, v and the analytical solution, u_ana as arrays and fill them with zero's 
time = np.zeros((tmax), dtype='d')
time_axis = np.zeros((tmax), dtype='d')
u_kt = np.zeros((tmax), dtype='d')# x-component velocity
v_kt = np.zeros((tmax), dtype='d')# y-component velocity

# Eimuiden Conditions Import
grad=pd.read_excel('GradientData.xlsx')
dpdx=grad['dpdx'].values
v_g=1/ro*(dpdx/1000)
#v_g=np.zeros(len(dpdx))
dpdy=grad['dpdy'].values
u_g=-1/ro*(dpdy/1000)
#u_g=np.zeros(len(dpdy))
u_Ei=grad['u_o'].values

# INITIAL CONDITION (t=0) : atmosphere in rest
t = 0
time[t] = 0
time_axis[t] = 0
u_kt[t]=u_Ei[t]
v_kt[t]=0

# TIME LOOP EULER FORWARD SCHEME
for t in range(len(time)-1): 
 time[t+1] = time[t]+dt 
 #Now we try out our 4th order Runge-Kutta Scheme
 #for u[t+1]
 u_kt1= u_kt[t]+(-((A/ro)* M.cos(omega*time[t]+phase))+f*(v_kt[t]-v_g[t])-lam*u_kt[t])*(dt/2)
 v_kt1= v_kt[t]+(-f*(u_kt[t]-u_g[t])-lam*v_kt[t])*(dt/2)
 u_kt2= u_kt[t]+(-((A/ro)* M.cos(omega*time[t]+phase))+f*(v_kt1-v_g[t])-lam*u_kt1)*(dt/2)
 v_kt2= v_kt[t]+(-f*(u_kt1-u_g[t])-lam*v_kt1)*(dt/2)
 u_kt3= u_kt[t]+(-((A/ro)* M.cos(omega*time[t]+phase))+f*(v_kt2-v_g[t])-lam*u_kt2)*dt
 v_kt3= v_kt[t]+(-f*(u_kt1-u_g[t])-lam*v_kt2)*dt
 u_kt4= u_kt[t]-(-((A/ro)* M.cos(omega*time[t]+phase))+f*(v_kt3-v_g[t])-lam*u_kt3)*(dt/2)
 v_kt4= v_kt[t]-(-f*(u_kt3-u_g[t])-lam*v_kt3)*(dt/2)
 u_kt[t+1]=(u_kt1+2*u_kt2+u_kt3-u_kt4)/3
 v_kt[t+1]=(v_kt1+2*v_kt2+v_kt3-v_kt4)/3
 
for t in range(len(time)):
 time_axis[t] = time[t] / 3600.
 
# MAKE PLOT of evolution in time of u and u_ana
plt.plot(time_axis,u_kt, color='black')
plt.plot(time_axis[::10],u_Ei[::10],'.', color='red')
plt.axis([0,time_axis[tmax-1],-25.0,25.0])  # define axes 
plt.xticks(np.arange(0,time_axis[tmax-1],6), fontsize=12) 
plt.yticks(np.arange(-25.0,25.0,5), fontsize=12) 
plt.xlabel('time [hours]', fontsize=14) # label along x-axes
plt.ylabel('u [m/s]', fontsize=14) # label along x-axes
plt.title('Eimuiden Sea Breeze - 7 and 8 May, 1976') # Title at top of plot
plt.grid(True)
#plt.savefig("SeabreezeSimulation.png") # save plot as png-file
plt.show() # show plot on screen