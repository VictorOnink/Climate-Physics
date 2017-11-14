"""
Code for the modelling of the sea-breeze in Eimuiden 
on the 7th and 8th of May, 1976

"""
import numpy as np  # http://www.numpy.org
import matplotlib.pyplot as plt   # http://matplotlib.org
import math as M  # https://docs.python.org/2/library/math.html
import pandas as pd
from statsmodels.tools import eval_measures

#Functions
def gradient(time,A,phase,omega,B):
    grad=np.zeros(len(time))
    for i in range(len(time)):
        grad[i]+=A*np.cos(omega*time[i]+phase)+B
    return grad

# PARAMETER VALUES
dt = 360.0 # time step in s
tmax = 480 # max. number of time steps in an integration 1440/60; 360/240
omega = 0.000072792  # angular velocity Earth [s^-1]
lat = 52 # latitude of IJmuiden in degrees
pi = np.pi 
ro = 1.25 # density
A = 0.0035 # Pa/m
phase = 4*3600# phase of surface pressure gradient in time
f = 2 * omega *  M.sin(lat * pi/180)  # Coriolis parameter
C1= - (A/(f*ro)) * ( (M.pow(omega,2) /(M.pow(f,2) - M.pow(omega,2))) + 1)
C3 = A * omega /(ro*(M.pow(f,2) - M.pow(omega,2)))
lam = 0.0016 # Rayleigh damping coefficient

#  DEFINE time, u, v and the analytical solution, u_ana as arrays and fill them with zero's 
time = np.zeros((tmax), dtype='d')
time_axis = np.zeros((tmax), dtype='d')
u_kt = np.zeros((tmax), dtype='d')# x-component velocity
v_kt = np.zeros((tmax), dtype='d')# y-component velocity

# Eimuiden Conditions Import
grad=pd.read_excel('GradientData.xlsx')
meandpdx=grad['mean dpdx'].values/1000
dpdx=grad['dpdx'].values/1000
v_g=1/ro*(dpdx)
#v_g=np.zeros(len(dpdx))
dpdy=grad['dpdy'].values/1000
u_g=-1/ro*(dpdy)
#u_g=np.zeros(len(dpdy))
u_Ei=grad['u_o'].values
    
# INITIAL CONDITION (t=0) : atmosphere in rest
t = 0
time[t] = 0
time_axis[t] = 0
u_kt[t]=u_Ei[t]
v_kt[t]=-7.5849 #initial v_0 measured in Eimuiden

# Time loop 4th order Runge-Kutta Scheme
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

##Determination of the proper phase and amplitude
#gradien=gradient(time[:len(time)/2],0.0035,1.25*pi,omega,0)
#
#plt.figure(1)
#plt.plot(meandpdx[0:24],color='black')
#plt.plot(time_axis[:len(time_axis)/2],gradien,color='red')

## MAKE PLOT of evolution in time of u and u_ana
#plt.figure(2)
#plt.plot(time_axis,u_kt, color='black')
#plt.plot(time_axis[::10],u_Ei[::10],'.', color='red')
#plt.axis([0,time_axis[tmax-1],-10.0,10.0])  # define axes 
#plt.xticks(np.arange(0,time_axis[tmax-1],6), fontsize=12) 
#plt.yticks(np.arange(-10.0,10.0,2), fontsize=12) 
#plt.xlabel('time [hours]', fontsize=14) # label along x-axes
#plt.ylabel('u [m/s]', fontsize=14) # label along x-axes
#plt.title('IJmuiden Sea Breeze - 7 and 8 May, 1976') # Title at top of plot
#plt.grid(True)
##plt.savefig("EimuidenSeabreezelambda0_005.png") # save plot as png-file
#plt.show() # show plot on screen

cor=np.corrcoef(u_Ei[::10],u_kt[::10])
print('correlation is '+str(cor[0,1]))
u_EiCom=u_Ei[::10]-np.mean(u_Ei[::10])
u_ktCom=u_kt[::10]-np.mean(u_kt[::10])
CPRMS=eval_measures.rmse(u_EiCom, u_ktCom, axis=0)
print('CP-RMS is '+str(CPRMS))
dev=np.std(u_kt)
print('STD is '+str(dev))