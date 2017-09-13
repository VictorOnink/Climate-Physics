    # Exercise 1.14 Atmospheric Dynamics: 
# http://www.staff.science.uu.nl/~delde102/AtmosphericDynamics2017aCh1.pdf

# https://docs.python.org/2.7/
# IMPORT MODULES
import numpy as np  # http://www.numpy.org
import matplotlib.pyplot as plt   # http://matplotlib.org
import math as M  # https://docs.python.org/2/library/math.html

# PARAMETER VALUES
dt = 360.0 # time step in s
tmax = 480 # max. number of time steps in an integration 1440/60; 360/240
omega = 0.000072792  # angular velocity Earth [s^-1]
lat = 52 # latitude of IJmuiden in degrees
pi = np.pi
phi = 0.0 # phase of the pressure gradient 
ro = 1.25 # density
A = 0.001 # Pa/m
phase = 0.0 # phase of surface pressure gradient in time
fcor = 2 * omega *  M.sin(lat * pi/180)  # Coriolis parameter
C1= - (A/(fcor*ro)) * ( (M.pow(omega,2) /(M.pow(fcor,2) - M.pow(omega,2))) + 1)
C3 = A * omega /(ro*(M.pow(fcor,2) - M.pow(omega,2)))

#  DEFINE time, u, v and the analytical solution, u_ana as arrays and fill them with zero's 
time = np.zeros((tmax), dtype='d')
time_axis = np.zeros((tmax), dtype='d')
u = np.zeros((tmax), dtype='d')# x-component velocity
v = np.zeros((tmax), dtype='d')# y-component velocity
u_ana = np.zeros((tmax), dtype='d')# analytical solution x-component velocity
u_kt = np.zeros((tmax), dtype='d')# x-component velocity
v_kt = np.zeros((tmax), dtype='d')# y-component velocity

# INITIAL CONDITION (t=0) : atmosphere in rest
t = 0
time[t] = 0
time_axis[t] = 0
u[t] = 0
v[t] = 0
u_ana[t] = 0
u_kt[t]=0
v_kt[t]=0

# TIME LOOP EULER FORWARD SCHEME
for t in range(len(time)-1): 
 du = dt * ((fcor*v[t]) - ((A/ro)* M.cos((omega*time[t])+phase)))
 dv = -dt * fcor * u[t]
 time[t+1] = time[t]+dt 
 u[t+1] = u[t] + du
 v[t+1] = v[t] + dv	
 u_ana[t+1] = (C1 * M.sin(fcor * time[t+1])) + ( C3* M.sin((omega * time[t+1]) + phase) )
 #Now we try out our 4th order Runge-Kutta Scheme
 #for u[t+1]
 u_kt1= u_kt[t]+(-((A/ro)* M.cos(omega*time[t]+phase))+fcor*v_kt[t])*(dt/2)
 v_kt1= v_kt[t]+(-fcor*u_kt[t])*(dt/2)
 u_kt2= u_kt[t]+(-((A/ro)* M.cos(omega*time[t]+phase))+fcor*v_kt1)*(dt/2)
 v_kt2= v_kt[t]+(-fcor*u_kt1)*(dt/2)
 u_kt3= u_kt[t]+(-((A/ro)* M.cos(omega*time[t]+phase))+fcor*v_kt2)*dt
 v_kt3= v_kt[t]+(-fcor*u_kt1)*dt
 u_kt4= u_kt[t]-(-((A/ro)* M.cos(omega*time[t]+phase))+fcor*v_kt3)*(dt/2)
 v_kt4= v_kt[t]-(-fcor*u_kt3)*(dt/2)
 u_kt[t+1]=(u_kt1+2*u_kt2+u_kt3-u_kt4)/3
 v_kt[t+1]=(v_kt1+2*v_kt2+v_kt3-v_kt4)/3
 
for t in range(len(time)):
 time_axis[t] = time[t] / 3600.
 
# MAKE PLOT of evolution in time of u and u_ana
plt.plot(time_axis, u_ana, color='black')
plt.plot(time_axis, u, color='red')
plt.plot(time_axis,u_kt,'-.', color='orange')
plt.axis([0,time_axis[tmax-1],-25.0,25.0])  # define axes 
plt.xticks(np.arange(0,time_axis[tmax-1],6), fontsize=12) 
plt.yticks(np.arange(-25.0,25.0,5), fontsize=12) 
plt.xlabel('time [hours]', fontsize=14) # label along x-axes
plt.ylabel('u [m/s]', fontsize=14) # label along x-axes
plt.title('Exercise 1.14 DYME') # Title at top of plot
plt.text(0.2, 23, 'u (Analytic): black line', fontsize=10, color='black')
plt.text(0.2, 21, 'u (Forward Euler): red line', fontsize=10, color='red')
plt.text(0.2, 19, 'u (4th order Runge-Kutta): orange line', fontsize=10, color='orange')
plt.grid(True)
plt.savefig("SeabreezeSimulation.png") # save plot as png-file
plt.show() # show plot on screen