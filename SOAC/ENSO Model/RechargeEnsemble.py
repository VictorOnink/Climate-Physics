# -*- coding: utf-8 -*-
"""
Created on Thu Oct  5 09:43:47 2017

Ensembles over the Recharge oscillator

@author: Paranoid Android
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import signal

days=1 # days per month for sampling rate. 1= 1sample per month, 30= 1 sample per day

# Define model parameters
b_o=2.5     # high-end value for coupling parameter
mu=2./3      # coupling coefficient. Jin(1997) showed this gave stable oscillation
b=b_o*mu    # relates thermocline gradient to easterly wind stress
gamma=0.75  # specifies feedback of  thermocline gradient on sst gradient
c=1.         # Damping rate of sst anomalies
R=gamma*b-c # describes Bjerknes positive feedback process
r=0.25      # damping of upper ocean  heat content
alpha=0.125 # relates easterly wind stress to recharge of ocean heat content
epsilon=0 # varies  degree of nonlinerarity

# Wind forcing parameters
f_ann=0     # annual forcing parameter
f_ran=0     # random forcing parameter
tau=12. # Dimensional. Divide by tscale when implementing
taucor=1./30 # 1 day, and we been running in months till now

frequency=np.sqrt(3./32) # oscillation frequency w_c
period=2*np.pi/frequency # oscillation period tau_c

# Scaling factors for non-dimensionalisation
T=7.5 # K, SST anomaly scale
h=150. # m, thermocline depth scale
tscale=2. # months, time scale

# Runge Kutta fourth order numerical approximation of SST and Thermocline Depth. Re-dimenionalises the output variables
def RungeKutta(R,gamma,epsilon,b,r,alpha,f_ann,f_ran):
    w=np.random.uniform(-1,1,len(time)) # white noise generation
    for t in range(len(time)-1):
        instability=(Therm[t]+b*SST[t]) # parameter determining stability of solution
        E=f_ann*np.cos(2*np.pi*t/(tau/tscale))+f_ran*w[t]*(taucor/dt)
        x1t=SST[t]+(R*SST[t]+gamma*Therm[t]-epsilon*(Therm[t]+b*SST[t])**3+gamma*E)*(dt/2)
        x1h=Therm[t]+(-r*Therm[t]-alpha*b*SST[t]-alpha*E)*(dt/2)
        x2t=SST[t]+(R*x1t+gamma*x1h-epsilon*(x1h+b*x1t)**3+gamma*E)*(dt/2)
        x2h=Therm[t]+(-r*x1h-alpha*b*x1t-alpha*E)*(dt/2)
        x3t=SST[t]+(R*x2t+gamma*x2h-epsilon*(x2h+b*x2t)**3+gamma*E)*dt
        x3h=Therm[t]+(-r*x2h-alpha*b*x2t-alpha*E)*dt
        x4t=SST[t]-(R*x3t+gamma*x3h-epsilon*(x3h+b*x3t)**3+gamma*E)*(dt/2)
        x4h=Therm[t]-(-r*x3h-alpha*b*x3t-alpha*E)*(dt/2)
        SST[t+1]=(1/3)*(x1t+2*x2t+x3t-x4t)
        Therm[t+1]=(1/3)*(x1h+2*x2h+x3h-x4h)
    # Re-dimensionalification and output 
#    return [SST*T,Therm*h,instability]
    # Or not
    return [SST,Therm,instability]

#%% Task F
# Ensemble forecast. Essentially, loop over initial SST and Thermocline depths

# Arrangements of SST's and depths, need to be realistic
Temparray=np.arange(-2,2,.4) # semi-open interval [) (start,stop,step)
Deptharray=np.arange(-.5,.5,0.1)
#Temparray=[1.125]
#Deptharray=[0]

# Define a time regime
dt=1./days # months, now changed to a day
tmax=40*int(round(period)) # 5 periods
time=np.zeros((tmax), dtype='d')
timeaxis=np.linspace(0,tmax,tmax) # (start,stop,number)  

# Empty arrays for prognostic variables SST and Thermocline Depth
SST=np.zeros((tmax), dtype='d')
Therm=np.zeros((tmax,), dtype='d') 

temp3d=np.zeros((tmax,len(Temparray),len(Deptharray)),dtype='d')
dep3d=np.zeros((tmax,len(Temparray),len(Deptharray)),dtype='d')
temp2d=np.zeros((tmax,len(Temparray)),dtype='d')
dep2d=np.zeros((tmax,len(Temparray)),dtype='d')
det=[]

for kk in range(len(Deptharray)):
    for jj in range(len(Temparray)):

        t=0
        SST[t]=Temparray[jj]
        Therm[t]=Deptharray[kk]
        
        [SST,Therm,instability]=RungeKutta(R,gamma,epsilon,b,r,alpha,f_ann,f_ran)
        
        # Find unstable solutions, save only stables
        if np.isnan(SST).any()==True or np.isinf(SST).any()==True \
         or np.isnan(Therm).any()==True or np.isinf(Therm).any()==True:
            print('Solution unstable at initial Temp=%s, Depth=%s'%(Temparray[jj],Deptharray[kk]))
            print('(h+bT)=%s'%(instability))
        else:        
            temp2d[:,jj]=SST
            dep2d[:,jj]=Therm
            det=np.append(det,instability)
            
    temp3d[:,:,kk]=temp2d
    dep3d[:,:,kk]=dep2d

#plt.figure(1)
#plt.plot(det)
#plt.title('h+bT of stable solutions, dt = 1 week')
#%%
#for kk in range(len(Deptharray)):
    
#    if kk % 2 == 0:
#plt.figure()
#plt.plot(timeaxis,temp3d[:,0,:])
#%%
# Data on ENSO 3.4 monthly SST anomaly
data=pd.read_csv('ensosstanom.csv')
#print(list(data))
data=data[' YR   MON  TOTAL ClimAdjust ANOM '].apply(lambda x: pd.Series(x.split('  ')))
data.columns=['year','month','total','climadjust','anomaly','nothing'] # nothing column is unfortunate consequence of the split spacer
data['anomaly']=data[['anomaly','nothing']].fillna('').sum(axis=1) # because spacer was uneven, anomaly data split over 'anomaly' and 'nothing'
del data['nothing']

def freq_to_year(datafreqs):
    V = (1./datafreqs)/(60*60*24*365)
    return ["%.2f" % z for z in V]

# power spectrum to compare
fs=1./(60*60*24*30) # sampling frequency (Hz). Of data, once per month
datafreqs,datapowerspec=signal.periodogram(data['anomaly'],fs)
peakfreq=datafreqs[int(np.argmax(datapowerspec))]
peakperiod=(1./peakfreq)/(60*60*24*365) # years

plt.figure()
plt.plot(data['anomaly'],label='NINO 3.4 SST Anomaly')
plt.ylabel('SST anomaly (C)')
plt.xlabel('Months since Jan 1950')
plt.title('NINO 3.4 SST Anomalies since 1950')
plt.plot(timeaxis,temp3d[:,0,0],label='Simulation')
plt.legend()

print('data peak period=',peakperiod,'years')

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()
# Add some extra space for the second axis at the bottom
fig.subplots_adjust(bottom=0.2)
ax1.plot(datafreqs[0:len(datafreqs)/2]/(1*10**-7),datapowerspec[0:len(datapowerspec)/2]/(1*10**8),alpha=1)
ax1.set_xlabel('Frequency ($10^{-7}$ Hz)')
ax1.set_ylabel('Spectral Density')
ax1.set_title('Power Spectrum, NINO3.4 SST Anomalies')
new_tick_locations = np.array([0.1, 0.2, 0.3, 0.4, 0.5,0.6,0.7,0.8,0.9])
# Move twinned axis ticks and label from top to bottom
ax2.xaxis.set_ticks_position("bottom")
ax2.xaxis.set_label_position("bottom")
# Offset the twin axis below the host
ax2.spines["bottom"].set_position(("axes", -0.2))
# Turn on the frame for the twin axis, but then hide all 
# but the bottom spine
ax2.set_frame_on(True)
ax2.patch.set_visible(False)
for sp in ax2.spines.values():
    sp.set_visible(False)
ax2.spines["bottom"].set_visible(True)
ax2.set_xticks(new_tick_locations)
ax2.set_xticklabels(freq_to_year(new_tick_locations*10**-7))
ax2.set_xlabel('Period (years)')
#%%
plt.close('all')
simfs=1./(60*60*24*(30/days)) # sampling frequency (Hz). 
for i in range(len(Deptharray)):
    for j in range(len(Temparray)):               
        simfreqs,simpowerspec=signal.periodogram(temp3d[:,i,j],simfs)       

#        simpeakfreq=simfreqs[int(np.argmax(simpowerspec))]
#        simpeakperiod=(1/simpeakfreq)/(60*60*24*365) # years
#        print('Simulation peak period=',simpeakperiod,'years. T(i)=',Temparray[j],', h(i)=',Deptharray[i],' dt=',dt,' month')
        if i in (1,2):
            if j in (3,7):
                fig = plt.figure()
                ax1 = fig.add_subplot(111)
                ax2 = ax1.twiny()
                # Add some extra space for the second axis at the bottom
                fig.subplots_adjust(bottom=0.2)
                ax1.plot(simfreqs[0:len(datafreqs)/2]/(1*10**-7),simpowerspec[0:len(datapowerspec)/2]/(1*10**8),alpha=1)
                #ax1.plot(simfreqs[0:len(datafreqs)/2],simpowerspec[0:len(datapowerspec)/2]/(1*10**8),alpha=1)
                ax1.set_xlabel('Frequency ($10^{-7}$ Hz)')
                ax1.set_ylabel('Spectral Density')
                ax1.set_title('Power Spectrum, Simulation: depth=%s, T=%s'%(Deptharray[i],Temparray[j]))
#                plt.tick_params(axis='both', which='major', labelsize=20)
#                plt.tick_params(axis='both', which='minor', labelsize=20)
                #ax1.set_xlim([0, 1])
                new_tick_locations = np.array([0.1, 0.2, 0.3, 0.4, 0.5,0.6,0.7,0.8,0.9])
                # Move twinned axis ticks and label from top to bottom
                ax2.xaxis.set_ticks_position("bottom")
                ax2.xaxis.set_label_position("bottom")
                # Offset the twin axis below the host
                ax2.spines["bottom"].set_position(("axes", -0.2))
                # Turn on the frame for the twin axis, but then hide all 
                # but the bottom spine
                ax2.set_frame_on(True)
                ax2.patch.set_visible(False)
                for sp in ax2.spines.values():
                    sp.set_visible(False)
                ax2.spines["bottom"].set_visible(True)
                ax2.set_xticks(new_tick_locations)
                ax2.set_xticklabels(freq_to_year(new_tick_locations*10**-7))
                ax2.set_xlabel('Period (years)')
#                plt.tick_params(axis='both', which='major', labelsize=20)
#                plt.tick_params(axis='both', which='minor', labelsize=20)
#                ax1.xaxis.set_tick_params(labelsize=20)
#                ax1.yaxis.set_tick_params(labelsize=20)
                
                plt.figure()
                plt.plot(data['anomaly'],label='NINO 3.4 SST Anomaly')
                plt.ylabel('SST anomaly (C)')
                plt.xlabel('Months since Jan 1950')
                plt.title('SST Anomalies and simulation:depth=%s, T=%s'%(Deptharray[i],Temparray[j]))
                plt.plot(timeaxis,temp3d[:,i,j],label='Simulation')
                plt.legend()
#                plt.tick_params(axis='both', which='major', labelsize=20)
#                plt.tick_params(axis='both', which='minor', labelsize=20)

plt.figure()
plt.plot(SST[0],Therm[0],'rx') # show the starting point
plt.plot(SST,Therm)
plt.xlabel('SST (C)')
plt.ylabel('Thermocline Depth (m)')
plt.title('Solution Trajectory')

#%%

#plt.figure()
#plt.ylabel('SST anomaly (C)',fontsize=20)
#plt.xlabel('Months since Jan 1950',fontsize=20)
#plt.title('Initial SST Anomaly Ensemble',fontsize=20)
#plt.plot(timeaxis,temp3d[:,:,0]*T)
##plt.legend()
#plt.ylim([-10,10])
#plt.tick_params(axis='both', which='major', labelsize=20)
#plt.tick_params(axis='both', which='minor', labelsize=20)
#%%
#sst_totalminmean=pd.to_numeric(data['total'])-np.mean(pd.to_numeric(data['total']))
#
#plt.figure()
#plt.plot(sst_totalminmean,'r',label='Total-mean')
#plt.plot(data['anomaly'],'b',label='NINO 3.4 SST Anomaly')
#plt.ylabel('SST anomaly (C)')
#plt.xlabel('Months since Jan 1950')
#plt.title('NINO 3.4 SST Anomalies since 1950')
##plt.plot(timeaxis,temp3d[:,0,0],label='Simulation')
#plt.legend()