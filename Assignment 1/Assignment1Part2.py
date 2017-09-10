#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  7 14:19:53 2017

@author: victoronink
Analysis and coding coefficient determination data sets
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from statsmodels.graphics.tsaplots import plot_acf
from statsmodels.graphics.tsaplots import plot_pacf

import statsmodels.formula.api as smf
import statsmodels.tsa.api as smt
import statsmodels.api as sm
import scipy.stats as scs

# import all the datasets
data_a=pd.read_fwf('data_a.txt')
data_b=pd.read_fwf('data_b.txt')
data_c=pd.read_fwf('data_c.txt')
data_d=pd.read_fwf('data_d.txt')
data_e=pd.read_fwf('data_e.txt')

#determine which of the data sets you want to analyse at each point
# 1= a, 2= b, 3= c, 4= d, 5= e
# also, we detrend all the data to get rid of linear trends
dataset=5
if dataset==1:
    data=data_a.values
    data=statsmodels.tsa.tsatools.detrend(data, order=1,axis=0)
elif dataset==2:
    data=data_b.values
    data=statsmodels.tsa.tsatools.detrend(data, order=1,axis=0)
elif dataset==3:
    data=data_c.values
    data=statsmodels.tsa.tsatools.detrend(data, order=1,axis=0)
elif dataset==4:
    data=data_d.values
    data=statsmodels.tsa.tsatools.detrend(data, order=1,axis=0)
elif dataset==5:
    data=data_e.values
    data=statsmodels.tsa.tsatools.detrend(data, order=1,axis=0)
    
#Plot the autocorrelation functions for the selected dataset
plt.figure(1)
plt.plot(data)
plt.title('data_e')
plt.savefig('data_e.pdf')

plt.figure(2)
acor_data=plot_acf(data)
plt.title('Autocorrelation Function for data_e')
plt.xlabel('Lag')
plt.ylabel('Correlation')
plt.xlim([-0.5,10])
plt.savefig('data_e_autocor.pdf')
plt.show(acor_data)

plt.figure(3)
pacor_data=plot_pacf(data, lags=100)
plt.title('Partial Autocorrelation Function for data_e')
plt.xlabel('Lag')
plt.ylabel('Correlation')
plt.xlim([-0.5,20])
plt.savefig('data_e_partialcor.pdf')
plt.show(pacor_data)

#Best fit section, determine the coefficient values
max_lag = 30
mdl = smt.ARMA(data, order=(2,0)).fit(
        maxlag=max_lag, method='css', trend='nc')
print mdl.summary() 