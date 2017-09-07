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

# import all the datasets
data_a=pd.read_fwf('data_a.txt')
data_b=pd.read_fwf('data_b.txt')
data_c=pd.read_fwf('data_c.txt')
data_d=pd.read_fwf('data_d.txt')
data_e=pd.read_fwf('data_e.txt')

plt.plot(data_e)