#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  8 13:42:31 2023

@author: aforte
"""

import sys
sys.path.insert(0,'/Users/aforte/Documents/GitHub/snowmelt-tectonics')
import stimpy as st
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from cmcrameri import cm
from matplotlib import colors

master_location='/Volumes/Choruh/Data/snowmelt_project/model_outputs/'
model_list1U=['gc025u','gc05u','gc1u','gc2u','gc4u','gc8u']
title_list=['GC025U','GC05U','GC1U','GC2U','GC4U','GC8U']
dict_list1U=[]

for i in range(len(model_list1U)):
    mObj=st.Stim1D(os.path.join(master_location,model_list1U[i]))
    d=mObj.parse_results(0,np.inf,True,0.1)
    dict_list1U.append(d)

xb=np.linspace(0,12,50)
yb=np.linspace(0,2,50)

# Generate supplemental figure
SMALL_SIZE = 8
MEDIUM_SIZE = 10
BIGGER_SIZE = 12
plt.rc('font', size=SMALL_SIZE,family='Futura')          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


f1=plt.figure(figsize=(8,5))
f1.set_dpi(250)

for i in range(len(model_list1U)):
    d_out=dict_list1U[i]
    r=d_out['mrout'].ravel()
    cr=d_out['crout'].ravel()
    
    ax1=plt.subplot(2,3,i+1)
    sc1=plt.hist2d(r,cr,[xb,yb],norm=colors.LogNorm(vmin=1,vmax=2500),cmap=cm.lajolla)
    if i>2:
        ax1.set_xlabel('Runoff [mm/day]')
    if (i==0) | (i==3):
        ax1.set_ylabel('Shape Parameter')
    ax1.set_xlim((0,12))
    ax1.set_ylim((0,2.5))
    cbar1=plt.colorbar(sc1[3],ax=ax1)
    if (i==2) | (i==5):
        cbar1.ax.set_ylabel('Density')
    ax1.set_title(title_list[i])


plt.tight_layout()
plt.rcdefaults()
