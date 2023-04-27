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
yb=np.linspace(0,2.5,50)


emp_master_location='/Volumes/Choruh/Data/snowmelt_project/'
# emp_master_location='/Users/aforte/Documents/Python/snowmelt/'

## Load global
df_global=pd.read_csv(emp_master_location+'wrr2_derived_data_v3.csv')
df_global=df_global.drop(index=df_global.index[np.isnan(df_global['mean_z'])])
df_global=df_global.reset_index(drop=True)

# Calc percents
global_perc_base=df_global['qsb']/df_global['mean_runoff']

# Calculate indices
grlf=df_global['max_z']-df_global['min_z']

# Set cutoffs
percb_cutoff=0.25

rem_idx=(grlf<=500) | (df_global['mean_z']<=250) | (global_perc_base>percb_cutoff) | (df_global['mean_rlf']<250)
df_global_s=df_global.drop(df_global.index[rem_idx])
df_global_s=df_global_s.reset_index(drop=True)
global_perc_snow=df_global_s['qsm']/df_global_s['mean_runoff']

plt.figure(2)
idx1=(df_global_s['r_c1']>0)

[h,xed,yed,im]=plt.hist2d(df_global_s.loc[idx1,'mean_runoff'],df_global_s.loc[idx1,'r_c1'],[xb,yb],norm=colors.LogNorm(vmin=1,vmax=1000),cmap=cm.grayC)

X=xed[0:-1]+np.diff(xed)[0]
Y=yed[0:-1]+np.diff(yed)[0]



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


f1=plt.figure(1,figsize=(8,5))
f1.set_dpi(250)

letters=['A','B','C','D','E','F']

for i in range(len(model_list1U)):
    d_out=dict_list1U[i]
    r=d_out['mrout'].ravel()
    cr=d_out['crout'].ravel()
    
    ax1=plt.subplot(2,3,i+1)
    sc1=plt.hist2d(r,cr,[xb,yb],norm=colors.LogNorm(vmin=1,vmax=2500),cmap=cm.lajolla)
    ax1.contour(X,Y,h.transpose(),10,colors='k',linewidths=0.5,linestyles=':')
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
    ax1.text(0.01, 0.99, letters[i],
            horizontalalignment='left',
            verticalalignment='top',
            transform=ax1.transAxes,
            fontsize=12,fontweight='extra bold')


plt.tight_layout()
plt.rcdefaults()

f1.savefig('P2_SUP_figureS2.pdf',dpi="figure")
