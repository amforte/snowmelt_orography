#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 17 14:10:43 2023

@author: aforte
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from cmcrameri import cm
from scipy.stats import linregress
from matplotlib import colors
from scipy import odr
import rasterio
import matplotlib.gridspec as gridspec
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from scipy import optimize
from scipy.special import gamma


def linear(B,x):
    return B[0]*x + B[1]

def rmsef(observed,predicted):
    return np.sqrt(np.sum((predicted-observed)**2)/len(predicted))

def odr_fit_ref(x,y):
    # Filter 0 values
    lx=np.log10(x[(x>0) & (y>0)])
    ly=np.log10(y[(x>0) & (y>0)])
    linmod=odr.Model(linear)
    
    fd=odr.Data(x,y)
    odrlin=odr.ODR(fd,linmod,beta0=[0.1,0.1])
    outlin=odrlin.run()
    slp=outlin.beta[0]
    yint=outlin.beta[1]
    lin_pred=slp*x+yint
    lin_rmse=rmsef(y,lin_pred)
    lin_chi2=outlin.res_var
    lin_sd=outlin.sd_beta
    
    fdlog=odr.Data(lx,ly)
    odrlog=odr.ODR(fdlog,linmod,beta0=[0.1,10])
    outlog=odrlog.run()
    logexp=outlog.beta[0]
    logcoeff=10**outlog.beta[1]
    nlog_pred=logcoeff*x**logexp
    log_rmse=rmsef(y,nlog_pred)
    log_chi2=outlog.res_var
    log_sd=outlog.sd_beta

    return slp,yint,lin_sd,lin_rmse,lin_chi2,logcoeff,logexp,log_sd,log_rmse,log_chi2 

# master_location='/Volumes/Choruh/Data/snowmelt_project/'
master_location='/Users/aforte/Documents/Python/snowmelt/'
repo_location='/Users/aforte/Documents/GitHub/snowmelt-tectonics/stimpy/'

## Load global
df_global=pd.read_csv(master_location+'wrr2_derived_data_v4.csv')
df_global=df_global.drop(index=df_global.index[np.isnan(df_global['mean_z'])])
df_global=df_global.reset_index(drop=True)

# Calc percents

global_perc_base=df_global['qsb']/df_global['mean_runoff']

# Calculate indices
grlf=df_global['max_z']-df_global['min_z']

# Set color
temp_cutoff=2
vmin=temp_cutoff-34
vmax=temp_cutoff+34

# Set cutoffs
perc_cutoff=0.275
percb_cutoff=0.25


rem_idx=(grlf<=500) | (df_global['mean_z']<=250) | (global_perc_base>percb_cutoff) | (df_global['mean_rlf25']<250)
df_global_s=df_global.drop(df_global.index[rem_idx])
df_global_s=df_global_s.reset_index(drop=True)
global_perc_snow=df_global_s['qsm']/df_global_s['mean_runoff']
s_f=df_global_s['mean_runoff']/gamma(1+1/df_global_s['r_c1'])
# Set cnorms
cnorm=colors.Normalize(vmin=vmin,vmax=vmax)

# Set vector
r=np.linspace(0.1,25,100)

# Define model
linmod=odr.Model(linear)

## Calculate singular rainfall dominated relationship between mean runoff and cR
gidx=(global_perc_snow<=0.35) & (df_global_s['r_c1']>0)

param1=[]
param2=[]
r_type=[]
rel=[]

RFcr=odr.Data(np.log10(df_global_s.loc[gidx,'mean_runoff']),np.log10(df_global_s.loc[gidx,'r_c1']))
RFcr_odr=odr.ODR(RFcr,linmod,beta0=[0.1,10])
RFcr_outlog=RFcr_odr.run()
RFcr_logexp=RFcr_outlog.beta[0]
RFcr_logcoeff=10**RFcr_outlog.beta[1]
param1.append(RFcr_logcoeff)
param2.append(RFcr_logexp)
r_type.append('power')
rel.append('rbar_to_cR')

RFsr=odr.Data(s_f[gidx],df_global_s.loc[gidx,'r_s1'])
RFsr_odr=odr.ODR(RFsr,linmod,beta0=[0.1,10])
RFsr_out=RFsr_odr.run()
RFsr_m=RFsr_out.beta[0]
RFsr_b=RFsr_out.beta[1]
param1.append(RFsr_m)
param2.append(RFsr_b)
r_type.append('linear')
rel.append('sRm_to_sRf')

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

f1=plt.figure(1,figsize=(6,3.5),layout='constrained',dpi=300)
xb=np.linspace(0,10,50)
yb=np.linspace(0,1.5,50)

lab1=r'$c_{R}$ = '+str(np.round(RFcr_logcoeff,2))+r'$\bar{R}$ $^{'+str(np.round(RFcr_logexp,2))+'}$'
lab2=r'$X_{0f}$ = '+str(np.round(RFsr_m,2))+r'$X_{0m}$ + '+str(np.round(RFsr_b,2))

ax1=plt.subplot(1,2,1)
sc1=plt.hist2d(df_global_s.loc[gidx,'mean_runoff'],df_global_s.loc[gidx,'r_c1'],[xb,yb],norm=colors.LogNorm(vmin=1,vmax=500),cmap=cm.lajolla)
cbar1=plt.colorbar(sc1[3],ax=ax1)
cbar1.ax.set_ylabel('Density')
plt.plot(r,RFcr_logcoeff*r**RFcr_logexp,c='k',linewidth=2,label=lab1)
plt.xlabel('WaterGAP3 Runoff [mm/day]')
plt.ylabel('WaterGAP3 Shape Parameter')
plt.xlim((0,10))
plt.ylim((0,1.5))
plt.legend(loc='lower left',bbox_to_anchor=[0.1,-0.3])
ax1.text(0.01, 0.99, 'A',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax1.transAxes,
        fontsize=12,fontweight='extra bold')

ax2=plt.subplot(1,2,2)
sc1=plt.hist2d(s_f[gidx],df_global_s.loc[gidx,'r_s1'],[xb,xb],norm=colors.LogNorm(vmin=1,vmax=500),cmap=cm.berlin)
cbar1=plt.colorbar(sc1[3],ax=ax2)
cbar1.ax.set_ylabel('Density')
plt.plot(r,RFsr_m*r + RFsr_b,c='k',linewidth=2,label=lab2)
plt.xlabel('Scale Estimated from WaterGAP3 Mean')
plt.ylabel('WaterGAP3 Scale Parameter')
plt.xlim((0,10))
plt.ylim((0,10))
plt.legend(loc='lower left',bbox_to_anchor=[0.1,-0.3])
ax2.text(0.01, 0.99, 'B',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax2.transAxes,
        fontsize=12,fontweight='extra bold')

plt.rcdefaults()

f1.savefig('P2_SUP_figureS7.pdf')