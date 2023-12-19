#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 14:32:04 2023

@author: aforte
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from cmcrameri import cm
from scipy import odr


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

def plot_and_fit_ref(ax2,ax4,dfs,bl,br,bb,bt,ttl_text):
    
    locat=[]
    descrip=[]
    param1=[]
    param2=[]
    
    spidx=(dfs['latitude']>=bb) & (dfs['latitude']<=bt) & (dfs['longitude']>=bl) & (dfs['longitude']<=br)
    
    
    max_z=np.max(np.ceil(dfs.loc[spidx,'max_z'].to_numpy()))
    max_r=np.max(np.ceil(dfs.loc[spidx,'mean_runoff'].to_numpy()))
    max_rlf=np.max(np.ceil(dfs.loc[spidx,'mean_rlf25'].to_numpy()))
    
    z=np.linspace(1,max_z,100)
    
    x=dfs.loc[spidx,'mean_rlf25'].to_numpy()
    y=dfs.loc[spidx,'mean_runoff'].to_numpy()
    s,yi,sd,lirmse,lichi2,co,ex,losd,lormse,lochi2=odr_fit_ref(x,y)
    bin_edges=np.histogram_bin_edges(dfs.loc[spidx,'mean_runoff'],'doane')
    bix=np.digitize(dfs.loc[spidx,'mean_runoff'],bin_edges)
    for i in range(len(bin_edges)-1):
        mx=np.mean(x[bix==i])
        stdx=np.std(x[bix==i])
        my=np.mean(y[bix==i])
        stdy=np.std(y[bix==i])
        ax2.scatter(mx,my,c='k',s=len(x[bix==i]))
        ax2.errorbar(mx,my,xerr=stdx,yerr=stdy,ecolor='k',elinewidth=0.5,zorder=0)
    # s,yi,sd,lirmse,lichi2,co,ex,losd,lormse,lochi2=odr_fit_ref(mx,my)
    # ax2.scatter(x,y,c='k',alpha=0.5,s=2,zorder=0)
    ax2.plot(z,co*z**ex,c='k',label='Local Relief - RMSE='+str(np.round(lormse,3)))   
    locat.append(ttl_text)
    descrip.append('rlf to mR')
    param1.append(co)
    param2.append(ex)
    
    x=dfs.loc[spidx,'max_z'].to_numpy()
    y=dfs.loc[spidx,'mean_runoff'].to_numpy()
    s,yi,sd,lirmse,lichi2,co,ex,losd,lormse,lochi2=odr_fit_ref(x,y)
    bin_edges=np.histogram_bin_edges(dfs.loc[spidx,'mean_runoff'],'doane')
    bix=np.digitize(dfs.loc[spidx,'mean_runoff'],bin_edges)
    for i in range(len(bin_edges)-1):
        mx=np.mean(x[bix==i])
        stdx=np.std(x[bix==i])
        my=np.mean(y[bix==i])
        stdy=np.std(y[bix==i])
        ax2.scatter(mx,my,c='r',s=len(x[bix==i]))
        ax2.errorbar(mx,my,xerr=stdx,yerr=stdy,ecolor='r',elinewidth=0.5,zorder=0)
    # s,yi,sd,lirmse,lichi2,co,ex,losd,lormse,lochi2=odr_fit_ref(mx,my)
    # ax2.scatter(x,y,c='r',alpha=0.5,s=2,zorder=0)
    ax2.plot(z,co*z**ex,c='r',label='Max Z - RMSE='+str(np.round(lormse,3)))
    
    x=dfs.loc[spidx,'mean_z'].to_numpy()
    y=dfs.loc[spidx,'mean_runoff'].to_numpy()
    s,yi,sd,lirmse,lichi2,co,ex,losd,lormse,lochi2=odr_fit_ref(x,y)
    bin_edges=np.histogram_bin_edges(dfs.loc[spidx,'mean_runoff'],'doane')
    bix=np.digitize(dfs.loc[spidx,'mean_runoff'],bin_edges)
    for i in range(len(bin_edges)-1):
        mx=np.mean(x[bix==i])
        stdx=np.std(x[bix==i])
        my=np.mean(y[bix==i])
        stdy=np.std(y[bix==i])
        ax2.scatter(mx,my,c='b',s=len(x[bix==i]))
        ax2.errorbar(mx,my,xerr=stdx,yerr=stdy,ecolor='b',elinewidth=0.5,zorder=0)
    # s,yi,sd,lirmse,lichi2,co,ex,losd,lormse,lochi2=odr_fit_ref(mx,my)
    # ax2.scatter(x,y,c='b',alpha=0.5,s=2,zorder=0)
    ax2.plot(z,co*z**ex,c='b',label='Mean Z - RMSE='+str(np.round(lormse,3)))    
    
    ax2.set_xlabel('Elevation [m]')
    ax2.set_ylabel('Mean Runoff [mm/day]')
    ax2.legend(loc='best')
    plt.ylim((0,max_r))
    plt.xlim((0,max_z))
    ax2.set_title(ttl_text)
    
    
    snmp=dfs.loc[spidx,'qsm']/dfs.loc[spidx,'mean_runoff']
    
    x=dfs.loc[spidx,'mean_rlf25'].to_numpy()
    y=snmp.to_numpy()
    s,yi,sd,lirmse,lichi2,co,ex,losd,lormse,lochi2=odr_fit_ref(x,y)
    bin_edges=np.histogram_bin_edges(snmp,'doane')
    bix=np.digitize(snmp,bin_edges)
    for i in range(len(bin_edges)-1):
        mx=np.mean(x[bix==i])
        stdx=np.std(x[bix==i])
        my=np.mean(y[bix==i])
        stdy=np.std(y[bix==i])
        ax4.scatter(mx,my,c='k',s=len(x[bix==i]))
        ax4.errorbar(mx,my,xerr=stdx,yerr=stdy,ecolor='k',elinewidth=0.5,zorder=0)
    # s,yi,sd,lirmse,lichi2,co,ex,losd,lormse,lochi2=odr_fit_ref(mx,my)
    # ax4.scatter(x,y,c='k',alpha=0.5,s=2,zorder=0) 
    ax4.plot(z,co*z**ex,c='k',label='Local Relief - RMSE='+str(np.round(lormse,3)))  

    x=dfs.loc[spidx,'max_z'].to_numpy()
    y=snmp.to_numpy()
    s,yi,sd,lirmse,lichi2,co,ex,losd,lormse,lochi2=odr_fit_ref(x,y)
    bin_edges=np.histogram_bin_edges(snmp,'doane')
    bix=np.digitize(snmp,bin_edges)
    for i in range(len(bin_edges)-1):
        mx=np.mean(x[bix==i])
        stdx=np.std(x[bix==i])
        my=np.mean(y[bix==i])
        stdy=np.std(y[bix==i])
        ax4.scatter(mx,my,c='r',s=len(x[bix==i]))
        ax4.errorbar(mx,my,xerr=stdx,yerr=stdy,ecolor='r',elinewidth=0.5,zorder=0)
    # s,yi,sd,lirmse,lichi2,co,ex,losd,lormse,lochi2=odr_fit_ref(mx,my)
    # ax4.scatter(x,y,c='r',alpha=0.5,s=2,zorder=0) 
    ax4.plot(z,co*z**ex,c='r',label='Max Z - RMSE='+str(np.round(lormse,3)))
    locat.append(ttl_text)
    descrip.append('maxZ to snowP')
    param1.append(co)
    param2.append(ex)

    x=dfs.loc[spidx,'mean_z'].to_numpy()
    y=snmp.to_numpy()
    s,yi,sd,lirmse,lichi2,co,ex,losd,lormse,lochi2=odr_fit_ref(x,y)
    bin_edges=np.histogram_bin_edges(snmp,'doane')
    bix=np.digitize(snmp,bin_edges)
    for i in range(len(bin_edges)-1):
        mx=np.mean(x[bix==i])
        stdx=np.std(x[bix==i])
        my=np.mean(y[bix==i])
        stdy=np.std(y[bix==i])
        ax4.scatter(mx,my,c='b',s=len(x[bix==i]))
        ax4.errorbar(mx,my,xerr=stdx,yerr=stdy,ecolor='b',elinewidth=0.5,zorder=0)
    # s,yi,sd,lirmse,lichi2,co,ex,losd,lormse,lochi2=odr_fit_ref(mx,my)
    # ax4.scatter(x,y,c='b',alpha=0.5,s=2,zorder=0) 
    ax4.plot(z,co*z**ex,c='b',label='Mean Z - RMSE='+str(np.round(lormse,3)))
    
    ax4.set_xlabel('Elevation [m]')
    ax4.set_ylabel('% Runoff from Snowmelt')
    ax4.legend(loc='best')
    plt.ylim((0,1)) 
    plt.tight_layout
    
    return locat,descrip,param1,param2


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


# Set cutoffs
percb_cutoff=0.25

rem_idx=(grlf<=500) | (df_global['mean_z']<=250) | (global_perc_base>percb_cutoff) | (df_global['mean_rlf25']<250)
df_global_s=df_global.drop(df_global.index[rem_idx])
df_global_s=df_global_s.reset_index(drop=True)
global_perc_snow=df_global_s['qsm']/df_global_s['mean_runoff']


f1=plt.figure(1,figsize=(6,8))
ax1=plt.subplot(3,2,1)
ax2=plt.subplot(3,2,2)
ax3=plt.subplot(3,2,3)
ax4=plt.subplot(3,2,4)
ax5=plt.subplot(3,2,5)
ax6=plt.subplot(3,2,6)

ax1.set_ylim((0,11))
ax3.set_ylim((0,11))
ax5.set_ylim((0,11))

ax2.set_ylim((0,1))
ax4.set_ylim((0,1))
ax6.set_ylim((0,1))

## Generate regional relationships
# Greater Caucasus
bl1=38
br1=51
bb1=39.5
bt1=45
[l1,d1,pp11,pp12]=plot_and_fit_ref(ax1,ax2,df_global_s,bl1,br1,bb1,bt1,'Greater Caucasus')

# Alps
bl2=5
br2=16
bb2=43
bt2=50
[l2,d2,pp21,pp22]=plot_and_fit_ref(ax3,ax4,df_global_s,bl2,br2,bb2,bt2,'Alps')

# British Columbia
bl3=-131
br3=-120
bb3=48
bt3=54
[l3,d3,pp31,pp32]=plot_and_fit_ref(ax5,ax6,df_global_s,bl3,br3,bb3,bt3,'British Columbia')

plt.tight_layout()

plt.rcdefaults()

f1.savefig('P1_SUP_FigureS2.pdf')