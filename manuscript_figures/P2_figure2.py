#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 24 12:49:54 2023

@author: aforte
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from cmcrameri import cm
from scipy.stats import linregress
from matplotlib import colors
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

def odr_fit_ref2(x,y):
    # Filter 0 values
    lx=np.log10(x[(x>0) & (y>0)])
    ly=np.log10(y[(x>0) & (y>0)])
    linmod=odr.Model(linear)
    
    fdlog=odr.Data(lx,ly)
    odrlog=odr.ODR(fdlog,linmod,beta0=[0.1,10])
    outlog=odrlog.run()
    logexp=outlog.beta[0]
    logcoeff=10**outlog.beta[1]

    return logcoeff,logexp


def plot_and_fit_ref(axn,axnt,dfs,bl,br,bb,bt,col,lab):
    spidx=(dfs['latitude']>=bb) & (dfs['latitude']<=bt) & (dfs['longitude']>=bl) & (dfs['longitude']<=br)
    
    max_z=4500
    max_r=10
    
    z=np.linspace(1,max_z,100)
    
    x=dfs.loc[spidx,'mean_rlf'].to_numpy()
    y=dfs.loc[spidx,'mean_runoff'].to_numpy()
    co,ex=odr_fit_ref2(x,y)
    bin_edges=np.histogram_bin_edges(dfs.loc[spidx,'mean_runoff'],'doane')
    bix=np.digitize(dfs.loc[spidx,'mean_runoff'],bin_edges)
    for i in range(len(bin_edges)-1):
        mx=np.mean(x[bix==i])
        stdx=np.std(x[bix==i])
        my=np.mean(y[bix==i])
        stdy=np.std(y[bix==i])
        axn.scatter(mx,my,c=col,s=len(x[bix==i]))
        axn.errorbar(mx,my,xerr=stdx,yerr=stdy,ecolor=col,elinewidth=0.5,zorder=0)
    axn.plot(z,co*z**ex,c=col,label=lab)  
    axn.set_xlabel('Mean Local Relief [m]')
    axn.set_ylabel('Mean Runoff [mm/day]')
    axn.set_xlim((0,2500))
    axn.set_ylim((0,max_r))

    snmp=dfs.loc[spidx,'qsm']/dfs.loc[spidx,'mean_runoff']
    x=dfs.loc[spidx,'max_z'].to_numpy()
    y=snmp.to_numpy()
    co,ex=odr_fit_ref2(x,y)
    bin_edges=np.histogram_bin_edges(snmp,'doane')
    bix=np.digitize(snmp,bin_edges)
    for i in range(len(bin_edges)-1):
        mx=np.mean(x[bix==i])
        stdx=np.std(x[bix==i])
        my=np.mean(y[bix==i])
        stdy=np.std(y[bix==i])
        axnt.scatter(mx,my,c=col,s=len(x[bix==i]),marker='s')
        axnt.errorbar(mx,my,xerr=stdx,yerr=stdy,ecolor=col,elinewidth=0.5,zorder=0)
    axnt.plot(z,co*z**ex,c=col,linestyle='--',label=lab)
    axnt.set_ylim((0,1))
    axnt.set_ylabel('Snowmelt Fraction')
    axnt.set_xlabel('Maximum Elevation [m]')
    
    axn.legend(loc='lower left',bbox_to_anchor= (0.05,-0.5))
    axnt.legend(loc='lower left',bbox_to_anchor= (0.05,-0.5))
        

# master_location='/Volumes/Choruh/Data/snowmelt_project/'
master_location='/Users/aforte/Documents/Python/snowmelt/'
repo_location='/Users/aforte/Documents/GitHub/snowmelt_orography/geospatial_codes/'

## Load global
df_global=pd.read_csv(master_location+'wrr2_derived_data_v3.csv')
df_global=df_global.drop(index=df_global.index[np.isnan(df_global['mean_z'])])
df_global=df_global.reset_index(drop=True)

# Calc percents
global_perc_base=df_global['qsb']/df_global['mean_runoff']

# Calculate indices
grlf=df_global['max_z']-df_global['min_z']

# Set cutoffs
perc_cutoff=0.275
percb_cutoff=0.25

rem_idx=(grlf<=500) | (df_global['mean_z']<=250) | (global_perc_base>percb_cutoff) | (df_global['mean_rlf']<250)
df_global_s=df_global.drop(df_global.index[rem_idx])
df_global_s=df_global_s.reset_index(drop=True)
global_perc_snow=df_global_s['qsm']/df_global_s['mean_runoff']


# Set vector
r=np.linspace(0.1,25,100)

# Define model
linmod=odr.Model(linear)



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

bv=np.concatenate((np.arange(0,0.80,0.05),[1.0]),axis=0)
xb=np.linspace(0,10,50)
yb=np.linspace(0,2,50)


f1=plt.figure(1,figsize=(8,4))
f1.set_dpi(250)

ax11=plt.subplot(1,3,1)
idx1=(df_global_s['r_c1']>0)

sc11=plt.hist2d(df_global_s.loc[idx1,'mean_runoff'],df_global_s.loc[idx1,'r_c1'],[xb,yb],norm=colors.LogNorm(vmin=1,vmax=1000),cmap=cm.grayC)
ax11.set_xlabel('WaterGAP3 Runoff [mm/day]')
ax11.set_ylabel('WaterGAP3 Shape Parameter')
ax11.set_xlim((0,10))
ax11.set_ylim((0,2.5))
cbar11=plt.colorbar(sc11[3],ax=ax11)
cbar11.ax.set_ylabel('Density')
ax11.text(0.01, 0.99, 'A',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax11.transAxes,
        fontsize=12,fontweight='extra bold')


cnt1=0
cnt2=0

for i in range(len(bv)-1):
    gidx=(global_perc_snow>=bv[i]) & (global_perc_snow<bv[i+1]) & (df_global_s['r_c1']>0)
        
    fdlog=odr.Data(np.log10(df_global_s.loc[gidx,'mean_runoff']),np.log10(df_global_s.loc[gidx,'r_c1']))
    odrlog=odr.ODR(fdlog,linmod,beta0=[0.1,10])
    outlog=odrlog.run()
    logexp=outlog.beta[0]
    logcoeff=10**outlog.beta[1]
    logchisquare=outlog.res_var
    cr1_pred=logcoeff*df_global_s.loc[gidx,'mean_runoff']**logexp
    rmse1=np.sqrt(np.sum((cr1_pred-df_global_s.loc[gidx,'r_c1'])**2)/len(cr1_pred))
    lab1='RMSE = '+str(np.round(rmse1,4))
    
    fdlin=odr.Data(df_global_s.loc[gidx,'mean_runoff'],df_global_s.loc[gidx,'r_c1'])
    odrlin=odr.ODR(fdlin,linmod,beta0=[1,0.1])
    outlin=odrlin.run()
    linslp=outlin.beta[0]
    linint=outlin.beta[1]
    linchisquare=outlin.res_var
    cr2_pred=linslp*df_global_s.loc[gidx,'mean_runoff']+linint
    rmse2=np.sqrt(np.sum((cr2_pred-df_global_s.loc[gidx,'r_c1'])**2)/len(cr2_pred))
    lab2='RMSE = '+str(np.round(rmse2,4))
    
    if bv[i+1]<=0.35:
        if cnt1==0:
            ax11.plot(r,logcoeff*r**logexp,c='indianred',linewidth=bv[i+1]*3,linestyle='-',label='Snowmelt Fraction < 0.35')
        else:
            ax11.plot(r,logcoeff*r**logexp,c='indianred',linewidth=bv[i+1]*3,linestyle='-') 
        cnt1=cnt1+1
    else:
        if cnt2==0:
            ax11.plot(r,linslp*r+linint,c='dodgerblue',linewidth=bv[i+1],linestyle='--',label='Snowmelt Fraction > 0.35')
        else:
            ax11.plot(r,linslp*r+linint,c='dodgerblue',linewidth=bv[i+1],linestyle='--')
        cnt2=cnt2+1
        

ax11.legend(loc='lower left',bbox_to_anchor= (-0.15,-0.5))

# Greater Caucasus
bl1=38
br1=51
bb1=39.5
bt1=45
gc_col='black'

# Alps
bl2=5
br2=16
bb2=43
bt2=50
alps_col='royalblue'

# British Columbia
bl3=-131
br3=-120
bb3=48
bt3=54
bc_col='orange'

ax13=plt.subplot(1,3,2)
ax14=plt.subplot(1,3,3)
plot_and_fit_ref(ax13,ax14,df_global_s,bl3,br3,bb3,bt3,bc_col,'British Columbia')
plot_and_fit_ref(ax13,ax14,df_global_s,bl2,br2,bb2,bt2,alps_col,'Alps')
plot_and_fit_ref(ax13,ax14,df_global_s,bl1,br1,bb1,bt1,gc_col,'Greater Caucasus')

ax13.text(0.01, 0.99, 'B',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax13.transAxes,
        fontsize=12,fontweight='extra bold')

ax14.text(0.01, 0.99, 'C',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax14.transAxes,
        fontsize=12,fontweight='extra bold')

plt.tight_layout()
f1.savefig('P2_figure2.pdf',dpi='figure')
plt.rcdefaults()