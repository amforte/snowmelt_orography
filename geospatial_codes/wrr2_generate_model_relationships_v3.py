#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script written by Adam M. Forte
aforte8@lsu.edu

This script reads the output of 'wrr2_derived_make_table_v4.py' and generates the
relationships between relief and mean runoff and maximum elevation and snowmelt 
used within the stimpy model to link growing topography to runoff.
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

def plot_and_fit_ref(fn,dfs,bl,br,bb,bt,ttl_text):
    
    locat=[]
    descrip=[]
    param1=[]
    param2=[]
    
    spidx=(dfs['latitude']>=bb) & (dfs['latitude']<=bt) & (dfs['longitude']>=bl) & (dfs['longitude']<=br)
    
    f=plt.figure(fn,figsize=(25,10)) 
    gs=f.add_gridspec(1,3) 
    
    max_z=np.max(np.ceil(dfs.loc[spidx,'max_z'].to_numpy()))
    max_r=np.max(np.ceil(dfs.loc[spidx,'mean_runoff'].to_numpy()))
    max_rlf=np.max(np.ceil(dfs.loc[spidx,'mean_rlf25'].to_numpy()))
    
    ax2=f.add_subplot(gs[1])
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
    plt.legend(loc='best')
    plt.ylim((0,max_r))
    plt.xlim((0,max_z))
    plt.title(ttl_text)
    
    
    ax4=f.add_subplot(gs[2])
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
    plt.legend(loc='best')
    plt.ylim((0,1)) 
    
    return locat,descrip,param1,param2

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

plt.figure(1,figsize=(30,30))
bv=np.concatenate((np.arange(0,0.80,0.05),[1.0]),axis=0)
xb=np.linspace(0,10,50)
yb=np.linspace(0,2,50)
rmse=np.zeros((len(bv)-1,4))


for i in range(len(bv)-1):
    ax1=plt.subplot(4,4,i+1)
    gidx=(global_perc_snow>=bv[i]) & (global_perc_snow<bv[i+1]) & (df_global_s['r_c1']>0)
    sc1=plt.hist2d(df_global_s.loc[gidx,'mean_runoff'],df_global_s.loc[gidx,'r_c1'],[xb,yb],norm=colors.LogNorm(vmin=1,vmax=500),cmap=cm.lajolla)
    
    
    fdlog=odr.Data(np.log10(df_global_s.loc[gidx,'mean_runoff']),np.log10(df_global_s.loc[gidx,'r_c1']))
    odrlog=odr.ODR(fdlog,linmod,beta0=[0.1,10])
    outlog=odrlog.run()
    logexp=outlog.beta[0]
    logcoeff=10**outlog.beta[1]
    logchisquare=outlog.res_var
    cr1_pred=logcoeff*df_global_s.loc[gidx,'mean_runoff']**logexp
    rmse1=np.sqrt(np.sum((cr1_pred-df_global_s.loc[gidx,'r_c1'])**2)/len(cr1_pred))
    lab1=r'$c_{R}$ = '+str(np.round(logcoeff,2))+' * $R^{'+str(np.round(logexp,2))+'}$; RMSE = '+str(np.round(rmse1,4))
    
    fdlin=odr.Data(df_global_s.loc[gidx,'mean_runoff'],df_global_s.loc[gidx,'r_c1'])
    odrlin=odr.ODR(fdlin,linmod,beta0=[1,0.1])
    outlin=odrlin.run()
    linslp=outlin.beta[0]
    linint=outlin.beta[1]
    linchisquare=outlin.res_var
    cr2_pred=linslp*df_global_s.loc[gidx,'mean_runoff']+linint
    rmse2=np.sqrt(np.sum((cr2_pred-df_global_s.loc[gidx,'r_c1'])**2)/len(cr2_pred))
    lab2=r'$c_{R}$ = '+str(np.round(linslp,2))+' * R + '+str(np.round(linint,2))+'; RMSE = '+str(np.round(rmse2,4))
    
    rmse_array=np.array([rmse1,rmse2])
    if np.argmin(rmse_array)==0:
        plt.plot(r,logcoeff*r**logexp,c='k',linewidth=2,label=lab1)
        plt.plot(r,linslp*r+linint,c='k',linestyle=':',label=lab2)
    elif np.argmin(rmse_array)==1:
        plt.plot(r,linslp*r+linint,c='k',linewidth=2,label=lab2)
        plt.plot(r,logcoeff*r**logexp,c='k',linestyle=':',label=lab1)

    cbar1=plt.colorbar(sc1[3],ax=ax1)
    cbar1.ax.set_ylabel('Density')
    plt.xlabel('Global Raster Runoff [mm/day]')
    plt.ylabel('Global Raster Variability')
    plt.xlim((0,10))
    plt.ylim((0,2.5))
    plt.legend(loc='best')
    
    plt.title(str(np.round(bv[i]*100,0))+' < % of Snowmelt Runoff <'+str(np.round(bv[i+1]*100,0))) 
    
## Calculate relationship between mean runoff and shape
# Initialize containers
# For mean runoff to shape (cr)
r_type=[]
p1=np.zeros((len(bv)-1))
p2=np.zeros((len(bv)-1)) 
sbl=np.zeros((len(bv)-1))
sbr=np.zeros((len(bv)-1))  

plt.figure(2,figsize=(30,30))
for i in range(len(bv)-1):
    ax1=plt.subplot(4,4,i+1)
    gidx=(global_perc_snow>=bv[i]) & (global_perc_snow<bv[i+1]) & (df_global_s['r_c1']>0)
    sc1=plt.scatter(df_global_s.loc[gidx,'mean_runoff'],df_global_s.loc[gidx,'r_c1'],c=df_global_s.loc[gidx,'matup'],cmap=cm.vik,norm=cnorm,s=3)

    sbl[i]=bv[i]
    sbr[i]=bv[i+1]    

    fdlog=odr.Data(np.log10(df_global_s.loc[gidx,'mean_runoff']),np.log10(df_global_s.loc[gidx,'r_c1']))
    odrlog=odr.ODR(fdlog,linmod,beta0=[0.1,10])
    outlog=odrlog.run()
    logexp=outlog.beta[0]
    logcoeff=10**outlog.beta[1]
    logchisquare=outlog.res_var
    cr1_pred=logcoeff*df_global_s.loc[gidx,'mean_runoff']**logexp
    rmse1=np.sqrt(np.sum((cr1_pred-df_global_s.loc[gidx,'r_c1'])**2)/len(cr1_pred))
    lab1=r'$c_{R}$ = '+str(np.round(logcoeff,2))+' * $R^{'+str(np.round(logexp,2))+'}$; RMSE = '+str(np.round(rmse1,4))
    
    fdlin=odr.Data(df_global_s.loc[gidx,'mean_runoff'],df_global_s.loc[gidx,'r_c1'])
    odrlin=odr.ODR(fdlin,linmod,beta0=[1,0.1])
    outlin=odrlin.run()
    linslp=outlin.beta[0]
    linint=outlin.beta[1]
    linchisquare=outlin.res_var
    cr2_pred=linslp*df_global_s.loc[gidx,'mean_runoff']+linint
    rmse2=np.sqrt(np.sum((cr2_pred-df_global_s.loc[gidx,'r_c1'])**2)/len(cr2_pred))
    lab2=r'$c_{R}$ = '+str(np.round(linslp,2))+' * R + '+str(np.round(linint,2))+'; RMSE = '+str(np.round(rmse2,4))
    
    rmse_array=np.array([rmse1,rmse2])
    if np.argmin(rmse_array)==0:
        plt.plot(r,logcoeff*r**logexp,c='k',linewidth=2,label=lab1)
        plt.plot(r,linslp*r+linint,c='k',linestyle=':',label=lab2)
        r_type.append('power')
        p1[i]=logcoeff
        p2[i]=logexp
    elif np.argmin(rmse_array)==1:
        plt.plot(r,linslp*r+linint,c='k',linewidth=2,label=lab2)
        plt.plot(r,logcoeff*r**logexp,c='k',linestyle=':',label=lab1)
        r_type.append('linear')
        p1[i]=linslp
        p2[i]=linint
    
    cbar1=plt.colorbar(sc1,ax=ax1)
    cbar1.ax.set_ylabel('Mat [C]')
    plt.xlabel('Global Raster Runoff [mm/day]')
    plt.ylabel('Global Raster Variability')
    plt.xscale('log')
    # plt.yscale('log')
    plt.xlim((0.1,10))
    plt.ylim((0.1,2.5))
    plt.legend(loc='best')
    
    plt.title(str(np.round(bv[i]*100,0))+' < % of Snowmelt Runoff <'+str(np.round(bv[i+1]*100,0)))
    
df_out1=pd.DataFrame(data={'snowP_left':sbl,
                           'snowP_right':sbr,
                           'relationship':r_type,
                           'param1':p1,
                           'param2':p2})
df_out1.to_csv(repo_location+'snowmelt_meanR_cR_relationships.csv',index=False)

## Calculate relationship between mean approximated scale and fit scale
# Initialize containers
# For approx scale to fit scale
r_type=[]
p1=np.zeros((len(bv)-1))
p2=np.zeros((len(bv)-1)) 
sbl=np.zeros((len(bv)-1))
sbr=np.zeros((len(bv)-1))  

plt.figure(3,figsize=(30,30))
for i in range(len(bv)-1):
    ax1=plt.subplot(4,4,i+1)
    gidx=(global_perc_snow>=bv[i]) & (global_perc_snow<bv[i+1]) & (df_global_s['r_c1']>0)
    sc1=plt.scatter(s_f[gidx],df_global_s.loc[gidx,'r_s1'],c=df_global_s.loc[gidx,'matup'],cmap=cm.vik,norm=cnorm,s=3)

    sbl[i]=bv[i]
    sbr[i]=bv[i+1]    

    
    fdlin=odr.Data(s_f[gidx],df_global_s.loc[gidx,'r_s1'])
    odrlin=odr.ODR(fdlin,linmod,beta0=[1,0.1])
    outlin=odrlin.run()
    linslp=outlin.beta[0]
    linint=outlin.beta[1]
    linchisquare=outlin.res_var
    
    plt.plot(r,linslp*r+linint,c='k',linewidth=2,label=lab2)

    r_type.append('linear')
    p1[i]=linslp
    p2[i]=linint
    
    cbar1=plt.colorbar(sc1,ax=ax1)
    cbar1.ax.set_ylabel('Mat [C]')
    plt.xlabel('Approximated Scale from Mean')
    plt.ylabel('Fit Scale')
    plt.xscale('log')
    # plt.yscale('log')
    plt.xlim((0.1,15))
    plt.ylim((0.1,15))
    plt.legend(loc='best')
    
    plt.title(str(np.round(bv[i]*100,0))+' < % of Snowmelt Runoff <'+str(np.round(bv[i+1]*100,0)))
    
df_out2=pd.DataFrame(data={'snowP_left':sbl,
                           'snowP_right':sbr,
                           'relationship':r_type,
                           'param1':p1,
                           'param2':p2})
df_out2.to_csv(repo_location+'snowmelt_approx_sR_fit_sR_relationships.csv',index=False)

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

f4=plt.figure(4,figsize=(10,5))

plt.subplot(1,2,1)
sc1=plt.hist2d(df_global_s.loc[gidx,'mean_runoff'],df_global_s.loc[gidx,'r_c1'],[xb,yb],norm=colors.LogNorm(vmin=1,vmax=500),cmap=cm.lajolla)
cbar1=plt.colorbar(sc1[3],ax=ax1)
cbar1.ax.set_ylabel('Density')
plt.plot(r,RFcr_logcoeff*r**RFcr_logexp,c='k',linewidth=2)
plt.xlabel('Global Raster Runoff [mm/day]')
plt.ylabel('Global Raster Variability')
plt.xlim((0,10))
plt.ylim((0,2.5))

plt.subplot(1,2,2)
sc1=plt.hist2d(s_f[gidx],df_global_s.loc[gidx,'r_s1'],[xb,xb],norm=colors.LogNorm(vmin=1,vmax=500),cmap=cm.berlin)
cbar1=plt.colorbar(sc1[3],ax=ax1)
cbar1.ax.set_ylabel('Density')
plt.plot(r,RFsr_m*r + RFsr_b,c='k',linewidth=2)
plt.xlabel('Scale Approximated From Mean')
plt.ylabel('Global Raster Scale')
plt.xlim((0,10))
plt.ylim((0,10))

df_out3=pd.DataFrame(data={'function':rel,
                           'relationshiip':r_type,
                           'param1':param1,
                           'param2':param2})
df_out3.to_csv(repo_location+'rainfall_relationships.csv',index=False)

plt.figure(5,figsize=(8,8))

## Generate regional relationships
# Greater Caucasus
bl1=38
br1=51
bb1=39.5
bt1=45
[l1,d1,pp11,pp12]=plot_and_fit_ref(5,df_global_s,bl1,br1,bb1,bt1,'Greater Caucasus')

# Alps
bl2=5
br2=16
bb2=43
bt2=50
[l2,d2,pp21,pp22]=plot_and_fit_ref(6,df_global_s,bl2,br2,bb2,bt2,'Alps')

# British Columbia
bl3=-131
br3=-120
bb3=48
bt3=54
[l3,d3,pp31,pp32]=plot_and_fit_ref(7,df_global_s,bl3,br3,bb3,bt3,'British Columbia')

df_out4=pd.DataFrame(data={'location':l1+l2+l3,
                           'relation_between':d1+d2+d3,
                           'param1':pp11+pp21+pp31,
                           'param2':pp12+pp22+pp32})
df_out4.to_csv(repo_location+'topo_snow_mean_relationships.csv',index=False)
