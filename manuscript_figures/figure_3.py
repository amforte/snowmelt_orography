#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  6 09:12:49 2023

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

master_location='/Volumes/Choruh/Data/snowmelt_project/'
repo_location='/Users/aforte/Documents/GitHub/snowmelt_orography/geospatial_codes/'

# Load original
gages2raster=pd.read_csv(repo_location+'gages2_wrr2_raster_values.csv')
gages2real=pd.read_csv(repo_location+'gages2_real_ts.csv')
gages2stats=pd.read_csv(repo_location+'gages2_station_stats.csv')

gages2ids1=pd.read_csv(master_location+'gagesII/basinchar_and_report_sept_2011/spreadsheets-in-csv-format/conterm_basinid.txt')
gages2ids2=pd.read_csv(master_location+'gagesII/basinchar_and_report_sept_2011/spreadsheets-in-csv-format/AKHIPR_basinid.txt')
gages2ids=pd.concat([gages2ids1,gages2ids2],axis=0)
gages2hcdn=gages2ids[['STAID','HCDN-2009']]

# Extract lat-lon of center
gages2morph1=pd.read_csv(master_location+'gagesII/basinchar_and_report_sept_2011/spreadsheets-in-csv-format/conterm_bas_morph.txt')
gages2morph2=pd.read_csv(master_location+'gagesII/basinchar_and_report_sept_2011/spreadsheets-in-csv-format/AKHIPR_bas_morph.txt')
gages2morph=pd.concat([gages2morph1,gages2morph2],axis=0)

# Merge
df=pd.merge(gages2raster,gages2hcdn,on='STAID',how='inner')
df=pd.merge(gages2real,df,on='STAID',how='inner')
df=pd.merge(gages2morph,df,on='STAID',how='inner')
df=pd.merge(gages2stats,df,on='STAID',how='inner')
percb_cutoff=0.25
perc_base=df['QSB']/df['R']
rlf=df['MAX_Z']-df['MIN_Z']

# Establish index and produce new dataset
hcdn_idx=(rlf>500) & (df['HCDN-2009']=='yes') & (df['MEAN_Z']>250) & (df['SlicedComp']>0.95) & (perc_base<percb_cutoff)
df_hcdn=df.loc[hcdn_idx,:]
perc_snow=df_hcdn['QSM']/df_hcdn['R']

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

ax11=plt.subplot(1,2,1)
idx1=(global_perc_snow<=0.35) & (df_global_s['r_c1']>0)
g2idx1=perc_snow<=0.35
# ax1.scatter(df_global_s.loc[idx1,'mean_runoff'],df_global_s.loc[idx1,'r_c1'],c='k',s=1)
ax11.scatter(df_hcdn.loc[g2idx1,'SlicedMeanR'],df_hcdn.loc[g2idx1,'SlicedRC1'],zorder=2,s=5,c='k',label='HCDN-2009')
sc11=plt.hist2d(df_global_s.loc[idx1,'mean_runoff'],df_global_s.loc[idx1,'r_c1'],[xb,yb],norm=colors.LogNorm(vmin=1,vmax=1000),cmap=cm.lajolla)
ax11.set_xlabel('WaterGAP3 Runoff [mm/day]')
ax11.set_ylabel('WaterGAP3 Shape Parameter')
ax11.set_xlim((0,10))
ax11.set_ylim((0,2.5))
ax11.set_title('Snowmelt Fraction < 0.35')
ax11.legend(loc=4)
cbar11=plt.colorbar(sc11[3],ax=ax11)
ax11.text(0.01, 0.99, 'A',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax11.transAxes,
        fontsize=12,fontweight='extra bold')

ax12=plt.subplot(1,2,2)
idx2=(global_perc_snow>0.35) & (df_global_s['r_c1']>0)
g2idx2=perc_snow>0.35
# ax2.scatter(df_global_s.loc[idx2,'mean_runoff'],df_global_s.loc[idx2,'r_c1'],c='k',s=1)
ax12.scatter(df_hcdn.loc[g2idx2,'SlicedMeanR'],df_hcdn.loc[g2idx2,'SlicedRC1'],zorder=2,s=5,c='k',label='HCDN-2009')
sc12=plt.hist2d(df_global_s.loc[idx2,'mean_runoff'],df_global_s.loc[idx2,'r_c1'],[xb,yb],norm=colors.LogNorm(vmin=1,vmax=1000),cmap=cm.lajolla)
ax12.set_xlabel('WaterGAP3 Runoff [mm/day]')
ax12.set_ylabel('WaterGAP3 Shape Parameter')
ax12.set_xlim((0,10))
ax12.set_ylim((0,2.5))
ax12.set_title('Snowmelt Fraction > 0.35')
ax12.legend(loc=4)
cbar12=plt.colorbar(sc12[3],ax=ax12)
cbar12.ax.set_ylabel('Density')
ax12.text(0.01, 0.99, 'B',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax12.transAxes,
        fontsize=12,fontweight='extra bold')

plt.tight_layout()

f2=plt.figure(2,figsize=(11,8))
f2.set_dpi(250)
plt.rc('axes', titlesize=SMALL_SIZE) 
plt.rc('axes', labelsize=SMALL_SIZE) 

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
    
    rmse_array=np.array([rmse1,rmse2])
    if np.argmin(rmse_array)==0:
        plt.plot(r,logcoeff*r**logexp,c='k',linewidth=1,label=lab1)
        plt.plot(r,linslp*r+linint,c='k',linestyle=':',label=lab2)
    elif np.argmin(rmse_array)==1:
        plt.plot(r,linslp*r+linint,c='k',linewidth=1,label=lab2)
        plt.plot(r,logcoeff*r**logexp,c='k',linestyle=':',label=lab1)
        
    if bv[i+1]<=0.35:
        ax11.plot(r,logcoeff*r**logexp,c='k',linewidth=bv[i+1]*3,linestyle='--')
    else:
        ax12.plot(r,linslp*r+linint,c='k',linewidth=bv[i+1],linestyle='--')

    
    cbar1=plt.colorbar(sc1[3],ax=ax1)

    if (i==0) | (i==4) | (i==8) | (i==12):
        plt.ylabel('WG3 Shape')
    if i>=12:
        plt.xlabel('WG3 Runoff [mm/day]')
    if (i==3) | (i==7) | (i==11) | (i==15):
        cbar1.ax.set_ylabel('Density')
        
    plt.xlim((0,10))
    plt.ylim((0,2.5))
    plt.legend(loc='best')
    
    plt.title(str(np.round(bv[i]*100,0))+' < % of Snowmelt Runoff <'+str(np.round(bv[i+1]*100,0))) 
    
plt.tight_layout()
plt.rcdefaults()

f1.savefig('figure_3.pdf',dpi="figure")
f2.savefig('figure_s1.pdf',dpi="figure")