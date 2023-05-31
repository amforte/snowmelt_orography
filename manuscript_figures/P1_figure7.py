#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 14 19:53:53 2022

@author: aforte
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from cmcrameri import cm
from matplotlib import colors
from scipy.stats import pearsonr
from scipy.stats import spearmanr
from scipy import stats
from scipy import odr
from scipy.stats import weibull_min
from scipy.special import gamma

import matplotlib.gridspec as gridspec
import cartopy.crs as ccrs
import cartopy.feature as cfeature

from netCDF4 import Dataset
import rasterio
import rasterio.transform as rt


def bin2d_corr(fn,dfs,x,y,x1,y1,min_num,p_sig,corr_type,bin_type='doane'):
    
    # Test overwriting bins
    x_idx=x<=np.percentile(x,99.9)
    y_idx=y<=np.percentile(y,99.9)
    xbin=np.histogram_bin_edges(x[x_idx],bin_type)
    ybin=np.histogram_bin_edges(y[y_idx],bin_type)    
    
    x_dim=len(xbin)-1
    y_dim=len(ybin)-1
    
    # Allocate
    X=np.zeros((y_dim,x_dim))
    P=np.ones((y_dim,x_dim))*0.25
    RHO=np.zeros(x.shape)*np.nan
    NUM=np.zeros(x.shape)*np.nan
    BN=np.zeros(x.shape)*np.nan
    
    x_ix=np.digitize(x,xbin)
    y_ix=np.digitize(y,ybin)

    for i in range(x_dim):
        for j in range(y_dim):
            idx_in=np.logical_and(x_ix==i+1,y_ix==j+1)
            NUM[idx_in]=np.sum(idx_in)
            if np.sum(idx_in)>min_num:
                if corr_type=='pearson':
                    [rho,p]=pearsonr(x1[idx_in],y1[idx_in])
                elif corr_type=='spearman':
                    [rho,p]=spearmanr(x1[idx_in],y1[idx_in])
                    
                if p<p_sig:
                    X[j,i]=rho
                    P[j,i]=np.nan
                    RHO[idx_in]=rho
                    BN[idx_in]=j+i
                else:
                    X[j,i]=np.nan
                    P[j,i]=1

            else:
                X[j,i]=np.nan
    return RHO,NUM,BN,X,P,xbin,ybin

## Load global
# master_location='/Volumes/Choruh/Data/snowmelt_project/'
master_location='/Users/aforte/Documents/Python/snowmelt/'

df_global=pd.read_csv(master_location+'wrr2_derived_data_v3.csv')
df_global=df_global.drop(index=df_global.index[np.isnan(df_global['mean_z'])])
df_global=df_global.reset_index(drop=True)

percb_cutoff=0.25
grlf=df_global['max_z']-df_global['min_z']
global_perc_base=df_global['qsb']/df_global['mean_runoff']

rem_idx=(grlf<=500) | (df_global['mean_z']<=250) | (global_perc_base>percb_cutoff)
df_global_s=df_global.drop(df_global.index[rem_idx])
df_global_s=df_global_s.reset_index(drop=True)

[RHO1a,NUM1a,BN1a,X1a,P1a,xbin1a,ybin1a]=bin2d_corr(1,df_global_s,df_global_s['mean_precip'],
                                            df_global_s['mat'],df_global_s['mean_rlf'],df_global_s['mean_runoff'],10,0.01,'spearman')

[RHO1b,NUM1b,BN1b,X1b,P1b,xbin1b,ybin1b]=bin2d_corr(1,df_global_s,df_global_s['mean_precip'],
                                            df_global_s['mat'],df_global_s['mean_z'],df_global_s['mean_runoff'],10,0.01,'spearman')

[RHO1c,NUM1c,BN1c,X1c,P1c,xbin1c,ybin1c]=bin2d_corr(1,df_global_s,df_global_s['mean_precip'],
                                            df_global_s['mat'],df_global_s['max_z'],df_global_s['mean_runoff'],10,0.01,'spearman')

[RHO2a,NUM2a,BN2a,X2a,P2a,xbin2a,ybin2a]=bin2d_corr(2,df_global_s,df_global_s['mean_precip'],
                                            df_global_s['mat'],df_global_s['max_z'],df_global_s['qsm']/df_global_s['mean_runoff'],10,0.01,'spearman')

[RHO2b,NUM2b,BN2b,X2b,P2b,xbin2b,ybin2b]=bin2d_corr(2,df_global_s,df_global_s['mean_precip'],
                                            df_global_s['mat'],df_global_s['mean_z'],df_global_s['qsm']/df_global_s['mean_runoff'],10,0.01,'spearman')

[RHO2c,NUM2c,BN2c,X2c,P2c,xbin2c,ybin2c]=bin2d_corr(2,df_global_s,df_global_s['mean_precip'],
                                            df_global_s['mat'],df_global_s['mean_rlf'],df_global_s['qsm']/df_global_s['mean_runoff'],10,0.01,'spearman')

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

f1=plt.figure(1,figsize=(8,2))
f1.set_dpi(250)
gs=f1.add_gridspec(1,3)   
 
ax1=f1.add_subplot(gs[0,0])
xb=np.histogram_bin_edges(df_global_s['mean_precip'],'auto')
yb=np.histogram_bin_edges(df_global_s['mat'],'auto')

sc1=plt.hist2d(df_global_s['mean_precip'],df_global_s['mat'],[xb,yb],norm=colors.LogNorm(),cmap=cm.lajolla)
for i in range(len(xbin1a)):
    ax1.axvline(xbin1a[i],c='k',linestyle=':',linewidth=0.25,alpha=0.50)
for i in range(len(ybin1a)):
    ax1.axhline(ybin1a[i],c='k',linestyle=':',linewidth=0.25,alpha=0.50)
cbar1=plt.colorbar(sc1[3],ax=ax1)
cbar1.ax.set_ylabel('Density')
plt.xlim((np.min(xbin1a),np.max(xbin1a)))
plt.ylim((np.min(ybin1a),np.max(ybin1a)))
plt.xlabel('Mean Precip [mm/day]')
plt.ylabel('MAT [C]')
ax1.text(0.90, 0.99, 'A',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax1.transAxes,
        fontsize=12,fontweight='extra bold')

ax2=f1.add_subplot(gs[0,1])
sc2=plt.imshow(X1a,origin='lower',
                extent=[np.min(xbin1a),np.max(xbin1a),np.min(ybin1a),np.max(ybin1a)],
                cmap=cm.bam,norm=colors.Normalize(vmin=-1,vmax=1),aspect='auto')
cbar2=plt.colorbar(sc2,ax=ax2)
cbar2.ax.set_ylabel(r'Local Relief to $\bar{R}$ Correlation')
plt.imshow(P1a,origin='lower',
                extent=[np.min(xbin1a),np.max(xbin1a),np.min(ybin1a),np.max(ybin1a)],
                cmap=cm.grayC,norm=colors.Normalize(vmin=0,vmax=1),aspect='auto')
plt.xlabel('Mean Precip [mm/day]')
plt.ylabel('MAT [C]')
ax2.text(0.90, 0.99, 'B',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax2.transAxes,
        fontsize=12,fontweight='extra bold')

ax3=f1.add_subplot(gs[0,2])
sc3=plt.imshow(X2a,origin='lower',
                extent=[np.min(xbin2a),np.max(xbin2a),np.min(ybin2a),np.max(ybin2a)],
                cmap=cm.bam,norm=colors.Normalize(vmin=-1,vmax=1),aspect='auto')
cbar3=plt.colorbar(sc3,ax=ax3)
cbar3.ax.set_ylabel('Max Z to SF Correlation')
plt.imshow(P2a,origin='lower',
                extent=[np.min(xbin2a),np.max(xbin2a),np.min(ybin2a),np.max(ybin2a)],
                cmap=cm.grayC,norm=colors.Normalize(vmin=0,vmax=1),aspect='auto')
plt.xlabel('Mean Precip [mm/day]')
plt.ylabel('MAT [C]')
ax3.text(0.90, 0.99, 'C',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax3.transAxes,
        fontsize=12,fontweight='extra bold')

plt.tight_layout()

#############################
f2=plt.figure(2,figsize=(8,4))
f2.set_dpi(250)
gs=f2.add_gridspec(2,3)

ax1=f2.add_subplot(gs[0,0])
im1=plt.imshow(X1b,origin='lower',
                extent=[np.min(xbin1b),np.max(xbin1b),np.min(ybin1b),np.max(ybin1b)],
                cmap=cm.bam,norm=colors.Normalize(vmin=-1,vmax=1),aspect='auto')
cbar1=plt.colorbar(im1,ax=ax1)
cbar1.ax.set_ylabel(r'Mean Z to $\bar{R}$')
plt.imshow(P1b,origin='lower',
                extent=[np.min(xbin1b),np.max(xbin1b),np.min(ybin1b),np.max(ybin1b)],
                cmap=cm.grayC,norm=colors.Normalize(vmin=0,vmax=1),aspect='auto')
# plt.xlabel('Mean Precip [mm/day]')
plt.ylabel('MAT [C]')
ax1.text(0.90, 0.99, 'A',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax1.transAxes,
        fontsize=12,fontweight='extra bold')

ax2=f2.add_subplot(gs[0,1])
im2=plt.imshow(X1c,origin='lower',
                extent=[np.min(xbin1c),np.max(xbin1c),np.min(ybin1c),np.max(ybin1c)],
                cmap=cm.bam,norm=colors.Normalize(vmin=-1,vmax=1),aspect='auto')
cbar2=plt.colorbar(im2,ax=ax2)
cbar2.ax.set_ylabel(r'Max Z to $\bar{R}$')
plt.imshow(P1c,origin='lower',
                extent=[np.min(xbin1c),np.max(xbin1c),np.min(ybin1c),np.max(ybin1c)],
                cmap=cm.grayC,norm=colors.Normalize(vmin=0,vmax=1),aspect='auto')
# plt.xlabel('Mean Precip [mm/day]')
# plt.ylabel('MAT [C]')
ax2.text(0.90, 0.99, 'B',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax2.transAxes,
        fontsize=12,fontweight='extra bold')

ax3=f2.add_subplot(gs[0,2])
im3=plt.imshow(X1a,origin='lower',
                extent=[np.min(xbin1a),np.max(xbin1a),np.min(ybin1a),np.max(ybin1a)],
                cmap=cm.bam,norm=colors.Normalize(vmin=-1,vmax=1),aspect='auto')
cbar3=plt.colorbar(im3,ax=ax3)
cbar3.ax.set_ylabel(r'Local Relief to $\bar{R}$')
plt.imshow(P1a,origin='lower',
                extent=[np.min(xbin1a),np.max(xbin1a),np.min(ybin1a),np.max(ybin1a)],
                cmap=cm.grayC,norm=colors.Normalize(vmin=0,vmax=1),aspect='auto')
# plt.xlabel('Mean Precip [mm/day]')
# plt.ylabel('MAT [C]')
ax3.text(0.90, 0.99, 'C',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax3.transAxes,
        fontsize=12,fontweight='extra bold')

ax4=f2.add_subplot(gs[1,0])
im4=plt.imshow(X2b,origin='lower',
                extent=[np.min(xbin2b),np.max(xbin2b),np.min(ybin2b),np.max(ybin2b)],
                cmap=cm.bam,norm=colors.Normalize(vmin=-1,vmax=1),aspect='auto')
cbar4=plt.colorbar(im4,ax=ax4)
cbar4.ax.set_ylabel('Mean Z to SF')
plt.imshow(P2b,origin='lower',
                extent=[np.min(xbin2b),np.max(xbin2b),np.min(ybin2b),np.max(ybin2b)],
                cmap=cm.grayC,norm=colors.Normalize(vmin=0,vmax=1),aspect='auto')
plt.xlabel('Mean Precip [mm/day]')
plt.ylabel('MAT [C]')
ax4.text(0.90, 0.99, 'D',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax4.transAxes,
        fontsize=12,fontweight='extra bold')

ax5=f2.add_subplot(gs[1,1])
im5=plt.imshow(X2a,origin='lower',
                extent=[np.min(xbin2a),np.max(xbin2a),np.min(ybin2a),np.max(ybin2a)],
                cmap=cm.bam,norm=colors.Normalize(vmin=-1,vmax=1),aspect='auto')
cbar5=plt.colorbar(im5,ax=ax5)
cbar5.ax.set_ylabel('Max Z to SF')
plt.imshow(P2a,origin='lower',
                extent=[np.min(xbin2a),np.max(xbin2a),np.min(ybin2a),np.max(ybin2a)],
                cmap=cm.grayC,norm=colors.Normalize(vmin=0,vmax=1),aspect='auto')
plt.xlabel('Mean Precip [mm/day]')
# plt.ylabel('MAT [C]')
ax5.text(0.90, 0.99, 'E',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax5.transAxes,
        fontsize=12,fontweight='extra bold')

ax6=f2.add_subplot(gs[1,2])
im6=plt.imshow(X2c,origin='lower',
                extent=[np.min(xbin2c),np.max(xbin2c),np.min(ybin2c),np.max(ybin2c)],
                cmap=cm.bam,norm=colors.Normalize(vmin=-1,vmax=1),aspect='auto')
cbar6=plt.colorbar(im6,ax=ax6)
cbar6.ax.set_ylabel('Local Relief to SF')
plt.imshow(P2c,origin='lower',
                extent=[np.min(xbin2c),np.max(xbin2c),np.min(ybin2c),np.max(ybin2c)],
                cmap=cm.grayC,norm=colors.Normalize(vmin=0,vmax=1),aspect='auto')
plt.xlabel('Mean Precip [mm/day]')
# plt.ylabel('MAT [C]')
ax6.text(0.90, 0.99, 'F',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax6.transAxes,
        fontsize=12,fontweight='extra bold')

plt.tight_layout()


f3=plt.figure(3,figsize=(4,3))
f3.set_dpi(250)
ax31=f3.add_subplot(1,1,1)
xb=np.histogram_bin_edges(df_global_s['mean_precip'],'auto')
yb=np.histogram_bin_edges(df_global_s['mat'],'auto')

sc1=plt.hist2d(df_global_s['mean_precip'],df_global_s['mat'],[xb,yb],norm=colors.LogNorm(),cmap=cm.lajolla)
for i in range(len(xbin1a)):
    ax31.axvline(xbin1a[i],c='k',linestyle=':',linewidth=0.25,alpha=0.50)
for i in range(len(ybin1a)):
    ax31.axhline(ybin1a[i],c='k',linestyle=':',linewidth=0.25,alpha=0.50)
cbar1=plt.colorbar(sc1[3],ax=ax31)
cbar1.ax.set_ylabel('Density')
plt.xlim((np.min(xbin1a),np.max(xbin1a)))
plt.ylim((np.min(ybin1a),np.max(ybin1a)))
plt.xlabel('Mean Precip [mm/day]')
plt.ylabel('MAT [C]')

plt.rcdefaults()
# f1.savefig('figure_x2a.pdf',dpi='figure')
f2.savefig('P1_figure7.pdf',dpi='figure')
# f3.savefig('figure_x2c.pdf',dpi='figure')



