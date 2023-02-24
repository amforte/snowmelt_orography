#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script written by Adam M. Forte
aforte8@lsu.edu

This script reads the output of 'wrr2_derived_make_table_v3.py' and 'wrr2_block_proc_v2.py' 
to calculate global binned correlations between relief and mean runoff and maximum elevation and
snowmelt.
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


def bin2d_rolling(dfs,x,y,x1,y1,min_num,p_sig,corr_type,xbin,ybin,bin_size,shift):
    
    # Define shift vector
    shift_vec=np.arange(0,bin_size,shift)
    # Pad bins
    xbin=np.concatenate(([xbin[0]-bin_size],xbin),axis=0)
    ybin=np.concatenate(([ybin[0]-bin_size],ybin),axis=0)
    # Define types
    x_dim=len(xbin)-1
    y_dim=len(ybin)-1
    # Num shifts
    num_shift=len(shift_vec)
    # Allocate container
    RHO=np.zeros((x.shape[0],num_shift**2))*np.nan
    NUM=np.zeros((x.shape[0],num_shift**2))*np.nan
    # Start counter
    cnt=0
    
    # Begin rolling
    for k in range(num_shift):
        for l in range(num_shift):
            print(f'Rolling - shifting x bins by {shift_vec[k]:1.3g} - shifting y bins by {shift_vec[l]:1.3g}')
            # Shift bins
            xbint=xbin+shift_vec[k]
            ybint=ybin+shift_vec[l]
            # Digitize
            x_ix=np.digitize(x,xbint)
            y_ix=np.digitize(y,ybint)
            
            for i in range(x_dim):
                for j in range(y_dim):
                    idx_in=np.logical_and(x_ix==i+1,y_ix==j+1)
                    NUM[idx_in,cnt]=np.sum(idx_in)
                    if np.sum(idx_in)>min_num:
                        if corr_type=='pearson':
                            [rho,p]=pearsonr(x1[idx_in],y1[idx_in])
                        elif corr_type=='spearman':
                            [rho,p]=spearmanr(x1[idx_in],y1[idx_in])
                            
                        if p<p_sig:
                            RHO[idx_in,cnt]=rho
                        # else:
                        #     RHO[idx_in,cnt]=np.nan
            cnt=+1
            
    
            
    return RHO,NUM

## Load global
master_location='/Volumes/Choruh/Data/snowmelt_project/'
df_global=pd.read_csv(master_location+'wrr2_derived_data_v3.csv')
df_global=df_global.drop(index=df_global.index[np.isnan(df_global['mean_z'])])
df_global=df_global.reset_index(drop=True)

# Parse global data
percb_cutoff=0.25
grlf=df_global['max_z']-df_global['min_z']
global_perc_base=df_global['qsb']/df_global['mean_runoff']

rem_idx=(grlf<=500) | (df_global['mean_z']<=250) | (global_perc_base>percb_cutoff) | (df_global['mean_rlf']<250)
df_global_s=df_global.drop(df_global.index[rem_idx])
df_global_s=df_global_s.reset_index(drop=True)

xbin=np.arange(-180,180,2)
ybin=np.arange(-90,90,2)

[RHO1,NUM1]=bin2d_rolling(df_global_s,df_global_s['longitude'],df_global_s['latitude'],df_global_s['mean_rlf'],df_global_s['mean_runoff'],
              10,0.01,'spearman',xbin,ybin,2,0.25)


[RHO2,NUM2]=bin2d_rolling(df_global_s,df_global_s['longitude'],df_global_s['latitude'],df_global_s['max_z'],df_global_s['qsm']/df_global_s['mean_runoff'],
              10,0.01,'spearman',xbin,ybin,2,0.25)


perc_nan1=np.sum(np.isnan(RHO1),axis=1)/RHO1.shape[1]
perc_nan2=np.sum(np.isnan(RHO2),axis=1)/RHO2.shape[1]

mn_rho1=np.nanmean(RHO1,axis=1)
std_rho1=np.nanstd(RHO1,axis=1)
mn_rho2=np.nanmean(RHO2,axis=1)
std_rho2=np.nanstd(RHO2,axis=1)

idx1=np.all(np.isnan(RHO1),axis=1)
idx2=np.all(np.isnan(RHO2),axis=1)
mn_rho1[idx1]=0
mn_rho2[idx2]=0

plt.figure(figsize=(10,10))
plt.hist2d(mn_rho1,mn_rho2,[np.linspace(-1,1,50),np.linspace(-1,1,50)],norm=colors.LogNorm(),cmap=cm.lajolla)
plt.xlabel('Relief to Runoff')
plt.ylabel('Max Elevation to Snowmelt')
plt.xlim((-1,1))
plt.ylim((-1,1))

f1=plt.figure(figsize=(15,20))
gs=f1.add_gridspec(2,1)
 
ax1=f1.add_subplot(gs[0],projection=ccrs.PlateCarree())
ax1.set_extent([-180, 180, -90, 90], crs=ccrs.PlateCarree())
ax1.gridlines(draw_labels=True)
ax1.coastlines(resolution='auto',color='k')  
plt.scatter(df_global_s['longitude'],df_global_s['latitude'],s=2,c=mn_rho1,cmap=cm.bam,norm=colors.Normalize(vmin=-1,vmax=1))

ax1=f1.add_subplot(gs[1],projection=ccrs.PlateCarree())
ax1.set_extent([-180, 180, -90, 90], crs=ccrs.PlateCarree())
ax1.gridlines(draw_labels=True)
ax1.coastlines(resolution='auto',color='k')  
plt.scatter(df_global_s['longitude'],df_global_s['latitude'],s=2,c=mn_rho2,cmap=cm.bam,norm=colors.Normalize(vmin=-1,vmax=1))


### OUTPUT RASTERS
# Load dimensions of underlying WRR2 dataset to generate transform
fh=Dataset(master_location+'wrr2_derived_new.nc')
# Lat - Lon
lat=fh.variables['lat'][:]
lon=fh.variables['lon'][:]
fh.close()
# Create a transform
res=lat[1]-lat[0]
[x,y]=np.meshgrid(lon,lat)
transform=rt.Affine.translation(x[0][0]-res/2,y[0][0]-res/2)*rt.Affine.scale(res,res)
# Generate row and column
row_col=rt.rowcol(transform,df_global_s['longitude'],df_global_s['latitude'])
rows=row_col[0]
cols=row_col[1]
# Generate blank images and fill
RHO1=np.zeros(x.shape)
RHO2=np.zeros(x.shape)
RHO1[:,:]=-9999
RHO2[:,:]=-9999
RHO1[rows,cols]=mn_rho1
RHO2[rows,cols]=mn_rho2

with rasterio.open(
        master_location+'wrr2_raster_outputs/rlf_runoff_corr.tif',
        'w',
        driver='GTiff',
        height=RHO1.shape[0],
        width=RHO1.shape[1],
        count=1,
        dtype=RHO1.dtype,
        crs='+proj=latlong',
        transform=transform,
        nodata=-9999,
    ) as dst: 
        dst.write(RHO1,1)

with rasterio.open(
        master_location+'wrr2_raster_outputs/maxz_snow_corr.tif',
        'w',
        driver='GTiff',
        height=RHO2.shape[0],
        width=RHO2.shape[1],
        count=1,
        dtype=RHO2.dtype,
        crs='+proj=latlong',
        transform=transform,
        nodata=-9999,
    ) as dst: 
        dst.write(RHO2,1)
