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
master_location='/Volumes/Choruh/Data/snowmelt_project/'
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

[RHO2a,NUM2a,BN2a,X2a,P2a,xbin2a,ybin2a]=bin2d_corr(2,df_global_s,df_global_s['mean_precip'],
                                            df_global_s['mat'],df_global_s['max_z'],df_global_s['qsm']/df_global_s['mean_runoff'],10,0.01,'spearman')

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

# ### OUTPUT RASTERS
# # Load dimensions of underlying WRR2 dataset to generate transform
# fh=Dataset(master_location+'wrr2_derived_new.nc')
# # Lat - Lon
# lat=fh.variables['lat'][:]
# lon=fh.variables['lon'][:]
# fh.close()
# # Create a transform
# res=lat[1]-lat[0]
# [x,y]=np.meshgrid(lon,lat)
# transform=rt.Affine.translation(x[0][0]-res/2,y[0][0]-res/2)*rt.Affine.scale(res,res)
# # Generate row and column
# row_col=rt.rowcol(transform,df_global_s['longitude'],df_global_s['latitude'])
# rows=row_col[0]
# cols=row_col[1]
# # Generate blank images and fill
# RHO1=np.zeros(x.shape)
# RHO2=np.zeros(x.shape)
# RHO1[:,:]=-9999
# RHO2[:,:]=-9999
# RHO1[rows,cols]=RHO1a
# RHO2[rows,cols]=RHO2a

# with rasterio.open(
#         master_location+'wrr2_raster_outputs/pmat_rlf_runoff_corr.tif',
#         'w',
#         driver='GTiff',
#         height=RHO1.shape[0],
#         width=RHO1.shape[1],
#         count=1,
#         dtype=RHO1.dtype,
#         crs='+proj=latlong',
#         transform=transform,
#         nodata=-9999,
#     ) as dst: 
#         dst.write(RHO1,1)

# with rasterio.open(
#         master_location+'wrr2_raster_outputs/pmat_maxz_snow_corr.tif',
#         'w',
#         driver='GTiff',
#         height=RHO2.shape[0],
#         width=RHO2.shape[1],
#         count=1,
#         dtype=RHO2.dtype,
#         crs='+proj=latlong',
#         transform=transform,
#         nodata=-9999,
#     ) as dst: 
#         dst.write(RHO2,1)
        
raster=master_location+'srtm30_plus/topo30_wrap.tif'
with rasterio.open(raster,'r') as src:
    raster_crs=src.crs
    left1, bottom1, right1, top1 = src.bounds
    dem=src.read()[0,:,:]
    
raster=master_location+'wrr2_raster_outputs/pmat_rlf_runoff_corr.tif'
with rasterio.open(raster,'r') as src:
    raster_crs=src.crs
    left, bottom, right, top = src.bounds
    rlfrun=src.read()[0,:,:]    

raster=master_location+'wrr2_raster_outputs/pmat_maxz_snow_corr.tif'
with rasterio.open(raster,'r') as src:
    raster_crs=src.crs
    left, bottom, right, top = src.bounds
    zsnow=src.read()[0,:,:]
    
rlfrun[rlfrun==-9999]=np.nan
zsnow[zsnow==-9999]=np.nan    
    
f2=plt.figure(2,figsize=(8,6)) 
f2.set_dpi(250)
gs=f2.add_gridspec(2,1)

ax1=f2.add_subplot(gs[0,0],projection=ccrs.PlateCarree())
ax1.set_extent([-180, 180, -90, 90], crs=ccrs.PlateCarree())
gl=ax1.gridlines(draw_labels=True,alpha=0.5,linewidth=0.25)
gl.right_labels=False
gl.bottom_labels=False
ax1.coastlines(resolution='auto',color='k',linewidth=0.25)
ax1.imshow(dem,extent=(left1,right1,bottom1,top1),cmap=cm.grayC,vmin=-4000,vmax=4000,alpha=0.5)
im1=ax1.imshow(rlfrun,extent=(left,right,bottom,top),cmap=cm.bam,vmin=-1,vmax=1)    
cbar1=plt.colorbar(im1,ax=ax1)
cbar1.ax.set_ylabel(r'Local Relief to $\bar{R}$ Correlation')
ax1.text(0.01, 0.99, 'A',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax1.transAxes,
        fontsize=12,fontweight='extra bold')

ax2=f2.add_subplot(gs[1,0],projection=ccrs.PlateCarree())
ax2.set_extent([-180, 180, -90, 90], crs=ccrs.PlateCarree())
gl=ax2.gridlines(draw_labels=True,alpha=0.5,linewidth=0.25)
gl.right_labels=False
gl.bottom_labels=False
ax2.coastlines(resolution='auto',color='k',linewidth=0.25)
ax2.imshow(dem,extent=(left1,right1,bottom1,top1),cmap=cm.grayC,vmin=-4000,vmax=4000,alpha=0.5)
im2=ax2.imshow(zsnow,extent=(left,right,bottom,top),cmap=cm.bam,vmin=-1,vmax=1)
cbar2=plt.colorbar(im2,ax=ax2)
cbar2.ax.set_ylabel('Max Z to SF Correlation')
ax2.text(0.01, 0.99, 'B',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax2.transAxes,
        fontsize=12,fontweight='extra bold')

plt.tight_layout()
plt.rcdefaults()

f1.savefig('figure_x2a.pdf',dpi='figure')
f2.savefig('figure_x2b.pdf',dpi='figure')

# def bin2d_corr_fig(fn,dfs,x,y,x1,y1,xname,yname,corr_title,min_num,p_sig,corr_type,asp='auto',one_to_one=False,save_fig=True,bin_type='doane'):
    
#     # Test overwriting bins
#     x_idx=x<=np.percentile(x,99.9)
#     y_idx=y<=np.percentile(y,99.9)
#     xbin=np.histogram_bin_edges(x[x_idx],bin_type)
#     ybin=np.histogram_bin_edges(y[y_idx],bin_type)    
    
#     x_dim=len(xbin)-1
#     y_dim=len(ybin)-1
    
#     # Allocate
#     X=np.zeros((y_dim,x_dim))
#     P=np.ones((y_dim,x_dim))*0.25
#     RHO=np.zeros(x.shape)*np.nan
#     NUM=np.zeros(x.shape)*np.nan
#     BN=np.zeros(x.shape)*np.nan
    
#     x_ix=np.digitize(x,xbin)
#     y_ix=np.digitize(y,ybin)

#     for i in range(x_dim):
#         for j in range(y_dim):
#             idx_in=np.logical_and(x_ix==i+1,y_ix==j+1)
#             NUM[idx_in]=np.sum(idx_in)
#             if np.sum(idx_in)>min_num:
#                 if corr_type=='pearson':
#                     [rho,p]=pearsonr(x1[idx_in],y1[idx_in])
#                 elif corr_type=='spearman':
#                     [rho,p]=spearmanr(x1[idx_in],y1[idx_in])
                    
#                 if p<p_sig:
#                     X[j,i]=rho
#                     P[j,i]=np.nan
#                     RHO[idx_in]=rho
#                     BN[idx_in]=j+i
#                 else:
#                     X[j,i]=np.nan
#                     P[j,i]=1

#             else:
#                 X[j,i]=np.nan
                    
#     SMALL_SIZE = 8
#     MEDIUM_SIZE = 10
#     BIGGER_SIZE = 12
#     plt.rc('font', size=SMALL_SIZE,family='Futura')          # controls default text sizes
#     plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
#     plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
#     plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
#     plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
#     plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
#     plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
                          
#     f=plt.figure(fn,figsize=(8,8))
 
#     gs=f.add_gridspec(2,2)   
 
#     ax1=f.add_subplot(gs[0,0])
#     xb=np.histogram_bin_edges(x,'auto')
#     yb=np.histogram_bin_edges(y,'auto')
    
#     sc1=plt.hist2d(x,y,[xb,yb],norm=colors.LogNorm(),cmap=cm.lajolla)
#     for i in range(len(xbin)):
#         ax1.axvline(xbin[i],c='k',linestyle=':',linewidth=0.25,alpha=0.50)
#     for i in range(len(ybin)):
#         ax1.axhline(ybin[i],c='k',linestyle=':',linewidth=0.25,alpha=0.50)
#     cbar1=plt.colorbar(sc1[3],ax=ax1)
#     cbar1.ax.set_ylabel('Density')
#     plt.xlim((np.min(xbin),np.max(xbin)))
#     plt.ylim((np.min(ybin),np.max(ybin)))
#     plt.xlabel(xname)
#     plt.ylabel(yname)

#     ax2=f.add_subplot(gs[0,1])
#     sc2=plt.imshow(X,origin='lower',
#                    extent=[np.min(xbin),np.max(xbin),np.min(ybin),np.max(ybin)],
#                    cmap=cm.bam,norm=colors.Normalize(vmin=-1,vmax=1),aspect=asp)
#     cbar3=plt.colorbar(sc2,ax=ax2)
#     plt.imshow(P,origin='lower',
#                    extent=[np.min(xbin),np.max(xbin),np.min(ybin),np.max(ybin)],
#                    cmap=cm.grayC,norm=colors.Normalize(vmin=0,vmax=1),aspect=asp)
    
#     if one_to_one:
#         plt.plot(ybin,ybin,c='k',linestyle=':')
    
#     if corr_type=='pearson':
#         cbar3.ax.set_ylabel('Pearson Correlation Coeffecient') 
#     elif corr_type=='spearman':
#         cbar3.ax.set_ylabel('Spearman Correlation Coeffecient')        
#     plt.xlabel(xname)
#     plt.ylabel(yname)
#     plt.title(corr_title)

#     ax3=f.add_subplot(gs[1,:],projection=ccrs.PlateCarree())
#     ax3.set_extent([-180, 180, -90, 90], crs=ccrs.PlateCarree())
#     ax3.gridlines(draw_labels=True)
#     ax3.coastlines(resolution='auto',color='k')  
#     plt.scatter(dfs['longitude'],dfs['latitude'],s=0.25,marker='s',c=RHO,cmap=cm.bam,norm=colors.Normalize(vmin=-1,vmax=1))
    
#     plt.tight_layout()
#     plt.rcdefaults()
  
    
#     if save_fig:
#         # gs.tight_layout(f)
#         f.savefig('/Users/aforte/Desktop/bin_figures/bin_correlate_'+str(fn)+'.png',bbox_inches='tight')
        
        
#     return RHO,NUM,BN
    
# ## Load global
# master_location='/Volumes/Choruh/Data/snowmelt_project/'
# df_global=pd.read_csv(master_location+'wrr2_derived_data_v3.csv')
# df_global=df_global.drop(index=df_global.index[np.isnan(df_global['mean_z'])])
# df_global=df_global.reset_index(drop=True)

# percb_cutoff=0.25
# grlf=df_global['max_z']-df_global['min_z']
# global_perc_base=df_global['qsb']/df_global['mean_runoff']

# rem_idx=(grlf<=500) | (df_global['mean_z']<=250) | (global_perc_base>percb_cutoff)
# df_global_s=df_global.drop(df_global.index[rem_idx])
# df_global_s=df_global_s.reset_index(drop=True)

# [RHO1a,NUM1a,BN1a]=bin2d_corr_fig(1,df_global_s,df_global_s['mean_precip'],df_global_s['mat'],df_global_s['mean_rlf'],df_global_s['mean_runoff'],
#                 'Mean Precip [mm/day]','MAT [C]','Mean Runoff to Mean Relief',10,0.01,'spearman',save_fig=False)


# [RHO2a,NUM2a,BN2a]=bin2d_corr_fig(2,df_global_s,df_global_s['mean_precip'],df_global_s['mat'],df_global_s['max_z'],df_global_s['qsm']/df_global_s['mean_runoff'],
#                 'Mean Precip [mm/day]','MAT [C]','Snowmelt Fraction to Max Elevation',10,0.01,'spearman',save_fig=False)

