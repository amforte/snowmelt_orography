#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  6 10:26:56 2023

@author: aforte
"""

import numpy as np
import pandas as pd
import rasterio
import matplotlib.gridspec as gridspec
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
from cmcrameri import cm
from matplotlib import colors
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy import odr

def linear(B,x):
    return B[0]*x + B[1]

def odr_fit_ref(x,y):
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


def plot_and_fit_ref(axn,dfs,bl,br,bb,bt,col,lbl_left,lbl_rgt):
    spidx=(dfs['latitude']>=bb) & (dfs['latitude']<=bt) & (dfs['longitude']>=bl) & (dfs['longitude']<=br)
    
    # max_z=np.max(np.ceil(dfs.loc[spidx,'max_z'].to_numpy()))
    # max_r=np.max(np.ceil(dfs.loc[spidx,'mean_runoff'].to_numpy()))
    # max_rlf=np.max(np.ceil(dfs.loc[spidx,'mean_rlf'].to_numpy()))
    max_z=4500
    max_r=10
    
    z=np.linspace(1,max_z,100)
    
    x=dfs.loc[spidx,'mean_rlf'].to_numpy()
    y=dfs.loc[spidx,'mean_runoff'].to_numpy()
    co,ex=odr_fit_ref(x,y)
    bin_edges=np.histogram_bin_edges(dfs.loc[spidx,'mean_runoff'],'doane')
    bix=np.digitize(dfs.loc[spidx,'mean_runoff'],bin_edges)
    for i in range(len(bin_edges)-1):
        mx=np.mean(x[bix==i])
        stdx=np.std(x[bix==i])
        my=np.mean(y[bix==i])
        stdy=np.std(y[bix==i])
        axn.scatter(mx,my,c=col,s=len(x[bix==i]))
        axn.errorbar(mx,my,xerr=stdx,yerr=stdy,ecolor=col,elinewidth=0.5,zorder=0)
    axn.plot(z,co*z**ex,c=col,label='Local Relief to Mean Runoff')  
    axn.set_xlabel('Topography [m]')
    if lbl_left:
        axn.set_ylabel('Mean Runoff [mm/day]')
    axn.set_xlim((0,max_z))
    axn.set_ylim((0,max_r))

    axnt=axn.twinx()
    snmp=dfs.loc[spidx,'qsm']/dfs.loc[spidx,'mean_runoff']
    x=dfs.loc[spidx,'max_z'].to_numpy()
    y=snmp.to_numpy()
    co,ex=odr_fit_ref(x,y)
    bin_edges=np.histogram_bin_edges(snmp,'doane')
    bix=np.digitize(snmp,bin_edges)
    for i in range(len(bin_edges)-1):
        mx=np.mean(x[bix==i])
        stdx=np.std(x[bix==i])
        my=np.mean(y[bix==i])
        stdy=np.std(y[bix==i])
        axnt.scatter(mx,my,c=col,s=len(x[bix==i]),marker='s')
        axnt.errorbar(mx,my,xerr=stdx,yerr=stdy,ecolor=col,elinewidth=0.5,zorder=0)
    axnt.plot(z,co*z**ex,c=col,linestyle='--',label='Max. Elev. to Snowmelt Fraction')
    axnt.set_ylim((0,1))
    if lbl_rgt:
        axnt.set_ylabel('Snowmelt Fraction')
        
    h1,l1=axn.get_legend_handles_labels()
    h2,l2=axnt.get_legend_handles_labels()
    axnt.legend(h1+h2,l1+l2,bbox_to_anchor= (-0.1,-0.3),loc='lower left')


# Set master location
master_location='/Volumes/Choruh/Data/snowmelt_project/'

raster=master_location+'srtm30_plus/topo30_wrap.tif'
with rasterio.open(raster,'r') as src:
    raster_crs=src.crs
    left1, bottom1, right1, top1 = src.bounds
    dem=src.read()[0,:,:]
    
raster=master_location+'wrr2_raster_outputs/rlf_runoff_corr.tif'
with rasterio.open(raster,'r') as src:
    raster_crs=src.crs
    left, bottom, right, top = src.bounds
    rlfrun=src.read()[0,:,:]    

raster=master_location+'wrr2_raster_outputs/maxz_snow_corr.tif'
with rasterio.open(raster,'r') as src:
    raster_crs=src.crs
    left, bottom, right, top = src.bounds
    zsnow=src.read()[0,:,:]

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
    
rlfrun[rlfrun==-9999]=np.nan
zsnow[zsnow==-9999]=np.nan    
    
f1=plt.figure(1,figsize=(8,11)) 
f1.set_dpi(250)
gs=f1.add_gridspec(3,3)

ax1=f1.add_subplot(gs[0,:],projection=ccrs.PlateCarree())
ax1.set_extent([-180, 180, -90, 90], crs=ccrs.PlateCarree())
gl=ax1.gridlines(draw_labels=True,alpha=0.5,linewidth=0.25)
gl.right_labels=False
gl.bottom_labels=False
ax1.coastlines(resolution='auto',color='k',linewidth=0.25)
ax1.imshow(dem,extent=(left1,right1,bottom1,top1),cmap=cm.grayC,vmin=-4000,vmax=4000,alpha=0.5)
im1=ax1.imshow(rlfrun,extent=(left,right,bottom,top),cmap=cm.vik_r,vmin=-1,vmax=1)

## GC INSET
# Plot Box
ax1.plot([bl1,bl1,br1,br1,bl1],[bb1,bt1,bt1,bb1,bb1],c=gc_col,linewidth=0.5)
# Set Insent
axins1 = inset_axes(ax1, width="50%", height="50%", loc='lower left',
                    bbox_to_anchor=(0.65, 0.1, 0.5, 0.5),
                    bbox_transform=ax1.transAxes,
                    axes_class=cartopy.mpl.geoaxes.GeoAxes,
                    axes_kwargs=dict(map_projection=ccrs.PlateCarree()))
axins1.set_extent([bl1,br1,bb1,bt1])
axins1.imshow(dem,extent=(left1,right1,bottom1,top1),cmap=cm.grayC,vmin=-4000,vmax=4000,alpha=0.5)
axins1.imshow(rlfrun,extent=(left,right,bottom,top),cmap=cm.vik_r,vmin=-1,vmax=1)
axins1.spines['geo'].set(edgecolor=gc_col)
# Draw Connectors
rect, connectors1 = ax1.indicate_inset_zoom(axins1, edgecolor=gc_col, alpha=0.5, transform=ax1.transAxes)
connectors1[0].set_visible(False)
connectors1[1].set_visible(True)
connectors1[3].set_visible(True)
connectors1[2].set_visible(False)
## Alps Inset
# Plot Box
ax1.plot([bl2,bl2,br2,br2,bl2],[bb2,bt2,bt2,bb2,bb2],c=alps_col,linewidth=0.5)
# Set Insent
axins2 = inset_axes(ax1, width="50%", height="50%", loc='lower left',
                    bbox_to_anchor=(0.30, 0.1, 0.5, 0.5),
                    bbox_transform=ax1.transAxes,
                    axes_class=cartopy.mpl.geoaxes.GeoAxes,
                    axes_kwargs=dict(map_projection=ccrs.PlateCarree()))
axins2.set_extent([bl2,br2,bb2,bt2])
axins2.imshow(dem,extent=(left1,right1,bottom1,top1),cmap=cm.grayC,vmin=-4000,vmax=4000,alpha=0.5)
axins2.imshow(rlfrun,extent=(left,right,bottom,top),cmap=cm.vik_r,vmin=-1,vmax=1)
axins2.spines['geo'].set(edgecolor=alps_col)
# Draw Connectors
rect, connectors2 = ax1.indicate_inset_zoom(axins2, edgecolor=alps_col, alpha=0.5, transform=ax1.transAxes)
connectors2[0].set_visible(False)
connectors2[1].set_visible(True)
connectors2[3].set_visible(True)
connectors2[2].set_visible(False)
## BC Inset
# Plot Box
ax1.plot([bl3,bl3,br3,br3,bl3],[bb3,bt3,bt3,bb3,bb3],c=bc_col,linewidth=0.5)
# Set Insent
axins3 = inset_axes(ax1, width="50%", height="50%", loc='lower left',
                    bbox_to_anchor=(0.01, 0.1, 0.5, 0.5),
                    bbox_transform=ax1.transAxes,
                    axes_class=cartopy.mpl.geoaxes.GeoAxes,
                    axes_kwargs=dict(map_projection=ccrs.PlateCarree()))
axins3.set_extent([bl3,br3,bb3,bt3])
axins3.imshow(dem,extent=(left1,right1,bottom1,top1),cmap=cm.grayC,vmin=-4000,vmax=4000,alpha=0.5)
axins3.imshow(rlfrun,extent=(left,right,bottom,top),cmap=cm.vik_r,vmin=-1,vmax=1)
axins3.spines['geo'].set(edgecolor=bc_col)
# Draw Connectors
rect, connectors3 = ax1.indicate_inset_zoom(axins3, edgecolor=bc_col, alpha=0.5, transform=ax1.transAxes)
connectors3[0].set_visible(False)
connectors3[1].set_visible(True)
connectors3[3].set_visible(True)
connectors3[2].set_visible(False)

cbar1=plt.colorbar(im1,ax=ax1)
cbar1.ax.set_ylabel('Local Relief to Runoff Correlation')
ax1.text(0.01, 0.99, 'A',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax1.transAxes,
        fontsize=12,fontweight='extra bold')

ax2=f1.add_subplot(gs[1,:],projection=ccrs.PlateCarree())
ax2.set_extent([-180, 180, -90, 90], crs=ccrs.PlateCarree())
gl=ax2.gridlines(draw_labels=True,alpha=0.5,linewidth=0.25)
gl.right_labels=False
gl.bottom_labels=False
ax2.coastlines(resolution='auto',color='k',linewidth=0.25)
ax2.imshow(dem,extent=(left1,right1,bottom1,top1),cmap=cm.grayC,vmin=-4000,vmax=4000,alpha=0.5)
im2=ax2.imshow(zsnow,extent=(left,right,bottom,top),cmap=cm.vik_r,vmin=-1,vmax=1)

## GC INSET
# Plot Box
ax2.plot([bl1,bl1,br1,br1,bl1],[bb1,bt1,bt1,bb1,bb1],c=gc_col,linewidth=0.5)
# Set Insent
axins1 = inset_axes(ax2, width="50%", height="50%", loc='lower left',
                    bbox_to_anchor=(0.65, 0.1, 0.5, 0.5),
                    bbox_transform=ax2.transAxes,
                    axes_class=cartopy.mpl.geoaxes.GeoAxes,
                    axes_kwargs=dict(map_projection=ccrs.PlateCarree()))
axins1.set_extent([bl1,br1,bb1,bt1])
axins1.imshow(dem,extent=(left1,right1,bottom1,top1),cmap=cm.grayC,vmin=-4000,vmax=4000,alpha=0.5)
axins1.imshow(zsnow,extent=(left,right,bottom,top),cmap=cm.vik_r,vmin=-1,vmax=1)
axins1.spines['geo'].set(edgecolor=gc_col)
# Draw Connectors
rect, connectors1 = ax2.indicate_inset_zoom(axins1, edgecolor=gc_col, alpha=0.5, transform=ax2.transAxes)
connectors1[0].set_visible(False)
connectors1[1].set_visible(True)
connectors1[3].set_visible(True)
connectors1[2].set_visible(False)
## Alps Inset
# Plot Box
ax2.plot([bl2,bl2,br2,br2,bl2],[bb2,bt2,bt2,bb2,bb2],c=alps_col,linewidth=0.5)
# Set Insent
axins2 = inset_axes(ax2, width="50%", height="50%", loc='lower left',
                    bbox_to_anchor=(0.30, 0.1, 0.5, 0.5),
                    bbox_transform=ax2.transAxes,
                    axes_class=cartopy.mpl.geoaxes.GeoAxes,
                    axes_kwargs=dict(map_projection=ccrs.PlateCarree()))
axins2.set_extent([bl2,br2,bb2,bt2])
axins2.imshow(dem,extent=(left1,right1,bottom1,top1),cmap=cm.grayC,vmin=-4000,vmax=4000,alpha=0.5)
axins2.imshow(zsnow,extent=(left,right,bottom,top),cmap=cm.vik_r,vmin=-1,vmax=1)
axins2.spines['geo'].set(edgecolor=alps_col)
# Draw Connectors
rect, connectors2 = ax2.indicate_inset_zoom(axins2, edgecolor=alps_col, alpha=0.5, transform=ax2.transAxes)
connectors2[0].set_visible(False)
connectors2[1].set_visible(True)
connectors2[3].set_visible(True)
connectors2[2].set_visible(False)
## BC Inset
# Plot Box
ax2.plot([bl3,bl3,br3,br3,bl3],[bb3,bt3,bt3,bb3,bb3],c=bc_col,linewidth=0.5)
# Set Insent
axins3 = inset_axes(ax2, width="50%", height="50%", loc='lower left',
                    bbox_to_anchor=(0.01, 0.1, 0.5, 0.5),
                    bbox_transform=ax2.transAxes,
                    axes_class=cartopy.mpl.geoaxes.GeoAxes,
                    axes_kwargs=dict(map_projection=ccrs.PlateCarree()))
axins3.set_extent([bl3,br3,bb3,bt3])
axins3.imshow(dem,extent=(left1,right1,bottom1,top1),cmap=cm.grayC,vmin=-4000,vmax=4000,alpha=0.5)
axins3.imshow(zsnow,extent=(left,right,bottom,top),cmap=cm.vik_r,vmin=-1,vmax=1)
axins3.spines['geo'].set(edgecolor=bc_col)
# Draw Connectors
rect, connectors3 = ax2.indicate_inset_zoom(axins3, edgecolor=bc_col, alpha=0.5, transform=ax2.transAxes)
connectors3[0].set_visible(False)
connectors3[1].set_visible(True)
connectors3[3].set_visible(True)
connectors3[2].set_visible(False)

cbar2=plt.colorbar(im2,ax=ax2)
cbar2.ax.set_ylabel('Maximum Elevation to Snowmelt Fraction Correlation')
ax2.text(0.01, 0.99, 'B',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax2.transAxes,
        fontsize=12,fontweight='extra bold') 

## Load global
df_global=pd.read_csv(master_location+'wrr2_derived_data_v3.csv')
df_global=df_global.drop(index=df_global.index[np.isnan(df_global['mean_z'])])
df_global=df_global.reset_index(drop=True)
# Calc percents
global_perc_base=df_global['qsb']/df_global['mean_runoff']
# Calculate indices
grlf=df_global['max_z']-df_global['min_z']
# Set cutoffs
percb_cutoff=0.25
# Index
rem_idx=(grlf<=500) | (df_global['mean_z']<=250) | (global_perc_base>percb_cutoff) | (df_global['mean_rlf']<250)
df_global_s=df_global.drop(df_global.index[rem_idx])
df_global_s=df_global_s.reset_index(drop=True)
global_perc_snow=df_global_s['qsm']/df_global_s['mean_runoff']


ax3=f1.add_subplot(gs[2,0])
plot_and_fit_ref(ax3,df_global_s,bl3,br3,bb3,bt3,bc_col,True,False)
ax3.set_title('British Columbia')
ax3.text(0.01, 0.99, 'C',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax3.transAxes,
        fontsize=12,fontweight='extra bold') 

ax4=f1.add_subplot(gs[2,1])
plot_and_fit_ref(ax4,df_global_s,bl2,br2,bb2,bt2,alps_col,False,False)
ax4.set_title('Alps')
ax4.text(0.01, 0.99, 'D',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax4.transAxes,
        fontsize=12,fontweight='extra bold')

ax5=f1.add_subplot(gs[2,2])
plot_and_fit_ref(ax5,df_global_s,bl1,br1,bb1,bt1,gc_col,False,True)
ax5.set_title('Greater Caucasus')
ax5.text(0.01, 0.99, 'E',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax5.transAxes,
        fontsize=12,fontweight='extra bold')

plt.tight_layout()
plt.rcdefaults()

