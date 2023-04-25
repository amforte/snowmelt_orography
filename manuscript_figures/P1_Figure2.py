#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 24 20:35:37 2023

@author: aforte
"""

import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
import matplotlib.gridspec as gridspec
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
from cmcrameri import cm
from matplotlib import colors
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

master_location='/Users/aforte/Documents/Python/snowmelt/'
# master_location='/Volumes/Choruh/Data/snowmelt_project/'
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

idx=(rlf>500) & (df['HCDN-2009']=='yes') & (df['MEAN_Z']>250) & (df['SlicedComp']>0.95) & (perc_base<percb_cutoff)


shape_name=master_location+'gagesII/boundaries_shapefiles_by_aggeco/bas_ref_all_wgs84.shp'
gdf=gpd.read_file(shape_name)

ids=df.loc[idx,'STAID'].to_numpy()
gagesall=gdf['GAGE_ID'].to_numpy().astype(int)
idx2=np.isin(gagesall,ids)
gagesdf=gdf.loc[idx2,:]


raster=master_location+'srtm30_plus/topo30_wrap.tif'
with rasterio.open(raster,'r') as src:
    raster_crs=src.crs
    left1, bottom1, right1, top1 = src.bounds
    dem=src.read()[0,:,:]
    
raster=master_location+'wrr2_raster_outputs/wrr2_R_unfiltered.tif'
with rasterio.open(raster,'r') as src:
    raster_crs=src.crs
    left, bottom, right, top = src.bounds
    runoff=src.read()[0,:,:]  

nodata=runoff[0,0]
runoff[runoff==nodata]=np.nan

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

f1=plt.figure(1,figsize=(6,4)) 
f1.set_dpi(250)
gs=f1.add_gridspec(1,1)

ax1=f1.add_subplot(gs[0],projection=ccrs.PlateCarree())
ax1.set_extent([-180, 180, -90, 90], crs=ccrs.PlateCarree())
gl=ax1.gridlines(draw_labels=True,alpha=0.5,linewidth=0.25)
gl.right_labels=False
gl.bottom_labels=False
ax1.coastlines(resolution='auto',color='k',linewidth=0.25)
ax1.imshow(dem,extent=(left1,right1,bottom1,top1),cmap=cm.grayC,vmin=-4000,vmax=4000,alpha=1)
im1=ax1.imshow(runoff,extent=(left,right,bottom,top),cmap=cm.roma,vmin=0,vmax=10,alpha=0.75)
cbar1=plt.colorbar(im1,ax=ax1)
cbar1.ax.set_ylabel('Mean Runoff [mm/day]')
gagesdf.plot(ax=ax1,color='k')

bl1=-130
br1=-70
bb1=30
bt1=50

# Plot Box
ax1.plot([bl1,bl1,br1,br1,bl1],[bb1,bt1,bt1,bb1,bb1],c='k',linewidth=0.5)
# Set Insent
axins1 = inset_axes(ax1, width="75%", height="75%", loc='lower left',
                    bbox_to_anchor=(0.01, -0.25, 1, 0.5),
                    bbox_transform=ax1.transAxes,
                    axes_class=cartopy.mpl.geoaxes.GeoAxes,
                    axes_kwargs=dict(map_projection=ccrs.PlateCarree()))
axins1.set_extent([bl1,br1,bb1,bt1])
axins1.imshow(dem,extent=(left1,right1,bottom1,top1),cmap=cm.grayC,vmin=-4000,vmax=4000)
axins1.imshow(runoff,extent=(left,right,bottom,top),cmap=cm.roma,vmin=0,vmax=10,alpha=0.75)
axins1.spines['geo'].set(edgecolor='k')
# Draw Connectors
rect, connectors1 = ax1.indicate_inset_zoom(axins1, edgecolor='k', alpha=0.5, transform=ax1.transAxes)
connectors1[0].set_visible(False)
connectors1[1].set_visible(True)
connectors1[3].set_visible(True)
connectors1[2].set_visible(False)
gagesdf.plot(ax=axins1,color='k')

plt.rcdefaults()

f1.savefig('P1_figure2.png',dpi="figure")

