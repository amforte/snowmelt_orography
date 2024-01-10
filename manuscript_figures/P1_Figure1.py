#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 24 17:16:58 2023

@author: aforte
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from cmcrameri import cm
from scipy.stats import linregress
from matplotlib import colors
import geopandas as gpd
import rasterio
import matplotlib.gridspec as gridspec
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.patches as mpatches

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

idx1=(df['Set']=='ref')
idx2=(rlf>500) & (df['HCDN-2009']=='yes') & (df['MEAN_Z']>250) & (df['SlicedComp']>0.95) & (perc_base<percb_cutoff)


shape_name=master_location+'gagesII/boundaries_shapefiles_by_aggeco/bas_ref_all_wgs84.shp'
gdf=gpd.read_file(shape_name)

ref_ids=df.loc[idx1,'STAID'].to_numpy()
hcdn_ids=df.loc[idx2,'STAID'].to_numpy()
gagesall=gdf['GAGE_ID'].to_numpy().astype(int)
gidx1=np.isin(gagesall,ref_ids)
gidx2=np.isin(gagesall,hcdn_ids)
gages_ref_df=gdf.loc[gidx1,:]
gages_hcdn_df=gdf.loc[gidx2,:]


R=np.linspace(0,10,500)
MAR=R*365.25

raster=master_location+'srtm30_plus/wc2.1_30s_elev.tif'
with rasterio.open(raster,'r') as src:
    raster_crs=src.crs
    left1, bottom1, right1, top1 = src.bounds
    dem=src.read()[0,:,:]

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

f1=plt.figure(1,figsize=(6,8))
f1.set_dpi(350)

gs=f1.add_gridspec(2,1)

ax1=f1.add_subplot(gs[0,0],projection=ccrs.PlateCarree())
ax1.set_extent([-128, -65, 20, 60], crs=ccrs.PlateCarree())
gl=ax1.gridlines(draw_labels=True,alpha=0.5,linewidth=0.25)
# gl.right_labels=False
gl.bottom_labels=False
ax1.coastlines(resolution='auto',color='k',linewidth=0.25)
ax1.imshow(dem,extent=(left1,right1,bottom1,top1),cmap=cm.grayC,vmin=0,vmax=4000,alpha=1)
gages_ref_df.plot(ax=ax1,color='steelblue')
gages_hcdn_df.plot(ax=ax1,color='skyblue')
# Dummy artists because matplotlib is really stupid and doesn't create legends automatically
# for polygons, because who would be crazy enough to want a LEGEND ENTRY FOR A POLYGON
ref_patch=mpatches.Patch(color='steelblue',label='Gages-II Reference')
hcdn_patch=mpatches.Patch(color='skyblue',label='Filtered HCDN-2009')
# plt.legend(handles=[ref_patch,hcdn_patch],bbox_to_anchor= (-0.1,-0.5),loc='lower left')
plt.legend(handles=[ref_patch,hcdn_patch],loc='upper right')
ax1.text(0.01, 0.99, 'A',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax1.transAxes,
        fontsize=12,fontweight='extra bold')

ax2=f1.add_subplot(gs[1,0])
sc1=plt.scatter(df.loc[idx1,'FullMeanR'],df.loc[idx1,'FullRC1'],s=15,c=df.loc[idx1,'MAT'],cmap=cm.vik,label='Gages-II Reference',vmin=-10,vmax=25)
plt.scatter(df.loc[idx2,'FullMeanR'],df.loc[idx2,'FullRC1'],s=20,c=df.loc[idx2,'MAT'],cmap=cm.vik,marker='s',edgecolors='k',label='Filtered HCDN-2009',vmin=-10,vmax=25)
plt.plot(R,0.109*MAR**0.21,c='k',label='Rainfall-dominated runoff (Rossi et al., 2016)')
cbar1=plt.colorbar(sc1,ax=ax2)
cbar1.ax.set_ylabel('MAT [$^\circ$C]')
plt.legend(loc='best')
# plt.legend(bbox_to_anchor= (-0.1,-0.5),loc='lower left')
ax2.text(0.01, 0.99, 'B',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax2.transAxes,
        fontsize=12,fontweight='extra bold')

plt.xlim((-0.1,10))
plt.ylim((0,2))
plt.xlabel('Mean Runoff [mm/day]')
plt.ylabel('Shape Parameter')

plt.rcdefaults()
plt.tight_layout()
f1.savefig('P1_figure1.pdf',dpi='figure')
# f1.savefig('P1_figure1.png',dpi='figure')




