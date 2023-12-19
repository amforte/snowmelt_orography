#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 19 10:53:05 2023

@author: aforte
"""

import pandas as pd
import numpy as np
import glob
import matplotlib.pyplot as plt
import geopandas as gpd
import matplotlib.gridspec as gridspec
import cartopy
import cartopy.crs as ccrs

grdc_gages=pd.read_csv('grdc_real_ts_v2.csv')
grdc_stats=pd.read_csv('grdc_station_stats.csv')
grdc_wrr2=pd.read_csv('grdc_wrr2_raster_outputs.csv')

df=pd.merge(grdc_stats,grdc_gages,on='GRDC_NO',how='inner')
df=pd.merge(df,grdc_wrr2,on='GRDC_NO',how='inner')

grdc_list=pd.read_csv('grdc_station_list.csv')

grdc_loc=pd.DataFrame({'GRDC_NO':grdc_list['grdc_no'],
                     'LAT':grdc_list['lat'],
                     'LON':grdc_list['long']}
                      )
df=pd.merge(df,grdc_loc,on='GRDC_NO',how='inner')

master_location='/Users/aforte/Documents/Python/snowmelt/'

# Drop incomplete
# df=df.drop(df.index[df['sliced_completeness']<0.90])
# df=df.reset_index(drop=True)

# Read in GRDC polygons
poly1=gpd.GeoDataFrame.from_file(master_location+'grdc_ts/chunk1_stationbasins.geojson')
poly2=gpd.GeoDataFrame.from_file(master_location+'grdc_ts/chunk2_stationbasins.geojson')
poly3=gpd.GeoDataFrame.from_file(master_location+'grdc_ts/chunk3_stationbasins.geojson')
poly4=gpd.GeoDataFrame.from_file(master_location+'grdc_ts/chunk4_stationbasins.geojson')
poly5=gpd.GeoDataFrame.from_file(master_location+'grdc_ts/chunk5_stationbasins.geojson')
poly6=gpd.GeoDataFrame.from_file(master_location+'grdc_ts/chunk6_stationbasins.geojson')
poly7=gpd.GeoDataFrame.from_file(master_location+'grdc_ts/chunk7_stationbasins.geojson')
poly8=gpd.GeoDataFrame.from_file(master_location+'grdc_ts/chunk8_stationbasins.geojson')
poly9=gpd.GeoDataFrame.from_file(master_location+'grdc_ts/chunk9_stationbasins.geojson')
poly10=gpd.GeoDataFrame.from_file(master_location+'grdc_ts/chunk10_stationbasins.geojson')
poly11=gpd.GeoDataFrame.from_file(master_location+'grdc_ts/chunk11_stationbasins.geojson')
poly12=gpd.GeoDataFrame.from_file(master_location+'grdc_ts/chunk12_stationbasins.geojson')

# Merge into one
grdc_poly=gpd.GeoDataFrame(pd.concat([poly1,poly2,poly3,poly4,poly5,poly6,poly7,poly8,poly9,poly10,poly11,poly12],ignore_index=True),crs=poly1.crs)

# Read in dam locations and index based on time and size
all_dams=gpd.GeoDataFrame.from_file(master_location+'GRanD_Version_1_3/GRanD_dams_v1_3.shp')
sel_dams=np.logical_and(all_dams.YEAR<2000,all_dams.AREA_POLY>5)
dams=all_dams.loc[sel_dams,:]

# Find GRDC polygons which contain dams
grdc_poly_with_dams=gpd.sjoin(grdc_poly,dams,how='left')
grouped=grdc_poly_with_dams.groupby('index_right')

# Get list of indices of GRDC stations which contain dams
ix=list(grouped.groups.items())
idx=[]
for i in range(len(ix)):
    ixi=ix[i][1]
    if len(ixi)==1:
        idx.append(ixi[0])
    else:
        for j in range(len(ixi)):
            idx.append(ixi[j])
            
idx=np.array(idx).astype(int)
idx=np.unique(idx)

# Extract GRDC ID numbers which contain dams
grdc_no_with_dams=grdc_poly.grdc_no[idx].to_numpy().astype(int)

# Find these within grdc_database
dam_idx=df['GRDC_NO'].isin(grdc_no_with_dams)



percb_cutoff=0.25
perc_base=df['QSB_MN_UPRASTER']/df['R_MN_UPRASTER']
rlf=df['MAX_Z']-df['MIN_Z']

idx=(rlf>500) & (df['MN_Z']>250) & (df['sliced_completeness']>0.95) & (perc_base<percb_cutoff) & (~dam_idx)
idx_dams=(rlf>500) & (df['MN_Z']>250) & (df['sliced_completeness']>0.95) & (perc_base<percb_cutoff) & (dam_idx)

dfm=df.loc[idx,:]
dfm_dam=df.loc[idx_dams,:]

f1=plt.figure(figsize=(8,8))
f1.set_dpi(250)

gs=f1.add_gridspec(2,2)

ax1=f1.add_subplot(gs[0,0])
plt.scatter(dfm['R_MN_UPRASTER'],dfm['sliced_reported_runoff'],c='k',s=10,label='Dam Free - GRanD')
plt.scatter(dfm_dam['R_MN_UPRASTER'],dfm_dam['sliced_reported_runoff'],c='r',s=10,label='Contains Dams - GRanD')
x=np.linspace(0,10)
plt.plot(x,x,c='k',linestyle='--',label='1:1 Reference Line')
ax1.set_xlabel('WaterGAP3 Mean Runoff [mm/day]')
ax1.set_ylabel('GRDC Mean Runoff [mm/day]')
ax1.legend(loc='best')

ax2=f1.add_subplot(gs[0,1])
plt.scatter(dfm['RC1_MN_UPRASTER'],dfm['sliced_reported_c1'],c='k',s=10,label='No Dams From GRaND')
plt.scatter(dfm_dam['RC1_MN_UPRASTER'],dfm_dam['sliced_reported_c1'],c='r',s=10,label='Dams From GRaND')
x=np.linspace(0,2)
plt.plot(x,x,c='k',linestyle='--',label='1:1 Reference Line')
ax2.set_xlabel('WaterGAP3 Shape ')
ax2.set_ylabel('GRDC Shape')

ax3=f1.add_subplot(gs[1,:],projection=ccrs.PlateCarree())
ax3.set_extent([-180, 180, -90, 90], crs=ccrs.PlateCarree())
gl=ax3.gridlines(draw_labels=True,alpha=0.5,linewidth=0.25)
ax3.coastlines(resolution='auto',color='k',linewidth=0.25)
ax3.scatter(dfm['LON'],dfm['LAT'],c='k',s=5)
ax3.scatter(dfm_dam['LON'],dfm_dam['LAT'],c='r',s=5)

plt.tight_layout()
