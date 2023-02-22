#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 11 14:46:17 2022

@author: aforte
"""

import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
import rasterio.mask
import matplotlib.pyplot as plt
import warnings

warnings.filterwarnings('ignore',message='Mean of empty slice')
warnings.filterwarnings('ignore',message='All-NaN slice encountered')
warnings.filterwarnings('ignore',message='Degrees of freedom <= 0 for slice.')

def raster_clip_stat_gages(shape_name,raster_name):
    # Creates and iterable object for a shapefile
    gdf=gpd.read_file(shape_name)
    # GAGES shapes have multiple polygons, so must iterate through them    
    num_poly=len(gdf)
    mn_val=np.zeros(num_poly)
    max_val=np.zeros(num_poly)
    min_val=np.zeros(num_poly)
    for i in range(num_poly):
        # Masks out the area within the polygon defined by the shapefile
        with rasterio.open(raster_name) as src:
            out_im,out_trn=rasterio.mask.mask(src,[gdf.geometry[i]],crop=True)
        # Set no data values to nan and find mean and standard deviation    
        out_array=out_im[0,:,:].astype('float')
        out_array[out_array==src.nodata]=np.nan
        mn_val[i]=np.nanmean(out_array)
        max_val[i]=np.nanmax(out_array)
        min_val[i]=np.nanmin(out_array)
    return gdf,mn_val,max_val,min_val


#################
## Define rasters
master_location='/Volumes/Choruh/Data/snowmelt_project/'
dem_raster=master_location+'srtm30_plus/topo30_wrap.tif'
et_an_raster=master_location+'Global-AI_ET0_v3_annual/et0_v3_yr.tif'

#################
## Start processing

## REF
shape_name=master_location+'gagesII/boundaries_shapefiles_by_aggeco/bas_ref_all_wgs84.shp'
[gdf,z_mn,z_max,z_min]=raster_clip_stat_gages(shape_name,dem_raster)
[gdf,et_mn,et_max,et_min]=raster_clip_stat_gages(shape_name,et_an_raster)

# Construct output dataframe
str_id=[]
for i in range(len(z_mn)):
    str_id.append('ref')
df_ref=pd.DataFrame({'STAID':gdf['GAGE_ID'],
                  'MEAN_Z':z_mn,
                  'MAX_Z':z_max,
                  'MIN_Z':z_min,
                  'MEAN_ET':et_mn/365.25,
                  'MAX_ET':et_max/365.25,
                  'MIN_ET':et_min/365.25})
print('Completed REF')

## AKHIPR
shape_name=master_location+'gagesII/boundaries_shapefiles_by_aggeco/bas_nonref_AKHIPR_wgs84.shp'
[gdf,z_mn,z_max,z_min]=raster_clip_stat_gages(shape_name,dem_raster)
[gdf,et_mn,et_max,et_min]=raster_clip_stat_gages(shape_name,et_an_raster)

# Construct output dataframe
str_id=[]
for i in range(len(z_mn)):
    str_id.append('AKHIPR')
df_akhipr=pd.DataFrame({'STAID':gdf['GAGE_ID'],
                  'MEAN_Z':z_mn,
                  'MAX_Z':z_max,
                  'MIN_Z':z_min,
                  'MEAN_ET':et_mn/365.25,
                  'MAX_ET':et_max/365.25,
                  'MIN_ET':et_min/365.25})
print('Completed AKHIPR')

## CntlPlains
shape_name=master_location+'gagesII/boundaries_shapefiles_by_aggeco/bas_nonref_CntlPlains_wgs84.shp'
[gdf,z_mn,z_max,z_min]=raster_clip_stat_gages(shape_name,dem_raster)
[gdf,et_mn,et_max,et_min]=raster_clip_stat_gages(shape_name,et_an_raster)

# Construct output dataframe
str_id=[]
for i in range(len(z_mn)):
    str_id.append('CntlPlains')
df_cntlplains=pd.DataFrame({'STAID':gdf['GAGE_ID'],
                  'MEAN_Z':z_mn,
                  'MAX_Z':z_max,
                  'MIN_Z':z_min,
                  'MEAN_ET':et_mn/365.25,
                  'MAX_ET':et_max/365.25,
                  'MIN_ET':et_min/365.25})
print('Completed CntlPlains')

## EastHghlnds
shape_name=master_location+'gagesII/boundaries_shapefiles_by_aggeco/bas_nonref_EastHghlnds_wgs84.shp'
[gdf,z_mn,z_max,z_min]=raster_clip_stat_gages(shape_name,dem_raster)
[gdf,et_mn,et_max,et_min]=raster_clip_stat_gages(shape_name,et_an_raster)

# Construct output dataframe
str_id=[]
for i in range(len(z_mn)):
    str_id.append('EastHghlnds')
df_easthghlnds=pd.DataFrame({'STAID':gdf['GAGE_ID'],
                  'MEAN_Z':z_mn,
                  'MAX_Z':z_max,
                  'MIN_Z':z_min,
                  'MEAN_ET':et_mn/365.25,
                  'MAX_ET':et_max/365.25,
                  'MIN_ET':et_min/365.25})
print('Completed EastHghlands')

## MxWdShld
shape_name=master_location+'gagesII/boundaries_shapefiles_by_aggeco/bas_nonref_MxWdShld_wgs84.shp'
[gdf,z_mn,z_max,z_min]=raster_clip_stat_gages(shape_name,dem_raster)
[gdf,et_mn,et_max,et_min]=raster_clip_stat_gages(shape_name,et_an_raster)

# Construct output dataframe
str_id=[]
for i in range(len(z_mn)):
    str_id.append('MxWdShld')
df_mxwdshld=pd.DataFrame({'STAID':gdf['GAGE_ID'],
                  'MEAN_Z':z_mn,
                  'MAX_Z':z_max,
                  'MIN_Z':z_min,
                  'MEAN_ET':et_mn/365.25,
                  'MAX_ET':et_max/365.25,
                  'MIN_ET':et_min/365.25})
print('Completed MxWdShld')

## NorthEast
shape_name=master_location+'gagesII/boundaries_shapefiles_by_aggeco/bas_nonref_NorthEast_wgs84.shp'
[gdf,z_mn,z_max,z_min]=raster_clip_stat_gages(shape_name,dem_raster)
[gdf,et_mn,et_max,et_min]=raster_clip_stat_gages(shape_name,et_an_raster)

# Construct output dataframe
str_id=[]
for i in range(len(z_mn)):
    str_id.append('NorthEast')
df_northeast=pd.DataFrame({'STAID':gdf['GAGE_ID'],
                  'MEAN_Z':z_mn,
                  'MAX_Z':z_max,
                  'MIN_Z':z_min,
                  'MEAN_ET':et_mn/365.25,
                  'MAX_ET':et_max/365.25,
                  'MIN_ET':et_min/365.25})
print('Completed NorthEast')

## SECstPlain
shape_name=master_location+'gagesII/boundaries_shapefiles_by_aggeco/bas_nonref_SECstPlain_wgs84.shp'
[gdf,z_mn,z_max,z_min]=raster_clip_stat_gages(shape_name,dem_raster)
[gdf,et_mn,et_max,et_min]=raster_clip_stat_gages(shape_name,et_an_raster)

# Construct output dataframe
str_id=[]
for i in range(len(z_mn)):
    str_id.append('SECstPlain')
df_secstplain=pd.DataFrame({'STAID':gdf['GAGE_ID'],
                  'MEAN_Z':z_mn,
                  'MAX_Z':z_max,
                  'MIN_Z':z_min,
                  'MEAN_ET':et_mn/365.25,
                  'MAX_ET':et_max/365.25,
                  'MIN_ET':et_min/365.25})
print('Completed SECstPlain')

## SEPlains
shape_name=master_location+'gagesII/boundaries_shapefiles_by_aggeco/bas_nonref_SEPlains_wgs84.shp'
[gdf,z_mn,z_max,z_min]=raster_clip_stat_gages(shape_name,dem_raster)
[gdf,et_mn,et_max,et_min]=raster_clip_stat_gages(shape_name,et_an_raster)

# Construct output dataframe
str_id=[]
for i in range(len(z_mn)):
    str_id.append('SEPlains')
df_seplains=pd.DataFrame({'STAID':gdf['GAGE_ID'],
                  'MEAN_Z':z_mn,
                  'MAX_Z':z_max,
                  'MIN_Z':z_min,
                  'MEAN_ET':et_mn/365.25,
                  'MAX_ET':et_max/365.25,
                  'MIN_ET':et_min/365.25})
print('Completed SEPlains')

## WestMnts
shape_name=master_location+'gagesII/boundaries_shapefiles_by_aggeco/bas_nonref_WestMnts_wgs84.shp'
[gdf,z_mn,z_max,z_min]=raster_clip_stat_gages(shape_name,dem_raster)
[gdf,et_mn,et_max,et_min]=raster_clip_stat_gages(shape_name,et_an_raster)

# Construct output dataframe
str_id=[]
for i in range(len(z_mn)):
    str_id.append('WestMnts')
df_westmnts=pd.DataFrame({'STAID':gdf['GAGE_ID'],
                  'MEAN_Z':z_mn,
                  'MAX_Z':z_max,
                  'MIN_Z':z_min,
                  'MEAN_ET':et_mn/365.25,
                  'MAX_ET':et_max/365.25,
                  'MIN_ET':et_min/365.25})
print('Completed WestMnts')

## WestPlains
shape_name=master_location+'gagesII/boundaries_shapefiles_by_aggeco/bas_nonref_WestPlains_wgs84.shp'
[gdf,z_mn,z_max,z_min]=raster_clip_stat_gages(shape_name,dem_raster)
[gdf,et_mn,et_max,et_min]=raster_clip_stat_gages(shape_name,et_an_raster)

# Construct output dataframe
str_id=[]
for i in range(len(z_mn)):
    str_id.append('WestPlains')
df_westplains=pd.DataFrame({'STAID':gdf['GAGE_ID'],
                  'MEAN_Z':z_mn,
                  'MAX_Z':z_max,
                  'MIN_Z':z_min,
                  'MEAN_ET':et_mn/365.25,
                  'MAX_ET':et_max/365.25,
                  'MIN_ET':et_min/365.25})
print('Completed WestPlains')

## WestXeric
shape_name=master_location+'gagesII/boundaries_shapefiles_by_aggeco/bas_nonref_WestXeric_wgs84.shp'
[gdf,z_mn,z_max,z_min]=raster_clip_stat_gages(shape_name,dem_raster)
[gdf,et_mn,et_max,et_min]=raster_clip_stat_gages(shape_name,et_an_raster)

# Construct output dataframe
str_id=[]
for i in range(len(z_mn)):
    str_id.append('WestXeric')
df_westxeric=pd.DataFrame({'STAID':gdf['GAGE_ID'],
                  'MEAN_Z':z_mn,
                  'MAX_Z':z_max,
                  'MIN_Z':z_min,
                  'MEAN_ET':et_mn/365.25,
                  'MAX_ET':et_max/365.25,
                  'MIN_ET':et_min/365.25})
print('Completed WestXeric')

## Concatenate
df=pd.concat([df_ref,df_akhipr,df_cntlplains,df_easthghlnds,df_mxwdshld,
              df_northeast,df_secstplain,df_seplains,df_westmnts,df_westplains,
              df_westxeric],axis=0)

df.to_csv('gages2_station_stats.csv')