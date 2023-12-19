#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  1 08:04:13 2022

@author: aforte
"""

import pandas as pd
import numpy as np
import geopandas as gpd
import fiona
import rasterio
import rasterio.mask
import glob
import warnings
from scipy.stats import linregress


# Suppress warning about centroid being innacurate
warnings.filterwarnings('ignore',message='Geometry is in a geographic CRS.')
# Suppres warnings related to nan sliced
warnings.filterwarnings('ignore',message='Mean of empty slice')
warnings.filterwarnings('ignore',message='All-NaN slice encountered')
warnings.filterwarnings('ignore',message='Degrees of freedom <= 0 for slice.')

def find_stations(record_length,area_cutoff,master_location):
    # Load in complete list of stations
    all_stat=pd.read_csv('grdc_station_list.csv')
    # Generate list for those stations with daily data and below area cutoff
    rec_len=all_stat['d_yrs'].to_numpy()
    idx1_1=~np.isnan(rec_len)
    idx1_2=rec_len>=record_length
    idx1=np.logical_and(idx1_1,idx1_2)
    area=all_stat['area'].to_numpy()
    idx2=area<=area_cutoff
    idx=np.logical_and(idx1,idx2)
    # Create list of station numbers
    stat_w_daily=all_stat['grdc_no'].to_numpy()[idx]
    # Extract length of record details
    d_start=all_stat['d_start'].to_numpy()[idx]
    d_end=all_stat['d_end'].to_numpy()[idx]
    d_len=rec_len[idx]
    # Generate list of stations with shapefiles
    shp_list=glob.glob(master_location+'grdc_basin_shp/*.shp')
    stat_w_shape=np.zeros(len(shp_list))
    for i in range(len(shp_list)):
        s=shp_list[i]
        num=''
        for ch in s:
            if ch.isdigit():
                num += ch
        stat_w_shape[i]=int(num)
    
    # Intersect arrays of station numbers
    [stat_w_both,in1,in2]=np.intersect1d(stat_w_daily,stat_w_shape.astype(int),
                                         return_indices=True)
    # Store in a pandas dataframe
    # Include zero values to be filled in later
    dummy=np.zeros(len(stat_w_both))
    d={'GRDC_NO':stat_w_both,
       'R_AREA':area[idx][in1],
       'H_AREA':dummy,
       'START':d_start[in1].astype(int),
       'END':d_end[in1].astype(int),
       'LENGTH':d_len[in1].astype(int),
       'MIN_Z':dummy,
       'MN_Z':dummy,
       'MAX_Z':dummy,
       'E_NORM_CRLF':dummy,
       'ERR_HI':dummy,
       'BW':dummy,
       'BE':dummy,
       'BN':dummy,
       'BS':dummy,
       'CX':dummy,
       'CY':dummy,
       'R_Q':dummy,
       'H_Q':dummy,
       'R_R':dummy,
       'H_R':dummy,
       }
    df=pd.DataFrame(data=d)
    return df,len(stat_w_both)

def raster_clip_stat(shape_name,raster_name):
    # Creates and iterable object for a shapefile
    with fiona.open(shape_name) as shapefile:
        shapes=[feature['geometry'] for feature in shapefile]
    # Masks out the area within the polygon defined by the shapefile
    with rasterio.open(raster_name) as src:
        out_im,out_trn=rasterio.mask.mask(src,shapes,crop=True)
    # Set no data values to nan and find mean and standard deviation    
    out_array=out_im[0,:,:].astype('float')
    out_array[out_array==src.nodata]=np.nan
    mn_val=np.nanmean(out_array)
    std_val=np.nanstd(out_array)
    max_val=np.nanmax(out_array)
    min_val=np.nanmin(out_array)
    return mn_val, std_val, max_val, min_val

def relate_rasters(rx,ry):
    # Assumes two raster arrays of equal dimensions and grid sizes
    rxf=rx.flatten()
    ryf=ry.flatten()
    idx=np.logical_and(~np.isnan(rxf),~np.isnan(ryf))
    res=linregress(rxf[idx],ryf[idx])
    slp=res.slope
    inter=res.intercept
    rsq=res.rvalue**2
    return slp,inter,rsq

############# BEGIN COMPILATION #####################

# Generate list of stations that meet criteria
area_cutoff=5000
master_location='/Users/aforte/Documents/Python/snowmelt/'
# master_location='/Volumes/Choruh/Data/snowmelt_project/'
[station_list,num_stations]=find_stations(10,area_cutoff,master_location)


dem_raster=master_location+'hyd_glo_dem_15s/hyd_glo_dem_15s.tif'
rlf_raster=master_location+'hyd_glo_dem_15s/hsheds_2500rlf_wgs84_rsmp.tif'

#### DEBUG
# num_stations=100

cnt=0

for i in range(num_stations):
    if np.mod(i,10)==0:
        print('Process is '+str(np.round((i/num_stations)*100,1))+'% complete')
        
    shp_name=master_location+'grdc_basin_shp/grdc_basins_smoothed_md_no_'+str(station_list['GRDC_NO'][i])+'.shp'
    df=gpd.read_file(shp_name)
    
    # There are isolated shapefiles that are missing parameters and so will throw
    # errors. The try except catches these and flags them as problematic
    try:
        area=df['AREA'].to_numpy()
        continue_flag=True
    except:
        station_list.loc[i,'R_AREA']=np.nan
        continue_flag=False
    
    if continue_flag:
        areah=df['AREA_HYS'].to_numpy() # Area of hydrosheds polygon
        # Populate the area
        station_list.loc[i,'R_AREA']=area[0]
        station_list.loc[i,'H_AREA']=areah[0]
        # Some of the basins have areas above the cutoff but did not have
        # reported areas, so this checks again and only proceeds if the hydrosheds
        # area is below the cutoff        
        if areah[0]<=area_cutoff:
            # Get bounds
            bnds=df.bounds # as pandas dataframe
            station_list.loc[i,'BW']=bnds.loc[0,'minx']
            station_list.loc[i,'BE']=bnds.loc[0,'maxx']
            station_list.loc[i,'BN']=bnds.loc[0,'maxy']
            station_list.loc[i,'BS']=bnds.loc[0,'miny']
        
            # Get centroid as point GeoSeries - it will warn you because this in lat-long
            center=df.centroid
            cx=center.x[0]
            cy=center.y[0]
            station_list.loc[i,'CX']=cx
            station_list.loc[i,'CY']=cy
        
            grdc_dis=df['LTA_DISCHA'].to_numpy() # Average discharge m^3/sec
            hydr_dis=df['DISC_HYS'].to_numpy() # Average modeled discharge m^3/sec
            

            station_list.loc[i,'R_Q']=grdc_dis[0]
            station_list.loc[i,'H_Q']=hydr_dis[0]
            if area<0:
                station_list.loc[i,'R_R']=(grdc_dis/(areah*1000*1000))*(100*10*60*60*24)            
            else:
                station_list.loc[i,'R_R']=(grdc_dis/(area*1000*1000))*(100*10*60*60*24)
            station_list.loc[i,'H_R']=(hydr_dis/(areah*1000*1000))*(100*10*60*60*24)
            
            # Begin processing rasters
            
            # Do DEM Stuff
            [z_mn,z_std,z_max,z_min]=raster_clip_stat(shp_name,dem_raster)
            station_list.loc[i,'MIN_Z']=z_min
            station_list.loc[i,'MN_Z']=z_mn
            station_list.loc[i,'MAX_Z']=z_max
            
            [r_mn,r_std,r_max,r_min]=raster_clip_stat(shp_name,rlf_raster)
            station_list.loc[i,'MIN_R']=r_min
            station_list.loc[i,'MN_R']=r_mn
            station_list.loc[i,'MAX_R']=r_max          
            

        else:
            station_list.loc[i,'R_Q']=np.nan
            station_list.loc[i,'H_Q']=np.nan
            station_list.loc[i,'R_R']=np.nan
            station_list.loc[i,'H_R']=np.nan
            station_list.loc[i,'BW']=np.nan
            station_list.loc[i,'BE']=np.nan
            station_list.loc[i,'BN']=np.nan
            station_list.loc[i,'BS']=np.nan         
            station_list.loc[i,'CX']=np.nan
            station_list.loc[i,'CY']=np.nan
            station_list.loc[i,'MIN_Z']=np.nan
            station_list.loc[i,'MN_Z']=np.nan
            station_list.loc[i,'MAX_Z']=np.nan
            station_list.loc[i,'MIN_R']=np.nan
            station_list.loc[i,'MN_R']=np.nan
            station_list.loc[i,'MAX_R']=np.nan            
    else:
        station_list.loc[i,'R_AREA']=np.nan
        station_list.loc[i,'H_AREA']=np.nan
        station_list.loc[i,'R_Q']=np.nan
        station_list.loc[i,'H_Q']=np.nan
        station_list.loc[i,'R_R']=np.nan
        station_list.loc[i,'H_R']=np.nan
        station_list.loc[i,'BW']=np.nan
        station_list.loc[i,'BE']=np.nan
        station_list.loc[i,'BN']=np.nan
        station_list.loc[i,'BS']=np.nan         
        station_list.loc[i,'CX']=np.nan
        station_list.loc[i,'CY']=np.nan
    cnt+=1
# Write out
station_list.to_csv('grdc_station_stats.csv')