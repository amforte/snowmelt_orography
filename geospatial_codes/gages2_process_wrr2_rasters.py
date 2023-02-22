#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 22 15:01:44 2022

@author: aforte
"""

import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
import rasterio.mask
import matplotlib.pyplot as plt
import warnings

warnings.filterwarnings(action='ignore', message='Mean of empty slice')
warnings.filterwarnings(action='ignore', message='Degrees of freedom <= 0 for slice')


def raster_clip_stat_gages(shape_name,raster_name):
    # Creates and iterable object for a shapefile
    gdf=gpd.read_file(shape_name)
    # GAGES shapes have multiple polygons, so must iterate through them    
    num_poly=len(gdf)
    mn_val=np.zeros(num_poly)
    std_val=np.zeros(num_poly)
    for i in range(num_poly):
        # Masks out the area within the polygon defined by the shapefile
        with rasterio.open(raster_name) as src:
            out_im,out_trn=rasterio.mask.mask(src,[gdf.geometry[i]],crop=True)
        # Set no data values to nan and find mean and standard deviation    
        out_array=out_im[0,:,:].astype('float')
        out_array[out_array==src.nodata]=np.nan
        mn_val[i]=np.nanmean(out_array)
        std_val[i]=np.nanstd(out_array)
    return gdf,mn_val,std_val


#################
## Define rasters
master_location='/Volumes/Choruh/Data/snowmelt_project/'
P_raster=master_location+'wrr2_raster_outputs/wrr2_P'
PC_raster=master_location+'wrr2_raster_outputs/wrr2_PC'
R_raster=master_location+'wrr2_raster_outputs/wrr2_R'
RC1_raster=master_location+'wrr2_raster_outputs/wrr2_RC1'
RC5_raster=master_location+'wrr2_raster_outputs/wrr2_RC5'
R2_5_raster=master_location+'wrr2_raster_outputs/wrr2_R2_5'
R5_raster=master_location+'wrr2_raster_outputs/wrr2_R5'
R10_raster=master_location+'wrr2_raster_outputs/wrr2_R10'
PS_raster=master_location+'wrr2_raster_outputs/wrr2_PS'
RS1_raster=master_location+'wrr2_raster_outputs/wrr2_RS1'
RS5_raster=master_location+'wrr2_raster_outputs/wrr2_RS5'
QS_raster=master_location+'wrr2_raster_outputs/wrr2_QS'
QSM_raster=master_location+'wrr2_raster_outputs/wrr2_QSM'
QSB_raster=master_location+'wrr2_raster_outputs/wrr2_QSB'
MAT_raster=master_location+'wrr2_raster_outputs/wrr2_MAT'

#################
## Start processing

## REF
shape_name=master_location+'gagesII/boundaries_shapefiles_by_aggeco/bas_ref_all_wgs84.shp'
[gdf,p_mn,p_std]=raster_clip_stat_gages(shape_name,P_raster+'_upsample.tif')
[gdf,p_c,p_c_std]=raster_clip_stat_gages(shape_name,PC_raster+'_upsample.tif')
[gdf,r_mn,r_std]=raster_clip_stat_gages(shape_name,R_raster+'_upsample.tif')
[gdf,r_c1,r_c1_std]=raster_clip_stat_gages(shape_name,RC1_raster+'_upsample.tif')
[gdf,r_c5,r_c5_std]=raster_clip_stat_gages(shape_name,RC5_raster+'_upsample.tif')
[gdf,r2_5,r2_5_std]=raster_clip_stat_gages(shape_name,R2_5_raster+'_upsample.tif')
[gdf,r5,r5_std]=raster_clip_stat_gages(shape_name,R5_raster+'_upsample.tif')
[gdf,r10,r10_std]=raster_clip_stat_gages(shape_name,R10_raster+'_upsample.tif')
[gdf,r_s1,r_s1_std]=raster_clip_stat_gages(shape_name,RS1_raster+'_upsample.tif')
[gdf,r_s5,r_s5_std]=raster_clip_stat_gages(shape_name,RS5_raster+'_upsample.tif')
[gdf,p_s,p_s_std]=raster_clip_stat_gages(shape_name,PS_raster+'_upsample.tif')
[gdf,qs,qs_std]=raster_clip_stat_gages(shape_name,QS_raster+'_upsample.tif')
[gdf,qsm,qsm_std]=raster_clip_stat_gages(shape_name,QSM_raster+'_upsample.tif')
[gdf,qsb,qsb_std]=raster_clip_stat_gages(shape_name,QSB_raster+'_upsample.tif')
[gdf,mat,mat_std]=raster_clip_stat_gages(shape_name,MAT_raster+'_upsample.tif')
# Construct output dataframe
str_id=[]
for i in range(len(p_mn)):
    str_id.append('ref')
df_ref=pd.DataFrame({'STAID':gdf['GAGE_ID'],
                  'AREA':gdf['AREA'],
                  'P':p_mn,
                  'PC':p_c,
                  'PS':p_s,
                  'R':r_mn,
                  'RC1':r_c1,
                  'RC5':r_c5,
                  'RS1':r_s1,
                  'RS5':r_s5,
                  'Rf2_5':r2_5,
                  'Rf5':r5,
                  'Rf10':r10,
                  'QS':qs,
                  'QSM':qsm,
                  'QSB':qsb,
                  'MAT':mat,
                  'Set':str_id})
print('Completed REF')

## AKHIPR
shape_name=master_location+'gagesII/boundaries_shapefiles_by_aggeco/bas_nonref_AKHIPR_wgs84.shp'
[gdf,p_mn,p_std]=raster_clip_stat_gages(shape_name,P_raster+'_upsample.tif')
[gdf,p_c,p_c_std]=raster_clip_stat_gages(shape_name,PC_raster+'_upsample.tif')
[gdf,r_mn,r_std]=raster_clip_stat_gages(shape_name,R_raster+'_upsample.tif')
[gdf,r_c1,r_c1_std]=raster_clip_stat_gages(shape_name,RC1_raster+'_upsample.tif')
[gdf,r_c5,r_c5_std]=raster_clip_stat_gages(shape_name,RC5_raster+'_upsample.tif')
[gdf,r2_5,r2_5_std]=raster_clip_stat_gages(shape_name,R2_5_raster+'_upsample.tif')
[gdf,r5,r5_std]=raster_clip_stat_gages(shape_name,R5_raster+'_upsample.tif')
[gdf,r10,r10_std]=raster_clip_stat_gages(shape_name,R10_raster+'_upsample.tif')
[gdf,r_s1,r_s1_std]=raster_clip_stat_gages(shape_name,RS1_raster+'_upsample.tif')
[gdf,r_s5,r_s5_std]=raster_clip_stat_gages(shape_name,RS5_raster+'_upsample.tif')
[gdf,p_s,p_s_std]=raster_clip_stat_gages(shape_name,PS_raster+'_upsample.tif')
[gdf,qs,qs_std]=raster_clip_stat_gages(shape_name,QS_raster+'_upsample.tif')
[gdf,qsm,qsm_std]=raster_clip_stat_gages(shape_name,QSM_raster+'_upsample.tif')
[gdf,qsb,qsb_std]=raster_clip_stat_gages(shape_name,QSB_raster+'_upsample.tif')
[gdf,mat,mat_std]=raster_clip_stat_gages(shape_name,MAT_raster+'_upsample.tif')
# Construct output dataframe
str_id=[]
for i in range(len(p_mn)):
    str_id.append('AKHIPR')
df_akhipr=pd.DataFrame({'STAID':gdf['GAGE_ID'],
                  'AREA':gdf['AREA'],
                  'P':p_mn,
                  'PC':p_c,
                  'PS':p_s,
                  'R':r_mn,
                  'RC1':r_c1,
                  'RC5':r_c5,
                  'RS1':r_s1,
                  'RS5':r_s5,
                  'Rf2_5':r2_5,
                  'Rf5':r5,
                  'Rf10':r10,
                  'QS':qs,
                  'QSM':qsm,
                  'QSB':qsb,
                  'MAT':mat,
                  'Set':str_id})
print('Completed AKHIPR')

## CntlPlains
shape_name=master_location+'gagesII/boundaries_shapefiles_by_aggeco/bas_nonref_CntlPlains_wgs84.shp'
[gdf,p_mn,p_std]=raster_clip_stat_gages(shape_name,P_raster+'_upsample.tif')
[gdf,p_c,p_c_std]=raster_clip_stat_gages(shape_name,PC_raster+'_upsample.tif')
[gdf,r_mn,r_std]=raster_clip_stat_gages(shape_name,R_raster+'_upsample.tif')
[gdf,r_c1,r_c1_std]=raster_clip_stat_gages(shape_name,RC1_raster+'_upsample.tif')
[gdf,r_c5,r_c5_std]=raster_clip_stat_gages(shape_name,RC5_raster+'_upsample.tif')
[gdf,r2_5,r2_5_std]=raster_clip_stat_gages(shape_name,R2_5_raster+'_upsample.tif')
[gdf,r5,r5_std]=raster_clip_stat_gages(shape_name,R5_raster+'_upsample.tif')
[gdf,r10,r10_std]=raster_clip_stat_gages(shape_name,R10_raster+'_upsample.tif')
[gdf,r_s1,r_s1_std]=raster_clip_stat_gages(shape_name,RS1_raster+'_upsample.tif')
[gdf,r_s5,r_s5_std]=raster_clip_stat_gages(shape_name,RS5_raster+'_upsample.tif')
[gdf,p_s,p_s_std]=raster_clip_stat_gages(shape_name,PS_raster+'_upsample.tif')
[gdf,qs,qs_std]=raster_clip_stat_gages(shape_name,QS_raster+'_upsample.tif')
[gdf,qsm,qsm_std]=raster_clip_stat_gages(shape_name,QSM_raster+'_upsample.tif')
[gdf,qsb,qsb_std]=raster_clip_stat_gages(shape_name,QSB_raster+'_upsample.tif')
[gdf,mat,mat_std]=raster_clip_stat_gages(shape_name,MAT_raster+'_upsample.tif')
# Construct output dataframe
str_id=[]
for i in range(len(p_mn)):
    str_id.append('CntlPlains')
df_cntlplains=pd.DataFrame({'STAID':gdf['GAGE_ID'],
                  'AREA':gdf['AREA'],
                  'P':p_mn,
                  'PC':p_c,
                  'PS':p_s,
                  'R':r_mn,
                  'RC1':r_c1,
                  'RC5':r_c5,
                  'RS1':r_s1,
                  'RS5':r_s5,
                  'Rf2_5':r2_5,
                  'Rf5':r5,
                  'Rf10':r10,
                  'QS':qs,
                  'QSM':qsm,
                  'QSB':qsb,
                  'MAT':mat,
                  'Set':str_id})
print('Completed CntlPlains')

## EastHghlnds
shape_name=master_location+'gagesII/boundaries_shapefiles_by_aggeco/bas_nonref_EastHghlnds_wgs84.shp'
[gdf,p_mn,p_std]=raster_clip_stat_gages(shape_name,P_raster+'_upsample.tif')
[gdf,p_c,p_c_std]=raster_clip_stat_gages(shape_name,PC_raster+'_upsample.tif')
[gdf,r_mn,r_std]=raster_clip_stat_gages(shape_name,R_raster+'_upsample.tif')
[gdf,r_c1,r_c1_std]=raster_clip_stat_gages(shape_name,RC1_raster+'_upsample.tif')
[gdf,r_c5,r_c5_std]=raster_clip_stat_gages(shape_name,RC5_raster+'_upsample.tif')
[gdf,r2_5,r2_5_std]=raster_clip_stat_gages(shape_name,R2_5_raster+'_upsample.tif')
[gdf,r5,r5_std]=raster_clip_stat_gages(shape_name,R5_raster+'_upsample.tif')
[gdf,r10,r10_std]=raster_clip_stat_gages(shape_name,R10_raster+'_upsample.tif')
[gdf,r_s1,r_s1_std]=raster_clip_stat_gages(shape_name,RS1_raster+'_upsample.tif')
[gdf,r_s5,r_s5_std]=raster_clip_stat_gages(shape_name,RS5_raster+'_upsample.tif')
[gdf,p_s,p_s_std]=raster_clip_stat_gages(shape_name,PS_raster+'_upsample.tif')
[gdf,qs,qs_std]=raster_clip_stat_gages(shape_name,QS_raster+'_upsample.tif')
[gdf,qsm,qsm_std]=raster_clip_stat_gages(shape_name,QSM_raster+'_upsample.tif')
[gdf,qsb,qsb_std]=raster_clip_stat_gages(shape_name,QSB_raster+'_upsample.tif')
[gdf,mat,mat_std]=raster_clip_stat_gages(shape_name,MAT_raster+'_upsample.tif')
# Construct output dataframe
str_id=[]
for i in range(len(p_mn)):
    str_id.append('EastHghlnds')
df_easthghlnds=pd.DataFrame({'STAID':gdf['GAGE_ID'],
                  'AREA':gdf['AREA'],
                  'P':p_mn,
                  'PC':p_c,
                  'PS':p_s,
                  'R':r_mn,
                  'RC1':r_c1,
                  'RC5':r_c5,
                  'RS1':r_s1,
                  'RS5':r_s5,
                  'Rf2_5':r2_5,
                  'Rf5':r5,
                  'Rf10':r10,
                  'QS':qs,
                  'QSM':qsm,
                  'QSB':qsb,
                  'MAT':mat,
                  'Set':str_id})
print('Completed EastHghlands')

## MxWdShld
shape_name=master_location+'gagesII/boundaries_shapefiles_by_aggeco/bas_nonref_MxWdShld_wgs84.shp'
[gdf,p_mn,p_std]=raster_clip_stat_gages(shape_name,P_raster+'_upsample.tif')
[gdf,p_c,p_c_std]=raster_clip_stat_gages(shape_name,PC_raster+'_upsample.tif')
[gdf,r_mn,r_std]=raster_clip_stat_gages(shape_name,R_raster+'_upsample.tif')
[gdf,r_c1,r_c1_std]=raster_clip_stat_gages(shape_name,RC1_raster+'_upsample.tif')
[gdf,r_c5,r_c5_std]=raster_clip_stat_gages(shape_name,RC5_raster+'_upsample.tif')
[gdf,r2_5,r2_5_std]=raster_clip_stat_gages(shape_name,R2_5_raster+'_upsample.tif')
[gdf,r5,r5_std]=raster_clip_stat_gages(shape_name,R5_raster+'_upsample.tif')
[gdf,r10,r10_std]=raster_clip_stat_gages(shape_name,R10_raster+'_upsample.tif')
[gdf,r_s1,r_s1_std]=raster_clip_stat_gages(shape_name,RS1_raster+'_upsample.tif')
[gdf,r_s5,r_s5_std]=raster_clip_stat_gages(shape_name,RS5_raster+'_upsample.tif')
[gdf,p_s,p_s_std]=raster_clip_stat_gages(shape_name,PS_raster+'_upsample.tif')
[gdf,qs,qs_std]=raster_clip_stat_gages(shape_name,QS_raster+'_upsample.tif')
[gdf,qsm,qsm_std]=raster_clip_stat_gages(shape_name,QSM_raster+'_upsample.tif')
[gdf,qsb,qsb_std]=raster_clip_stat_gages(shape_name,QSB_raster+'_upsample.tif')
[gdf,mat,mat_std]=raster_clip_stat_gages(shape_name,MAT_raster+'_upsample.tif')
# Construct output dataframe
str_id=[]
for i in range(len(p_mn)):
    str_id.append('MxWdShld')
df_mxwdshld=pd.DataFrame({'STAID':gdf['GAGE_ID'],
                  'AREA':gdf['AREA'],
                  'P':p_mn,
                  'PC':p_c,
                  'PS':p_s,
                  'R':r_mn,
                  'RC1':r_c1,
                  'RC5':r_c5,
                  'RS1':r_s1,
                  'RS5':r_s5,
                  'Rf2_5':r2_5,
                  'Rf5':r5,
                  'Rf10':r10,
                  'QS':qs,
                  'QSM':qsm,
                  'QSB':qsb,
                  'MAT':mat,
                  'Set':str_id})
print('Completed MxWdShld')

## NorthEast
shape_name=master_location+'gagesII/boundaries_shapefiles_by_aggeco/bas_nonref_NorthEast_wgs84.shp'
[gdf,p_mn,p_std]=raster_clip_stat_gages(shape_name,P_raster+'_upsample.tif')
[gdf,p_c,p_c_std]=raster_clip_stat_gages(shape_name,PC_raster+'_upsample.tif')
[gdf,r_mn,r_std]=raster_clip_stat_gages(shape_name,R_raster+'_upsample.tif')
[gdf,r_c1,r_c1_std]=raster_clip_stat_gages(shape_name,RC1_raster+'_upsample.tif')
[gdf,r_c5,r_c5_std]=raster_clip_stat_gages(shape_name,RC5_raster+'_upsample.tif')
[gdf,r2_5,r2_5_std]=raster_clip_stat_gages(shape_name,R2_5_raster+'_upsample.tif')
[gdf,r5,r5_std]=raster_clip_stat_gages(shape_name,R5_raster+'_upsample.tif')
[gdf,r10,r10_std]=raster_clip_stat_gages(shape_name,R10_raster+'_upsample.tif')
[gdf,r_s1,r_s1_std]=raster_clip_stat_gages(shape_name,RS1_raster+'_upsample.tif')
[gdf,r_s5,r_s5_std]=raster_clip_stat_gages(shape_name,RS5_raster+'_upsample.tif')
[gdf,p_s,p_s_std]=raster_clip_stat_gages(shape_name,PS_raster+'_upsample.tif')
[gdf,qs,qs_std]=raster_clip_stat_gages(shape_name,QS_raster+'_upsample.tif')
[gdf,qsm,qsm_std]=raster_clip_stat_gages(shape_name,QSM_raster+'_upsample.tif')
[gdf,qsb,qsb_std]=raster_clip_stat_gages(shape_name,QSB_raster+'_upsample.tif')
[gdf,mat,mat_std]=raster_clip_stat_gages(shape_name,MAT_raster+'_upsample.tif')
# Construct output dataframe
str_id=[]
for i in range(len(p_mn)):
    str_id.append('NorthEast')
df_northeast=pd.DataFrame({'STAID':gdf['GAGE_ID'],
                  'AREA':gdf['AREA'],
                  'P':p_mn,
                  'PC':p_c,
                  'PS':p_s,
                  'R':r_mn,
                  'RC1':r_c1,
                  'RC5':r_c5,
                  'RS1':r_s1,
                  'RS5':r_s5,
                  'Rf2_5':r2_5,
                  'Rf5':r5,
                  'Rf10':r10,
                  'QS':qs,
                  'QSM':qsm,
                  'QSB':qsb,
                  'MAT':mat,
                  'Set':str_id})
print('Completed NorthEast')

## SECstPlain
shape_name=master_location+'gagesII/boundaries_shapefiles_by_aggeco/bas_nonref_SECstPlain_wgs84.shp'
[gdf,p_mn,p_std]=raster_clip_stat_gages(shape_name,P_raster+'_upsample.tif')
[gdf,p_c,p_c_std]=raster_clip_stat_gages(shape_name,PC_raster+'_upsample.tif')
[gdf,r_mn,r_std]=raster_clip_stat_gages(shape_name,R_raster+'_upsample.tif')
[gdf,r_c1,r_c1_std]=raster_clip_stat_gages(shape_name,RC1_raster+'_upsample.tif')
[gdf,r_c5,r_c5_std]=raster_clip_stat_gages(shape_name,RC5_raster+'_upsample.tif')
[gdf,r2_5,r2_5_std]=raster_clip_stat_gages(shape_name,R2_5_raster+'_upsample.tif')
[gdf,r5,r5_std]=raster_clip_stat_gages(shape_name,R5_raster+'_upsample.tif')
[gdf,r10,r10_std]=raster_clip_stat_gages(shape_name,R10_raster+'_upsample.tif')
[gdf,r_s1,r_s1_std]=raster_clip_stat_gages(shape_name,RS1_raster+'_upsample.tif')
[gdf,r_s5,r_s5_std]=raster_clip_stat_gages(shape_name,RS5_raster+'_upsample.tif')
[gdf,p_s,p_s_std]=raster_clip_stat_gages(shape_name,PS_raster+'_upsample.tif')
[gdf,qs,qs_std]=raster_clip_stat_gages(shape_name,QS_raster+'_upsample.tif')
[gdf,qsm,qsm_std]=raster_clip_stat_gages(shape_name,QSM_raster+'_upsample.tif')
[gdf,qsb,qsb_std]=raster_clip_stat_gages(shape_name,QSB_raster+'_upsample.tif')
[gdf,mat,mat_std]=raster_clip_stat_gages(shape_name,MAT_raster+'_upsample.tif')
# Construct output dataframe
str_id=[]
for i in range(len(p_mn)):
    str_id.append('SECstPlain')
df_secstplain=pd.DataFrame({'STAID':gdf['GAGE_ID'],
                  'AREA':gdf['AREA'],
                  'P':p_mn,
                  'PC':p_c,
                  'PS':p_s,
                  'R':r_mn,
                  'RC1':r_c1,
                  'RC5':r_c5,
                  'RS1':r_s1,
                  'RS5':r_s5,
                  'Rf2_5':r2_5,
                  'Rf5':r5,
                  'Rf10':r10,
                  'QS':qs,
                  'QSM':qsm,
                  'QSB':qsb,
                  'MAT':mat,
                  'Set':str_id})
print('Completed SECstPlain')

## SEPlains
shape_name=master_location+'gagesII/boundaries_shapefiles_by_aggeco/bas_nonref_SEPlains_wgs84.shp'
[gdf,p_mn,p_std]=raster_clip_stat_gages(shape_name,P_raster+'_upsample.tif')
[gdf,p_c,p_c_std]=raster_clip_stat_gages(shape_name,PC_raster+'_upsample.tif')
[gdf,r_mn,r_std]=raster_clip_stat_gages(shape_name,R_raster+'_upsample.tif')
[gdf,r_c1,r_c1_std]=raster_clip_stat_gages(shape_name,RC1_raster+'_upsample.tif')
[gdf,r_c5,r_c5_std]=raster_clip_stat_gages(shape_name,RC5_raster+'_upsample.tif')
[gdf,r2_5,r2_5_std]=raster_clip_stat_gages(shape_name,R2_5_raster+'_upsample.tif')
[gdf,r5,r5_std]=raster_clip_stat_gages(shape_name,R5_raster+'_upsample.tif')
[gdf,r10,r10_std]=raster_clip_stat_gages(shape_name,R10_raster+'_upsample.tif')
[gdf,r_s1,r_s1_std]=raster_clip_stat_gages(shape_name,RS1_raster+'_upsample.tif')
[gdf,r_s5,r_s5_std]=raster_clip_stat_gages(shape_name,RS5_raster+'_upsample.tif')
[gdf,p_s,p_s_std]=raster_clip_stat_gages(shape_name,PS_raster+'_upsample.tif')
[gdf,qs,qs_std]=raster_clip_stat_gages(shape_name,QS_raster+'_upsample.tif')
[gdf,qsm,qsm_std]=raster_clip_stat_gages(shape_name,QSM_raster+'_upsample.tif')
[gdf,qsb,qsb_std]=raster_clip_stat_gages(shape_name,QSB_raster+'_upsample.tif')
[gdf,mat,mat_std]=raster_clip_stat_gages(shape_name,MAT_raster+'_upsample.tif')
# Construct output dataframe
str_id=[]
for i in range(len(p_mn)):
    str_id.append('SEPlains')
df_seplains=pd.DataFrame({'STAID':gdf['GAGE_ID'],
                  'AREA':gdf['AREA'],
                  'P':p_mn,
                  'PC':p_c,
                  'PS':p_s,
                  'R':r_mn,
                  'RC1':r_c1,
                  'RC5':r_c5,
                  'RS1':r_s1,
                  'RS5':r_s5,
                  'Rf2_5':r2_5,
                  'Rf5':r5,
                  'Rf10':r10,
                  'QS':qs,
                  'QSM':qsm,
                  'QSB':qsb,
                  'MAT':mat,
                  'Set':str_id})
print('Completed SEPlains')

## WestMnts
shape_name=master_location+'gagesII/boundaries_shapefiles_by_aggeco/bas_nonref_WestMnts_wgs84.shp'
[gdf,p_mn,p_std]=raster_clip_stat_gages(shape_name,P_raster+'_upsample.tif')
[gdf,p_c,p_c_std]=raster_clip_stat_gages(shape_name,PC_raster+'_upsample.tif')
[gdf,r_mn,r_std]=raster_clip_stat_gages(shape_name,R_raster+'_upsample.tif')
[gdf,r_c1,r_c1_std]=raster_clip_stat_gages(shape_name,RC1_raster+'_upsample.tif')
[gdf,r_c5,r_c5_std]=raster_clip_stat_gages(shape_name,RC5_raster+'_upsample.tif')
[gdf,r2_5,r2_5_std]=raster_clip_stat_gages(shape_name,R2_5_raster+'_upsample.tif')
[gdf,r5,r5_std]=raster_clip_stat_gages(shape_name,R5_raster+'_upsample.tif')
[gdf,r10,r10_std]=raster_clip_stat_gages(shape_name,R10_raster+'_upsample.tif')
[gdf,r_s1,r_s1_std]=raster_clip_stat_gages(shape_name,RS1_raster+'_upsample.tif')
[gdf,r_s5,r_s5_std]=raster_clip_stat_gages(shape_name,RS5_raster+'_upsample.tif')
[gdf,p_s,p_s_std]=raster_clip_stat_gages(shape_name,PS_raster+'_upsample.tif')
[gdf,qs,qs_std]=raster_clip_stat_gages(shape_name,QS_raster+'_upsample.tif')
[gdf,qsm,qsm_std]=raster_clip_stat_gages(shape_name,QSM_raster+'_upsample.tif')
[gdf,qsb,qsb_std]=raster_clip_stat_gages(shape_name,QSB_raster+'_upsample.tif')
[gdf,mat,mat_std]=raster_clip_stat_gages(shape_name,MAT_raster+'_upsample.tif')
# Construct output dataframe
str_id=[]
for i in range(len(p_mn)):
    str_id.append('WestMnts')
df_westmnts=pd.DataFrame({'STAID':gdf['GAGE_ID'],
                  'AREA':gdf['AREA'],
                  'P':p_mn,
                  'PC':p_c,
                  'PS':p_s,
                  'R':r_mn,
                  'RC1':r_c1,
                  'RC5':r_c5,
                  'RS1':r_s1,
                  'RS5':r_s5,
                  'Rf2_5':r2_5,
                  'Rf5':r5,
                  'Rf10':r10,
                  'QS':qs,
                  'QSM':qsm,
                  'QSB':qsb,
                  'MAT':mat,
                  'Set':str_id})
print('Completed WestMnts')

## WestPlains
shape_name=master_location+'gagesII/boundaries_shapefiles_by_aggeco/bas_nonref_WestPlains_wgs84.shp'
[gdf,p_mn,p_std]=raster_clip_stat_gages(shape_name,P_raster+'_upsample.tif')
[gdf,p_c,p_c_std]=raster_clip_stat_gages(shape_name,PC_raster+'_upsample.tif')
[gdf,r_mn,r_std]=raster_clip_stat_gages(shape_name,R_raster+'_upsample.tif')
[gdf,r_c1,r_c1_std]=raster_clip_stat_gages(shape_name,RC1_raster+'_upsample.tif')
[gdf,r_c5,r_c5_std]=raster_clip_stat_gages(shape_name,RC5_raster+'_upsample.tif')
[gdf,r2_5,r2_5_std]=raster_clip_stat_gages(shape_name,R2_5_raster+'_upsample.tif')
[gdf,r5,r5_std]=raster_clip_stat_gages(shape_name,R5_raster+'_upsample.tif')
[gdf,r10,r10_std]=raster_clip_stat_gages(shape_name,R10_raster+'_upsample.tif')
[gdf,r_s1,r_s1_std]=raster_clip_stat_gages(shape_name,RS1_raster+'_upsample.tif')
[gdf,r_s5,r_s5_std]=raster_clip_stat_gages(shape_name,RS5_raster+'_upsample.tif')
[gdf,p_s,p_s_std]=raster_clip_stat_gages(shape_name,PS_raster+'_upsample.tif')
[gdf,qs,qs_std]=raster_clip_stat_gages(shape_name,QS_raster+'_upsample.tif')
[gdf,qsm,qsm_std]=raster_clip_stat_gages(shape_name,QSM_raster+'_upsample.tif')
[gdf,qsb,qsb_std]=raster_clip_stat_gages(shape_name,QSB_raster+'_upsample.tif')
[gdf,mat,mat_std]=raster_clip_stat_gages(shape_name,MAT_raster+'_upsample.tif')
# Construct output dataframe
str_id=[]
for i in range(len(p_mn)):
    str_id.append('WestPlains')
df_westplains=pd.DataFrame({'STAID':gdf['GAGE_ID'],
                  'AREA':gdf['AREA'],
                  'P':p_mn,
                  'PC':p_c,
                  'PS':p_s,
                  'R':r_mn,
                  'RC1':r_c1,
                  'RC5':r_c5,
                  'RS1':r_s1,
                  'RS5':r_s5,
                  'Rf2_5':r2_5,
                  'Rf5':r5,
                  'Rf10':r10,
                  'QS':qs,
                  'QSM':qsm,
                  'QSB':qsb,
                  'MAT':mat,
                  'Set':str_id})
print('Completed WestPlains')

## WestXeric
shape_name=master_location+'gagesII/boundaries_shapefiles_by_aggeco/bas_nonref_WestXeric_wgs84.shp'
[gdf,p_mn,p_std]=raster_clip_stat_gages(shape_name,P_raster+'_upsample.tif')
[gdf,p_c,p_c_std]=raster_clip_stat_gages(shape_name,PC_raster+'_upsample.tif')
[gdf,r_mn,r_std]=raster_clip_stat_gages(shape_name,R_raster+'_upsample.tif')
[gdf,r_c1,r_c1_std]=raster_clip_stat_gages(shape_name,RC1_raster+'_upsample.tif')
[gdf,r_c5,r_c5_std]=raster_clip_stat_gages(shape_name,RC5_raster+'_upsample.tif')
[gdf,r2_5,r2_5_std]=raster_clip_stat_gages(shape_name,R2_5_raster+'_upsample.tif')
[gdf,r5,r5_std]=raster_clip_stat_gages(shape_name,R5_raster+'_upsample.tif')
[gdf,r10,r10_std]=raster_clip_stat_gages(shape_name,R10_raster+'_upsample.tif')
[gdf,r_s1,r_s1_std]=raster_clip_stat_gages(shape_name,RS1_raster+'_upsample.tif')
[gdf,r_s5,r_s5_std]=raster_clip_stat_gages(shape_name,RS5_raster+'_upsample.tif')
[gdf,p_s,p_s_std]=raster_clip_stat_gages(shape_name,PS_raster+'_upsample.tif')
[gdf,qs,qs_std]=raster_clip_stat_gages(shape_name,QS_raster+'_upsample.tif')
[gdf,qsm,qsm_std]=raster_clip_stat_gages(shape_name,QSM_raster+'_upsample.tif')
[gdf,qsb,qsb_std]=raster_clip_stat_gages(shape_name,QSB_raster+'_upsample.tif')
[gdf,mat,mat_std]=raster_clip_stat_gages(shape_name,MAT_raster+'_upsample.tif')
# Construct output dataframe
str_id=[]
for i in range(len(p_mn)):
    str_id.append('WestXeric')
df_westxeric=pd.DataFrame({'STAID':gdf['GAGE_ID'],
                  'AREA':gdf['AREA'],
                  'P':p_mn,
                  'PC':p_c,
                  'PS':p_s,
                  'R':r_mn,
                  'RC1':r_c1,
                  'RC5':r_c5,
                  'RS1':r_s1,
                  'RS5':r_s5,
                  'Rf2_5':r2_5,
                  'Rf5':r5,
                  'Rf10':r10,
                  'QS':qs,
                  'QSM':qsm,
                  'QSB':qsb,
                  'MAT':mat,
                  'Set':str_id})
print('Completed WestXeric')

## Concatenate and save
df=pd.concat([df_ref,df_akhipr,df_cntlplains,df_easthghlnds,df_mxwdshld,
              df_northeast,df_secstplain,df_seplains,df_westmnts,df_westplains,
              df_westxeric],axis=0)
df.to_csv('gages2_wrr2_raster_values.csv',index=False)
