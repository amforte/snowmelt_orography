#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  9 20:24:38 2022

@author: aforte
"""

import numpy as np
import pandas as pd
import rasterio
import fiona
import rasterio.mask


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
    return mn_val,std_val


master_location='/Users/aforte/Documents/Python/snowmelt/'

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


# Read GRDC stations stats
df=pd.read_csv('grdc_station_stats.csv')
df=df.drop(df.index[np.isnan(df['MIN_Z'])])
df=df.reset_index(drop=True)

# Set arrays
grdc_id=np.zeros(len(df)).astype('int')

mn_p=np.zeros(len(df))
mn_p_up=np.zeros(len(df))
std_p=np.zeros(len(df))
std_p_up=np.zeros(len(df))

mn_pc=np.zeros(len(df))
mn_pc_up=np.zeros(len(df))
std_pc=np.zeros(len(df))
std_pc_up=np.zeros(len(df))

mn_r=np.zeros(len(df))
mn_r_up=np.zeros(len(df))
std_r=np.zeros(len(df))
std_r_up=np.zeros(len(df))

mn_rc1=np.zeros(len(df))
mn_rc1_up=np.zeros(len(df))
std_rc1=np.zeros(len(df))
std_rc1_up=np.zeros(len(df))

mn_rc5=np.zeros(len(df))
mn_rc5_up=np.zeros(len(df))
std_rc5=np.zeros(len(df))
std_rc5_up=np.zeros(len(df))

mn_r2_5=np.zeros(len(df))
mn_r2_5_up=np.zeros(len(df))
std_r2_5=np.zeros(len(df))
std_r2_5_up=np.zeros(len(df))

mn_r5=np.zeros(len(df))
mn_r5_up=np.zeros(len(df))
std_r5=np.zeros(len(df))
std_r5_up=np.zeros(len(df))

mn_r10=np.zeros(len(df))
mn_r10_up=np.zeros(len(df))
std_r10=np.zeros(len(df))
std_r10_up=np.zeros(len(df))


mn_ps=np.zeros(len(df))
mn_ps_up=np.zeros(len(df))
std_ps=np.zeros(len(df))
std_ps_up=np.zeros(len(df))

mn_rs1=np.zeros(len(df))
mn_rs1_up=np.zeros(len(df))
std_rs1=np.zeros(len(df))
std_rs1_up=np.zeros(len(df))

mn_rs5=np.zeros(len(df))
mn_rs5_up=np.zeros(len(df))
std_rs5=np.zeros(len(df))
std_rs5_up=np.zeros(len(df))

mn_qs=np.zeros(len(df))
mn_qs_up=np.zeros(len(df))
std_qs=np.zeros(len(df))
std_qs_up=np.zeros(len(df))

mn_qsm=np.zeros(len(df))
mn_qsm_up=np.zeros(len(df))
std_qsm=np.zeros(len(df))
std_qsm_up=np.zeros(len(df))

mn_qsb=np.zeros(len(df))
mn_qsb_up=np.zeros(len(df))
std_qsb=np.zeros(len(df))
std_qsb_up=np.zeros(len(df))

mn_mat=np.zeros(len(df))
mn_mat_up=np.zeros(len(df))
std_mat=np.zeros(len(df))
std_mat_up=np.zeros(len(df))


for i in range(len(df)):
    if np.mod(i,100)==0:
        print(i)
    shp_name=master_location+'grdc_basin_shp/grdc_basins_smoothed_md_no_'+str(df.loc[i,'GRDC_NO'])+'.shp'
    grdc_id[i]=df.loc[i,'GRDC_NO']
    
    [mn_p[i],std_p[i]]=raster_clip_stat(shp_name,P_raster+'.tif')
    [mn_p_up[i],std_p_up[i]]=raster_clip_stat(shp_name,P_raster+'_upsample.tif')
    
    [mn_pc[i],std_pc[i]]=raster_clip_stat(shp_name,PC_raster+'.tif')
    [mn_pc_up[i],std_pc_up[i]]=raster_clip_stat(shp_name,PC_raster+'_upsample.tif')

    [mn_ps[i],std_ps[i]]=raster_clip_stat(shp_name,PS_raster+'.tif')
    [mn_ps_up[i],std_ps_up[i]]=raster_clip_stat(shp_name,PS_raster+'_upsample.tif')
    
    [mn_r[i],std_r[i]]=raster_clip_stat(shp_name,R_raster+'.tif')
    [mn_r_up[i],std_r_up[i]]=raster_clip_stat(shp_name,R_raster+'_upsample.tif')
    
    [mn_rc1[i],std_rc1[i]]=raster_clip_stat(shp_name,RC1_raster+'.tif')
    [mn_rc1_up[i],std_rc1_up[i]]=raster_clip_stat(shp_name,RC1_raster+'_upsample.tif')

    [mn_rs1[i],std_rs1[i]]=raster_clip_stat(shp_name,RS1_raster+'.tif')
    [mn_rs1_up[i],std_rs1_up[i]]=raster_clip_stat(shp_name,RS1_raster+'_upsample.tif')
    
    [mn_rc5[i],std_rc5[i]]=raster_clip_stat(shp_name,RC5_raster+'.tif')
    [mn_rc5_up[i],std_rc5_up[i]]=raster_clip_stat(shp_name,RC5_raster+'_upsample.tif')
    
    [mn_rs5[i],std_rs5[i]]=raster_clip_stat(shp_name,RS5_raster+'.tif')
    [mn_rs5_up[i],std_rs5_up[i]]=raster_clip_stat(shp_name,RS5_raster+'_upsample.tif')
    
    [mn_r2_5[i],std_r2_5[i]]=raster_clip_stat(shp_name,R2_5_raster+'.tif')
    [mn_r2_5_up[i],std_r2_5_up[i]]=raster_clip_stat(shp_name,R2_5_raster+'_upsample.tif')
    
    [mn_r5[i],std_r5[i]]=raster_clip_stat(shp_name,R5_raster+'.tif')
    [mn_r5_up[i],std_r5_up[i]]=raster_clip_stat(shp_name,R5_raster+'_upsample.tif')
    
    [mn_r10[i],std_r10[i]]=raster_clip_stat(shp_name,R10_raster+'.tif')
    [mn_r10_up[i],std_r10_up[i]]=raster_clip_stat(shp_name,R10_raster+'_upsample.tif')

    [mn_qs[i],std_qs[i]]=raster_clip_stat(shp_name,QS_raster+'.tif')
    [mn_qs_up[i],std_qs_up[i]]=raster_clip_stat(shp_name,QS_raster+'_upsample.tif')

    [mn_qsm[i],std_qsm[i]]=raster_clip_stat(shp_name,QSM_raster+'.tif')
    [mn_qsm_up[i],std_qsm_up[i]]=raster_clip_stat(shp_name,QSM_raster+'_upsample.tif')
    
    [mn_qsb[i],std_qsb[i]]=raster_clip_stat(shp_name,QSB_raster+'.tif')
    [mn_qsb_up[i],std_qsb_up[i]]=raster_clip_stat(shp_name,QSB_raster+'_upsample.tif')    

    [mn_mat[i],std_mat[i]]=raster_clip_stat(shp_name,MAT_raster+'.tif')
    [mn_mat_up[i],std_mat_up[i]]=raster_clip_stat(shp_name,MAT_raster+'_upsample.tif')


df_out=pd.DataFrame({'GRDC_NO':grdc_id,
                     'P_MN_RASTER':mn_p,
                     'P_MN_UPRASTER':mn_p_up,
                     'P_STD_RASTER':std_p,
                     'P_STD_UPRASTER':std_p_up,
                     'PC_MN_RASTER':mn_pc,
                     'PC_MN_UPRASTER':mn_pc_up,
                     'PC_STD_RASTER':std_pc,
                     'PC_STD_UPRASTER':std_pc_up,
                     'PS_MN_RASTER':mn_ps,
                     'PS_MN_UPRASTER':mn_ps_up,
                     'PS_STD_RASTER':std_ps,
                     'PS_STD_UPRASTER':std_ps_up,
                     'R_MN_RASTER':mn_r,
                     'R_MN_UPRASTER':mn_r_up,
                     'R_STD_RASTER':std_r,
                     'R_STD_UPRASTER':std_r_up,
                     'RC1_MN_RASTER':mn_rc1,
                     'RC1_MN_UPRASTER':mn_rc1_up,
                     'RC1_STD_RASTER':std_rc1,
                     'RC1_STD_UPRASTER':std_rc1_up,
                     'RC5_MN_RASTER':mn_rc5,
                     'RC5_MN_UPRASTER':mn_rc5_up,
                     'RC5_STD_RASTER':std_rc5,
                     'RC5_STD_UPRASTER':std_rc5_up,
                     'RS1_MN_RASTER':mn_rs1,
                     'RS1_MN_UPRASTER':mn_rs1_up,
                     'RS1_STD_RASTER':std_rs1,
                     'RS1_STD_UPRASTER':std_rs1_up,
                     'RS5_MN_RASTER':mn_rs5,
                     'RS5_MN_UPRASTER':mn_rs5_up,
                     'RS5_STD_RASTER':std_rs5,
                     'RS5_STD_UPRASTER':std_rs5_up,
                     'R2_5_MN_RASTER':mn_r2_5,
                     'R2_5_MN_UPRASTER':mn_r2_5_up,
                     'R2_5_STD_RASTER':std_r2_5,
                     'R2_5_STD_UPRASTER':std_r2_5_up,
                     'R5_MN_RASTER':mn_r5,
                     'R5_MN_UPRASTER':mn_r5_up,
                     'R5_STD_RASTER':std_r5,
                     'R5_STD_UPRASTER':std_r5_up,
                     'R10_MN_RASTER':mn_r10,
                     'R10_MN_UPRASTER':mn_r10_up,
                     'R10_STD_RASTER':std_r10,
                     'R10_STD_UPRASTER':std_r10_up,
                     'QS_MN_RASTER':mn_qs,
                     'QS_MN_UPRASTER':mn_qs_up,
                     'QS_STD_RASTER':std_qs,
                     'QS_STD_UPRASTER':std_qs_up,
                     'QSM_MN_RASTER':mn_qsm,
                     'QSM_MN_UPRASTER':mn_qsm_up,
                     'QSM_STD_RASTER':std_qsm,
                     'QSM_STD_UPRASTER':std_qsm_up,
                     'QSB_MN_RASTER':mn_qsb,
                     'QSB_MN_UPRASTER':mn_qsb_up,
                     'QSB_STD_RASTER':std_qsb,
                     'QSB_STD_UPRASTER':std_qsb_up,
                     'MAT_MN_RASTER':mn_mat,
                     'MAT_MN_UPRASTER':mn_mat_up,
                     'MAT_STD_RASTER':std_mat,
                     'MAT_STD_UPRASTER':std_mat_up,})
df_out.to_csv('grdc_wrr2_raster_outputs.csv')