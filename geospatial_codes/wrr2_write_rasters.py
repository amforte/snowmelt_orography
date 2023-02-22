# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

from netCDF4 import Dataset
import numpy as np
import pandas as pd
import rasterio
import rasterio.transform as rt
from rasterio.enums import Resampling

def write_tiffs(Z,name,upscale_factor):
    # Write rasters
    with rasterio.open(
            master_location+'wrr2_raster_outputs/'+name+'.tif',
            'w',
            driver='GTiff',
            height=Z.shape[0],
            width=Z.shape[1],
            count=1,
            dtype=Z.dtype,
            crs='+proj=latlong',
            transform=transform,
            nodata=-9999,
        ) as dst: 
            dst.write(Z,1)
        
    with rasterio.open(master_location+'wrr2_raster_outputs/'+name+'.tif') as dataset:
    
        # resample data to target shape
        data = dataset.read(
            out_shape=(
                dataset.count,
                int(dataset.height * upscale_factor),
                int(dataset.width * upscale_factor)
            ),
            resampling=Resampling.nearest
        )
    
        # scale image transform
        transform_up = dataset.transform * dataset.transform.scale(
            (dataset.width / data.shape[-1]),
            (dataset.height / data.shape[-2])
        )    
    
    Z_up=data[0,:,:]
    with rasterio.open(
            master_location+'wrr2_raster_outputs/'+name+'_upsample.tif',
            'w',
            driver='GTiff',
            height=Z_up.shape[0],
            width=Z_up.shape[1],
            count=1,
            dtype=Z_up.dtype,
            crs='+proj=latlong',
            transform=transform_up,
            nodata=-9999,
        ) as dst: 
            dst.write(Z_up,1)    

master_location='/Volumes/Choruh/Data/snowmelt_project/'

# Load global filtered dataset
df=pd.read_csv(master_location+'wrr2_derived_data_v3.csv')
df=df.drop(['Unnamed: 0'],axis=1)
df=df.drop(index=df.index[df['p_s']>35])
df=df.drop(index=df.index[np.isnan(df['mean_z'])])

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
row_col=rt.rowcol(transform,df['longitude'],df['latitude'])
rows=row_col[0]
cols=row_col[1]
# Generate blank images
P=np.zeros(x.shape)
PC=np.zeros(x.shape)
R=np.zeros(x.shape)
RC1=np.zeros(x.shape)
RC5=np.zeros(x.shape)
R2_5=np.zeros(x.shape)
R5=np.zeros(x.shape)
R10=np.zeros(x.shape)
RS1=np.zeros(x.shape)
RS5=np.zeros(x.shape)
PS=np.zeros(x.shape)
QS=np.zeros(x.shape)
QSM=np.zeros(x.shape)
QSB=np.zeros(x.shape)
MAT=np.zeros(x.shape)

P[:,:]=-9999
PC[:,:]=-9999
R[:,:]=-9999
RC1[:,:]=-9999
RC5[:,:]=-9999
R2_5[:,:]=-9999
R5[:,:]=-9999
R10[:,:]=-9999
RS1[:,:]=-9999
RS5[:,:]=-9999
PS[:,:]=-9999
QS[:,:]=-9999
QSM[:,:]=-9999
QSB[:,:]=-9999
MAT[:,:]=-9999

# Populate rasters
P[rows,cols]=df['mean_precip'] 
PC[rows,cols]=df['p_c']
R[rows,cols]=df['mean_runoff']
RC1[rows,cols]=df['r_c1']
RC5[rows,cols]=df['r_c5']
R2_5[rows,cols]=df['r2_5']
R5[rows,cols]=df['r5']
R10[rows,cols]=df['r10']
RS1[rows,cols]=df['r_s1']
RS5[rows,cols]=df['r_s5']
PS[rows,cols]=df['p_s']
QS[rows,cols]=df['qs']
QSM[rows,cols]=df['qsm']
QSB[rows,cols]=df['qsb']
MAT[rows,cols]=df['mat']

# Write TIFFs
write_tiffs(P,'wrr2_P',10)
write_tiffs(PC,'wrr2_PC',10)
write_tiffs(R,'wrr2_R',10)
write_tiffs(RC1,'wrr2_RC1',10)
write_tiffs(RC5,'wrr2_RC5',10)
write_tiffs(R2_5,'wrr2_R2_5',10)
write_tiffs(R5,'wrr2_R5',10)
write_tiffs(R10,'wrr2_R10',10)
write_tiffs(RS1,'wrr2_RS1',10)
write_tiffs(RS5,'wrr2_RS5',10)
write_tiffs(PS,'wrr2_PS',10)
write_tiffs(QS,'wrr2_QS',10)
write_tiffs(QSM,'wrr2_QSM',10)
write_tiffs(QSB,'wrr2_QSB',10)
write_tiffs(MAT,'wrr2_MAT',10)
