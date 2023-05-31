#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 23 09:29:50 2023

@author: aforte
"""

from netCDF4 import Dataset,num2date
import numpy as np
import pandas as pd
import time
import matplotlib.pyplot as plt
from scipy import ndimage
import rasterio.transform as rt
import rasterio


# Lots of random warnings for specific cells
import warnings
warnings.filterwarnings('ignore')

def estimate_area(clat,dlat,dlon,area):
    # Set constants
    # Equatorial radius
    a=6378137 #m
    # Eccentricity squared
    e2=0.00669437999014
    # 1 degree of latitude in km
    delta_lat=((np.pi * a * (1-e2))/(180 * (1-e2*(np.sin(np.radians(clat))**2))**(3/2)))/1000
    delta_lon=((np.pi * a * np.cos(np.radians(clat)))/(180 * np.sqrt(1 - e2*(np.sin(np.radians(clat))**2))))/1000
    # Estimate rectangular area of bounding box
    rect_area=(delta_lat*dlat)*(delta_lon*dlon)
    # Estimate eliptical area
    radius=np.sqrt(area/np.pi) # Radius in degrees
    rad_lat=radius*delta_lat # Radius in latitude direction
    rad_lon=radius*delta_lon # Radisu in longitude direction
    ellip_area=rad_lat*rad_lon*np.pi
    return rect_area,ellip_area    
    
def identify_events(tr_rstr,qs,qsm,qsb,threshold,transform,date,rstrd,dem):
    a=tr_rstr>=threshold
    # Generate labels
    lbl,nlbl=ndimage.label(a)
    
    # Control for empty array
    if nlbl>0:    
        # Index of all labels
        lbls=np.arange(1,nlbl+1)
        # Find sum of values within labeled areas
        sums=ndimage.labeled_comprehension(a,lbl,lbls,np.sum,int,0)
        # Find sums and means of different components
        qs_sums=ndimage.labeled_comprehension(qs,lbl,lbls,np.sum,float,0)
        qsm_sums=ndimage.labeled_comprehension(qsm,lbl,lbls,np.sum,float,0)
        qsb_sums=ndimage.labeled_comprehension(qsb,lbl,lbls,np.sum,float,0)
        qs_mean=ndimage.mean(qs,lbl,lbls)
        qsm_mean=ndimage.mean(qsm,lbl,lbls)
        qsb_mean=ndimage.mean(qsb,lbl,lbls)
        # Find centers
        centers=ndimage.center_of_mass(a,lbl,lbls)
        # Find objects
        objs=ndimage.find_objects(lbl)
        # Find sizes of bounding box
        bbox_size=np.array([[s.stop-s.start for s in object_slice] for object_slice in objs])
        
        center_x=np.zeros((nlbl))
        center_y=np.zeros((nlbl))
        dim_x=np.zeros((nlbl))
        dim_y=np.zeros((nlbl))
        area_pixels=np.zeros((nlbl))
        area=np.zeros((nlbl))
        rect_area=np.zeros((nlbl))
        ellip_area=np.zeros((nlbl))
        mean_z=np.zeros((nlbl))
        qs_sum=np.zeros((nlbl))
        qsm_sum=np.zeros((nlbl))
        qsb_sum=np.zeros((nlbl))
        qs_mn=np.zeros((nlbl))
        qsm_mn=np.zeros((nlbl))
        qsb_mn=np.zeros((nlbl))
        
        date_list=[]
        
        for i in range(nlbl):
            lonlat=rt.xy(transform,centers[i][0],centers[i][1])
            
            center_y[i]=lonlat[1]
            center_x[i]=lonlat[0]
            dim_y[i]=bbox_size[i,0] * 0.25 # account for size of pixel
            dim_x[i]=bbox_size[i,1] * 0.25 # account for size of pixel 
            area_pixels[i]=sums[i]
            area[i]=sums[i] * (0.25**2) # account for size of pixel
            date_list.append(date)
            qs_sum[i]=qs_sums[i]
            qsm_sum[i]=qsm_sums[i]
            qsb_sum[i]=qsb_sums[i]
            qs_mn[i]=qs_mean[i]
            qsm_mn[i]=qsm_mean[i]
            qsb_mn[i]=qsb_mean[i]
            
            # Estimate areas
            rect_area[i],ellip_area[i]=estimate_area(np.abs(center_y[i]),dim_y[i],dim_x[i],area[i])
            
            # Determine bounds
            ur=rstrd.index(center_x[i]+dim_x[i]/2,center_y[i]+dim_y[i]/2)
            ul=rstrd.index(center_x[i]-dim_x[i]/2,center_y[i]+dim_y[i]/2)
            lr=rstrd.index(center_x[i]+dim_x[i]/2,center_y[i]-dim_y[i]/2)
            ll=rstrd.index(center_x[i]-dim_x[i]/2,center_y[i]-dim_y[i]/2)
    
            z=dem[ul[0]:ll[0],ul[1]:ur[1]]
            idx=z==rstrd.nodata
            mean_z[i]=np.mean(z[~idx])
    else:
        date_list=[date]
        center_x=np.nan
        center_y=np.nan
        dim_x=np.nan 
        dim_y=np.nan 
        area=np.nan
        area_pixels=np.nan
        rect_area=np.nan 
        ellip_area=np.nan 
        mean_z=np.nan
        qs_sum=np.nan 
        qsm_sum=np.nan 
        qsb_sum=np.nan 
        qs_mn=np.nan 
        qsm_mn=np.nan 
        qsb_mn=np.nan
        
    df=pd.DataFrame(data={'date':date_list,
                          'center_lon':center_x,
                          'center_lat':center_y,
                          'dim_lon':dim_x,
                          'dim_lat':dim_y,
                          'area_deg':area,
                          'area_pixels':area_pixels,
                          'rect_area_km':rect_area,
                          'elliptical_area_km':ellip_area,
                          'mean_z':mean_z,
                          'qs_sum':qs_sum,
                          'qsm_sum':qsm_sum,
                          'qsb_sum':qsb_sum,
                          'qs_mn':qs_mn,
                          'qsm_mn':qsm_mn,
                          'qsb_mn':qsb_mn})
    return df


##### THRESHOLD IN mm/day
threshold_list=[25,30]

master_location='/Volumes/Choruh/Data/snowmelt_project/'

# Precipitation
P_layer='Precip'
P_mfile=[master_location+'e2o_univk_wrr2/e2o_univk_wrr2_glob15_day_Precip_1980-1989.nc',
         master_location+'e2o_univk_wrr2/e2o_univk_wrr2_glob15_day_Precip_1990-1999.nc']
# Surface runoff
Qs_layer='Qs'
Qs_mfile=[master_location+'e2o_univk_wrr2/e2o_univk_wrr2_glob15_day_Qs_1980-1989.nc',
          master_location+'e2o_univk_wrr2/e2o_univk_wrr2_glob15_day_Qs_1990-1999.nc']
# Snowmelt
Qsm_layer='Qsm'
Qsm_mfile=[master_location+'e2o_univk_wrr2/e2o_univk_wrr2_glob15_day_Qsm_1980-1989.nc',
          master_location+'e2o_univk_wrr2/e2o_univk_wrr2_glob15_day_Qsm_1990-1999.nc']
# Baseflow
Qsb_layer='Qsb'
Qsb_mfile=[master_location+'e2o_univk_wrr2/e2o_univk_wrr2_glob15_day_Qsb_1980-1989.nc',
          master_location+'e2o_univk_wrr2/e2o_univk_wrr2_glob15_day_Qsb_1990-1999.nc']


rstr=rasterio.open(master_location+'srtm30_plus/wc2.1_30s_elev.tif')
dem=rstr.read(1)


print('Rasters Loaded')

# Open Datasets
P_fh0=Dataset(P_mfile[0])
P_fh1=Dataset(P_mfile[1])
Qs_fh0=Dataset(Qs_mfile[0])
Qs_fh1=Dataset(Qs_mfile[1])
Qsm_fh0=Dataset(Qsm_mfile[0])
Qsm_fh1=Dataset(Qsm_mfile[1])
Qsb_fh0=Dataset(Qsb_mfile[0])
Qsb_fh1=Dataset(Qsb_mfile[1])

# Open Datasets
P_fh0=Dataset(P_mfile[0])
P_fh1=Dataset(P_mfile[1])
Qs_fh0=Dataset(Qs_mfile[0])
Qs_fh1=Dataset(Qs_mfile[1])
Qsm_fh0=Dataset(Qsm_mfile[0])
Qsm_fh1=Dataset(Qsm_mfile[1])
Qsb_fh0=Dataset(Qsb_mfile[0])
Qsb_fh1=Dataset(Qsb_mfile[1])

# Create indices to not select leap days
T0=P_fh0.variables['time'][:]
T0obj=P_fh0.variables['time']
dt0=num2date(T0obj[:],T0obj.units)
DT0=pd.DataFrame(data={'date':dt0})
DT0['date']=pd.to_datetime(DT0['date'].astype('str'))
day=pd.DatetimeIndex(DT0['date']).day
mnth=pd.DatetimeIndex(DT0['date']).month
DT0=DT0.drop(index=DT0.index[np.logical_and(day==29,mnth==2)])
idx0=DT0.index

T1=P_fh1.variables['time'][:]
T1obj=P_fh1.variables['time']
dt1=num2date(T1obj[:],T1obj.units)
DT1=pd.DataFrame(data={'date':dt1})
DT1['date']=pd.to_datetime(DT1['date'].astype('str'))
day=pd.DatetimeIndex(DT1['date']).day
mnth=pd.DatetimeIndex(DT1['date']).month
DT1=DT1.drop(index=DT1.index[np.logical_and(day==29,mnth==2)])
idx1=DT1.index

# Create a transform
lat=P_fh1.variables['lat'][:]
lon=P_fh1.variables['lon'][:]
res=lat[1]-lat[0]
[x,y]=np.meshgrid(lon,lat)
transform=rt.Affine.translation(x[0][0]-res/2,y[0][0]-res/2)*rt.Affine.scale(res,res)

# Define year list    
Y0=np.arange(1980,1990,1)
Y1=np.arange(1990,2000,1)
  

# Define block indices
bix=np.arange(0,3651,365)
for l in range(len(threshold_list)):
    threshold=threshold_list[l]
    k=0  
    # Begin loop
    for i in range(len(bix)-1):
        print('Working on block '+str(i+1)+' of '+str(len(bix)-1))
        
        P0=[]
        P1=[]
        R0=[]
        R1=[]
    
        for j in range(bix[i],bix[i+1]):
            if np.mod(j+1,25)==0:
                print('      day '+str(j+1))
            # Grab index
            ix0=idx0[j]
            ix1=idx1[j]
            
            date0=DT0.iloc[j][0]
            date1=DT1.iloc[j][0]
            
            # 1980-1989
            p0=P_fh0.variables[P_layer][ix0,:,:]
            qs0=Qs_fh0.variables[Qs_layer][ix0,:,:]
            qsm0=Qsm_fh0.variables[Qsm_layer][ix0,:,:]
            qsb0=Qsb_fh0.variables[Qsb_layer][ix0,:,:]
            # 1990-1999
            p1=P_fh1.variables[P_layer][ix1,:,:]
            qs1=Qs_fh1.variables[Qs_layer][ix1,:,:]
            qsm1=Qsm_fh1.variables[Qsm_layer][ix1,:,:]
            qsb1=Qsb_fh1.variables[Qsb_layer][ix1,:,:]  
            
            # Flatten masked arrays
            pf0=p0.filled(np.nan)
            pf1=p1.filled(np.nan)
            qsf0=qs0.filled(np.nan)
            qsf1=qs1.filled(np.nan)
            qsmf0=qsm0.filled(np.nan)
            qsmf1=qsm1.filled(np.nan)
            qsbf0=qsb0.filled(np.nan)
            qsbf1=qsb1.filled(np.nan)
            
            # Transform
            pf0=(pf0/1000)*(60*60*24*100*10)
            qsf0=(qsf0/1000)*(60*60*24*100*10)*-1
            qsmf0=(qsmf0/1000)*(60*60*24*100*10)
            qsbf0=(qsbf0/1000)*(60*60*24*100*10)*-1
            rf0=qsf0+qsmf0+qsbf0
            pf0[pf0<0]=0
            rf0[rf0<0]=0
            pf1=(pf1/1000)*(60*60*24*100*10)
            qsf1=(qsf1/1000)*(60*60*24*100*10)*-1
            qsmf1=(qsmf1/1000)*(60*60*24*100*10)
            qsbf1=(qsbf1/1000)*(60*60*24*100*10)*-1 
            rf1=qsf1+qsmf1+qsbf1
            pf1[pf1<0]=0
            rf1[rf1<0]=0
            
            pdf0=identify_events(pf0,pf0,pf0,pf0,threshold,transform,date0,rstr,dem)
            rdf0=identify_events(rf0,qsf0,qsmf0,qsbf0,threshold,transform,date0,rstr,dem)
            pdf1=identify_events(pf1,pf1,pf1,pf1,threshold,transform,date1,rstr,dem)
            rdf1=identify_events(rf1,qsf1,qsmf1,qsbf1,threshold,transform,date1,rstr,dem)
            
            P0.append(pdf0)
            P1.append(pdf1)
            R0.append(rdf0)
            R1.append(rdf1)
        
        P0out=pd.concat(P0,ignore_index=True)
        P1out=pd.concat(P1,ignore_index=True)
        R0out=pd.concat(R0,ignore_index=True)
        R1out=pd.concat(R1,ignore_index=True)
        
        P0out.to_csv(master_location+'wrr2_events/Precip_Threshold_'+str(threshold)+'_'+str(Y0[k])+'.csv',index=False)
        P1out.to_csv(master_location+'wrr2_events/Precip_Threshold_'+str(threshold)+'_'+str(Y1[k])+'.csv',index=False)
        R0out.to_csv(master_location+'wrr2_events/Runoff_Threshold_'+str(threshold)+'_'+str(Y0[k])+'.csv',index=False)
        R1out.to_csv(master_location+'wrr2_events/Runoff_Threshold_'+str(threshold)+'_'+str(Y1[k])+'.csv',index=False)
        
        k=k+1
    
    
        
P_fh0.close()
P_fh1.close()
Qs_fh0.close()
Qs_fh1.close()
Qsm_fh0.close()
Qsm_fh1.close()
Qsb_fh0.close()
Qsb_fh1.close()