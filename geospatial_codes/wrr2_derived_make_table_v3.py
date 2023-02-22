#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 12 14:36:12 2022

@author: aforte
"""

from netCDF4 import Dataset
import numpy as np
import pandas as pd
import rasterio

master_location='/Users/aforte/Documents/Python/snowmelt/'
# master_location='/Volumes/Choruh/Data/snowmelt_project/'

fh=Dataset(master_location+'wrr2_derived_v2.nc')

p_mean=fh.variables['pmean'][:,:]
p_std=fh.variables['pstd'][:,:]
r_mean=fh.variables['rmean'][:,:]
r_std=fh.variables['rstd'][:,:]
et_mean=fh.variables['etmean'][:,:]
et_std=fh.variables['etstd'][:,:]
pet_mean=fh.variables['petmean'][:,:]
pet_std=fh.variables['petstd'][:,:]
p_c=fh.variables['pshp'][:,:]
p_s=fh.variables['pscl'][:,:]
r_c1=fh.variables['rshp1'][:,:]
r_s1=fh.variables['rscl1'][:,:]
r_s5=fh.variables['rscl5'][:,:]
r_c5=fh.variables['rshp5'][:,:]
rnu=fh.variables['rnu'][:,:]
r2_5=fh.variables['r25'][:,:]
r5=fh.variables['r5'][:,:]
r10=fh.variables['r10'][:,:]
t_mean=fh.variables['tmean'][:,:]
pe_sn=fh.variables['ps'][:,:]
qs=fh.variables['qs'][:,:]
qsm=fh.variables['qsm'][:,:]
qsb=fh.variables['qsb'][:,:]
rcom=fh.variables['rcom'][:,:]
pcom=fh.variables['pcom'][:,:]
latv=fh.variables['lat'][:]
lonv=fh.variables['lon'][:]

fh2=Dataset(master_location+'wrr2_derived_Q.nc')

q_mean=fh2.variables['qmean'][:,:]
q_std=fh2.variables['qstd'][:,:]
q_c1=fh2.variables['qshp1'][:,:]
q_s1=fh2.variables['qscl1'][:,:]
q_s5=fh2.variables['qscl5'][:,:]
q_c5=fh2.variables['qshp5'][:,:]
qnu=fh2.variables['qnu'][:,:]
q2_5=fh2.variables['q25'][:,:]
q5=fh2.variables['q5'][:,:]
q10=fh2.variables['q10'][:,:]
qcom=fh2.variables['qcom'][:,:]
swe_mean=fh2.variables['swemean'][:,:]
swe_std=fh2.variables['swestd'][:,:]

[LON,LAT]=np.meshgrid(lonv,latv)

rc1=r_c1[~r_c1.mask]
rc5=r_c5[~r_c5.mask]
rs1=r_s1[~r_s1.mask]
rs5=r_s5[~r_s5.mask]
r2_5=r2_5[~r2_5.mask]
r5=r5[~r5.mask]
r10=r10[~r10.mask]
rm=r_mean[~r_mean.mask]
rst=r_std[~r_std.mask]
p=p_mean[~p_mean.mask]
pst=p_std[~p_std.mask]
rcom=rcom[~rcom.mask]
pcom=pcom[~pcom.mask]
lat=LAT[~r_c1.mask]
lon=LON[~r_c1.mask]
pe_sn=pe_sn[~pe_sn.mask]
t=t_mean[~t_mean.mask]
pc=p_c[~p_c.mask]
ps=p_s[~p_s.mask]
rnu=rnu[~rnu.mask]
qsb=qsb[~qsb.mask]
qs=qs[~qs.mask]
qsm=qsm[~qsm.mask]
et=et_mean[~et_mean.mask]
etst=et_std[~et_std.mask]
pet=pet_mean[~pet_mean.mask]
petst=pet_std[~pet_std.mask]

qc1=q_c1[~q_c1.mask]
qc5=q_c5[~q_c5.mask]
qs1=q_s1[~q_s1.mask]
qs5=q_s5[~q_s5.mask]
q2_5=q2_5[~q2_5.mask]
q5=q5[~q5.mask]
q10=q10[~q10.mask]
qm=q_mean[~q_mean.mask]
qst=q_std[~q_std.mask]
qcom=qcom[~qcom.mask]
qnu=qnu[~qnu.mask]
swem=swe_mean[~swe_mean.mask]
swest=swe_std[~swe_std.mask]


# Load one raster in early to deal with bounds that are not fully global
rstr=rasterio.open(master_location+'hyd_glo_dem_15s/hyd_glo_dem_15s.tif')
dem=rstr.read(1)
bnds=rstr.bounds


idx=[(~np.isnan(qc1)) & (~np.isnan(qm)) & (~np.isnan(p)) & (p>0) & (rm>0) & (qm>0) &
     (qcom==1) & (~np.isnan(t)) & (qnu>500)  & (~np.isclose(r2_5,r5,1e-2,1e-5)) 
     & (~np.isclose(r5,r10,1e-2,1e-5)) & (~np.isclose(q2_5,q5,1e-2,1e-5)) 
     & (~np.isclose(q5,q10,1e-2,1e-5)) & (lat>bnds[1]+(0.25/2)) & (lat<bnds[3]-(0.25/2))] 

rc1=rc1[idx]
rc5=rc5[idx]
rs1=rs1[idx]
rs5=rs5[idx]
r2_5=r2_5[idx]
r5=r5[idx]
r10=r10[idx]
rm=rm[idx]
rst=rst[idx]
p=p[idx]
pst=pst[idx]
pc=pc[idx]
ps=ps[idx]
rcom=rcom[idx]
pcom=pcom[idx]
lat=lat[idx]
lon=lon[idx]
pe_sn=pe_sn[idx]
t=t[idx]
rnu=rnu[idx]
qs=qs[idx]
qsm=qsm[idx]
qsb=qsb[idx]
et=et[idx]
etst=etst[idx]
pet=pet[idx]
petst=petst[idx]

qc1=qc1[idx]
qc5=qc5[idx]
qs1=qs1[idx]
qs5=qs5[idx]
q2_5=q2_5[idx]
q5=q5[idx]
q10=q10[idx]
qm=qm[idx]
qst=qst[idx]
qcom=qcom[idx]
qnu=qnu[idx]
swem=swem[idx]
swest=swest[idx]

fh.close()
fh2.close()

# TOPO DATA
rstr2=rasterio.open(master_location+'hyd_glo_aca_15s/hyd_glo_aca_15s.tif')
area=rstr2.read(1)

rstr1=rasterio.open(master_location+'hyd_glo_dem_15s/rlf_17nghbr.tif')
rlf=rstr1.read(1)

rstr3=rasterio.open(master_location+'hyd_glo_dem_15s/maxZ_upstream.tif')
zup=rstr3.read(1)

r_res=360/rstr.width

z_min=np.zeros(lon.shape)
z_mean=np.zeros(lon.shape)
z_max=np.zeros(lon.shape)

a_max=np.zeros(lon.shape)
zu_max=np.zeros(lon.shape)
zu_mean=np.zeros(lon.shape)

rlf_min=np.zeros(lon.shape)
rlf_mean=np.zeros(lon.shape)
rlf_max=np.zeros(lon.shape)

# Accumulated Upstream Values
p_rstr=rasterio.open(master_location+'wrr2_raster_outputs/upstream_means/wrr2_up_P.tif')
pup=p_rstr.read(1)
pc_rstr=rasterio.open(master_location+'wrr2_raster_outputs/upstream_means/wrr2_up_PC.tif')
pcup=pc_rstr.read(1)
ps_rstr=rasterio.open(master_location+'wrr2_raster_outputs/upstream_means/wrr2_up_PS.tif')
psup=ps_rstr.read(1)
rc1_rstr=rasterio.open(master_location+'wrr2_raster_outputs/upstream_means/wrr2_up_RC1.tif')
rc1up=rc1_rstr.read(1)
rs1_rstr=rasterio.open(master_location+'wrr2_raster_outputs/upstream_means/wrr2_up_RS1.tif')
rs1up=rs1_rstr.read(1)
rc5_rstr=rasterio.open(master_location+'wrr2_raster_outputs/upstream_means/wrr2_up_RC5.tif')
rc5up=rc5_rstr.read(1)
rs5_rstr=rasterio.open(master_location+'wrr2_raster_outputs/upstream_means/wrr2_up_RS5.tif')
rs5up=rs5_rstr.read(1)
psnow_rstr=rasterio.open(master_location+'wrr2_raster_outputs/upstream_means/wrr2_up_PSNOW.tif')
psnowup=psnow_rstr.read(1)
pbase_rstr=rasterio.open(master_location+'wrr2_raster_outputs/upstream_means/wrr2_up_PBASE.tif')
pbaseup=pbase_rstr.read(1)
mat_rstr=rasterio.open(master_location+'wrr2_raster_outputs/upstream_means/wrr2_up_MAT.tif')
matup=mat_rstr.read(1)
rcom_rstr=rasterio.open(master_location+'wrr2_raster_outputs/upstream_means/wrr2_up_RCOM.tif')
rcomup=rcom_rstr.read(1)
r_rstr=rasterio.open(master_location+'wrr2_raster_outputs/upstream_means/wrr2_up_R.tif')
rup=r_rstr.read(1)
up_res=360/rstr.width

pup_a=np.zeros(lon.shape)
pcup_a=np.zeros(lon.shape)
psup_a=np.zeros(lon.shape)
rup_a=np.zeros(lon.shape)
rc1up_a=np.zeros(lon.shape)
rs1up_a=np.zeros(lon.shape)
rc5up_a=np.zeros(lon.shape)
rs5up_a=np.zeros(lon.shape)
psnowup_a=np.zeros(lon.shape)
pbaseup_a=np.zeros(lon.shape)
matup_a=np.zeros(lon.shape)
rcomup_a=np.zeros(lon.shape)

for i in range(len(lon)):
    if np.mod(i,1000)==0:
        print(i)
    
    x=lon[i]
    y=lat[i]
    ## TOPO
    xv=np.arange(x-0.25/2,(x+0.25/2),r_res)
    yv=np.arange(y-0.25/2,(y+0.25/2),r_res)
    [X,Y]=np.meshgrid(xv,yv)
    [r,c]=rstr.index(X.ravel(),Y.ravel())
    z=dem[r,c]
    # Filter out no data values
    idx=z==rstr.nodata
    zf=z[~idx]
    if zf.shape[0]>0:
        z_min[i]=np.amin(zf)
        z_mean[i]=np.mean(zf)
        z_max[i]=np.amax(zf)
    else:
        z_min[i]=np.nan 
        z_mean[i]=np.nan
        z_max[i]=np.nan
        
    a=area[r,c]
    # Filter out no data values
    idx=a==rstr2.nodata
    af=a[~idx]
    if af.shape[0]>0:
        a_max[i]=np.amax(af/100) # Convert to km2 from hectares
    else:
        a_max[i]=np.nan 
        
    zu=zup[r,c]
    # Filter out no data values
    idx=zu==rstr3.nodata
    zupf=zu[~idx]
    if zupf.shape[0]>0:
        zu_max[i]=np.amax(zupf)
        zu_mean[i]=np.mean(zupf)
    else:
        zu_max[i]=np.nan 
        zu_mean[i]=np.nan
     
    rlf0=rlf[r,c]
    idx=rlf0==rstr1.nodata
    rlf0f=rlf0[~idx]
    if rlf0f.shape[0]>0:
        rlf_min[i]=np.nanmin(rlf0f)
        rlf_mean[i]=np.nanmean(rlf0f)
        rlf_max[i]=np.nanmax(rlf0f)
    else:
        rlf_min[i]=np.nan
        rlf_mean[i]=np.nan
        rlf_max[i]=np.nan
        
    ## UPSTREAM RASTERS
    xv=np.arange(x-0.25/2,(x+0.25/2),up_res)
    yv=np.arange(y-0.25/2,(y+0.25/2),up_res)
    [X,Y]=np.meshgrid(xv,yv)
    [r,c]=p_rstr.index(X.ravel(),Y.ravel())

    val=pup[r,c]
    idx=np.isnan(val)
    valf=val[~idx]
    if valf.shape[0]>0:
        pup_a[i]=np.nanmax(valf)
    else:
        pup_a[i]=np.nan
        
    val=pcup[r,c]
    idx=np.isnan(val)
    valf=val[~idx]
    if valf.shape[0]>0:
        pcup_a[i]=np.nanmax(valf)
    else:
        pcup_a[i]=np.nan
    
    val=psup[r,c]
    idx=np.isnan(val)
    valf=val[~idx]
    if valf.shape[0]>0:
        psup_a[i]=np.nanmax(valf)
    else:
        psup_a[i]=np.nan
    
    val=rup[r,c]
    idx=np.isnan(val)
    valf=val[~idx]
    if valf.shape[0]>0:
        rup_a[i]=np.nanmax(valf)
    else:
        rup_a[i]=np.nan
        
    val=rc1up[r,c]
    idx=np.isnan(val)
    valf=val[~idx]
    if valf.shape[0]>0:
        rc1up_a[i]=np.nanmax(valf)
    else:
        rc1up_a[i]=np.nan
        
    val=rs1up[r,c]
    idx=np.isnan(val)
    valf=val[~idx]
    if valf.shape[0]>0:
        rs1up_a[i]=np.nanmax(valf)
    else:
        rs1up_a[i]=np.nan 
        
    val=rc5up[r,c]
    idx=np.isnan(val)
    valf=val[~idx]
    if valf.shape[0]>0:
        rc5up_a[i]=np.nanmax(valf)
    else:
        rc5up_a[i]=np.nan
        
    val=rs5up[r,c]
    idx=np.isnan(val)
    valf=val[~idx]
    if valf.shape[0]>0:
        rs5up_a[i]=np.nanmax(valf)
    else:
        rs5up_a[i]=np.nan

    val=psnowup[r,c]
    idx=np.isnan(val)
    valf=val[~idx]
    if valf.shape[0]>0:
        psnowup_a[i]=np.nanmax(valf)
    else:
        psnowup_a[i]=np.nan
        
    val=pbaseup[r,c]
    idx=np.isnan(val)
    valf=val[~idx]
    if valf.shape[0]>0:
        pbaseup_a[i]=np.nanmax(valf)
    else:
        pbaseup_a[i]=np.nan
        
    val=matup[r,c]
    idx=np.isnan(val)
    valf=val[~idx]
    if valf.shape[0]>0:
        matup_a[i]=np.nanmean(valf)
    else:
        matup_a[i]=np.nan
        
    val=rcomup[r,c]
    idx=np.isnan(val)
    valf=val[~idx]
    if valf.shape[0]>0:
        rcomup_a[i]=np.nanmean(valf)
    else:
        rcomup_a[i]=np.nan
    
df=pd.DataFrame(data={'longitude':lon,
                      'latitude':lat,
                      'min_z':z_min,
                      'mean_z':z_mean,
                      'max_z':z_max,
                      'max_A':a_max,
                      'mean_zu':zu_mean,
                      'max_zu':zu_max,
                      'min_rlf':rlf_min,
                      'mean_rlf':rlf_mean,
                      'max_rlf':rlf_max,
                      'mean_runoff':rm,
                      'mean_precip':p,
                      'r_std':rst,
                      'r_c1':rc1,
                      'r_c5':rc5,
                      'r_s1':rs1,
                      'r_s5':rs5,
                      'p_c':pc,
                      'p_s':ps,
                      'p_std':pst,
                      'r2_5':r2_5,
                      'r5':r5,
                      'r10':r10,
                      'pup':pup_a,
                      'pcup':pcup_a,
                      'psup':psup_a,
                      'rup':rup_a,
                      'rc1up':rc1up_a,
                      'rs1up':rs1up_a,
                      'rc5up':rc5up_a,
                      'rs5up':rs5up_a,
                      'psnowup':psnowup_a,
                      'pbaseup':pbaseup_a,
                      'rcomup':rcomup_a,
                      'matup':matup_a,
                      'mean_discharge':qm,
                      'q_std':qst,
                      'q_c1':qc1,
                      'q_s1':qs1,
                      'q_c5':qc5,
                      'q_s5':qs5,
                      'q10':q10,
                      'q5':q5,
                      'q2_5':q2_5,
                      'mean_swe':swem,
                      'std_swe':swest,
                      'qs':qs,
                      'qsm':qsm,
                      'qsb':qsb,
                      'pe_sn':pe_sn,
                      'mat':t,
                      'et_mean':et,
                      'et_std':etst,
                      'pet_mean':pet,
                      'pet_std':petst})
df.to_csv('wrr2_derived_data_v3.csv')


