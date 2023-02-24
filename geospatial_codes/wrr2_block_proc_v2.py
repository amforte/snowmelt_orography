#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script written by Adam M. Forte
aforte8@lsu.edu

This script processes individual pixels within the global daily WaterGAP3 
timeseries and calculates various time series quantities (e.g., mean runoff,
variability, etc.) and saves them in a single netcdf4 file, 'wrr2_derived_v2.nc'
"""

from netCDF4 import Dataset,num2date
import numpy as np
import pandas as pd
import time
### xarray is stupid ####

def survive(ts):
    ts_sort=np.sort(ts)
    tsn=len(ts_sort)
    tsrank=np.arange(1,tsn+1,1)
    ts_freq_excd=(tsn+1-tsrank)/tsn
    return ts_sort,ts_freq_excd

def weibull_tail_fit(x,y,thresh):
    ix=np.nonzero(y<thresh)[0][:1][0]
    xtrim=x[ix:len(x)]
    ytrim=y[ix:len(x)]
    xts=np.log(xtrim)
    yts=np.log(-np.log(ytrim))       
    [lin,r,rnk,sng,V]=np.polyfit(xts,yts,1,full=True)
    c=lin[0]
    s=np.exp(-1*lin[1]/c)
    return c,s

## Requires:
## Water Resources Reanalysis v2 - WaterGAP3 Data, available here:
## http://www.earth2observe.eu/

## Define master location
master_location='/Volumes/Choruh/Data/snowmelt_project/'

## Define file locations for datasets of interest
# Surfex-Trip Temperature
T_layer='AvgSurfT'
T_mfile=[master_location+'e2o_metfr_wrr2/e2o_metfr_wrr2_glob15_day_AvgSurfT_1980-1989.nc',
         master_location+'e2o_metfr_wrr2/e2o_metfr_wrr2_glob15_day_AvgSurfT_1990-1999.nc']
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

# Potential Evapotranspiration
PET_layer='PotEvap'
PET_mfile=[master_location+'e2o_univk_wrr2/e2o_univk_wrr2_glob15_day_PotEvap_1980-1989.nc',
           master_location+'e2o_univk_wrr2/e2o_univk_wrr2_glob15_day_PotEvap_1990-1999.nc']

# Evapotranspiration
ET_layer='Evap'
ET_mfile=[master_location+'e2o_univk_wrr2/e2o_univk_wrr2_glob15_day_Evap_1980-1989.nc',
           master_location+'e2o_univk_wrr2/e2o_univk_wrr2_glob15_day_Evap_1990-1999.nc']
 
## Open Datasets
T_fh0=Dataset(T_mfile[0])
T_fh1=Dataset(T_mfile[1])
P_fh0=Dataset(P_mfile[0])
P_fh1=Dataset(P_mfile[1])
Qs_fh0=Dataset(Qs_mfile[0])
Qs_fh1=Dataset(Qs_mfile[1])
Qsm_fh0=Dataset(Qsm_mfile[0])
Qsm_fh1=Dataset(Qsm_mfile[1])
Qsb_fh0=Dataset(Qsb_mfile[0])
Qsb_fh1=Dataset(Qsb_mfile[1])
PET_fh0=Dataset(PET_mfile[0])
PET_fh1=Dataset(PET_mfile[1])
ET_fh0=Dataset(ET_mfile[0])
ET_fh1=Dataset(ET_mfile[1])


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


T0=T_fh0.variables['time'][:]
T0obj=T_fh0.variables['time']
dt0=num2date(T0obj[:],T0obj.units)
DT0=pd.DataFrame(data={'date':dt0})
DT0['date']=pd.to_datetime(DT0['date'].astype('str'))
day=pd.DatetimeIndex(DT0['date']).day
mnth=pd.DatetimeIndex(DT0['date']).month
DT0=DT0.drop(index=DT0.index[np.logical_and(day==29,mnth==2)])
Tidx0=DT0.index

T1=T_fh1.variables['time'][:]
T1obj=T_fh1.variables['time']
dt1=num2date(T1obj[:],T1obj.units)
DT1=pd.DataFrame(data={'date':dt1})
DT1['date']=pd.to_datetime(DT1['date'].astype('str'))
day=pd.DatetimeIndex(DT1['date']).day
mnth=pd.DatetimeIndex(DT1['date']).month
DT1=DT1.drop(index=DT1.index[np.logical_and(day==29,mnth==2)])
Tidx1=DT1.index

# Define dimensions
r=720
c=1440

# Generate index grids
[CIX,RIX]=np.meshgrid(np.arange(0,c,1),np.arange(0,r,1))

# Define blocks
blocks=5
rix=np.arange(0,r+1,r/blocks).astype(int)
cix=np.arange(0,c+1,c/blocks).astype(int)

# Track block process
block_count=1
total_blocks=blocks**2

# Allocate arrays
p_mean=np.zeros((r,c))
p_std=np.zeros((r,c))
p_c=np.zeros((r,c))
p_s=np.zeros((r,c))
p_comp=np.zeros((r,c))

r_mean=np.zeros((r,c))
r_std=np.zeros((r,c))
r_c1=np.zeros((r,c))
r_s1=np.zeros((r,c))
r_c5=np.zeros((r,c))
r_s5=np.zeros((r,c))
r_comp=np.zeros((r,c))
r_nu=np.zeros((r,c))
r_2_5=np.zeros((r,c))
r_5=np.zeros((r,c))
r_10=np.zeros((r,c))

qs_mean=np.zeros((r,c))
qsb_mean=np.zeros((r,c))
qsm_mean=np.zeros((r,c))

pet_mean=np.zeros((r,c))
pet_std=np.zeros((r,c))
et_mean=np.zeros((r,c))
et_std=np.zeros((r,c))

t_mean=np.zeros((r,c))
t_std=np.zeros((r,c))

perc_s=np.zeros((r,c))

# Start master nested loop
for i in range(blocks):
    for j in range(blocks):
        st_t=time.time()
        # Block index
        bcix=CIX[rix[i]:rix[i+1],cix[j]:cix[j+1]]
        brix=RIX[rix[i]:rix[i+1],cix[j]:cix[j+1]]
        
        # Load precip blocks
        P_block0=P_fh0.variables[P_layer][idx0,rix[i]:rix[i+1],cix[j]:cix[j+1]]
        P_block1=P_fh1.variables[P_layer][idx1,rix[i]:rix[i+1],cix[j]:cix[j+1]]
        # Load full temp blocks
        T_block0=T_fh0.variables[T_layer][:,rix[i]:rix[i+1],cix[j]:cix[j+1]]
        T_block1=T_fh1.variables[T_layer][:,rix[i]:rix[i+1],cix[j]:cix[j+1]]
        for k in range(P_block0.shape[1]):
            for l in range(P_block0.shape[2]):
                MASK=~P_block0[0,:,:].mask
                if MASK[k,l]:
                    # Concatenate
                    p_ts=np.concatenate((P_block0[:,k,l],P_block1[:,k,l]))
                    t_ts=np.concatenate((T_block0[Tidx0,k,l],T_block1[Tidx1,k,l]))
                    # Convert time series to mm/day
                    p_ts=(p_ts/1000)*(60*60*24*100*10)
                    t_ts=t_ts-273.15
                    # Remove bad values
                    pidx=p_ts<0
                    o_len=len(p_ts)
                    p_ts[pidx]=0
                    n_len=o_len-np.sum(pidx)
                    p_comp[brix[k,l],bcix[k,l]]=n_len/o_len
                    # Do calculations
                    if n_len/o_len>0.50:
                        p_mean[brix[k,l],bcix[k,l]]=np.nanmean(p_ts)
                        p_std[brix[k,l],bcix[k,l]]=np.nanstd(p_ts)
                        [p_sort,p_freq]=survive(p_ts)
                        [pc,ps]=weibull_tail_fit(p_sort,p_freq,0.01)
                        p_c[brix[k,l],bcix[k,l]]=pc
                        p_s[brix[k,l],bcix[k,l]]=ps
                        perc_s[brix[k,l],bcix[k,l]]=np.sum(p_ts[t_ts<0])/np.sum(p_ts)
                    else:
                        p_mean[brix[k,l],bcix[k,l]]=np.nan
                        p_std[brix[k,l],bcix[k,l]]=np.nan
                        p_c[brix[k,l],bcix[k,l]]=np.nan
                        p_s[brix[k,l],bcix[k,l]]=np.nan
                        perc_s[brix[k,l],bcix[k,l]]=np.nan                        
                else:
                    p_mean[brix[k,l],bcix[k,l]]=np.nan
                    p_std[brix[k,l],bcix[k,l]]=np.nan
                    p_c[brix[k,l],bcix[k,l]]=np.nan
                    p_s[brix[k,l],bcix[k,l]]=np.nan
                    perc_s[brix[k,l],bcix[k,l]]=np.nan
                    p_comp[brix[k,l],bcix[k,l]]=np.nan
        # Load runoff blocks and et blocks
        Qs_block0=Qs_fh0.variables[Qs_layer][idx0,rix[i]:rix[i+1],cix[j]:cix[j+1]]
        Qs_block1=Qs_fh1.variables[Qs_layer][idx1,rix[i]:rix[i+1],cix[j]:cix[j+1]]
        Qsm_block0=Qsm_fh0.variables[Qsm_layer][idx0,rix[i]:rix[i+1],cix[j]:cix[j+1]]
        Qsm_block1=Qsm_fh1.variables[Qsm_layer][idx1,rix[i]:rix[i+1],cix[j]:cix[j+1]]
        Qsb_block0=Qsb_fh0.variables[Qsb_layer][idx0,rix[i]:rix[i+1],cix[j]:cix[j+1]]
        Qsb_block1=Qsb_fh1.variables[Qsb_layer][idx1,rix[i]:rix[i+1],cix[j]:cix[j+1]]
        PET_block0=PET_fh0.variables[PET_layer][idx0,rix[i]:rix[i+1],cix[j]:cix[j+1]]
        PET_block1=PET_fh1.variables[PET_layer][idx1,rix[i]:rix[i+1],cix[j]:cix[j+1]]
        ET_block0=ET_fh0.variables[ET_layer][idx0,rix[i]:rix[i+1],cix[j]:cix[j+1]]
        ET_block1=ET_fh1.variables[ET_layer][idx1,rix[i]:rix[i+1],cix[j]:cix[j+1]]        
        
        for k in range(P_block0.shape[1]):
            for l in range(P_block0.shape[2]):
                MASK=~P_block0[0,:,:].mask
                if MASK[k,l]:
                    # Concatenate
                    qs_ts=np.concatenate((Qs_block0[:,k,l],Qs_block1[:,k,l]))
                    qsm_ts=np.concatenate((Qsm_block0[:,k,l],Qsm_block1[:,k,l]))
                    qsb_ts=np.concatenate((Qsb_block0[:,k,l],Qsb_block1[:,k,l]))
                    pet_ts=np.concatenate((PET_block0[:,k,l],PET_block1[:,k,l]))
                    et_ts=np.concatenate((ET_block0[:,k,l],ET_block1[:,k,l]))
                    # Convert time series to mm/day
                    qs_ts=(qs_ts/1000)*(60*60*24*100*10)*-1
                    qsm_ts=(qsm_ts/1000)*(60*60*24*100*10)
                    qsb_ts=(qsb_ts/1000)*(60*60*24*100*10)*-1   
                    pet_ts=(pet_ts/1000)*(60*60*24*100*10)*-1
                    et_ts=(et_ts/1000)*(60*60*24*100*10)*-1
                    qs_mean[brix[k,l],bcix[k,l]]=np.nanmean(qs_ts)
                    qsm_mean[brix[k,l],bcix[k,l]]=np.nanmean(qsm_ts)
                    qsb_mean[brix[k,l],bcix[k,l]]=np.nanmean(qsb_ts)                    
                    r_ts=qs_ts+qsm_ts+qsb_ts
                    # Remove bad values
                    ridx=r_ts<0
                    o_len=len(r_ts)
                    r_ts[ridx]=0
                    n_len=o_len-np.sum(ridx)
                    r_comp[brix[k,l],bcix[k,l]]=n_len/o_len
                    # Do calculations
                    if n_len/o_len>0.50:
                        # Calculate mean
                        r_mean[brix[k,l],bcix[k,l]]=np.nanmean(r_ts)
                        r_std[brix[k,l],bcix[k,l]]=np.nanstd(r_ts)
                        pet_mean[brix[k,l],bcix[k,l]]=np.nanmean(pet_ts)
                        pet_std[brix[k,l],bcix[k,l]]=np.nanstd(pet_ts)
                        et_mean[brix[k,l],bcix[k,l]]=np.nanmean(et_ts)
                        et_std[brix[k,l],bcix[k,l]]=np.nanstd(et_ts)                        
                        # Survival function
                        [r_sort,r_freq]=survive(r_ts)
                        # At 1%
                        [rc1,rs1]=weibull_tail_fit(r_sort,r_freq,0.01)
                        r_c1[brix[k,l],bcix[k,l]]=rc1
                        r_s1[brix[k,l],bcix[k,l]]=rs1
                        # At 5%
                        [rc5,rs5]=weibull_tail_fit(r_sort,r_freq,0.05)
                        r_c5[brix[k,l],bcix[k,l]]=rc5
                        r_s5[brix[k,l],bcix[k,l]]=rs5                        
                        # Calculate num unique runoffs
                        r_nu[brix[k,l],bcix[k,l]]=len(np.unique(r_ts))
                        # Calculate return floods
                        r_2_5[brix[k,l],bcix[k,l]]=r_sort[np.argmin(np.abs(r_freq-(1/(2.5*365))))]
                        r_5[brix[k,l],bcix[k,l]]=r_sort[np.argmin(np.abs(r_freq-(1/(5*365))))]
                        r_10[brix[k,l],bcix[k,l]]=r_sort[np.argmin(np.abs(r_freq-(1/(10*365))))]
                    else:
                        r_mean[brix[k,l],bcix[k,l]]=np.nan
                        r_std[brix[k,l],bcix[k,l]]=np.nan
                        pet_mean[brix[k,l],bcix[k,l]]=np.nan
                        pet_std[brix[k,l],bcix[k,l]]=np.nan
                        et_mean[brix[k,l],bcix[k,l]]=np.nan
                        et_std[brix[k,l],bcix[k,l]]=np.nan 
                        r_c1[brix[k,l],bcix[k,l]]=np.nan
                        r_s1[brix[k,l],bcix[k,l]]=np.nan
                        r_c5[brix[k,l],bcix[k,l]]=np.nan
                        r_s5[brix[k,l],bcix[k,l]]=np.nan 
                        r_nu[brix[k,l],bcix[k,l]]=np.nan 
                        r_2_5[brix[k,l],bcix[k,l]]=np.nan
                        r_5[brix[k,l],bcix[k,l]]=np.nan 
                        r_10[brix[k,l],bcix[k,l]]=np.nan 
                else:
                    r_mean[brix[k,l],bcix[k,l]]=np.nan
                    r_std[brix[k,l],bcix[k,l]]=np.nan
                    pet_mean[brix[k,l],bcix[k,l]]=np.nan
                    pet_std[brix[k,l],bcix[k,l]]=np.nan
                    et_mean[brix[k,l],bcix[k,l]]=np.nan
                    et_std[brix[k,l],bcix[k,l]]=np.nan 
                    r_c1[brix[k,l],bcix[k,l]]=np.nan
                    r_s1[brix[k,l],bcix[k,l]]=np.nan 
                    r_c5[brix[k,l],bcix[k,l]]=np.nan
                    r_s5[brix[k,l],bcix[k,l]]=np.nan 
                    r_comp[brix[k,l],bcix[k,l]]=np.nan
                    r_nu[brix[k,l],bcix[k,l]]=np.nan 
                    r_2_5[brix[k,l],bcix[k,l]]=np.nan
                    r_5[brix[k,l],bcix[k,l]]=np.nan 
                    r_10[brix[k,l],bcix[k,l]]=np.nan 
                    qs_mean[brix[k,l],bcix[k,l]]=np.nan
                    qsm_mean[brix[k,l],bcix[k,l]]=np.nan
                    qsb_mean[brix[k,l],bcix[k,l]]=np.nan
                    
                    
        # Use full temp block for mean temps
        for k in range(T_block0.shape[1]):
            for l in range(T_block0.shape[2]):
                MASK=~T_block0[0,:,:].mask
                if MASK[k,l]:
                    # Concatenate
                    t_ts=np.concatenate((T_block0[:,k,l],T_block1[:,k,l]))
                    t_ts=t_ts-273.15
                    # Do calculations
                    t_mean[brix[k,l],bcix[k,l]]=np.nanmean(t_ts)
                    t_std[brix[k,l],bcix[k,l]]=np.nanstd(t_ts)
                else:
                    t_mean[brix[k,l],bcix[k,l]]=np.nan
                    t_std[brix[k,l],bcix[k,l]]=np.nan 
        stp_t=time.time()
        print('Completed block '+str(block_count)+' of '+str(total_blocks)+' in '+str(np.round((stp_t-st_t)/60,2))+' minutes')
        block_count=block_count+1

# Start writing out dataset
nc_out=Dataset('wrr2_derived_v2.nc',mode='w',format='NETCDF4_CLASSIC')

# Create dimensions
lat_dim=nc_out.createDimension('lat',P_fh0.variables['lat'].shape[0])
lon_dim=nc_out.createDimension('lon',P_fh0.variables['lon'].shape[0])
time_dim=nc_out.createDimension('time',1)

# Create title
nc_out.title='WRR2 Derived Values'

# Create variables
lat=nc_out.createVariable('lat',np.float32,('lat',))
lat.units='degrees_north'
lat.long_name='latitude'

lon=nc_out.createVariable('lon',np.float32,('lon',))
lon.units='degrees_east'
lon.long_name='longitude'

pmean=nc_out.createVariable('pmean',np.float64,('lat','lon'))
pmean.units='mm'
pmean.long_name='mean daily precipitation'
pmean.source='WaterGAP3'
pmean.duration='Jan 1 1980 - Dec 31 1999'

pstd=nc_out.createVariable('pstd',np.float64,('lat','lon'))
pstd.units='mm'
pstd.long_name='stdev daily precipitation'
pstd.source='WaterGAP3'
pstd.duration='Jan 1 1980 - Dec 31 1999'

pshp=nc_out.createVariable('pshp',np.float64,('lat','lon'))
pshp.units='None'
pshp.long_name='precipitation weibull shape 1 percent threshold'
pshp.source='WaterGAP3'
pshp.duration='Jan 1 1980 - Dec 31 1999'

pscl=nc_out.createVariable('pscl',np.float64,('lat','lon'))
pscl.units='None'
pscl.long_name='precipitation weibull scale 1 percent threshold'
pscl.source='WaterGAP3'
pscl.duration='Jan 1 1980 - Dec 31 1999'

pcom=nc_out.createVariable('pcom',np.float64,('lat','lon'))
pcom.units='Percent'
pcom.long_name='precipitation completeness'
pcom.source='WaterGAP3'
pcom.duration='Jan 1 1980 - Dec 31 1999'

ps=nc_out.createVariable('ps',np.float64,('lat','lon'))
ps.units='Percent'
ps.long_name='percent precipitation as snow'
ps.source='WaterGAP3 and SuflexTrip'
ps.duration='Jan 1 1980 - Dec 31 1999'

qs=nc_out.createVariable('qs',np.float64,('lat','lon'))
qs.units='mm'
qs.long_name='mean qs'
qs.source='WaterGAP3'
qs.duration='Jan 1 1980 - Dec 31 1999'

qsm=nc_out.createVariable('qsm',np.float64,('lat','lon'))
qsm.units='mm'
qsm.long_name='mean qsm'
qsm.source='WaterGAP3'
qsm.duration='Jan 1 1980 - Dec 31 1999'

qsb=nc_out.createVariable('qsb',np.float64,('lat','lon'))
qsb.units='mm'
qsb.long_name='mean qsb'
qsb.source='WaterGAP3'
qsb.duration='Jan 1 1980 - Dec 31 1999'

rmean=nc_out.createVariable('rmean',np.float64,('lat','lon'))
rmean.units='mm'
rmean.long_name='mean daily total runoff'
rmean.source='WaterGAP3'
rmean.duration='Jan 1 1980 - Dec 31 1999'

rstd=nc_out.createVariable('rstd',np.float64,('lat','lon'))
rstd.units='mm'
rstd.long_name='stdev daily total runoff'
rstd.source='WaterGAP3'
rstd.duration='Jan 1 1980 - Dec 31 1999'

rshp1=nc_out.createVariable('rshp1',np.float64,('lat','lon'))
rshp1.units='None'
rshp1.long_name='total runoff weibull shape 1 percent threshold'
rshp1.source='WaterGAP3'
rshp1.duration='Jan 1 1980 - Dec 31 1999'

rscl1=nc_out.createVariable('rscl1',np.float64,('lat','lon'))
rscl1.units='None'
rscl1.long_name='total runoff weibull scale 1 percent threshold'
rscl1.source='WaterGAP3'
rscl1.duration='Jan 1 1980 - Dec 31 1999'

rshp5=nc_out.createVariable('rshp5',np.float64,('lat','lon'))
rshp5.units='None'
rshp5.long_name='total runoff weibull shape 5 percent threshold'
rshp5.source='WaterGAP3'
rshp5.duration='Jan 1 1980 - Dec 31 1999'

rscl5=nc_out.createVariable('rscl5',np.float64,('lat','lon'))
rscl5.units='None'
rscl5.long_name='total runoff weibull scale 5 percent threshold'
rscl5.source='WaterGAP3'
rscl5.duration='Jan 1 1980 - Dec 31 1999'

rnu=nc_out.createVariable('rnu',np.float64,('lat','lon'))
rnu.units='None'
rnu.long_name='number of unique runoff days'
rnu.source='WaterGAP3'
rnu.duration='Jan 1 1980 - Dec 31 1999'

r25=nc_out.createVariable('r25',np.float64,('lat','lon'))
r25.units='mm'
r25.long_name='magnitude of 2.5 year flood'
r25.source='WaterGAP3'
r25.duration='Jan 1 1980 - Dec 31 1999'

r5=nc_out.createVariable('r5',np.float64,('lat','lon'))
r5.units='mm'
r5.long_name='magnitude of 5 year flood'
r5.source='WaterGAP3'
r5.duration='Jan 1 1980 - Dec 31 1999'

r10=nc_out.createVariable('r10',np.float64,('lat','lon'))
r10.units='mm'
r10.long_name='magnitude of 10 year flood'
r10.source='WaterGAP3'
r10.duration='Jan 1 1980 - Dec 31 1999'

rcom=nc_out.createVariable('rcom',np.float64,('lat','lon'))
rcom.units='Percent'
rcom.long_name='total runoff completeness'
rcom.source='WaterGAP3'
rcom.duration='Jan 1 1980 - Dec 31 1999'

etmean=nc_out.createVariable('etmean',np.float64,('lat','lon'))
etmean.units='mm'
etmean.long_name='mean daily et'
etmean.source='WaterGAP3'
etmean.duration='Jan 1 1980 - Dec 31 1999'

etstd=nc_out.createVariable('etstd',np.float64,('lat','lon'))
etstd.units='mm'
etstd.long_name='stdev daily et'
etstd.source='WaterGAP3'
etstd.duration='Jan 1 1980 - Dec 31 1999'

petmean=nc_out.createVariable('petmean',np.float64,('lat','lon'))
petmean.units='mm'
petmean.long_name='mean daily pet'
petmean.source='WaterGAP3'
petmean.duration='Jan 1 1980 - Dec 31 1999'

petstd=nc_out.createVariable('petstd',np.float64,('lat','lon'))
petstd.units='mm'
petstd.long_name='stdev daily pet'
petstd.source='WaterGAP3'
petstd.duration='Jan 1 1980 - Dec 31 1999'

tmean=nc_out.createVariable('tmean',np.float64,('lat','lon'))
tmean.units='C'
tmean.long_name='mean daily temperature'
tmean.source='SuflexTrip'
tmean.duration='Jan 1 1980 - Dec 31 1999'

tstd=nc_out.createVariable('tstd',np.float64,('lat','lon'))
tstd.units='C'
tstd.long_name='standard deviation daily temperature'
tstd.source='SuflexTrip'
tstd.duration='Jan 1 1980 - Dec 31 1999'

# Write variables
lat[:]=P_fh0.variables['lat'][:]
lon[:]=P_fh0.variables['lon'][:]

# WaterGAP3 Mask
MASK=P_fh0.variables['Precip'][0,:,:].mask
pmean[:,:]=np.ma.array(p_mean,mask=MASK)
pstd[:,:]=np.ma.array(p_std,mask=MASK)
pcom[:,:]=np.ma.array(p_comp,mask=MASK)
pshp[:,:]=np.ma.array(p_c,mask=MASK)
pscl[:,:]=np.ma.array(p_s,mask=MASK)
ps[:,:]=np.ma.array(perc_s,mask=MASK)
rmean[:,:]=np.ma.array(r_mean,mask=MASK)
rstd[:,:]=np.ma.array(r_std,mask=MASK)
rcom[:,:]=np.ma.array(r_comp,mask=MASK)
rshp1[:,:]=np.ma.array(r_c1,mask=MASK)
rscl1[:,:]=np.ma.array(r_s1,mask=MASK)
rshp5[:,:]=np.ma.array(r_c5,mask=MASK)
rscl5[:,:]=np.ma.array(r_s5,mask=MASK)
rnu[:,:]=np.ma.array(r_nu,mask=MASK)
r25[:,:]=np.ma.array(r_2_5,mask=MASK)
r5[:,:]=np.ma.array(r_5,mask=MASK)
r10[:,:]=np.ma.array(r_10,mask=MASK)

qs[:,:]=np.ma.array(qs_mean,mask=MASK)
qsm[:,:]=np.ma.array(qsm_mean,mask=MASK)
qsb[:,:]=np.ma.array(qsb_mean,mask=MASK)

etmean[:,:]=np.ma.array(et_mean,mask=MASK)
etstd[:,:]=np.ma.array(et_std,mask=MASK)
petmean[:,:]=np.ma.array(pet_mean,mask=MASK)
petstd[:,:]=np.ma.array(pet_std,mask=MASK)

# SuflexTrip Mask
# MASK=T_fh0.variables['AvgSurfT'][0,:,:].mask
tmean[:,:]=np.ma.array(t_mean,mask=MASK)
tstd[:,:]=np.ma.array(t_std,mask=MASK)

# Close and write
nc_out.close()

# Close other datasets
T_fh0.close()
T_fh1.close()
P_fh0.close()
P_fh1.close()
Qs_fh0.close()
Qs_fh1.close()
Qsm_fh0.close()
Qsm_fh1.close()
Qsb_fh0.close()
Qsb_fh1.close()
ET_fh0.close()
ET_fh1.close()
PET_fh0.close()
PET_fh1.close()




