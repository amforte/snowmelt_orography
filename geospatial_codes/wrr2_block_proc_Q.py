#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 11 08:24:11 2022

@author: aforte
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

# Precipitation
Q_layer='RivOut'
Q_mfile=[master_location+'e2o_univk_wrr2/e2o_univk_wrr2_glob15_day_RivOut_1980-1989.nc',
         master_location+'e2o_univk_wrr2/e2o_univk_wrr2_glob15_day_RivOut_1990-1999.nc']
# Surface runoff
SWE_layer='SWE'
SWE_mfile=[master_location+'e2o_univk_wrr2/e2o_univk_wrr2_glob15_day_SWE_1980-1989.nc',
          master_location+'e2o_univk_wrr2/e2o_univk_wrr2_glob15_day_SWE_1990-1999.nc']

## Open Datasets
Q_fh0=Dataset(Q_mfile[0])
Q_fh1=Dataset(Q_mfile[1])
SWE_fh0=Dataset(SWE_mfile[0])
SWE_fh1=Dataset(SWE_mfile[1])


# Create indices to not select leap days
Q0=Q_fh0.variables['time'][:]
Q0obj=Q_fh0.variables['time']
dt0=num2date(Q0obj[:],Q0obj.units)
DT0=pd.DataFrame(data={'date':dt0})
DT0['date']=pd.to_datetime(DT0['date'].astype('str'))
day=pd.DatetimeIndex(DT0['date']).day
mnth=pd.DatetimeIndex(DT0['date']).month
DT0=DT0.drop(index=DT0.index[np.logical_and(day==29,mnth==2)])
idx0=DT0.index

Q1=Q_fh1.variables['time'][:]
Q1obj=Q_fh1.variables['time']
dt1=num2date(Q1obj[:],Q1obj.units)
DT1=pd.DataFrame(data={'date':dt1})
DT1['date']=pd.to_datetime(DT1['date'].astype('str'))
day=pd.DatetimeIndex(DT1['date']).day
mnth=pd.DatetimeIndex(DT1['date']).month
DT1=DT1.drop(index=DT1.index[np.logical_and(day==29,mnth==2)])
idx1=DT1.index


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
q_mean=np.zeros((r,c))
q_std=np.zeros((r,c))
q_c1=np.zeros((r,c))
q_s1=np.zeros((r,c))
q_c5=np.zeros((r,c))
q_s5=np.zeros((r,c))
q_comp=np.zeros((r,c))
q_nu=np.zeros((r,c))
q_2_5=np.zeros((r,c))
q_5=np.zeros((r,c))
q_10=np.zeros((r,c))

swe_mean=np.zeros((r,c))
swe_std=np.zeros((r,c))


# Start master nested loop
for i in range(blocks):
    for j in range(blocks):
        st_t=time.time()
        # Block index
        bcix=CIX[rix[i]:rix[i+1],cix[j]:cix[j+1]]
        brix=RIX[rix[i]:rix[i+1],cix[j]:cix[j+1]]
        
        # Load precip blocks
        Q_block0=Q_fh0.variables[Q_layer][idx0,rix[i]:rix[i+1],cix[j]:cix[j+1]]
        Q_block1=Q_fh1.variables[Q_layer][idx1,rix[i]:rix[i+1],cix[j]:cix[j+1]]
        # Load full temp blocks
        SWE_block0=SWE_fh0.variables[SWE_layer][:,rix[i]:rix[i+1],cix[j]:cix[j+1]]
        SWE_block1=SWE_fh1.variables[SWE_layer][:,rix[i]:rix[i+1],cix[j]:cix[j+1]]       

        for k in range(Q_block0.shape[1]):
            for l in range(Q_block0.shape[2]):
                MASK=~Q_block0[0,:,:].mask
                if MASK[k,l]:
                    # Concatenate
                    q_ts=np.concatenate((Q_block0[:,k,l],Q_block1[:,k,l]))
                    swe_ts=np.concatenate((SWE_block0[:,k,l],SWE_block1[:,k,l]))
                    
                    # Convert swe to mm
                    swe_ts=(swe_ts/1000)*(100*10)


                    # Remove bad values
                    qidx=q_ts<0
                    o_len=len(q_ts)
                    q_ts[qidx]=0
                    n_len=o_len-np.sum(qidx)
                    q_comp[brix[k,l],bcix[k,l]]=n_len/o_len
                    # Do calculations
                    if n_len/o_len>0.50:
                        # Calculate mean discharge in m3/sec
                        q_mean[brix[k,l],bcix[k,l]]=np.nanmean(q_ts)
                        q_std[brix[k,l],bcix[k,l]]=np.nanstd(q_ts)
                        # Calculate mean snowmelt water equiavlent
                        swe_mean[brix[k,l],bcix[k,l]]=np.nanmean(swe_ts)
                        swe_std[brix[k,l],bcix[k,l]]=np.nanstd(swe_ts)
                        # Survival function
                        [q_sort,q_freq]=survive(q_ts)
                        # At 1%
                        [qc1,qs1]=weibull_tail_fit(q_sort,q_freq,0.01)
                        q_c1[brix[k,l],bcix[k,l]]=qc1
                        q_s1[brix[k,l],bcix[k,l]]=qs1
                        # At 5%
                        [qc5,qs5]=weibull_tail_fit(q_sort,q_freq,0.05)
                        q_c5[brix[k,l],bcix[k,l]]=qc5
                        q_s5[brix[k,l],bcix[k,l]]=qs5                        
                        # Calculate num unique runoffs
                        q_nu[brix[k,l],bcix[k,l]]=len(np.unique(q_ts))
                        # Calculate return floods
                        q_2_5[brix[k,l],bcix[k,l]]=q_sort[np.argmin(np.abs(q_freq-(1/(2.5*365))))]
                        q_5[brix[k,l],bcix[k,l]]=q_sort[np.argmin(np.abs(q_freq-(1/(5*365))))]
                        q_10[brix[k,l],bcix[k,l]]=q_sort[np.argmin(np.abs(q_freq-(1/(10*365))))]
                    else:
                        q_mean[brix[k,l],bcix[k,l]]=np.nan
                        q_std[brix[k,l],bcix[k,l]]=np.nan
                        swe_mean[brix[k,l],bcix[k,l]]=np.nan
                        swe_std[brix[k,l],bcix[k,l]]=np.nan
                        q_c1[brix[k,l],bcix[k,l]]=np.nan
                        q_s1[brix[k,l],bcix[k,l]]=np.nan
                        q_c5[brix[k,l],bcix[k,l]]=np.nan
                        q_s5[brix[k,l],bcix[k,l]]=np.nan 
                        q_nu[brix[k,l],bcix[k,l]]=np.nan 
                        q_2_5[brix[k,l],bcix[k,l]]=np.nan
                        q_5[brix[k,l],bcix[k,l]]=np.nan 
                        q_10[brix[k,l],bcix[k,l]]=np.nan 
                else:
                    q_mean[brix[k,l],bcix[k,l]]=np.nan
                    q_std[brix[k,l],bcix[k,l]]=np.nan
                    swe_mean[brix[k,l],bcix[k,l]]=np.nan
                    swe_std[brix[k,l],bcix[k,l]]=np.nan
                    q_c1[brix[k,l],bcix[k,l]]=np.nan
                    q_s1[brix[k,l],bcix[k,l]]=np.nan
                    q_c5[brix[k,l],bcix[k,l]]=np.nan
                    q_s5[brix[k,l],bcix[k,l]]=np.nan 
                    q_nu[brix[k,l],bcix[k,l]]=np.nan 
                    q_2_5[brix[k,l],bcix[k,l]]=np.nan
                    q_5[brix[k,l],bcix[k,l]]=np.nan 
                    q_10[brix[k,l],bcix[k,l]]=np.nan

        stp_t=time.time()
        print('Completed block '+str(block_count)+' of '+str(total_blocks)+' in '+str(np.round((stp_t-st_t)/60,2))+' minutes')
        block_count=block_count+1

# Start writing out dataset
nc_out=Dataset('wrr2_derived_Q.nc',mode='w',format='NETCDF4_CLASSIC')

# Create dimensions
lat_dim=nc_out.createDimension('lat',Q_fh0.variables['lat'].shape[0])
lon_dim=nc_out.createDimension('lon',Q_fh0.variables['lon'].shape[0])
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

qmean=nc_out.createVariable('qmean',np.float64,('lat','lon'))
qmean.units='m3/sec'
qmean.long_name='mean daily discharge'
qmean.source='WaterGAP3'
qmean.duration='Jan 1 1980 - Dec 31 1999'

qstd=nc_out.createVariable('qstd',np.float64,('lat','lon'))
qstd.units='m3/sec'
qstd.long_name='stdev daily discharge'
qstd.source='WaterGAP3'
qstd.duration='Jan 1 1980 - Dec 31 1999'


qshp1=nc_out.createVariable('qshp1',np.float64,('lat','lon'))
qshp1.units='None'
qshp1.long_name='discharge weibull shape 1 percent threshold'
qshp1.source='WaterGAP3'
qshp1.duration='Jan 1 1980 - Dec 31 1999'

qscl1=nc_out.createVariable('qscl1',np.float64,('lat','lon'))
qscl1.units='None'
qscl1.long_name='discharge weibull scale 1 percent threshold'
qscl1.source='WaterGAP3'
qscl1.duration='Jan 1 1980 - Dec 31 1999'

qshp5=nc_out.createVariable('qshp5',np.float64,('lat','lon'))
qshp5.units='None'
qshp5.long_name='discharge weibull shape 5 percent threshold'
qshp5.source='WaterGAP3'
qshp5.duration='Jan 1 1980 - Dec 31 1999'

qscl5=nc_out.createVariable('qscl5',np.float64,('lat','lon'))
qscl5.units='None'
qscl5.long_name='discharge weibull scale 5 percent threshold'
qscl5.source='WaterGAP3'
qscl5.duration='Jan 1 1980 - Dec 31 1999'

qnu=nc_out.createVariable('qnu',np.float64,('lat','lon'))
qnu.units='None'
qnu.long_name='number of unique discharge days'
qnu.source='WaterGAP3'
qnu.duration='Jan 1 1980 - Dec 31 1999'

q25=nc_out.createVariable('q25',np.float64,('lat','lon'))
q25.units='m3/sec'
q25.long_name='magnitude of 2.5 year flood'
q25.source='WaterGAP3'
q25.duration='Jan 1 1980 - Dec 31 1999'

q5=nc_out.createVariable('q5',np.float64,('lat','lon'))
q5.units='m3/sec'
q5.long_name='magnitude of 5 year flood'
q5.source='WaterGAP3'
q5.duration='Jan 1 1980 - Dec 31 1999'

q10=nc_out.createVariable('q10',np.float64,('lat','lon'))
q10.units='m3/sec'
q10.long_name='magnitude of 10 year flood'
q10.source='WaterGAP3'
q10.duration='Jan 1 1980 - Dec 31 1999'

qcom=nc_out.createVariable('qcom',np.float64,('lat','lon'))
qcom.units='Percent'
qcom.long_name='total discharge completeness'
qcom.source='WaterGAP3'
qcom.duration='Jan 1 1980 - Dec 31 1999'

swemean=nc_out.createVariable('swemean',np.float64,('lat','lon'))
swemean.units='mm'
swemean.long_name='mean snowmelt water equivalent'
swemean.source='WaterGAP3'
swemean.duration='Jan 1 1980 - Dec 31 1999'

swestd=nc_out.createVariable('swestd',np.float64,('lat','lon'))
swestd.units='mm'
swestd.long_name='stdev snowmelt water equivalent'
swestd.source='WaterGAP3'
swestd.duration='Jan 1 1980 - Dec 31 1999'


# Write variables
lat[:]=Q_fh0.variables['lat'][:]
lon[:]=Q_fh0.variables['lon'][:]

# WaterGAP3 Mask
MASK=Q_fh0.variables['RivOut'][0,:,:].mask
qmean[:,:]=np.ma.array(q_mean,mask=MASK)
qstd[:,:]=np.ma.array(q_std,mask=MASK)
qcom[:,:]=np.ma.array(q_comp,mask=MASK)
qshp1[:,:]=np.ma.array(q_c1,mask=MASK)
qscl1[:,:]=np.ma.array(q_s1,mask=MASK)
qshp5[:,:]=np.ma.array(q_c5,mask=MASK)
qscl5[:,:]=np.ma.array(q_s5,mask=MASK)
qnu[:,:]=np.ma.array(q_nu,mask=MASK)
q25[:,:]=np.ma.array(q_2_5,mask=MASK)
q5[:,:]=np.ma.array(q_5,mask=MASK)
q10[:,:]=np.ma.array(q_10,mask=MASK)

swemean[:,:]=np.ma.array(swe_mean,mask=MASK)
swestd[:,:]=np.ma.array(swe_std,mask=MASK)

# Close and write
nc_out.close()

# Close other datasets
Q_fh0.close()
Q_fh1.close()
SWE_fh0.close()
SWE_fh1.close()

