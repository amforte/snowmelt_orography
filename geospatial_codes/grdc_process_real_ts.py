#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  9 12:20:40 2022

@author: aforte
"""

import pandas as pd
import numpy as np

def survive(ts):
    ts_sort=np.sort(ts)
    tsn=len(ts_sort)
    tsrank=np.arange(1,tsn+1,1)
    ts_freq_excd=(tsn+1-tsrank)/tsn
    return ts_sort,ts_freq_excd

def weibull_tail_fit(x,y,thresh):
    ix=np.nonzero(y<thresh)[0][:1][0]
    xtrim=x[ix:]
    ytrim=y[ix:]
    xts=np.log(xtrim)
    yts=np.log(-np.log(ytrim))      
    [lin,r,rnk,sng,V]=np.polyfit(xts,yts,1,full=True)
    c=lin[0]
    s=np.exp(-1*lin[1]/c)
    return c,s 


master_location='/Users/aforte/Documents/Python/snowmelt/'

# Read GRDC stations stats
df=pd.read_csv('grdc_station_stats.csv')
df=df.drop(df.index[np.isnan(df['MIN_Z'])])
df=df.reset_index(drop=True)

# Count stations
num_stat=len(df)

# Set temporal ranges
start_=pd.to_datetime('1980-01-01')
stop_=pd.to_datetime('1999-12-31')
data_dr=pd.date_range(start_,stop_)

# Generate arrays
grdc_id=np.zeros(num_stat).astype('int')

q=np.zeros(num_stat)
rr=np.zeros(num_stat)
hr=np.zeros(num_stat)
fc=np.zeros(num_stat)
rr_c1=np.zeros(num_stat)
hr_c1=np.zeros(num_stat)
rr_s1=np.zeros(num_stat)
hr_s1=np.zeros(num_stat)
rr_c5=np.zeros(num_stat)
hr_c5=np.zeros(num_stat)
rr_s5=np.zeros(num_stat)
hr_s5=np.zeros(num_stat)
rr_2_5=np.zeros(num_stat)
hr_2_5=np.zeros(num_stat)
rr_5=np.zeros(num_stat)
hr_5=np.zeros(num_stat)
rr_10=np.zeros(num_stat)
hr_10=np.zeros(num_stat)

sq=np.zeros(num_stat)
srr=np.zeros(num_stat)
shr=np.zeros(num_stat)
sc=np.zeros(num_stat)
srr_c1=np.zeros(num_stat)
shr_c1=np.zeros(num_stat)
srr_s1=np.zeros(num_stat)
shr_s1=np.zeros(num_stat)
srr_c5=np.zeros(num_stat)
shr_c5=np.zeros(num_stat)
srr_s5=np.zeros(num_stat)
shr_s5=np.zeros(num_stat)
srr_2_5=np.zeros(num_stat)
shr_2_5=np.zeros(num_stat)
srr_5=np.zeros(num_stat)
shr_5=np.zeros(num_stat)
srr_10=np.zeros(num_stat)
shr_10=np.zeros(num_stat)


for i in range(num_stat):
    if np.mod(i,100)==0:
        print(i)
    
    grdc_id[i]=df.loc[i,'GRDC_NO']
    
    # Extract areas
    r_area=df.loc[i,'R_AREA']
    h_area=df.loc[i,'H_AREA']
    
    # Check areas for missing reported areas
    if r_area<0:
        r_area=h_area

        
    try:
        # Read in GRDC timeseries
        fn=master_location+'grdc_ts/'+str(df.loc[i,'GRDC_NO'])+'_Q_Day.Cmd.txt'
        dfoi=pd.read_table(fn,header=36,sep=';',engine='python',encoding='unicode_escape',parse_dates=[0])
        
        # Drop missing values
        idx=dfoi[' Value']<0
        dfoi=dfoi.drop(dfoi.index[idx])
        dfoi=dfoi.reset_index(drop=True)
        
        # Determine true length and completeness
        ts_dr=pd.date_range(dfoi.iloc[0,0],dfoi.iloc[-1,0])
        fc[i]=len(dfoi)/len(ts_dr)
        
        # Calculate runoffs
        dfoi['R_R']=(dfoi[' Value']/(r_area*1000*1000))*(100*10*60*60*24) 
        dfoi['H_R']=(dfoi[' Value']/(h_area*1000*1000))*(100*10*60*60*24)
        
        # Calculate means
        q[i]=np.nanmean(dfoi[' Value'])
        rr[i]=np.nanmean(dfoi['R_R'])
        hr[i]=np.nanmean(dfoi['H_R'])
        
        # Calculate shape and scales
        rr_v=dfoi['R_R'].to_numpy()
        rr_v=rr_v[~np.isnan(rr_v)]
        if len(np.nonzero(rr_v)[0])>1000:
            [rr_sr,rr_f]=survive(rr_v)
            [rr_c1[i],rr_s1[i]]=weibull_tail_fit(rr_sr,rr_f,0.01)
            [rr_c5[i],rr_s5[i]]=weibull_tail_fit(rr_sr,rr_f,0.05)
            rr_2_5[i]=rr_sr[np.argmin(np.abs(rr_f-(1/(2.5*365))))]
            rr_5[i]=rr_sr[np.argmin(np.abs(rr_f-(1/(5*365))))]
            rr_10[i]=rr_sr[np.argmin(np.abs(rr_f-(1/(10*365))))]
        else:
            rr_c1[i]=np.nan 
            rr_s1[i]=np.nan
            rr_c5[i]=np.nan 
            rr_s5[i]=np.nan
            rr_2_5[i]=np.nan
            rr_5[i]=np.nan
            rr_10[i]=np.nan       
        
        hr_v=dfoi['H_R'].to_numpy()
        hr_v=hr_v[~np.isnan(hr_v)]
        if len(np.nonzero(hr_v)[0])>1000:
            [hr_sr,hr_f]=survive(hr_v)
            [hr_c1[i],hr_s1[i]]=weibull_tail_fit(hr_sr,hr_f,0.01)
            [hr_c5[i],hr_s5[i]]=weibull_tail_fit(hr_sr,hr_f,0.05)
            hr_2_5[i]=hr_sr[np.argmin(np.abs(hr_f-(1/(2.5*365))))]
            hr_5[i]=hr_sr[np.argmin(np.abs(hr_f-(1/(5*365))))]
            hr_10[i]=hr_sr[np.argmin(np.abs(hr_f-(1/(10*365))))]
        else:
            hr_c1[i]=np.nan 
            hr_s1[i]=np.nan
            hr_c5[i]=np.nan 
            hr_s5[i]=np.nan 
            hr_2_5[i]=np.nan
            hr_5[i]=np.nan
            hr_10[i]=np.nan
            
        # Build temporal index and create sliced time series
        idx=(dfoi['YYYY-MM-DD']>=start_) & (dfoi['YYYY-MM-DD']<=stop_)
        dfois=dfoi.loc[idx,:]
        dfois=dfois.reset_index(drop=True)
        
        # Calculate completion of sliced time series
        sc[i]=len(dfois)/len(data_dr)
        
        # Calculate means    
        if len(dfois)>0:
            sq[i]=np.nanmean(dfois[' Value'])
            srr[i]=np.nanmean(dfois['R_R'])
            shr[i]=np.nanmean(dfois['H_R'])
        else:
            sq[i]=np.nan
            srr[i]=np.nan
            shr[i]=np.nan
        
        
        # Calculate shape and scales
        rr_v=dfois['R_R'].to_numpy()
        rr_v=rr_v[~np.isnan(rr_v)]
        if len(np.nonzero(rr_v)[0])>1000:
            [rr_sr,rr_f]=survive(rr_v)
            [srr_c1[i],srr_s1[i]]=weibull_tail_fit(rr_sr,rr_f,0.01)
            [srr_c5[i],srr_s5[i]]=weibull_tail_fit(rr_sr,rr_f,0.05)
            srr_2_5[i]=rr_sr[np.argmin(np.abs(rr_f-(1/(2.5*365))))]
            srr_5[i]=rr_sr[np.argmin(np.abs(rr_f-(1/(5*365))))]
            srr_10[i]=rr_sr[np.argmin(np.abs(rr_f-(1/(10*365))))]
        else:
            srr_c1[i]=np.nan 
            srr_s1[i]=np.nan
            srr_c5[i]=np.nan 
            srr_s5[i]=np.nan 
            srr_2_5[i]=np.nan
            srr_5[i]=np.nan
            srr_10[i]=np.nan 
        
        hr_v=dfois['H_R'].to_numpy()
        hr_v=hr_v[~np.isnan(hr_v)]
        if len(np.nonzero(hr_v)[0])>1000:
            [hr_sr,hr_f]=survive(hr_v)
            [shr_c1[i],shr_s1[i]]=weibull_tail_fit(hr_sr,hr_f,0.01)
            [shr_c5[i],shr_s5[i]]=weibull_tail_fit(hr_sr,hr_f,0.05)
            shr_2_5[i]=hr_sr[np.argmin(np.abs(hr_f-(1/(2.5*365))))]
            shr_5[i]=hr_sr[np.argmin(np.abs(hr_f-(1/(5*365))))]
            shr_10[i]=hr_sr[np.argmin(np.abs(hr_f-(1/(10*365))))]
        else:
            shr_c1[i]=np.nan 
            shr_s1[i]=np.nan
            shr_c5[i]=np.nan 
            shr_s5[i]=np.nan 
            shr_2_5[i]=np.nan
            shr_5[i]=np.nan
            shr_10[i]=np.nan
    except:
        q[i]=np.nan
        rr[i]=np.nan
        hr[i]=np.nan
        fc[i]=np.nan
        rr_c1[i]=np.nan
        hr_c1[i]=np.nan
        rr_s1[i]=np.nan
        hr_s1[i]=np.nan
        rr_c5[i]=np.nan
        hr_c5[i]=np.nan
        rr_s5[i]=np.nan
        hr_s5[i]=np.nan
        rr_2_5[i]=np.nan
        hr_2_5[i]=np.nan
        rr_5[i]=np.nan
        hr_5[i]=np.nan
        rr_10[i]=np.nan
        hr_10[i]=np.nan

        sq[i]=np.nan
        srr[i]=np.nan
        shr[i]=np.nan
        sc[i]=np.nan
        srr_c1[i]=np.nan
        shr_c1[i]=np.nan
        srr_s1[i]=np.nan
        shr_s1[i]=np.nan
        srr_c5[i]=np.nan
        shr_c5[i]=np.nan
        srr_s5[i]=np.nan
        shr_s5[i]=np.nan
        srr_2_5[i]=np.nan
        shr_2_5[i]=np.nan
        srr_5[i]=np.nan
        shr_5[i]=np.nan
        srr_10[i]=np.nan
        shr_10[i]=np.nan

        
    
# Compile and write out    
df_out=pd.DataFrame({'GRDC_NO':grdc_id,
                     'full_completeness':fc,
                     'full_discharge':q,
                     'full_reported_runoff':rr,
                     'full_hydrosheds_runoff':hr,
                     'full_reported_c1':rr_c1,
                     'full_reported_s1':rr_s1,
                     'full_reported_c5':rr_c5,
                     'full_reported_s5':rr_s5,
                     'full_hydrosheds_c1':hr_c1,
                     'full_hydrosheds_s1':hr_s1,
                     'full_hydrosheds_c5':hr_c5,
                     'full_hydrosheds_s5':hr_s5,
                     'full_reported_2_5flood':rr_2_5,
                     'full_reported_5flood':rr_5,
                     'full_reported_10flood':rr_10,
                     'full_hydrosheds_2_5flood':hr_2_5,
                     'full_hydrosheds_5flood':hr_5,
                     'full_hydrosheds_10flood':hr_10,
                     'sliced_completeness':sc,
                     'sliced_discharge':sq,
                     'sliced_reported_runoff':srr,
                     'sliced_hydrosheds_runoff':shr,
                     'sliced_reported_c1':srr_c1,
                     'sliced_reported_s1':srr_s1,
                     'sliced_reported_c5':srr_c5,
                     'sliced_reported_s5':srr_s5,
                     'sliced_hydrosheds_c1':shr_c1,
                     'sliced_hydrosheds_s1':shr_s1,
                     'sliced_hydrosheds_c5':shr_c5,
                     'sliced_hydrosheds_s5':shr_s5,
                     'sliced_reported_2_5flood':srr_2_5,
                     'sliced_reported_5flood':srr_5,
                     'sliced_reported_10flood':srr_10,
                     'sliced_hydrosheds_2_5flood':shr_2_5,
                     'sliced_hydrosheds_5flood':shr_5,
                     'sliced_hydrosheds_10flood':shr_10,})    
df_out.to_csv('grdc_real_ts_v2.csv',index=False)