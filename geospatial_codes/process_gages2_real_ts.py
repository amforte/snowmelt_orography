#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  6 20:08:49 2022

@author: aforte
"""
import pandas as pd
import numpy as np
import glob

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

# Set master location
master_location='/Volumes/Choruh/Data/snowmelt_project/'
gages2=pd.read_csv('gages2_wrr2_raster_values.csv')
# Build file list
fL=glob.glob(master_location+'gagesII/gagesII_ts/*.txt')
# Number files
num_files=len(fL)

# Generate output arrays
stat=np.zeros(num_files)
fc=np.zeros(num_files)
flen=np.zeros(num_files)
fmr=np.zeros(num_files)
frc1=np.zeros(num_files)
frs1=np.zeros(num_files)
frc5=np.zeros(num_files)
frs5=np.zeros(num_files)
fr2_5=np.zeros(num_files)
fr5=np.zeros(num_files)
fr10=np.zeros(num_files)

fstartY=np.zeros(num_files)
fstartM=np.zeros(num_files)
fstartD=np.zeros(num_files)
fstopY=np.zeros(num_files)
fstopM=np.zeros(num_files)
fstopD=np.zeros(num_files)

sc=np.zeros(num_files)
slen=np.zeros(num_files)
smr=np.zeros(num_files)
src1=np.zeros(num_files)
srs1=np.zeros(num_files)
src5=np.zeros(num_files)
srs5=np.zeros(num_files)
sr2_5=np.zeros(num_files)
sr5=np.zeros(num_files)
sr10=np.zeros(num_files)

sstartY=np.zeros(num_files)
sstartM=np.zeros(num_files)
sstartD=np.zeros(num_files)
sstopY=np.zeros(num_files)
sstopM=np.zeros(num_files)
sstopD=np.zeros(num_files)

# Set temporal ranges
start_=pd.to_datetime('1980-01-01')
stop_=pd.to_datetime('1999-12-31')
data_dr=pd.date_range(start_,stop_)

# What follows is a kludgey and hacked together chain of conditionals
# because the USGS is an agent of chaos and created perhaps one of the
# most awful datsets to try to automate processing of, ever.

## NEED TO TRACK START AND END DATES OF BOTH

for i in range(num_files):
    if np.mod(i,100)==0:
        print(i)
    # Load
    df=pd.read_table(fL[i],comment='#',sep='\t')
    # Drop first line because skiprows doesn't seem to work with the comment command
    df=df.drop([0])
    df=df.reset_index(drop=True)
    
    # Control for empty datatables
    if len(df)>1:
        # Convert date to datetime format
        df[['datetime']]=df[['datetime']].apply(pd.to_datetime)
        # Extract the discharge component and it's quality indicators
        qdf=df.filter(like='_00060_00003')
        
        # Filter for missing mean value entirely
        if len(qdf.columns)>0:
            # Do initial filtering based on whether it includes any of the magic
            # codes that seem to be undocummented except in indvidual files
            # because the person who coded this data hated life itself
            # ix=qdf.index[(qdf.iloc[:,0]=='Ice') | (qdf.iloc[:,0]=='Eqp') | (qdf.iloc[:,0]=='Dis') | (qdf.iloc[:,0]=='Rat')]
            
            # JUST FILTER ALL PROVISIONAL DATA BECAUSE TRACKING ALL THE CODES IS SILLY
            ix=qdf.index[qdf.iloc[:,1]=='P']
            
            # Drop 
            df=df.drop(ix)
            df=df.reset_index(drop=True)
            qdf=qdf.drop(ix)
            qdf=qdf.reset_index(drop=True)        
    
            # Filter again based on empty mean values
            ix=qdf.index[np.isnan(qdf.iloc[:,0].to_numpy().astype(float))]
            # Drop 
            df=df.drop(ix)
            df=df.reset_index(drop=True)
            qdf=qdf.drop(ix)
            qdf=qdf.reset_index(drop=True)
            
            # Extract discharge in cfs
            q=qdf.iloc[:,0].to_numpy().astype(float)
            # Extract quality
            qual=qdf.iloc[:,1]        
            # Conver to m3/sec
            q=q*0.028316846592
    
            # Extract id number
            idnum=int(df.loc[0,'site_no'])
            stat[i]=idnum
            # Find area 
            ix=gages2.index[gages2['ID']==idnum]
            area=gages2.loc[ix,'AREA'].to_numpy()[0]
            # Convert to runoff in mm/day
            r=(q/area)*(24*60*60*100*10)
        
            # Determine completeness of full record
            ts_dr=pd.date_range(df.iloc[0,2],df.iloc[-1,2])
            fc[i]=len(df)/len(ts_dr) 
            flen[i]=len(df)
            
            # Extract start and stop details
            fyr=df['datetime'].dt.year
            fmnth=df['datetime'].dt.month
            fday=df['datetime'].dt.day
            fstartY[i]=fyr.iloc[0]
            fstopY[i]=fyr.iloc[-1]
            fstartM[i]=fmnth.iloc[0]
            fstopM[i]=fmnth.iloc[-1]
            fstartD[i]=fday.iloc[0]
            fstopD[i]=fday.iloc[-1]
            
            # Do calculations for full time series
            # Find mean runoff 
            fmr[i]=np.mean(r)
            
            # Control for records that are mostly zeros
            if len(np.nonzero(r)[0])>3650:
                # Find shape and scale
                [fr_sort,fr_freq]=survive(r)
                [frc1[i],frs1[i]]=weibull_tail_fit(fr_sort,fr_freq,0.01)
                [frc5[i],frs5[i]]=weibull_tail_fit(fr_sort,fr_freq,0.05)
                fr2_5[i]=fr_sort[np.argmin(np.abs(fr_freq-(1/(2.5*365))))]
                fr5[i]=fr_sort[np.argmin(np.abs(fr_freq-(1/(5*365))))]
                fr10[i]=fr_sort[np.argmin(np.abs(fr_freq-(1/(10*365))))]
            else:
                frc1[i]=np.nan
                frs1[i]=np.nan 
                frc5[i]=np.nan 
                frs5[i]=np.nan
                fr2_5[i]=np.nan
                fr5[i]=np.nan 
                fr10[i]=np.nan
            
            # Build temporal index and create sliced time series
            idx=(df['datetime']>=start_) & (df['datetime']<=stop_)
            dfs=df.loc[idx,:]
            dfs=dfs.reset_index(drop=True)
            # Calculate completion of sliced time series
            sc[i]=len(dfs)/len(data_dr)
            slen[i]=len(dfs)
            
            # Generate sliced runoff
            sr=r[idx]
            
            if len(sr)>0:
                smr[i]=np.mean(sr)
                # Extract start and stop details
                syr=dfs['datetime'].dt.year
                smnth=dfs['datetime'].dt.month
                sday=dfs['datetime'].dt.day
                sstartY[i]=syr.iloc[0]
                sstopY[i]=syr.iloc[-1]
                sstartM[i]=smnth.iloc[0]
                sstopM[i]=smnth.iloc[-1]
                sstartD[i]=sday.iloc[0]
                sstopD[i]=sday.iloc[-1]
            else:
                smr[i]=np.nan
                sstartY[i]=np.nan
                sstopY[i]=np.nan
                sstartM[i]=np.nan
                sstopM[i]=np.nan
                sstartD[i]=np.nan
                sstopD[i]=np.nan
            
            # Do calculations for sliced time series
            if len(np.nonzero(sr)[0])>3650:
                # Find shape and scale
                [sr_sort,sr_freq]=survive(sr)
                [src1[i],srs1[i]]=weibull_tail_fit(sr_sort,sr_freq,0.01)
                [src5[i],srs5[i]]=weibull_tail_fit(sr_sort,sr_freq,0.05)
                sr2_5[i]=sr_sort[np.argmin(np.abs(sr_freq-(1/(2.5*365))))]
                sr5[i]=sr_sort[np.argmin(np.abs(sr_freq-(1/(5*365))))]
                sr10[i]=sr_sort[np.argmin(np.abs(sr_freq-(1/(10*365))))]
            else: 
                src1[i]=np.nan 
                srs1[i]=np.nan
                src5[i]=np.nan
                srs5[i]=np.nan
                sr2_5[i]=np.nan
                sr5[i]=np.nan 
                sr10[i]=np.nan
        else:
            stat[i]=np.nan
            fc[i]=np.nan
            flen[i]=np.nan
            fmr[i]=np.nan
            frc1[i]=np.nan
            frs1[i]=np.nan 
            frc5[i]=np.nan 
            frs5[i]=np.nan
            fr2_5[i]=np.nan
            fr5[i]=np.nan 
            fr10[i]=np.nan
            fstartY[i]=np.nan
            fstopY[i]=np.nan
            fstartM[i]=np.nan
            fstopM[i]=np.nan
            fstartD[i]=np.nan
            fstopD[i]=np.nan
            sc[i]=np.nan
            slen[i]=np.nan
            smr[i]=np.nan 
            src1[i]=np.nan 
            srs1[i]=np.nan
            src5[i]=np.nan
            srs5[i]=np.nan
            sr2_5[i]=np.nan
            sr5[i]=np.nan 
            sr10[i]=np.nan
            sstartY[i]=np.nan
            sstopY[i]=np.nan
            sstartM[i]=np.nan
            sstopM[i]=np.nan
            sstartD[i]=np.nan
            sstopD[i]=np.nan
    else:
        stat[i]=np.nan
        fc[i]=np.nan
        flen[i]=np.nan
        fmr[i]=np.nan
        frc1[i]=np.nan
        frs1[i]=np.nan 
        frc5[i]=np.nan 
        frs5[i]=np.nan
        fr2_5[i]=np.nan
        fr5[i]=np.nan 
        fr10[i]=np.nan
        fstartY[i]=np.nan
        fstopY[i]=np.nan
        fstartM[i]=np.nan
        fstopM[i]=np.nan
        fstartD[i]=np.nan
        fstopD[i]=np.nan
        sc[i]=np.nan
        slen[i]=np.nan
        smr[i]=np.nan 
        src1[i]=np.nan 
        srs1[i]=np.nan
        src5[i]=np.nan
        srs5[i]=np.nan
        sr2_5[i]=np.nan
        sr5[i]=np.nan 
        sr10[i]=np.nan
        sstartY[i]=np.nan
        sstopY[i]=np.nan
        sstartM[i]=np.nan
        sstopM[i]=np.nan
        sstartD[i]=np.nan
        sstopD[i]=np.nan

dfout=pd.DataFrame({'STAID':stat,
                    'FullComp':fc,
                    'FullLen':flen,
                    'FullStartY':fstartY,
                    'FullStartM':fstartM,
                    'FullStartD':fstartD,
                    'FullStopY':fstopY,
                    'FullStopM':fstopM,
                    'FullStopD':fstopD,
                    'FullMeanR':fmr,
                    'FullRC1':frc1,
                    'FullRS1':frs1,
                    'FullRC5':frc5,
                    'FullRS5':frs5,
                    'Full2_5':fr2_5,
                    'Full5':fr5,
                    'Full10':fr10,
                    'SlicedComp':sc,
                    'SlicedLen':slen,
                    'SlicedStartY':sstartY,
                    'SlicedStartM':sstartM,
                    'SlicedStartD':sstartD,
                    'SlicedStopY':sstopY,
                    'SlicedStopM':sstopM,
                    'SlicedStopD':sstopD,
                    'SlicedMeanR':smr,
                    'SlicedRC1':src1,
                    'SlicedRS1':srs1,
                    'SlicedRC5':src5,
                    'SlicedRS5':srs5,
                    'Sliced2_5':sr2_5,
                    'Sliced5':sr5,
                    'Sliced10':sr10})

dfout.to_csv('gages2_real_ts.csv',index=False)
    
    
    
    