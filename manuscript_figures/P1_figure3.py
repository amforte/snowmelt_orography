#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 29 10:24:04 2022

@author: aforte
"""

import pandas as pd
import numpy as np
import glob
import matplotlib.pyplot as plt

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

repo_location='/Users/aforte/Documents/GitHub/snowmelt_orography/geospatial_codes/'
master_location='/Volumes/Choruh/Data/snowmelt_project/'

## Requires gagesII data, available here:
## https://cmerwebmap.cr.usgs.gov/catalog/item/5788f619e4b0d27deb389055

# Load original
gages2raster=pd.read_csv(repo_location+'gages2_wrr2_raster_values.csv')
gages2real=pd.read_csv(repo_location+'gages2_real_ts.csv')
gages2stats=pd.read_csv(repo_location+'gages2_station_stats.csv')

gages2ids1=pd.read_csv(master_location+'gagesII/basinchar_and_report_sept_2011/spreadsheets-in-csv-format/conterm_basinid.txt')
gages2ids2=pd.read_csv(master_location+'gagesII/basinchar_and_report_sept_2011/spreadsheets-in-csv-format/AKHIPR_basinid.txt')
gages2ids=pd.concat([gages2ids1,gages2ids2],axis=0)
gages2hcdn=gages2ids[['STAID','HCDN-2009']]

# Extract lat-lon of center
gages2morph1=pd.read_csv(master_location+'gagesII/basinchar_and_report_sept_2011/spreadsheets-in-csv-format/conterm_bas_morph.txt')
gages2morph2=pd.read_csv(master_location+'gagesII/basinchar_and_report_sept_2011/spreadsheets-in-csv-format/AKHIPR_bas_morph.txt')
gages2morph=pd.concat([gages2morph1,gages2morph2],axis=0)

# Merge
df=pd.merge(gages2raster,gages2hcdn,on='STAID',how='inner')
df=pd.merge(gages2real,df,on='STAID',how='inner')
df=pd.merge(gages2morph,df,on='STAID',how='inner')
df=pd.merge(gages2stats,df,on='STAID',how='inner')
percb_cutoff=0.25
perc_base=df['QSB']/df['R']
rlf=df['MAX_Z']-df['MIN_Z']

# Establish index and produce new dataset
hcdn_idx=(rlf>500) & (df['HCDN-2009']=='yes') & (df['MEAN_Z']>250) & (df['SlicedComp']>0.95) & (perc_base<percb_cutoff)
df_hcdn=df.loc[hcdn_idx,:]

# Build file list of routed gages
fL=glob.glob(master_location+'gagesII_wrr2_hcdn2009_routed/gage*.txt')
# Number files
num_files=len(fL)

staid=np.zeros(num_files).astype(int)
mnr=np.zeros(num_files)
rc1=np.zeros(num_files)
rs1=np.zeros(num_files)
for i in range(num_files):
    fN=fL[i]
    targ_str=master_location+'gagesII_wrr2_hcdn2009_routed/gage_'
    num=fN.replace(targ_str,'')
    num=num.replace('.txt','')
    staid[i]=int(num)
    # Load
    dfts=pd.read_csv(fL[i])
    # Conver to numpy and strip nans
    r=dfts['runoff'].to_numpy()
    r=r[~np.isnan(r)]
    # Calculate variability
    [ts_sort,ts_freq]=survive(r)
    [rc1[i],rs1[i]]=weibull_tail_fit(ts_sort,ts_freq,0.01)
    mnr[i]=np.mean(r)
    
# Dump to pandas dataframe and merge
df_route=pd.DataFrame(data={'STAID':staid,
                            'R_route':mnr,
                            'RC1_route':rc1,
                            'RS1_route':rs1})  
df_merge=pd.merge(df_route,df_hcdn,on='STAID',how='inner')  


SMALL_SIZE = 8
MEDIUM_SIZE = 10
BIGGER_SIZE = 12
plt.rc('font', size=SMALL_SIZE,family='Futura')          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

f1=plt.figure(figsize=(8,8))
f1.set_dpi(250)

ax1=plt.subplot(2,2,1)
plt.scatter(df_hcdn['R'],df_hcdn['SlicedMeanR'],c='w',marker='s',edgecolors='k',s=20,label='Unrouted')
plt.scatter(df_merge['R_route'],df_merge['SlicedMeanR'],c='k',s=25,label='Routed')

x=np.linspace(0,10)
plt.plot(x,x,c='k',linestyle='--')
plt.xlim((0,np.max(x)))
plt.ylim((0,np.max(x)))
plt.xlabel('WaterGAP3 Runoff [mm/day]')
plt.ylabel('Gages-II Runoff [mm/day]')
plt.legend(loc=4)
ax1.text(0.01, 0.99, 'A',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax1.transAxes,
        fontsize=12,fontweight='extra bold')

ax2=plt.subplot(2,2,2)
plt.scatter(df_hcdn['RC1'],df_hcdn['SlicedRC1'],c='w',marker='s',edgecolors='k',s=20,label='Unrouted')
plt.scatter(df_merge['RC1_route'],df_merge['SlicedRC1'],c='k',s=25,label='Routed')
x=np.linspace(0,2.25)
plt.plot(x,x,c='k',linestyle='--')
plt.xlim((0,np.max(x)))
plt.ylim((0,np.max(x)))
plt.xlabel('WaterGAP3 Shape Parameter')
plt.ylabel('Gages-II Shape Parameter')
plt.legend(loc=4)
ax2.text(0.01, 0.99, 'B',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax2.transAxes,
        fontsize=12,fontweight='extra bold')

ax3=plt.subplot(2,2,3)
plt.scatter(df_hcdn['SlicedMeanR']-df_hcdn['R'],df_hcdn['SlicedRC1']-df_hcdn['RC1'],c='w',marker='s',edgecolors='k',s=20,label='Unrouted')
plt.scatter(df_merge['SlicedMeanR']-df_merge['R_route'],df_merge['SlicedRC1']-df_merge['RC1_route'],c='k',s=25,label='Routed')
plt.axvline(0,c='k',linestyle='--')
plt.axhline(0,c='k',linestyle='--')
plt.xlabel('Mean Runoff Residual (GII-WG3)')
plt.ylabel('Shape Residual (GII-WG3)')
plt.legend(loc='best')
ax3.text(0.01, 0.99, 'C',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax3.transAxes,
        fontsize=12,fontweight='extra bold')

ax4=plt.subplot(2,2,4)
plt.scatter(df_hcdn['MAT'],df_hcdn['SlicedRC1']-df_hcdn['RC1'],c='w',marker='s',edgecolors='k',s=20,label='Unrouted')
plt.scatter(df_merge['MAT'],df_merge['SlicedRC1']-df_merge['RC1_route'],c='k',s=25,label='Routed')
plt.axhline(0,c='k',linestyle='--')
plt.xlabel('Mean Annual Temperature [C]')
plt.ylabel('Shape Residual (GII-WG3)')
ax4.text(0.01, 0.99, 'D',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax4.transAxes,
        fontsize=12,fontweight='extra bold')

plt.tight_layout()
plt.rcdefaults()

f1.savefig('P1_figure3.pdf',dpi="figure")