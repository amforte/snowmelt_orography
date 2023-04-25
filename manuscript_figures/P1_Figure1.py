#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 24 17:16:58 2023

@author: aforte
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from cmcrameri import cm
from scipy.stats import linregress
from matplotlib import colors

master_location='/Users/aforte/Documents/Python/snowmelt/'
# master_location='/Volumes/Choruh/Data/snowmelt_project/'
repo_location='/Users/aforte/Documents/GitHub/snowmelt_orography/geospatial_codes/'

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

idx1=(df['Set']=='ref')
idx2=(rlf>500) & (df['HCDN-2009']=='yes') & (df['MEAN_Z']>250) & (df['SlicedComp']>0.95) & (perc_base<percb_cutoff)

R=np.linspace(0,10,500)
MAR=R*365.25


SMALL_SIZE = 8
MEDIUM_SIZE = 10
BIGGER_SIZE = 12
plt.rc('font', size=SMALL_SIZE,family='Futura')          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

f1=plt.figure(1,figsize=(6,6))
f1.set_dpi(250)
ax1=plt.subplot(1,1,1)
sc1=plt.scatter(df.loc[idx1,'FullMeanR'],df.loc[idx1,'FullRC1'],s=15,c=df.loc[idx1,'MAT'],cmap=cm.vik,label='Gages-II Reference',vmin=-10,vmax=25)
plt.scatter(df.loc[idx2,'FullMeanR'],df.loc[idx2,'FullRC1'],s=20,c=df.loc[idx2,'MAT'],cmap=cm.vik,marker='s',edgecolors='k',label='Filtered HCDN-2009',vmin=-10,vmax=25)
plt.plot(R,0.065*MAR**0.23,c='k',linestyle='--',label='Puerto Rico (Rossi et al., 2016)')
plt.plot(R,0.109*MAR**0.21,c='k',label='CONUS (Rossi et al., 2016)')

cbar1=plt.colorbar(sc1,ax=ax1)
cbar1.ax.set_ylabel('MAT [C]')
plt.legend(loc='best')

plt.xlim((-0.1,10))
plt.ylim((0,2))
plt.xlabel('Mean Runoff [mm/day]')
plt.ylabel('Shape Parameter')

plt.rcdefaults()
f1.savefig('P1_figure1.pdf',dpi='figure')




