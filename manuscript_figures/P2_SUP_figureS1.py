#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 07:27:04 2022

@author: aforte
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import rasterio
from scipy.stats import linregress
import matplotlib.gridspec as gridspec
from scipy import stats

def cdf(data,bins):
    H,X1 = np.histogram(data,bins=bins,density=True )
    dx = X1[1] - X1[0]
    F1 = np.cumsum(H)*dx
    return F1

def bootstrap_conf(data,bins,sample_size,num_trials):
    X=np.zeros((num_trials,len(bins)-1))
    datas=np.random.choice(data,(num_trials,sample_size))
    for i in range(num_trials):
        X[i,:]=cdf(datas[i,:],bins)
    x5=np.percentile(X,5,axis=0)
    x95=np.percentile(X,95,axis=0)
    return x5,x95


repo_location='/Users/aforte/Documents/GitHub/snowmelt_orography/geospatial_codes/'
gc=pd.read_csv(repo_location+'gc_ksn_rlf.csv')
alps=pd.read_csv(repo_location+'alps_ksn_rlf.csv')
bc=pd.read_csv(repo_location+'bc_ksn_rlf.csv')

gc_col='black'
alps_col='royalblue'
bc_col='orange'

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


f1=plt.figure(figsize=(8,10))
# f1.set_dpi(250)

gs=f1.add_gridspec(3,6)

# ax1=plt.subplot(3,2,1)
ax1=f1.add_subplot(gs[0,0:3])
plt.scatter(gc['mean_ksn'],gc['mean_rlf2500'],s=2,c=gc_col,label='Greater Caucasus')
plt.scatter(alps['mean_ksn'],alps['mean_rlf2500'],s=2,c=alps_col,label='Alps')
plt.scatter(bc['mean_ksn'],bc['mean_rlf2500'],s=2,c=bc_col,label='BC')
plt.ylim((0,2500))
plt.xlabel(r'Mean $k_{sn}$ [m]')
plt.ylabel('Mean Local 2.5 km Relief [m]')
plt.legend(loc=4)
plt.title('All Basins')
ax1.text(0.01, 0.99, 'A',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax1.transAxes,
        fontsize=12,fontweight='extra bold')

chi_R2=0.90
gc_idx=gc['chi_R_squared']>=chi_R2
alps_idx=alps['chi_R_squared']>=chi_R2
bc_idx=bc['chi_R_squared']>=chi_R2

[gc_slp,gc_int,gc_r,_,_]=linregress(gc.loc[gc_idx,'mean_ksn'],gc.loc[gc_idx,'mean_rlf2500'])
[alps_slp,alps_int,alps_r,_,_]=linregress(alps.loc[alps_idx,'mean_ksn'],alps.loc[alps_idx,'mean_rlf2500'])
[bc_slp,bc_int,bc_r,_,_]=linregress(bc.loc[bc_idx,'mean_ksn'],bc.loc[bc_idx,'mean_rlf2500'])

gc_r2=str(np.round(gc_r**2,2))
alps_r2=str(np.round(alps_r**2,2))
bc_r2=str(np.round(bc_r**2,2))
gc_slp_s=str(np.round(gc_slp,2))
alps_slp_s=str(np.round(alps_slp,2))
bc_slp_s=str(np.round(bc_slp,2))
gc_int_s=str(np.round(gc_int,2))
alps_int_s=str(np.round(alps_int,2))
bc_int_s=str(np.round(bc_int,2))

gc_label=r'LR = '+gc_slp_s+'$k_{sn}$ + '+gc_int_s+'; $R^{2}$ = '+gc_r2
alps_label=r'LR = '+alps_slp_s+'$k_{sn}$ + '+alps_int_s+'; $R^{2}$ = '+alps_r2
bc_label=r'LR = '+bc_slp_s+'$k_{sn}$ + '+bc_int_s+'; $R^{2}$ = '+bc_r2

ax3=f1.add_subplot(gs[1,0:3])
ksn_vec=np.linspace(0,800,100)
plt.scatter(gc.loc[gc_idx,'mean_ksn'],gc.loc[gc_idx,'mean_rlf2500'],s=2,c=gc_col)
plt.scatter(alps.loc[alps_idx,'mean_ksn'],alps.loc[alps_idx,'mean_rlf2500'],s=2,c=alps_col)
plt.scatter(bc.loc[bc_idx,'mean_ksn'],bc.loc[bc_idx,'mean_rlf2500'],s=2,c=bc_col)
plt.plot(ksn_vec,gc_slp*ksn_vec + gc_int,c=gc_col,label=gc_label)
plt.plot(ksn_vec,alps_slp*ksn_vec + alps_int,c=alps_col,label=alps_label)
plt.plot(ksn_vec,bc_slp*ksn_vec + bc_int,c=bc_col,label=bc_label)
plt.ylim((0,2500))
plt.xlabel(r'Mean $k_{sn}$ [m]')
plt.ylabel('Mean Local 2.5 km Relief [m]')
plt.legend(loc=4)
plt.title(r'Filtered to $\chi$ $R^{2}$ > 0.90')
ax3.text(0.01, 0.99, 'C',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax3.transAxes,
        fontsize=12,fontweight='extra bold')

ax2=f1.add_subplot(gs[0,3:])
plt.scatter(gc['mean_gradient'],gc['mean_rlf2500'],s=2,c=gc_col,label='Greater Caucasus')
plt.scatter(alps['mean_gradient'],alps['mean_rlf2500'],s=2,c=alps_col,label='Alps')
plt.scatter(bc['mean_gradient'],bc['mean_rlf2500'],s=2,c=bc_col,label='BC')
plt.ylim((0,2500))
plt.xlabel('Mean Gradient [m/m]')
plt.ylabel('Mean Local 2.5 km Relief [m]')
# plt.legend(loc=4)
plt.title('All Basins')
ax2.text(0.01, 0.99, 'B',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax2.transAxes,
        fontsize=12,fontweight='extra bold')


ax4=f1.add_subplot(gs[1,3:])
plt.scatter(gc.loc[gc_idx,'mean_gradient'],gc.loc[gc_idx,'mean_rlf2500'],s=2,c=gc_col,label='Greater Caucasus')
plt.scatter(alps.loc[alps_idx,'mean_gradient'],alps.loc[alps_idx,'mean_rlf2500'],s=2,c=alps_col,label='Alps')
plt.scatter(bc.loc[bc_idx,'mean_gradient'],bc.loc[bc_idx,'mean_rlf2500'],s=2,c=bc_col,label='BC')
plt.ylim((0,2500))
plt.xlabel('Mean Gradient [m/m]')
plt.ylabel('Mean Local 2.5 km Relief [m]')
# plt.legend(loc=4)
plt.title(r'Filtered to $\chi$ $R^{2}$ > 0.90')
ax4.text(0.01, 0.99, 'D',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax4.transAxes,
        fontsize=12,fontweight='extra bold')


master_location='/Users/aforte/Documents/Python/snowmelt/'
# master_location='/Volumes/Choruh/Data/snowmelt_project'

rstr1=rasterio.open(master_location+'hyd_glo_dem_15s/hsheds_2500rlf_wgs84_rsmp.tif')
rlf=rstr1.read(1)

rstr2=rasterio.open(master_location+'regions/ALPS_SRTM90_RLF25.tif')
alps_rlf=rstr2.read(1)

rstr3=rasterio.open(master_location+'regions/GC_SRTM90_RLF25.tif')
gc_rlf=rstr3.read(1)

rstr4=rasterio.open(master_location+'regions/BC_SRTM90_RLF25.tif')
bc_rlf=rstr4.read(1)

bins=np.linspace(0,2500,250)

# Alps
bl2=5
br2=16
bb2=43
bt2=50
alps_col='royalblue'

[bl2r,bl2c]=rasterio.transform.rowcol(rstr1.transform,bl2,bb2)
[br2r,br2c]=rasterio.transform.rowcol(rstr1.transform,br2,bb2)
[bb2r,bb2c]=rasterio.transform.rowcol(rstr1.transform,bl2,bb2)
[bt2r,bt2c]=rasterio.transform.rowcol(rstr1.transform,bl2,bt2)

alps_rlf_clip=rlf[bt2r:bb2r+1,bl2c:br2c+1]

alps_rlf=alps_rlf.ravel()
alps_rlf=alps_rlf[(~np.isnan(alps_rlf)) & (alps_rlf>250)]

alps_rlf_clip=alps_rlf_clip.ravel()
alps_rlf_clip=alps_rlf_clip[(alps_rlf_clip!=rstr1.nodata) & (alps_rlf_clip>250)]

alps_cdf=cdf(alps_rlf,bins)
alps_clip_cdf=cdf(alps_rlf_clip,bins)

# Greater Caucasus
bl1=38
br1=51
bb1=39.5
bt1=45
gc_col='black'

[bl1r,bl1c]=rasterio.transform.rowcol(rstr1.transform,bl1,bb1)
[br1r,br1c]=rasterio.transform.rowcol(rstr1.transform,br1,bb1)
[bb1r,bb1c]=rasterio.transform.rowcol(rstr1.transform,bl1,bb1)
[bt1r,bt1c]=rasterio.transform.rowcol(rstr1.transform,bl1,bt1)

gc_rlf_clip=rlf[bt1r:bb1r+1,bl1c:br1c+1]

gc_rlf=gc_rlf.ravel()
gc_rlf=gc_rlf[(~np.isnan(gc_rlf)) & (gc_rlf>250)]

gc_rlf_clip=gc_rlf_clip.ravel()
gc_rlf_clip=gc_rlf_clip[(gc_rlf_clip!=rstr1.nodata) & (gc_rlf_clip>250)]

gc_cdf=cdf(gc_rlf,bins)
gc_clip_cdf=cdf(gc_rlf_clip,bins)

# British Columbia
bl3=-131
br3=-120
bb3=48
bt3=54
bc_col='orange'

[bl3r,bl3c]=rasterio.transform.rowcol(rstr1.transform,bl3,bb3)
[br3r,br3c]=rasterio.transform.rowcol(rstr1.transform,br3,bb3)
[bb3r,bb3c]=rasterio.transform.rowcol(rstr1.transform,bl3,bb3)
[bt3r,bt3c]=rasterio.transform.rowcol(rstr1.transform,bl3,bt3)

bc_rlf_clip=rlf[bt3r:bb3r+1,bl3c:br3c+1]

bc_rlf=bc_rlf.ravel()
bc_rlf=bc_rlf[(~np.isnan(bc_rlf)) & (bc_rlf>250)]

bc_rlf_clip=bc_rlf_clip.ravel()
bc_rlf_clip=bc_rlf_clip[(bc_rlf_clip!=rstr1.nodata) & (bc_rlf_clip>250)]

bc_cdf=cdf(bc_rlf,bins)
bc_clip_cdf=cdf(bc_rlf_clip,bins)


print('Starting Bootstrap')
[alps5,alps95]=bootstrap_conf(alps_rlf,bins,500,1000)
[alpsc5,alpsc95]=bootstrap_conf(alps_rlf_clip,bins,500,1000)
[gc5,gc95]=bootstrap_conf(gc_rlf,bins,500,1000)
[gcc5,gcc95]=bootstrap_conf(gc_rlf_clip,bins,500,1000)
[bc5,bc95]=bootstrap_conf(bc_rlf,bins,500,1000)
[bcc5,bcc95]=bootstrap_conf(bc_rlf_clip,bins,500,1000)

# ax5=f1.add_subplot(gs[2,0:])

ax5=f1.add_subplot(gs[2,2:4])
plt.plot(bins[1:],alps_cdf,color=alps_col,linewidth=2,label='SRTM-90')
plt.plot(bins[1:],alps_clip_cdf,color=alps_col,linewidth=2,linestyle='--',label='HS-15s')
plt.plot(bins[1:],alps5,color=alps_col,linewidth=0.5)
plt.plot(bins[1:],alps95,color=alps_col,linewidth=0.5)
plt.plot(bins[1:],alpsc5,color=alps_col,linewidth=0.5,linestyle='--')
plt.plot(bins[1:],alpsc95,color=alps_col,linewidth=0.5,linestyle='--')
plt.xlabel('Local Relief 2.5 km')
plt.ylabel('Cumulative Proability')
plt.legend(loc='lower right')
plt.title('Alps')
plt.xlim((250,2500))
ax5.text(0.01, 0.99, 'F',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax5.transAxes,
        fontsize=12,fontweight='extra bold')


ax6=f1.add_subplot(gs[2,4:])
plt.plot(bins[1:],gc_cdf,color=gc_col,linewidth=2,label='SRTM-90')
plt.plot(bins[1:],gc_clip_cdf,color=gc_col,linewidth=2,linestyle='--',label='HS-15s')
plt.plot(bins[1:],gc5,color=gc_col,linewidth=0.5)
plt.plot(bins[1:],gc95,color=gc_col,linewidth=0.5)
plt.plot(bins[1:],gcc5,color=gc_col,linewidth=0.5,linestyle='--')
plt.plot(bins[1:],gcc95,color=gc_col,linewidth=0.5,linestyle='--')
plt.xlabel('Local Relief 2.5 km')
plt.ylabel('Cumulative Proability')
plt.legend(loc='lower right')
plt.title('GC')
plt.xlim((250,2500))
ax6.text(0.01, 0.99, 'G',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax6.transAxes,
        fontsize=12,fontweight='extra bold')

ax7=f1.add_subplot(gs[2,0:2])
plt.plot(bins[1:],bc_cdf,color=bc_col,linewidth=2,label='SRTM-90')
plt.plot(bins[1:],bc_clip_cdf,color=bc_col,linewidth=2,linestyle='--',label='HS-15s')
plt.plot(bins[1:],bc5,color=bc_col,linewidth=0.5)
plt.plot(bins[1:],bc95,color=bc_col,linewidth=0.5)
plt.plot(bins[1:],bcc5,color=bc_col,linewidth=0.5,linestyle='--')
plt.plot(bins[1:],bcc95,color=bc_col,linewidth=0.5,linestyle='--')
plt.xlabel('Local Relief 2.5 km')
plt.ylabel('Cumulative Proability')
plt.legend(loc='lower right')
plt.title('BC')
plt.xlim((250,2500))
ax7.text(0.01, 0.99, 'E',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax7.transAxes,
        fontsize=12,fontweight='extra bold')



plt.tight_layout()
plt.rcdefaults()
f1.savefig('P2_SUP_figureS1.pdf')
