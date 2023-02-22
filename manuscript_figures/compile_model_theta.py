#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  7 07:43:18 2023

@author: aforte
"""

import sys
sys.path.insert(0,'/Users/aforte/Documents/GitHub/snowmelt-tectonics')
import stimpy as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from cmcrameri import cm
from scipy import odr

def linear(B,x):
    return B[0]*x + B[1]

def odr_fit_ref(x,y):
    # Filter 0 values
    lx=x[(x>0) & (y>0)]
    ly=y[(x>0) & (y>0)]
    linmod=odr.Model(linear)
    
    fdlog=odr.Data(lx,ly)
    odrlog=odr.ODR(fdlog,linmod,beta0=[0.1,10])
    outlog=odrlog.run()
    slp=outlog.beta[0]
    yint=outlog.beta[1]

    return slp,yint

def plot_and_fit_ref(axn,dfs,bl,br,bb,bt,col,lbl_left):
    spidx=(dfs['latitude']>=bb) & (dfs['latitude']<=bt) & (dfs['longitude']>=bl) & (dfs['longitude']<=br)

    max_r=12
    max_p=12
    
    r=np.linspace(1,max_r,100)
    
    x=dfs.loc[spidx,'mean_runoff'].to_numpy()
    y=dfs.loc[spidx,'mean_precip'].to_numpy()
    s,yi=odr_fit_ref(x,y)
    bin_edges=np.histogram_bin_edges(dfs.loc[spidx,'mean_precip'],'doane')
    bix=np.digitize(dfs.loc[spidx,'mean_precip'],bin_edges)
    for i in range(len(bin_edges)-1):
        mx=np.mean(x[bix==i])
        stdx=np.std(x[bix==i])
        my=np.mean(y[bix==i])
        stdy=np.std(y[bix==i])
        axn.scatter(x,y,c=col,s=1,zorder=1,alpha=0.25)
        axn.scatter(mx,my,c=col,s=len(x[bix==i]),zorder=2)
        axn.errorbar(mx,my,xerr=stdx,yerr=stdy,ecolor=col,elinewidth=0.5,zorder=1)
    axn.plot([0,max_r],[0,max_p],c='gray',linestyle=':',zorder=0)
    axn.plot(r,s*r+yi,c=col,label='Mean Runoff to Mean Precip')  
    axn.set_xlabel('Mean Runoff [m]')
    if lbl_left:
        axn.set_ylabel('Mean Precip [mm/day]')
    axn.set_xlim((0,max_r))
    axn.set_ylim((0,max_p))
    return s,yi




# 
# Greater Caucasus
bl1=38
br1=51
bb1=39.5
bt1=45
gc_col='black'

# Alps
bl2=5
br2=16
bb2=43
bt2=50
alps_col='royalblue'

# British Columbia
bl3=-131
br3=-120
bb3=48
bt3=54
bc_col='orange'

## Determine runoff to precip relationships
## Load global
df_global=pd.read_csv('/Volumes/Choruh/Data/snowmelt_project/wrr2_derived_data_v3.csv')
df_global=df_global.drop(index=df_global.index[np.isnan(df_global['mean_z'])])
df_global=df_global.reset_index(drop=True)
# Calc percents
global_perc_base=df_global['qsb']/df_global['mean_runoff']
# Calculate indices
grlf=df_global['max_z']-df_global['min_z']
# Set cutoffs
percb_cutoff=0.25
# Index
rem_idx=(grlf<=500) | (df_global['mean_z']<=250) | (global_perc_base>percb_cutoff) | (df_global['mean_rlf']<250)
df_global_s=df_global.drop(df_global.index[rem_idx])
df_global_s=df_global_s.reset_index(drop=True)

f1=plt.figure(figsize=(8,2.5))
f1.set_dpi(250)
ax1=plt.subplot(1,3,1)
ax2=plt.subplot(1,3,2)
ax3=plt.subplot(1,3,3)

[gcs,gcyi]=plot_and_fit_ref(ax3,df_global_s,bl1,br1,bb1,bt1,gc_col,False)
[als,alyi]=plot_and_fit_ref(ax2,df_global_s,bl2,br2,bb2,bt2,alps_col,False)
[bcs,bcyi]=plot_and_fit_ref(ax1,df_global_s,bl3,br3,bb3,bt3,bc_col,True)

master_location='/Volumes/Choruh/Data/snowmelt_project/model_outputs/'



gc_mlu=['gc025u','gc05u','gc1u','gc2u','gc4u','gc8u']
gc_plu=['ts_','ts_','ts_','ts_','ts_','ts_']
gc_dlu=['gcu','gcu','gcu','gcu','gcu','gcu']

gc_mll=['gc025l','gc05l','gc1l','gc2l','gc4l','gc8l']
gc_pll=['ts_','ts_','ts_','ts_','ts_','ts_']
gc_dll=['gcl','gcl','gcl','gcl','gcl','gcl']

gc_mlu10=['gc025u_10l','gc05u_10l','gc1u_10l','gc2u_10l','gc4u_10l','gc8u_10l']
gc_plu10=['ts_','ts_','ts_','ts_','ts_','ts_']
gc_dlu10=['gcu_10l','gcu_10l','gcu_10l','gcu_10l','gcu_10l','gcu_10l']

gc_mll10=['gc025l_10l','gc05l_10l','gc1l_10l','gc2l_10l','gc4l_10l','gc8l_10l']
gc_pll10=['ts_','ts_','ts_','ts_','ts_','ts_']
gc_dll10=['gcl_10l','gcl_10l','gcl_10l','gcl_10l','gcl_10l','gcl_10l']

rp_slp=[]
rp_yint=[]
for i in range(len(gc_mlu)):
    rp_slp.append(gcs)
    rp_yint.append(gcyi)

gc_mcu=st.ModelComparison(master_location,gc_mlu,gc_plu,gc_dlu)
[thc1,ths1,thcr1,thcp1]=gc_mcu.theta_final_ts(rp_slp,rp_yint)

gc_mcl=st.ModelComparison(master_location,gc_mll,gc_pll,gc_dll)
[thc2,ths2,thcr2,thcp2]=gc_mcl.theta_final_ts(rp_slp,rp_yint)

gc_mcu10=st.ModelComparison(master_location,gc_mlu10,gc_plu10,gc_dlu10)
[thc3,ths3,thcr3,thcp3]=gc_mcu10.theta_final_ts(rp_slp,rp_yint)

gc_mcl10=st.ModelComparison(master_location,gc_mll10,gc_pll10,gc_dll10)
[thc4,ths4,thcr4,thcp4]=gc_mcl10.theta_final_ts(rp_slp,rp_yint)


a_mlu=['a025u','a05u','a1u','a2u','a4u','a8u']
a_plu=['ts_','ts_','ts_','ts_','ts_','ts_']
a_dlu=['au','au','au','au','au','au']

a_mll=['a025l','a05l','a1l','a2l','a4l','a8l']
a_pll=['ts_','ts_','ts_','ts_','ts_','ts_']
a_dll=['al','al','al','al','al','al']

rp_slp=[]
rp_yint=[]
for i in range(len(a_mlu)):
    rp_slp.append(als)
    rp_yint.append(alyi)
    
a_mcu=st.ModelComparison(master_location,a_mlu,a_plu,a_dlu)
[thc5,ths5,thcr5,thcp5]=a_mcu.theta_final_ts(rp_slp,rp_yint)

a_mcl=st.ModelComparison(master_location,a_mll,a_pll,a_dll)
[thc6,ths6,thcr6,thcp6]=a_mcl.theta_final_ts(rp_slp,rp_yint)

bc_mlu=['bc025u','bc05u','bc1u','bc2u','bc4u','bc8u']
bc_plu=['ts_','ts_','ts_','ts_','ts_','ts_']
bc_dlu=['bcu','bcu','bcu','bcu','bcu','bcu']

bc_mll=['bc025l','bc05l','bc1l','bc2l','bc4l','bc8l']
bc_pll=['ts_','ts_','ts_','ts_','ts_','ts_']
bc_dll=['bcl','bcl','bcl','bcl','bcl','bcl']

rp_slp=[]
rp_yint=[]
for i in range(len(bc_mlu)):
    rp_slp.append(bcs)
    rp_yint.append(bcyi)
    
bc_mcu=st.ModelComparison(master_location,bc_mlu,bc_plu,bc_dlu)
[thc7,ths7,thcr7,thcp7]=bc_mcu.theta_final_ts(rp_slp,rp_yint)

bc_mcl=st.ModelComparison(master_location,bc_mll,bc_pll,bc_dll)
[thc8,ths8,thcr8,thcp8]=bc_mcl.theta_final_ts(rp_slp,rp_yint)

# Package
model=gc_mlu+gc_mll+gc_mlu10+gc_mll10+a_mlu+a_mll+bc_mlu+bc_mll
group=gc_dlu+gc_dll+gc_dlu10+gc_dll10+a_dlu+a_dll+bc_dlu+bc_dll
theta_chi=np.concatenate((thc1,thc2,thc3,thc4,thc5,thc6,thc7,thc8),axis=0)
theta_slope=np.concatenate((ths1,ths2,ths3,ths4,ths5,ths6,ths7,ths8),axis=0)
theta_chir=np.concatenate((thcr1,thcr2,thcr3,thcr4,thcr5,thcr6,thcr7,thcr8),axis=0)
theta_chip=np.concatenate((thcp1,thcp2,thcp3,thcp4,thcp5,thcp6,thcp7,thcp8),axis=0)

df=pd.DataFrame(data={'Model':model,
                      'Group':group,
                      'theta_chi':theta_chi,
                      'theta_slope':theta_slope,
                      'theta_chi_r':theta_chir,
                      'theta_chi_p':theta_chip})






