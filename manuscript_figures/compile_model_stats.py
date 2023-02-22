#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 20 07:22:07 2023

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
[mn_ksn1,std_ksn1,mn_ksn1qr,std_ksn1qr,mn_ksn1qp,std_ksn1qp,mn_E1,std_E1]=gc_mcu.ksn_final_ts(rp_slp,rp_yint)
[mn_R1,mn_P1,cr1]=gc_mcu.rpv_final_ts(rp_slp,rp_yint)

gc_mcl=st.ModelComparison(master_location,gc_mll,gc_pll,gc_dll)
[mn_ksn2,std_ksn2,mn_ksn2qr,std_ksn2qr,mn_ksn2qp,std_ksn2qp,mn_E2,std_E2]=gc_mcl.ksn_final_ts(rp_slp,rp_yint)
[mn_R2,mn_P2,cr2]=gc_mcl.rpv_final_ts(rp_slp,rp_yint)

gc_mcu10=st.ModelComparison(master_location,gc_mlu10,gc_plu10,gc_dlu10)
[mn_ksn3,std_ksn3,mn_ksn3qr,std_ksn3qr,mn_ksn3qp,std_ksn3qp,mn_E3,std_E3]=gc_mcu10.ksn_final_ts(rp_slp,rp_yint)
[mn_R3,mn_P3,cr3]=gc_mcu10.rpv_final_ts(rp_slp,rp_yint)

gc_mcl10=st.ModelComparison(master_location,gc_mll10,gc_pll10,gc_dll10)
[mn_ksn4,std_ksn4,mn_ksn4qr,std_ksn4qr,mn_ksn4qp,std_ksn4qp,mn_E4,std_E4]=gc_mcl10.ksn_final_ts(rp_slp,rp_yint)
[mn_R4,mn_P4,cr4]=gc_mcl10.rpv_final_ts(rp_slp,rp_yint)


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
[mn_ksn5,std_ksn5,mn_ksn5qr,std_ksn5qr,mn_ksn5qp,std_ksn5qp,mn_E5,std_E5]=a_mcu.ksn_final_ts(rp_slp,rp_yint)
[mn_R5,mn_P5,cr5]=a_mcu.rpv_final_ts(rp_slp,rp_yint)

a_mcl=st.ModelComparison(master_location,a_mll,a_pll,a_dll)
[mn_ksn6,std_ksn6,mn_ksn6qr,std_ksn6qr,mn_ksn6qp,std_ksn6qp,mn_E6,std_E6]=a_mcl.ksn_final_ts(rp_slp,rp_yint)
[mn_R6,mn_P6,cr6]=a_mcl.rpv_final_ts(rp_slp,rp_yint)

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
[mn_ksn7,std_ksn7,mn_ksn7qr,std_ksn7qr,mn_ksn7qp,std_ksn7qp,mn_E7,std_E7]=bc_mcu.ksn_final_ts(rp_slp,rp_yint)
[mn_R7,mn_P7,cr7]=bc_mcu.rpv_final_ts(rp_slp,rp_yint)

bc_mcl=st.ModelComparison(master_location,bc_mll,bc_pll,bc_dll)
[mn_ksn8,std_ksn8,mn_ksn8qr,std_ksn8qr,mn_ksn8qp,std_ksn8qp,mn_E8,std_E8]=bc_mcl.ksn_final_ts(rp_slp,rp_yint)
[mn_R8,mn_P8,cr8]=bc_mcl.rpv_final_ts(rp_slp,rp_yint)

# Package
model=gc_mlu+gc_mll+gc_mlu10+gc_mll10+a_mlu+a_mll+bc_mlu+bc_mll
group=gc_dlu+gc_dll+gc_dlu10+gc_dll10+a_dlu+a_dll+bc_dlu+bc_dll
mn_ksn=np.concatenate((mn_ksn1,mn_ksn2,mn_ksn3,mn_ksn4,mn_ksn5,mn_ksn6,mn_ksn7,mn_ksn8),axis=0)
std_ksn=np.concatenate((std_ksn1,std_ksn2,std_ksn3,std_ksn4,std_ksn5,std_ksn6,std_ksn7,std_ksn8),axis=0)
mn_ksnqr=np.concatenate((mn_ksn1qr,mn_ksn2qr,mn_ksn3qr,mn_ksn4qr,mn_ksn5qr,mn_ksn6qr,mn_ksn7qr,mn_ksn8qr),axis=0)
std_ksnqr=np.concatenate((std_ksn1qr,std_ksn2qr,std_ksn3qr,std_ksn4qr,std_ksn5qr,std_ksn6qr,std_ksn7qr,std_ksn8qr),axis=0)
mn_ksnqp=np.concatenate((mn_ksn1qp,mn_ksn2qp,mn_ksn3qp,mn_ksn4qp,mn_ksn5qp,mn_ksn6qp,mn_ksn7qp,mn_ksn8qp),axis=0)
std_ksnqp=np.concatenate((std_ksn1qp,std_ksn2qp,std_ksn3qp,std_ksn4qp,std_ksn5qp,std_ksn6qp,std_ksn7qp,std_ksn8qp),axis=0)
mn_E=np.concatenate((mn_E1,mn_E2,mn_E3,mn_E4,mn_E5,mn_E6,mn_E7,mn_E8),axis=0)
std_E=np.concatenate((std_E1,std_E2,std_E3,std_E4,std_E5,std_E6,std_E7,std_E8),axis=0)
mn_R=np.concatenate((mn_R1,mn_R2,mn_R3,mn_R4,mn_R5,mn_R6,mn_R7,mn_R8),axis=0)
mn_P=np.concatenate((mn_P1,mn_P2,mn_P3,mn_P4,mn_P5,mn_P6,mn_P7,mn_P8),axis=0)
cr=np.concatenate((cr1,cr2,cr3,cr4,cr5,cr6,cr7,cr8),axis=0)

df=pd.DataFrame(data={'Model':model,
                      'Group':group,
                      'mn_ksn':mn_ksn,
                      'std_ksn':std_ksn,
                      'mn_ksnqr':mn_ksnqr,
                      'std_ksnqr':std_ksnqr,
                      'mn_ksnqp':mn_ksnqp,
                      'std_ksnqp':std_ksnqp,
                      'mn_E':mn_E,
                      'std_E':std_E,
                      'mn_R':mn_R,
                      'mn_P':mn_P,
                      'cr':cr})




