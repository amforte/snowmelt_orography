#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 20 07:22:07 2023

@author: aforte
"""

import sys
sys.path.insert(0,'/Users/aforte/Documents/GitHub/snowmelt_orography')
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

# master_location='/Volumes/Choruh/Data/snowmelt_project/'
master_data_location='/Users/aforte/Documents/Python/snowmelt/'

## Determine runoff to precip relationships
## Load global
df_global=pd.read_csv(master_data_location+'wrr2_derived_data_v4.csv')
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

master_location=master_data_location+'model_outputs_v2/'

u_vec=np.array([250,500,1000,2000,4000,8000])

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

gc_mluA=['gc025u_Area','gc05u_Area','gc1u_Area','gc2u_Area','gc4u_Area','gc8u_Area']
gc_pluA=['ts_','ts_','ts_','ts_','ts_','ts_']
gc_dluA=['gcu_A','gcu_A','gcu_A','gcu_A','gcu_A','gcu_A']

gc_mllA=['gc025l_Area','gc05l_Area','gc1l_Area','gc2l_Area','gc4l_Area','gc8l_Area']
gc_pllA=['ts_','ts_','ts_','ts_','ts_','ts_']
gc_dllA=['gcl_A','gcl_A','gcl_A','gcl_A','gcl_A','gcl_A']

gc_mluR=['gc025u_RO','gc05u_RO','gc1u_RO','gc2u_RO','gc4u_RO','gc8u_RO']
gc_pluR=['ts_','ts_','ts_','ts_','ts_','ts_']
gc_dluR=['gcu_R','gcu_R','gcu_R','gcu_R','gcu_R','gcu_R']

gc_mllR=['gc025l_RO','gc05l_RO','gc1l_RO','gc2l_RO','gc4l_RO','gc8l_RO']
gc_pllR=['ts_','ts_','ts_','ts_','ts_','ts_']
gc_dllR=['gcl_R','gcl_R','gcl_R','gcl_R','gcl_R','gcl_R']

rp_slp=[]
rp_yint=[]
for i in range(len(gc_mlu)):
    rp_slp.append(gcs)
    rp_yint.append(gcyi)

gc_mcu=st.ModelComparison(master_location,gc_mlu,gc_plu,gc_dlu)
[mn_ksn1,std_ksn1,stde_ksn1,mn_ksn1qr,std_ksn1qr,stde_ksn1qr,mn_ksn1qp,std_ksn1qp,stde_ksn1qp,mn_E1,std_E1,stde_E1,p25_E1,p75_E1]=gc_mcu.ksn_final_ts(rp_slp,rp_yint)
[mn_R1,mn_P1,cr1,sr1]=gc_mcu.rpv_final_ts(rp_slp,rp_yint)

gc_mcl=st.ModelComparison(master_location,gc_mll,gc_pll,gc_dll)
[mn_ksn2,std_ksn2,stde_ksn2,mn_ksn2qr,std_ksn2qr,stde_ksn2qr,mn_ksn2qp,std_ksn2qp,stde_ksn2qp,mn_E2,std_E2,stde_E2,p25_E2,p75_E2]=gc_mcl.ksn_final_ts(rp_slp,rp_yint)
[mn_R2,mn_P2,cr2,sr2]=gc_mcl.rpv_final_ts(rp_slp,rp_yint)

gc_mcu10=st.ModelComparison(master_location,gc_mlu10,gc_plu10,gc_dlu10)
[mn_ksn3,std_ksn3,stde_ksn3,mn_ksn3qr,std_ksn3qr,stde_ksn3qr,mn_ksn3qp,std_ksn3qp,stde_ksn3qp,mn_E3,std_E3,stde_E3,p25_E3,p75_E3]=gc_mcu10.ksn_final_ts(rp_slp,rp_yint)
[mn_R3,mn_P3,cr3,sr3]=gc_mcu10.rpv_final_ts(rp_slp,rp_yint)

gc_mcl10=st.ModelComparison(master_location,gc_mll10,gc_pll10,gc_dll10)
[mn_ksn4,std_ksn4,stde_ksn4,mn_ksn4qr,std_ksn4qr,stde_ksn4qr,mn_ksn4qp,std_ksn4qp,stde_ksn4qp,mn_E4,std_E4,stde_E4,p25_E4,p75_E4]=gc_mcl10.ksn_final_ts(rp_slp,rp_yint)
[mn_R4,mn_P4,cr4,sr4]=gc_mcl10.rpv_final_ts(rp_slp,rp_yint)

gc_mcuA=st.ModelComparison(master_location,gc_mluA,gc_pluA,gc_dluA)
[mn_ksn9,std_ksn9,stde_ksn9,mn_ksn9qr,std_ksn9qr,stde_ksn9qr,mn_ksn9qp,std_ksn9qp,stde_ksn9qp,mn_E9,std_E9,stde_E9,p25_E9,p75_E9]=gc_mcuA.ksn_final_ts(rp_slp,rp_yint)
[mn_R9,mn_P9,cr9,sr9]=gc_mcuA.rpv_final_ts(rp_slp,rp_yint)

gc_mclA=st.ModelComparison(master_location,gc_mllA,gc_pllA,gc_dllA)
[mn_ksn10,std_ksn10,stde_ksn10,mn_ksn10qr,std_ksn10qr,stde_ksn10qr,mn_ksn10qp,std_ksn10qp,stde_ksn10qp,mn_E10,std_E10,stde_E10,p25_E10,p75_E10]=gc_mclA.ksn_final_ts(rp_slp,rp_yint)
[mn_R10,mn_P10,cr10,sr10]=gc_mclA.rpv_final_ts(rp_slp,rp_yint)

gc_mcuR=st.ModelComparison(master_location,gc_mluR,gc_pluR,gc_dluR)
[mn_ksn11,std_ksn11,stde_ksn11,mn_ksn11qr,std_ksn11qr,stde_ksn11qr,mn_ksn11qp,std_ksn11qp,stde_ksn11qp,mn_E11,std_E11,stde_E11,p25_E11,p75_E11]=gc_mcuR.ksn_final_ts(rp_slp,rp_yint)
[mn_R11,mn_P11,cr11,sr11]=gc_mcuR.rpv_final_ts(rp_slp,rp_yint)

gc_mclR=st.ModelComparison(master_location,gc_mllR,gc_pllR,gc_dllR)
[mn_ksn12,std_ksn12,stde_ksn12,mn_ksn12qr,std_ksn12qr,stde_ksn12qr,mn_ksn12qp,std_ksn12qp,stde_ksn12qp,mn_E12,std_E12,stde_E12,p25_E12,p75_E12]=gc_mclR.ksn_final_ts(rp_slp,rp_yint)
[mn_R12,mn_P12,cr12,sr12]=gc_mclR.rpv_final_ts(rp_slp,rp_yint)


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
[mn_ksn5,std_ksn5,stde_ksn5,mn_ksn5qr,std_ksn5qr,stde_ksn5qr,mn_ksn5qp,std_ksn5qp,stde_ksn5qp,mn_E5,std_E5,stde_E5,p25_E5,p75_E5]=a_mcu.ksn_final_ts(rp_slp,rp_yint)
[mn_R5,mn_P5,cr5,sr5]=a_mcu.rpv_final_ts(rp_slp,rp_yint)

a_mcl=st.ModelComparison(master_location,a_mll,a_pll,a_dll)
[mn_ksn6,std_ksn6,stde_ksn6,mn_ksn6qr,std_ksn6qr,stde_ksn6qr,mn_ksn6qp,std_ksn6qp,stde_ksn6qp,mn_E6,std_E6,stde_E6,p25_E6,p75_E6]=a_mcl.ksn_final_ts(rp_slp,rp_yint)
[mn_R6,mn_P6,cr6,sr6]=a_mcl.rpv_final_ts(rp_slp,rp_yint)

bc_mlu=['bc025u','bc05u','bc1u','bc2u','bc4u','bc8u']
bc_plu=['ts_','ts_','ts_','ts_','ts_','ts_']
bc_dlu=['bcu','bcu','bcu','bcu','bcu','bcu']

bc_mll=['bc025l','bc05l','bc1l','bc2l','bc4l','bc8l']
bc_pll=['ts_','ts_','ts_','ts_','ts_','ts_']
bc_dll=['bcl','bcl','bcl','bcl','bcl','bcl']

bc_mluR=['bc025u_RO','bc05u_RO','bc1u_RO','bc2u_RO','bc4u_RO','bc8u_RO']
bc_pluR=['ts_','ts_','ts_','ts_','ts_','ts_']
bc_dluR=['bcu_R','bcu_R','bcu_R','bcu_R','bcu_R','bcu_R']

bc_mllR=['bc025l_RO','bc05l_RO','bc1l_RO','bc2l_RO','bc4l_RO','bc8l_RO']
bc_pllR=['ts_','ts_','ts_','ts_','ts_','ts_']
bc_dllR=['bcl_R','bcl_R','bcl_R','bcl_R','bcl_R','bcl_R']

rp_slp=[]
rp_yint=[]
for i in range(len(bc_mlu)):
    rp_slp.append(bcs)
    rp_yint.append(bcyi)
    
bc_mcu=st.ModelComparison(master_location,bc_mlu,bc_plu,bc_dlu)
[mn_ksn7,std_ksn7,stde_ksn7,mn_ksn7qr,std_ksn7qr,stde_ksn7qr,mn_ksn7qp,std_ksn7qp,stde_ksn7qp,mn_E7,std_E7,stde_E7,p25_E7,p75_E7]=bc_mcu.ksn_final_ts(rp_slp,rp_yint)
[mn_R7,mn_P7,cr7,sr7]=bc_mcu.rpv_final_ts(rp_slp,rp_yint)

bc_mcl=st.ModelComparison(master_location,bc_mll,bc_pll,bc_dll)
[mn_ksn8,std_ksn8,stde_ksn8,mn_ksn8qr,std_ksn8qr,stde_ksn8qr,mn_ksn8qp,std_ksn8qp,stde_ksn8qp,mn_E8,std_E8,stde_E8,p25_E8,p75_E8]=bc_mcl.ksn_final_ts(rp_slp,rp_yint)
[mn_R8,mn_P8,cr8,sr8]=bc_mcl.rpv_final_ts(rp_slp,rp_yint)

bc_mcuR=st.ModelComparison(master_location,bc_mluR,bc_pluR,bc_dluR)
[mn_ksn13,std_ksn13,stde_ksn13,mn_ksn13qr,std_ksn13qr,stde_ksn13qr,mn_ksn13qp,std_ksn13qp,stde_ksn13qp,mn_E13,std_E13,stde_E13,p25_E13,p75_E13]=bc_mcuR.ksn_final_ts(rp_slp,rp_yint)
[mn_R13,mn_P13,cr13,sr13]=bc_mcuR.rpv_final_ts(rp_slp,rp_yint)

bc_mclR=st.ModelComparison(master_location,bc_mllR,bc_pllR,bc_dllR)
[mn_ksn14,std_ksn14,stde_ksn14,mn_ksn14qr,std_ksn14qr,stde_ksn14qr,mn_ksn14qp,std_ksn14qp,stde_ksn14qp,mn_E14,std_E14,stde_E14,p25_E14,p75_E14]=bc_mclR.ksn_final_ts(rp_slp,rp_yint)
[mn_R14,mn_P14,cr14,sr14]=bc_mclR.rpv_final_ts(rp_slp,rp_yint)

# Package
model=gc_mlu+gc_mll+gc_mlu10+gc_mll10+a_mlu+a_mll+bc_mlu+bc_mll+bc_mluR+bc_mllR+gc_mluA+gc_mllA+gc_mluR+gc_mllR
group=gc_dlu+gc_dll+gc_dlu10+gc_dll10+a_dlu+a_dll+bc_dlu+bc_dll+bc_dluR+bc_dllR+gc_dluA+gc_dllA+gc_dluR+gc_dllR
mn_ksn=np.concatenate((mn_ksn1,mn_ksn2,mn_ksn3,mn_ksn4,mn_ksn5,mn_ksn6,mn_ksn7,mn_ksn8,mn_ksn13,mn_ksn14,mn_ksn9,mn_ksn10,mn_ksn11,mn_ksn12),axis=0)
std_ksn=np.concatenate((std_ksn1,std_ksn2,std_ksn3,std_ksn4,std_ksn5,std_ksn6,std_ksn7,std_ksn8,std_ksn13,std_ksn14,std_ksn9,std_ksn10,std_ksn11,std_ksn12),axis=0)
stde_ksn=np.concatenate((stde_ksn1,stde_ksn2,stde_ksn3,stde_ksn4,stde_ksn5,stde_ksn6,stde_ksn7,stde_ksn8,stde_ksn13,stde_ksn14,stde_ksn9,stde_ksn10,stde_ksn11,stde_ksn12),axis=0)
mn_ksnqr=np.concatenate((mn_ksn1qr,mn_ksn2qr,mn_ksn3qr,mn_ksn4qr,mn_ksn5qr,mn_ksn6qr,mn_ksn7qr,mn_ksn8qr,mn_ksn13qr,mn_ksn14qr,mn_ksn9qr,mn_ksn10qr,mn_ksn11qr,mn_ksn12qr),axis=0)
std_ksnqr=np.concatenate((std_ksn1qr,std_ksn2qr,std_ksn3qr,std_ksn4qr,std_ksn5qr,std_ksn6qr,std_ksn7qr,std_ksn8qr,std_ksn13qr,std_ksn14qr,std_ksn9qr,std_ksn10qr,std_ksn11qr,std_ksn12qr),axis=0)
stde_ksnqr=np.concatenate((stde_ksn1qr,stde_ksn2qr,stde_ksn3qr,stde_ksn4qr,stde_ksn5qr,stde_ksn6qr,stde_ksn7qr,stde_ksn8qr,stde_ksn13qr,stde_ksn14qr,stde_ksn9qr,stde_ksn10qr,stde_ksn11qr,stde_ksn12qr),axis=0)
mn_ksnqp=np.concatenate((mn_ksn1qp,mn_ksn2qp,mn_ksn3qp,mn_ksn4qp,mn_ksn5qp,mn_ksn6qp,mn_ksn7qp,mn_ksn8qp,mn_ksn13qp,mn_ksn14qp,mn_ksn9qp,mn_ksn10qp,mn_ksn11qp,mn_ksn12qp),axis=0)
std_ksnqp=np.concatenate((std_ksn1qp,std_ksn2qp,std_ksn3qp,std_ksn4qp,std_ksn5qp,std_ksn6qp,std_ksn7qp,std_ksn8qp,std_ksn13qp,std_ksn14qp,std_ksn9qp,std_ksn10qp,std_ksn11qp,std_ksn12qp),axis=0)
stde_ksnqp=np.concatenate((stde_ksn1qp,stde_ksn2qp,stde_ksn3qp,stde_ksn4qp,stde_ksn5qp,stde_ksn6qp,stde_ksn7qp,stde_ksn8qp,stde_ksn13qp,stde_ksn14qp,stde_ksn9qp,stde_ksn10qp,stde_ksn11qp,stde_ksn12qp),axis=0)
mn_E=np.concatenate((mn_E1,mn_E2,mn_E3,mn_E4,mn_E5,mn_E6,mn_E7,mn_E8,mn_E13,mn_E14,mn_E9,mn_E10,mn_E11,mn_E12),axis=0)
std_E=np.concatenate((std_E1,std_E2,std_E3,std_E4,std_E5,std_E6,std_E7,std_E8,std_E13,std_E14,std_E9,std_E10,std_E11,std_E12),axis=0)
stde_E=np.concatenate((stde_E1,stde_E2,stde_E3,stde_E4,stde_E5,stde_E6,stde_E7,stde_E8,stde_E13,stde_E14,stde_E9,stde_E10,stde_E11,stde_E12),axis=0)
p25_E=np.concatenate((p25_E1,p25_E2,p25_E3,p25_E4,p25_E5,p25_E6,p25_E7,p25_E8,p25_E13,p25_E14,p25_E9,p25_E10,p25_E11,p25_E12),axis=0)
p75_E=np.concatenate((p75_E1,p75_E2,p75_E3,p75_E4,p75_E5,p75_E6,p75_E7,p75_E8,p75_E13,p75_E14,p75_E9,p75_E10,p75_E11,p75_E12),axis=0)
mn_R=np.concatenate((mn_R1,mn_R2,mn_R3,mn_R4,mn_R5,mn_R6,mn_R7,mn_R8,mn_R13,mn_R14,mn_R9,mn_R10,mn_R11,mn_R12),axis=0)
mn_P=np.concatenate((mn_P1,mn_P2,mn_P3,mn_P4,mn_P5,mn_P6,mn_P7,mn_P8,mn_P13,mn_P14,mn_P9,mn_P10,mn_P11,mn_P12),axis=0)
cr=np.concatenate((cr1,cr2,cr3,cr4,cr5,cr6,cr7,cr8,cr13,cr14,cr9,cr10,cr11,cr12),axis=0)
sr=np.concatenate((sr1,sr2,sr3,sr4,sr5,sr6,sr7,sr8,sr13,sr14,sr9,sr10,sr11,sr12),axis=0)
U=np.tile(u_vec,14)
perc_E_std=(std_E/mn_E)*100
perc_U=((mn_E-U)/U)*100
int_qrt=p75_E-p25_E

df=pd.DataFrame(data={'Model':model,
                      'Group':group,
                      'mn_ksn':mn_ksn,
                      'std_ksn':std_ksn,
                      'stde_ksn':stde_ksn,
                      'mn_ksnqr':mn_ksnqr,
                      'std_ksnqr':std_ksnqr,
                      'stde_ksnqr':stde_ksnqr,
                      'mn_ksnqp':mn_ksnqp,
                      'std_ksnqp':std_ksnqp,
                      'stde_ksnqp':stde_ksnqp,
                      'mn_E':mn_E,
                      'std_E':std_E,
                      'stde_E':stde_E,
                      'p25_E':p25_E,
                      'p75_E':p75_E,
                      'mn_R':mn_R,
                      'mn_P':mn_P,
                      'cr':cr,
                      'sr':sr,
                      'U':U,
                      'perc_E_std':perc_E_std,
                      'perc_U':perc_U,
                      'int_qrt_E':int_qrt})

df.to_csv('model_final_stats.csv',index=False)




