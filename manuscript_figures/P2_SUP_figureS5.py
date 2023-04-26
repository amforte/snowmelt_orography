#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  8 11:17:03 2023

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

model_list1=['gc025u','gc05u','gc1u','gc2u','gc4u','gc8u',
            'gc025l','gc05l','gc1l','gc2l','gc4l','gc8l']

prefix_list1=['ts_','ts_','ts_','ts_','ts_','ts_','ts_','ts_','ts_','ts_','ts_','ts_']

descript_list1=['GC STIM Unlinked; 0.25 mm/yr','GC STIM Unlinked; 0.5 mm/yr','GC STIM Unlinked; 1 mm/yr',
               'GC STIM Unlinked; 2 mm/yr','GC STIM Unlinked; 4 mm/yr','GC STIM Unlinked; 8 mm/yr',
               'GC STIM Linked; 0.25 mm/yr','GC STIM Linked; 0.5 mm/yr','GC STIM Linked; 1 mm/yr',
               'GC STIM Linked; 2 mm/yr','GC STIM Linked; 4 mm/yr','GC STIM Linked; 8 mm/yr']

group_list1=['GC Unlinked','GC Unlinked','GC Unlinked','GC Unlinked','GC Unlinked','GC Unlinked',
             'GC Linked','GC Linked','GC Linked','GC Linked','GC Linked','GC Linked']

model_list2=['a025u','a05u','a1u','a2u','a4u','a8u',
            'a025l','a05l','a1l','a2l','a4l','a8l']

prefix_list2=['ts_','ts_','ts_','ts_','ts_','ts_','ts_','ts_','ts_','ts_','ts_','ts_']

descript_list2=['Alps STIM Unlinked; 0.25 mm/yr','Alps STIM Unlinked; 0.5 mm/yr','Alps STIM Unlinked; 1 mm/yr',
               'Alps STIM Unlinked; 2 mm/yr','Alps STIM Unlinked; 4 mm/yr','Alps STIM Unlinked; 8 mm/yr',
               'Alps STIM Linked; 0.25 mm/yr','Alps STIM Linked; 0.5 mm/yr','Alps STIM Linked; 1 mm/yr',
               'Alps STIM Linked; 2 mm/yr','Alps STIM Linked; 4 mm/yr','Alps STIM Linked; 8 mm/yr']

group_list2=['Alps Unlinked','Alps Unlinked','Alps Unlinked','Alps Unlinked','Alps Unlinked','Alps Unlinked',
             'Alps Linked','Alps Linked','Alps Linked','Alps Linked','Alps Linked','Alps Linked']

model_list3=['bc025u','bc05u','bc1u','bc2u','bc4u','bc8u',
            'bc025l','bc05l','bc1l','bc2l','bc4l','bc8l']

prefix_list3=['ts_','ts_','ts_','ts_','ts_','ts_','ts_','ts_','ts_','ts_','ts_','ts_']

descript_list3=['BC STIM Unlinked; 0.25 mm/yr','BC STIM Unlinked; 0.5 mm/yr','BC STIM Unlinked; 1 mm/yr',
               'BC STIM Unlinked; 2 mm/yr','BC STIM Unlinked; 4 mm/yr','BC STIM Unlinked; 8 mm/yr',
               'BC STIM Linked; 0.25 mm/yr','BC STIM Linked; 0.5 mm/yr','BC STIM Linked; 1 mm/yr',
               'BC STIM Linked; 2 mm/yr','BC STIM Linked; 4 mm/yr','BC STIM Linked; 8 mm/yr']

group_list3=['BC Unlinked','BC Unlinked','BC Unlinked','BC Unlinked','BC Unlinked','BC Unlinked',
             'BC Linked','BC Linked','BC Linked','BC Linked','BC Linked','BC Linked']

model_list=model_list1+model_list2+model_list3
prefix_list=prefix_list1+prefix_list2+prefix_list3
descript_list=descript_list1+descript_list2+descript_list3
group_list=group_list1+group_list2+group_list3

mc=st.ModelComparison(master_location,model_list,prefix_list,descript_list)

col_list=['black','black','royalblue','royalblue','orange','orange']
shape_list=['o','s','o','s','o','s']
shape_filled_list=['filled','open','filled','open','filled','open']
line_list=['-','--','-','--','-','--']
grp=['GC Unlinked','GC Linked','Alps Unlinked','Alps Linked','BC Unlinked','BC Linked']

rp_slp1=[]
rp_yint1=[]
for i in range(len(model_list1)):
    rp_slp1.append(gcs)
    rp_yint1.append(gcyi)
rp_slp2=[]
rp_yint2=[]
for i in range(len(model_list2)):
    rp_slp2.append(als)
    rp_yint2.append(alyi)
rp_slp3=[]
rp_yint3=[]
for i in range(len(model_list3)):
    rp_slp3.append(bcs)
    rp_yint3.append(bcyi)
rp_slp=rp_slp1+rp_slp2+rp_slp3
rp_yint=rp_yint1+rp_yint2+rp_yint3

[fit_df,f1,f2]=mc.comp_final_ts(group_list,col_list,shape_list,shape_filled_list,line_list,grp,rp_slp,rp_yint)
fit_df.to_csv('ksn_e_fit.csv',index=False)

f2.savefig('P2_SUP_figureS5.pdf',dpi="figure")