#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 09:56:06 2023

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

[gcs,gcyi]=plot_and_fit_ref(ax1,df_global_s,bl1,br1,bb1,bt1,gc_col,False)

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


model_list2=['gc025u_10l','gc05u_10l','gc1u_10l','gc2u_10l','gc4u_10l','gc8u_10l']
prefix_list2=['ts_','ts_','ts_','ts_','ts_','ts_']
descript_list2=['GC STIM Unlinked; 0.25 mm/yr','GC STIM Unlinked; 0.5 mm/yr','GC STIM Unlinked; 1 mm/yr',
               'GC STIM Unlinked; 2 mm/yr','GC STIM Unlinked; 4 mm/yr','GC STIM Unlinked; 8 mm/yr']
group_list2=['GC Unlinked 10 km','GC Unlinked 10 km','GC Unlinked 10 km',
             'GC Unlinked 10 km','GC Unlinked 10 km','GC Unlinked 10 km']


model_list3=['gc025l_10l','gc05l_10l','gc1l_10l','gc2l_10l','gc4l_10l','gc8l_10l']
prefix_list3=['ts_','ts_','ts_','ts_','ts_','ts_']
descript_list3=['GC STIM Linked; 0.25 mm/yr','GC STIM Linked; 0.5 mm/yr','GC STIM Linked; 1 mm/yr',
               'GC STIM Linked; 2 mm/yr','GC STIM Linked; 4 mm/yr','GC STIM Linked; 8 mm/yr']
group_list3=['GC Linked 10 km','GC Linked 10 km','GC Linked 10 km',
             'GC Linked 10 km','GC Linked 10 km','GC Linked 10 km']

model_list=model_list1+model_list2+model_list3
prefix_list=prefix_list1+prefix_list2+prefix_list3
descript_list=descript_list1+descript_list2+descript_list3
group_list=group_list1+group_list2 +group_list3

mc=st.ModelComparison(master_location,model_list,prefix_list,descript_list)

col_list=['black','black','gray','gray']
shape_list=['o','s','o','s']
shape_filled_list=['filled','open','filled','open']
line_list=['-','--','-','--']
grp=['GC Unlinked','GC Linked','GC Unlinked 10 km','GC Linked 10 km']


rp_slp1=[]
rp_yint1=[]
for i in range(len(model_list1)):
    rp_slp1.append(gcs)
    rp_yint1.append(gcyi)
rp_slp2=[]
rp_yint2=[]
for i in range(len(model_list2)):
    rp_slp2.append(gcs)
    rp_yint2.append(gcyi)
rp_slp3=[]
rp_yint3=[]
for i in range(len(model_list3)):
    rp_slp3.append(gcs)
    rp_yint3.append(gcyi)
rp_slp=rp_slp1+rp_slp2+rp_slp3
rp_yint=rp_yint1+rp_yint2+rp_yint3
[fit_df,_,_]=mc.comp_final_ts(group_list,col_list,shape_list,shape_filled_list,line_list,grp,rp_slp,rp_yint)

# Generate supplemental figure
f1=plt.figure(figsize=(8,4))
f1.set_dpi(250)
e_vec=np.logspace(1,4,100)
gc_ero=pd.read_csv('/Users/aforte/Documents/GitHub/Caucasus_Erosion/data_tables/gc_ero_master_table.csv')

mn_ksn=gc_ero['mean_ksn'].to_numpy()
se_ksn=gc_ero['se_ksn'].to_numpy()
mn_E=gc_ero['St_E_rate_m_Myr'].to_numpy()
se_E=gc_ero['St_Ext_Unc'].to_numpy()
da=gc_ero['drainage_area'].to_numpy()

sc1=plt.scatter(mn_E,mn_ksn,s=30,c=np.log10(da),label='Observed Erosion Rates',zorder=1,cmap=cm.lajolla,edgecolors='k')
plt.errorbar(mn_E,mn_ksn,se_ksn,se_E,ecolor='gray',linestyle='',zorder=0,elinewidth=0.5)
plt.plot(e_vec,fit_df.loc[0,'C']*(e_vec/1e6)**fit_df.loc[0,'phi'],c='k',label=r'Unlinked 50 km (LogDA = 2.69 $km^{2}$)')
plt.plot(e_vec,fit_df.loc[1,'C']*(e_vec/1e6)**fit_df.loc[1,'phi'],c='k',linestyle='--',label=r'Linked 50 km (LogDA = 2.69 $km^{2}$)')
plt.plot(e_vec,fit_df.loc[2,'C']*(e_vec/1e6)**fit_df.loc[2,'phi'],c='gray',linestyle='-',label=r'Unlinked 10 km (LogDA = 1.41 $km^{2}$)')
plt.plot(e_vec,fit_df.loc[3,'C']*(e_vec/1e6)**fit_df.loc[3,'phi'],c='gray',linestyle='--',label=r'Linked 10 km (LogDA = 1.41 $km^{2}$)')
plt.xlabel('Erosion Rate [m/Myr]')
plt.ylabel(r'$k_{sn}$ [m]')
plt.legend(bbox_to_anchor= (1.3,0.99),loc='upper left')
plt.xscale('log')
ax1=plt.gca()
cbar1=plt.colorbar(sc1,ax=ax1)
cbar1.ax.set_ylabel(r'Log Drainage Area [$km^{2}$]')

plt.tight_layout()
plt.rcdefaults()

f1.savefig('figure_s7.pdf',dpi="figure")