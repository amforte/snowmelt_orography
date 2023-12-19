#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 11 12:31:41 2023

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

master_dir='/Users/aforte/Documents/Python/snowmelt'

## Determine runoff to precip relationships
## Load global
df_global=pd.read_csv(master_dir+'/wrr2_derived_data_v4.csv')
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
ax1=plt.subplot(1,1,1)

[gc_slp,gc_yint]=plot_and_fit_ref(ax1,df_global_s,bl1,br1,bb1,bt1,gc_col,False)

master_location=master_dir+'/model_outputs_v2/'

# Figure S3
model_list=['gc1u','gc1u_5b','gc1u_10b','gc1l']
prefix_list=['ts_','ts_','ts_','ts_']
descript_list=['GC1U','GC1U-5B','GC1U-10B','GC1L']

mc1=st.ModelComparison(master_location,model_list,prefix_list,descript_list)
[mn_ksn1,std_ksn1,stde_ksn1,mn_ksnqr1,std_ksnqr1,stde_ksnqr1,mn_ksnqp1,std_ksnqp1,stde_ksnqp1,mn_E1,std_E1,stde_E1,p25_E1,p75_E1]=mc1.ksn_final_ts([gc_slp,gc_slp,gc_slp,gc_slp],[gc_yint,gc_yint,gc_yint,gc_yint])

# Figure S4
model_list=['gc1u','gc1u_2000r','gc1u_1500r']
prefix_list=['ts_','ts_','ts_']
descript_list=['GC1U','GC1U-2000R','GC1U-1500R']

mc2=st.ModelComparison(master_location,model_list,prefix_list,descript_list)
[mn_ksn2,std_ksn2,stde_ksn2,mn_ksnqr2,std_ksnqr2,stde_ksnqr2,mn_ksnqp2,std_ksnqp2,stde_ksnqp2,mn_E2,std_E2,stde_E2,p25_E2,p75_E2]=mc2.ksn_final_ts([gc_slp,gc_slp,gc_slp],[gc_yint,gc_yint,gc_yint])

# Figure S5
model_list=['gc1u','gc1u_40l','gc1u_30l','gc1u_20l','gc1u_10l','gc1u_100l']
prefix_list=['ts_','ts_','ts_','ts_','ts_','ts_']
descript_list=['GC1U','GC1U-40L','GC1U-30L','GC1U-20L','GC1U-10L','GC1U-1OOL']

mc3=st.ModelComparison(master_location,model_list,prefix_list,descript_list)
[mn_ksn3,std_ksn3,stde_ksn3,mn_ksnqr3,std_ksnqr3,stde_ksnqr3,mn_ksnqp3,std_ksnqp3,stde_ksnqp3,mn_E3,std_E3,stde_E3,p25_E3,p75_E3]=mc3.ksn_final_ts([gc_slp,gc_slp,gc_slp,gc_slp,gc_slp,gc_slp],[gc_yint,gc_yint,gc_yint,gc_yint,gc_yint,gc_yint])

# Figure S5
model_list=['gc1u_5b','gc1u_10l','gc1u_10b','gc1u_10l_1b'] 
prefix_list=['ts_','ts_','ts_','ts_']
descript_list=['GC1U-5B','GC1U-10L','GC1U-10B','GC1U-10L-1B']

mc4=st.ModelComparison(master_location,model_list,prefix_list,descript_list)
[mn_ksn4,std_ksn4,stde_ksn4,mn_ksnqr4,std_ksnqr4,stde_ksnqr4,mn_ksnqp4,std_ksnqp4,stde_ksnqp4,mn_E4,std_E4,stde_E4,p25_E4,p75_E4]=mc4.ksn_final_ts([gc_slp,gc_slp,gc_slp,gc_slp,gc_slp],[gc_yint,gc_yint,gc_yint,gc_yint,gc_yint])


# Area bins to regular bins
model_list=['gc025u','gc05u','gc1u','gc2u','gc4u','gc8u'] 
prefix_list=['ts_','ts_','ts_','ts_','ts_','ts_']
descript_list=['GC025U','GC05U','GC1U','GC2U','GC4U','GC8U']

mc5=st.ModelComparison(master_location,model_list,prefix_list,descript_list)
[mn_ksn5,std_ksn5,stde_ksn5,mn_ksnqr5,std_ksnqr5,stde_ksnqr5,mn_ksnqp5,std_ksnqp5,stde_ksnqp5,mn_E5,std_E5,stde_E5,p25_E5,p75_E5]=mc5.ksn_final_ts([gc_slp,gc_slp,gc_slp,gc_slp,gc_slp,gc_slp,gc_slp],[gc_yint,gc_yint,gc_yint,gc_yint,gc_yint,gc_yint,gc_yint])

# Area bins to regular bins
model_list=['gc025l','gc05l','gc1l','gc2l','gc4l','gc8l'] 
prefix_list=['ts_','ts_','ts_','ts_','ts_','ts_']
descript_list=['GC025L','GC05L','GC1L','GC2L','GC4L','GC8L']

mc6=st.ModelComparison(master_location,model_list,prefix_list,descript_list)
[mn_ksn6,std_ksn6,stde_ksn6,mn_ksnqr6,std_ksnqr6,stde_ksnqr6,mn_ksnqp6,std_ksnqp6,stde_ksnqp6,mn_E6,std_E6,stde_E6,p25_E6,p75_E6]=mc6.ksn_final_ts([gc_slp,gc_slp,gc_slp,gc_slp,gc_slp,gc_slp,gc_slp],[gc_yint,gc_yint,gc_yint,gc_yint,gc_yint,gc_yint,gc_yint])

model_list=['gc025u_Area','gc05u_Area','gc1u_Area','gc2u_Area','gc4u_Area','gc8u_Area'] 
prefix_list=['ts_','ts_','ts_','ts_','ts_','ts_']
descript_list=['GC025U-AB','GC05U-AB','GC1U-AB','GC2U-AB','GC4U-AB','GC8U-AB']

mc7=st.ModelComparison(master_location,model_list,prefix_list,descript_list)
[mn_ksn7,std_ksn7,stde_ksn7,mn_ksnqr7,std_ksnqr7,stde_ksnqr7,mn_ksnqp7,std_ksnqp7,stde_ksnqp7,mn_E7,std_E7,stde_E7,p25_E7,p75_E7]=mc7.ksn_final_ts([gc_slp,gc_slp,gc_slp,gc_slp,gc_slp,gc_slp,gc_slp],[gc_yint,gc_yint,gc_yint,gc_yint,gc_yint,gc_yint,gc_yint])

# Area bins to regular bins
model_list=['gc025l_Area','gc05l_Area','gc1l_Area','gc2l_Area','gc4l_Area','gc8l_Area'] 
prefix_list=['ts_','ts_','ts_','ts_','ts_','ts_']
descript_list=['GC025L-AB','GC05L-AB','GC1L-AB','GC2L-AB','GC4L-AB','GC8L-AB']

mc8=st.ModelComparison(master_location,model_list,prefix_list,descript_list)
[mn_ksn8,std_ksn8,stde_ksn8,mn_ksnqr8,std_ksnqr8,stde_ksnqr8,mn_ksnqp8,std_ksnqp8,stde_ksnqp8,mn_E8,std_E8,stde_E8,p25_E8,p75_E8]=mc8.ksn_final_ts([gc_slp,gc_slp,gc_slp,gc_slp,gc_slp,gc_slp,gc_slp],[gc_yint,gc_yint,gc_yint,gc_yint,gc_yint,gc_yint,gc_yint])



f1=plt.figure(figsize=(6,8))
f1.set_dpi(250)
ax1=f1.add_subplot(3,2,1)
ax1.set_xlabel('Bin Size [km]')
ax1.set_ylabel(r'$k_{sn}$ [m]')
ax1.errorbar([2],mn_ksn1[0],yerr=std_ksn1[0],ecolor='k',elinewidth=0.5,linestyle='',zorder=0)
ax1.errorbar([5,10],mn_ksn1[1:3],yerr=std_ksn1[1:3],ecolor='gray',elinewidth=0.5,linestyle='',zorder=0)
ax1.errorbar([2],mn_ksn1[3],yerr=std_ksn1[3],ecolor='k',elinewidth=0.5,linestyle='',zorder=0)
ax1.scatter([2],mn_ksn1[0],c='k',s=50,label='GC1U',zorder=1)

# ax1.scatter([2+0.25],mn_ksn1[4],c='k',s=50,marker='^',zorder=1,label='GC1U-AB')
# ax1.scatter([2+0.25],mn_ksn1[5],c='k',s=50,marker='v',zorder=1,label='GC1L-AB')
# ax1.errorbar([2+0.25],mn_ksn1[4],yerr=std_ksn1[4],ecolor='k',elinewidth=0.5,linestyle='',zorder=0)
# ax1.errorbar([2_0.25],mn_ksn1[5],yerr=std_ksn1[5],ecolor='k',elinewidth=0.5,linestyle='',zorder=0)

ax1.scatter([5,10],mn_ksn1[1:3],c='gray',s=50,zorder=1)
ax1.scatter([2],mn_ksn1[3],c='k',marker='s',s=50,label='GC1L',zorder=1)
ax1.set_ylim((200,455))
ax1.set_xticks([2,5,10])
ax1.set_xlim((1.5,10.5))
ax1.legend(loc='upper right')
ax1.text(0.01, 0.99, 'A',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax1.transAxes,
        fontsize=12,fontweight='extra bold') 

ax1a=ax1.twiny()
ax1a.set_xlim((1.5,10.5))
ax1a.set_xticks([2,5,10])
ax1a.set_xticklabels(['25','10','5'])
ax1a.set_xlabel('Number of Bins')

ax2=f1.add_subplot(3,2,2)
ax2.set_xlabel('Imposed Maximum Relief [km]')
ax2.set_ylabel(r'$k_{sn}$ [m]')
ax2.errorbar([2.5],mn_ksn2[0],yerr=std_ksn2[0],ecolor='k',elinewidth=0.5,linestyle='',zorder=0)
ax2.errorbar([2,1.5],mn_ksn2[1:],yerr=std_ksn2[1:],ecolor='gray',elinewidth=0.5,linestyle='',zorder=0)
ax2.scatter([2.5],mn_ksn2[0],c='k',s=50,zorder=1)
ax2.scatter([2,1.5],mn_ksn2[1:],c='gray',s=50,zorder=1)
ax2.set_xticks([1.5,2,2.5])
ax2.set_ylim((200,455))
ax2.text(0.01, 0.99, 'B',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax2.transAxes,
        fontsize=12,fontweight='extra bold') 

ax3=f1.add_subplot(3,2,3)
ax3.set_xlabel('Profile Length [km]')
ax3.set_ylabel(r'$k_{sn}$ [m]')
ax3.errorbar([50],mn_ksn3[0],yerr=std_ksn3[0],ecolor='k',elinewidth=0.5,linestyle='')
ax3.errorbar([40,30,20,10,100],mn_ksn3[1:],yerr=std_ksn3[1:],ecolor='gray',elinewidth=0.5,linestyle='')
ax3.scatter([50],mn_ksn3[0],c='k',s=50)
ax3.scatter([40,30,20,10,100],mn_ksn3[1:],c='gray',s=50)
ax3.set_ylim((200,455))
ax3.set_xlim((5,105))
ax3.set_xticks([10,20,30,40,50,100])
ax3.text(0.01, 0.99, 'C',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax3.transAxes,
        fontsize=12,fontweight='extra bold')

ax3a=ax3.twiny()
ax3a.set_xlim((5/2,105/2))
ax3a.set_xticks([5,10,15,20,25,50])
ax3a.set_xlabel('Number of Bins')
 

# ['gc1u_5b','gc1u_10l','gc1u_10b','gc1u_10l_1b'] 

ax4=f1.add_subplot(3,2,4)
ax4.set_xlabel('Number of Bins')
ax4.set_ylabel(r'$k_{sn}$ [m]')
ax4.errorbar([10.2],mn_ksn4[0],yerr=std_ksn4[0],ecolor='gray',elinewidth=0.5,linestyle='')
ax4.scatter([10.2],mn_ksn4[0],c='gray',s=50,marker='o',label='GC1U-5B')
ax4.errorbar([4.8],mn_ksn4[1],yerr=std_ksn4[1],ecolor='k',elinewidth=0.5,linestyle='',zorder=0)
ax4.scatter([4.8],mn_ksn4[1],c='gray',s=50,marker='s',label='GC1U-10L',edgecolor='k')
ax4.errorbar([5.2],mn_ksn4[2],yerr=std_ksn4[2],ecolor='k',elinewidth=0.5,linestyle='',zorder=0)
ax4.scatter([5.2],mn_ksn4[2],c='gray',s=50,marker='o',label='GC1U-10B',edgecolor='k')
ax4.errorbar([9.8],mn_ksn4[3],yerr=std_ksn4[3],ecolor='gray',elinewidth=0.5,linestyle='')
ax4.scatter([9.8],mn_ksn4[3],c='gray',s=50,marker='s',label='GC1U-10L-1B')

ax4.set_ylim((200,455))
ax4.set_xlim((3,12))
ax4.set_xticks([5,10])
ax4.legend(loc='upper left')
ax4.text(0.92, 0.99, 'D',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax4.transAxes,
        fontsize=12,fontweight='extra bold') 


ax5=f1.add_subplot(3,2,5)
ax5.set_xlabel('Uplift Rate [m/Myr]')
ax5.set_ylabel(r'$k_{sn}$ [m]')
ax5.set_xscale('log')
U=np.array([250,500,1000,2000,4000,8000])
ax5.set_ylim((100,550))
ax5.set_xticks(ticks=U,labels=['250','500','1000','2000','4000','8000'])
ax5.errorbar(U,mn_ksn5,yerr=std_ksn5,ecolor='k',elinewidth=0.5,linestyle='')
ax5.scatter(U,mn_ksn5,c='k',s=50,marker='o',label='Unlinked')
ax5.errorbar(U,mn_ksn6,yerr=std_ksn6,ecolor='k',elinewidth=0.5,linestyle='')
ax5.scatter(U,mn_ksn6,c='k',s=50,marker='s',label='Linked')

ax5.errorbar(U,mn_ksn7,yerr=std_ksn7,ecolor='gray',elinewidth=0.5,linestyle='')
ax5.scatter(U,mn_ksn7,c='gray',s=50,marker='^',label='Unlinked-AB')
ax5.errorbar(U,mn_ksn8,yerr=std_ksn8,ecolor='gray',elinewidth=0.5,linestyle='')
ax5.scatter(U,mn_ksn8,c='gray',s=50,marker='v',label='Linked-AB')
ax5.legend(loc='upper right', bbox_to_anchor=(1.6, 1))
ax5.text(0.01, 0.99, 'E',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax5.transAxes,
        fontsize=12,fontweight='extra bold') 

plt.tight_layout()


# Generate representative streams
# s1=st.Stream(50000,25,dx=100,bin_size=2000) # GC1U
# s2=st.Stream(50000,25,dx=100,bin_size=1.95e7,bin_dim='A') #GC1U-Area
# s3=st.Stream(100000,25,dx=100,bin_size=2000) #GC1U-100L
# s4=st.Stream(40000,25,dx=100,bin_size=2000) #GC1U-40L
# s5=st.Stream(30000,25,dx=100,bin_size=2000) #GC1U-30L
# s6=st.Stream(20000,25,dx=100,bin_size=2000) #GC1U-20L
# s7=st.Stream(10000,25,dx=100,bin_size=2000) #GC1U-10L
# s8=st.Stream(10000,25,dx=100,bin_size=1000) #GC1U-10L-1B
# s9=st.Stream(50000,25,dx=100,bin_size=10000) #GC1U-10B
# s10=st.Stream(50000,25,dx=100,bin_size=5000) # GC1U-5B


# f2=plt.figure(figsize=(6,6))
# ax1f2=f2.add_subplot(1,1,1)
# sc1=ax1f2.scatter(np.mean(s1.A_cents),mn_ksn1[0],c=2000,s=50,label='25 bins',vmin=1000,vmax=10000,cmap=cm.vik)
# # ax1f2.errorbar(np.mean(s1.A_cents),mn_ksn1[0],xerr=np.std(s1.A_cents),yerr=std_ksn1[0],ecolor='k',elinewidth=0.5,linestyle='')

# ax1f2.scatter(np.mean(s2.A_cents),mn_ksn1[4],c='gray',s=50,marker='^',label='25 area bins')
# # ax1f2.errorbar(np.mean(s2.A_cents),mn_ksn1[4],xerr=np.std(s2.A_cents),yerr=std_ksn1[4],ecolor='gray',elinewidth=0.5,linestyle='',zorder=0)

# ax1f2.scatter(np.mean(s3.A_cents),mn_ksn3[5],c=2000,s=100,vmin=1000,vmax=10000,marker='p',label='50 bins',cmap=cm.vik)
# ax1f2.scatter(np.mean(s4.A_cents),mn_ksn3[1],c=2000,s=40,vmin=1000,vmax=10000,marker='d',label='20 bins',cmap=cm.vik)
# ax1f2.scatter(np.mean(s5.A_cents),mn_ksn3[2],c=2000,s=30,vmin=1000,vmax=10000,marker='D',label='15 bins',cmap=cm.vik)
# ax1f2.scatter(np.mean(s6.A_cents),mn_ksn3[3],c=2000,s=20,vmin=1000,vmax=10000,marker='s',label='10 bins',cmap=cm.vik)
# ax1f2.scatter(np.mean(s7.A_cents),mn_ksn3[4],c=2000,s=10,vmin=1000,vmax=10000,marker='*',label='5 bins',cmap=cm.vik)
# ax1f2.scatter(np.mean(s8.A_cents),mn_ksn4[3],c=1000,s=10,vmin=1000,vmax=10000,marker='s',cmap=cm.vik)
# ax1f2.scatter(np.mean(s9.A_cents),mn_ksn4[2],c=10000,s=50,vmin=1000,vmax=10000,marker='*',cmap=cm.vik)
# ax1f2.scatter(np.mean(s10.A_cents),mn_ksn4[0],c=5000,s=50,vmin=1000,vmax=10000,marker='s',cmap=cm.vik)

# ax1f2.plot([np.mean(s1.A_cents),np.mean(s1.A_cents)],[mn_ksn1[0],mn_ksn1[4]],c='k',linestyle=':',zorder=0)
# ax1f2.plot([np.mean(s7.A_cents),np.mean(s9.A_cents)],[mn_ksn3[4],mn_ksn4[2]],c='k',linestyle=':',zorder=0)
# ax1f2.plot([np.mean(s8.A_cents),np.mean(s6.A_cents),np.mean(s10.A_cents)],[mn_ksn4[3],mn_ksn3[3],mn_ksn4[0]],c='k',linestyle=':',zorder=0)


# cbar=plt.colorbar(sc1,ax=ax1f2)
# cbar.ax.set_ylabel('Bin Length [m]')

# ax1f2.set_xscale('log')
# ax1f2.set_xlabel('Bin Drainage Area [m$^{2}$]')
# ax1f2.set_ylabel('Mean k$_{sn}$ [m]')
# plt.legend(loc='upper left')


# f3=plt.figure(figsize=(6,6))
# ax1f3=f3.add_subplot(1,1,1)
# ax1f3.scatter(np.mean(s1.A_cents)/50000,mn_ksn1[0],c=2000,s=50,label='GC1U',vmin=1000,vmax=10000)
# # ax1f2.errorbar(np.mean(s1.A_cents),mn_ksn1[0],xerr=np.std(s1.A_cents),yerr=std_ksn1[0],ecolor='k',elinewidth=0.5,linestyle='')

# ax1f3.scatter(np.mean(s2.A_cents)/50000,mn_ksn1[4],c='gray',s=50,marker='^',label='GC1U-Area')
# # ax1f2.errorbar(np.mean(s2.A_cents),mn_ksn1[4],xerr=np.std(s2.A_cents),yerr=std_ksn1[4],ecolor='gray',elinewidth=0.5,linestyle='',zorder=0)

# ax1f3.scatter(np.mean(s3.A_cents)/100000,mn_ksn3[5],c=2000,s=100,label='GC1U-100L',vmin=1000,vmax=10000)
# ax1f3.scatter(np.mean(s4.A_cents)/40000,mn_ksn3[1],c=2000,s=40,label='GC1U-40L',vmin=1000,vmax=10000)
# ax1f3.scatter(np.mean(s5.A_cents)/30000,mn_ksn3[2],c=2000,s=30,label='GC1U-30L',vmin=1000,vmax=10000)
# ax1f3.scatter(np.mean(s6.A_cents)/20000,mn_ksn3[3],c=2000,s=20,label='GC1U-20L',vmin=1000,vmax=10000)
# ax1f3.scatter(np.mean(s7.A_cents)/10000,mn_ksn3[4],c=2000,s=10,label='GC1U-10L',vmin=1000,vmax=10000)
# ax1f3.scatter(np.mean(s8.A_cents)/10000,mn_ksn4[3],c=1000,s=10,label='GC1U-10L-1B',vmin=1000,vmax=10000)
# ax1f3.scatter(np.mean(s9.A_cents)/50000,mn_ksn4[2],c=10000,s=50,label='GC1U-10B',vmin=1000,vmax=10000)
# ax1f3.scatter(np.mean(s10.A_cents)/50000,mn_ksn4[0],c=5000,s=50,label='GC1U-5B',vmin=1000,vmax=10000)

# ax1f3.set_xscale('log')
# ax1f3.set_xlabel('Bin Drainage Area [m$^{2}$] / Stream Length [m]')
# ax1f3.set_ylabel('Mean k$_{sn}$ [m]')
# plt.tight_layout()

plt.rcdefaults()

f1.savefig('P2_figure7.pdf',dpi="figure")
# f2.savefig('P2_figure7_alt.pdf')





