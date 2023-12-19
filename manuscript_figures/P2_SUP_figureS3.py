#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 16 13:44:46 2023

@author: aforte
"""

import sys
sys.path.insert(0,'/Users/aforte/Documents/GitHub/snowmelt_orography')
import stimpy as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from cmcrameri import cm
from matplotlib import colors
from matplotlib import cm as cmm
import matplotlib.gridspec as gridspec


# master_location='/Volumes/Choruh/Data/snowmelt_project/'
master_data_location='/Users/aforte/Documents/Python/snowmelt/'
master_location=master_data_location+'model_outputs_v2/'

# Set number of timesteps to average across
num_avg=40


# Bins
ml1=['gc1u','gc1u_5b','gc1u_10b','gc1l']
pl1=['ts_','ts_','ts_','ts_']
dl1=['GC1U','GC1U-5B','GC1U-10B','GC1L']
sl1=[st.Stream(50000,25,dx=100,bin_size=2000),
     st.Stream(50000,25,dx=100,bin_size=5000),
     st.Stream(50000,25,dx=100,bin_size=10000),
     st.Stream(50000,25,dx=100,bin_size=2000)]
cl1=['black','olivedrab','limegreen','black']
st1=['-','--','-','--']

mc1=st.ModelComparison(master_location,ml1,pl1,dl1)
d1=mc1.output_result_dict(np.tile([np.inf],4),last_n_ts=num_avg)

# Relief
ml2=['gc1u','gc1u_2000r','gc1u_1500r']
pl2=['ts_','ts_','ts_',]
dl2=['GC1U','GC1U-2000R','GC1U-1500R']
sl2=[st.Stream(50000,25,dx=100,bin_size=2000),
     st.Stream(50000,25,dx=100,bin_size=2000),
     st.Stream(50000,25,dx=100,bin_size=2000)]
cl2=['black','tomato','darksalmon']
st2=['-','--','-']

mc2=st.ModelComparison(master_location,ml2,pl2,dl2)
d2=mc2.output_result_dict(np.tile([np.inf],3),last_n_ts=num_avg)

# Length
ml3=['gc1u_10l','gc1u_20l','gc1u_30l','gc1u_40l','gc1u','gc1u_100l']
pl3=['ts_','ts_','ts_','ts_','ts_','ts_']
dl3=['GC1U-10L','GC1U-20L','GC1U-30L','GC1U-40L','GC1U','GC1U-100L']
sl3=[st.Stream(10000,25,dx=100,bin_size=2000),
     st.Stream(20000,25,dx=100,bin_size=2000),
     st.Stream(30000,25,dx=100,bin_size=2000),
     st.Stream(40000,25,dx=100,bin_size=2000),
     st.Stream(50000,25,dx=100,bin_size=2000),
     st.Stream(100000,25,dx=100,bin_size=2000)]
cl3=['dodgerblue','lightsteelblue','deepskyblue','aqua','black','darkblue']
st3=['-','--','-','--','-','--']


mc3=st.ModelComparison(master_location,ml3,pl3,dl3)
d3=mc3.output_result_dict(np.tile([np.inf],6),last_n_ts=num_avg)

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

f1=plt.figure(1,figsize=(8,10),dpi=300,layout='constrained')
gs=gridspec.GridSpec(5,7,figure=f1)
# Elevation
ax1=f1.add_subplot(gs[0,0:2])
ax1.set_ylabel('Elevation [km]')
ax1.set_ylim((-0.1,3))
ax1.set_xlim((-1,51))
ax1.set_title('Bins')

ax2=f1.add_subplot(gs[0,2:4])
ax2.set_title('Relief')
ax2.set_ylim((-0.1,3))
ax2.set_xlim((-1,51))

ax3=f1.add_subplot(gs[0,4:6])
ax3.set_title('Length')
ax3.set_ylim((-0.1,3))
ax3.set_xlim((-1,101))


# ksn
ax4=f1.add_subplot(gs[1,0:2])
ax4.set_ylabel(r'$k_{sn}$ [m]')
ax4.set_xlim((-1,51))
ax4.set_ylim((100,550))

ax5=f1.add_subplot(gs[1,2:4])
ax5.set_xlim((-1,51))
ax5.set_ylim((100,550))

ax6=f1.add_subplot(gs[1,4:6])
ax6.set_xlim((-1,101))
ax6.set_ylim((100,550))


# Snow Fraction
ax7=f1.add_subplot(gs[2,0:2])
ax7.set_ylabel('Snow Fraction')
ax7.set_ylim((0,1))
ax7.set_xlim((-1,51))
ax7.axhline(0.35,c='gray',linestyle=':',zorder=0)

ax8=f1.add_subplot(gs[2,2:4])
ax8.set_ylim((0,1))
ax8.set_xlim((-1,51))
ax8.axhline(0.35,c='gray',linestyle=':',zorder=0)

ax9=f1.add_subplot(gs[2,4:6])
ax9.set_ylim((0,1))
ax9.set_xlim((-1,101))
ax9.axhline(0.35,c='gray',linestyle=':',zorder=0)

# Mean Runoff
ax10=f1.add_subplot(gs[3,0:2])
ax10.set_ylabel(r'$\bar{R}$ [mm/day]')
ax10.set_xlim((-1,51))
ax10.set_ylim((0,10))

ax11=f1.add_subplot(gs[3,2:4])
ax11.set_xlim((-1,51))
ax11.set_ylim((0,10))

ax12=f1.add_subplot(gs[3,4:6])
ax12.set_xlim((-1,101))
ax12.set_ylim((0,10))

# Variability
ax13=f1.add_subplot(gs[4,0:2])
ax13.set_ylabel(r'Variability $c_{R}$')
ax13.set_xlabel('Stream Distance [km]')
ax13.set_xlim((-1,51))
ax13.set_ylim((0.45,1.4))

ax14=f1.add_subplot(gs[4,2:4])
ax14.set_xlabel('Stream Distance [km]')
ax14.set_xlim((-1,51))
ax14.set_ylim((0.45,1.4))

ax15=f1.add_subplot(gs[4,4:6])
ax15.set_xlabel('Stream Distance [km]')
ax15.set_xlim((-1,101))
ax15.set_ylim((0.45,1.4))

for i in range(len(ml1)):
    dI=d1[i]
    sObj=sl1[i]

    # Extract shared x
    x=dI['x']/1000
    xc=dI['x_center']/1000
    
    dIz=np.mean(dI['zout']/1000,axis=0)
    dIz=np.mean(dI['zout']/1000,axis=0)
    dIsn=np.mean(dI['snp'],axis=0)
    dIsn=np.mean(dI['snp'],axis=0)  
    dImr=np.mean(dI['mrout'],axis=0)
    dImr=np.mean(dI['mrout'],axis=0)    
    dIcr=np.mean(dI['crout'],axis=0)
    dIcr=np.mean(dI['crout'],axis=0)

    dIksn=np.mean(np.diff(dI['zout'],axis=1)/np.diff(dI['chi']),axis=0)
    dIksn=np.concatenate(([dIksn[0]],dIksn),axis=0)
    dIksn=np.bincount(sObj.ix,dIksn,sObj.num_bins)[1:sObj.num_bins+1]/np.bincount(sObj.ix,None,sObj.num_bins)[1:sObj.num_bins+1]
    
    dIksn=np.mean(np.diff(dI['zout'],axis=1)/np.diff(dI['chi']),axis=0)
    dIksn=np.concatenate(([dIksn[0]],dIksn),axis=0)
    dIksn=np.bincount(sObj.ix,dIksn,sObj.num_bins)[1:sObj.num_bins+1]/np.bincount(sObj.ix,None,sObj.num_bins)[1:sObj.num_bins+1]
    

    ax1.plot(x,dIz,linestyle=st1[i],c=cl1[i])

    ax4.plot(xc,dIksn,linestyle=st1[i],c=cl1[i])
    
    ax7.plot(xc,dIsn,linestyle=st1[i],c=cl1[i])
    
    ax10.plot(xc,dImr,linestyle=st1[i],c=cl1[i])
    
    ax13.plot(xc,dIcr,linestyle=st1[i],c=cl1[i],label=dl1[i])
    ax13.legend(loc='upper left')

for i in range(len(ml2)):
    dI=d2[i]
    sObj=sl2[i]

    # Extract shared x
    x=dI['x']/1000
    xc=dI['x_center']/1000
    
    dIz=np.mean(dI['zout']/1000,axis=0)
    dIz=np.mean(dI['zout']/1000,axis=0)
    dIsn=np.mean(dI['snp'],axis=0)
    dIsn=np.mean(dI['snp'],axis=0)  
    dImr=np.mean(dI['mrout'],axis=0)
    dImr=np.mean(dI['mrout'],axis=0)    
    dIcr=np.mean(dI['crout'],axis=0)
    dIcr=np.mean(dI['crout'],axis=0)

    dIksn=np.mean(np.diff(dI['zout'],axis=1)/np.diff(dI['chi']),axis=0)
    dIksn=np.concatenate(([dIksn[0]],dIksn),axis=0)
    dIksn=np.bincount(sObj.ix,dIksn,sObj.num_bins)[1:sObj.num_bins+1]/np.bincount(sObj.ix,None,sObj.num_bins)[1:sObj.num_bins+1]
    
    dIksn=np.mean(np.diff(dI['zout'],axis=1)/np.diff(dI['chi']),axis=0)
    dIksn=np.concatenate(([dIksn[0]],dIksn),axis=0)
    dIksn=np.bincount(sObj.ix,dIksn,sObj.num_bins)[1:sObj.num_bins+1]/np.bincount(sObj.ix,None,sObj.num_bins)[1:sObj.num_bins+1]
    

    ax2.plot(x,dIz,linestyle=st2[i],c=cl2[i])

    ax5.plot(xc,dIksn,linestyle=st2[i],c=cl2[i])
    
    ax8.plot(xc,dIsn,linestyle=st2[i],c=cl2[i])
    
    ax11.plot(xc,dImr,linestyle=st2[i],c=cl2[i])
    
    ax14.plot(xc,dIcr,linestyle=st2[i],c=cl2[i],label=dl2[i])
    ax14.legend(loc='upper left')
    
for i in range(len(ml3)):
    dI=d3[i]
    sObj=sl3[i]

    # Extract shared x
    x=dI['x']/1000
    xc=dI['x_center']/1000
    
    dIz=np.mean(dI['zout']/1000,axis=0)
    dIz=np.mean(dI['zout']/1000,axis=0)
    dIsn=np.mean(dI['snp'],axis=0)
    dIsn=np.mean(dI['snp'],axis=0)  
    dImr=np.mean(dI['mrout'],axis=0)
    dImr=np.mean(dI['mrout'],axis=0)    
    dIcr=np.mean(dI['crout'],axis=0)
    dIcr=np.mean(dI['crout'],axis=0)

    dIksn=np.mean(np.diff(dI['zout'],axis=1)/np.diff(dI['chi']),axis=0)
    dIksn=np.concatenate(([dIksn[0]],dIksn),axis=0)
    dIksn=np.bincount(sObj.ix,dIksn,sObj.num_bins)[1:sObj.num_bins+1]/np.bincount(sObj.ix,None,sObj.num_bins)[1:sObj.num_bins+1]
    
    dIksn=np.mean(np.diff(dI['zout'],axis=1)/np.diff(dI['chi']),axis=0)
    dIksn=np.concatenate(([dIksn[0]],dIksn),axis=0)
    dIksn=np.bincount(sObj.ix,dIksn,sObj.num_bins)[1:sObj.num_bins+1]/np.bincount(sObj.ix,None,sObj.num_bins)[1:sObj.num_bins+1]
    

    ax3.plot(x,dIz,linestyle=st3[i],c=cl3[i])
    
    ax6.plot(xc,dIksn,linestyle=st3[i],c=cl3[i])
    
    ax9.plot(xc,dIsn,linestyle=st3[i],c=cl3[i])
    
    ax12.plot(xc,dImr,linestyle=st3[i],c=cl3[i])
    
    ax15.plot(xc,dIcr,linestyle=st3[i],c=cl3[i],label=dl3[i])
    ax15.legend(loc='upper left')
    
plt.rcdefaults()
f1.savefig('P2_SUP_figureS3.pdf',dpi='figure')