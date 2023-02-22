#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  3 09:31:27 2023

@author: aforte
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colors
from cmcrameri import cm

import sys
sys.path.insert(0,'/Users/aforte/Documents/GitHub/snowmelt_orography')
import stimpy as st

def hypo_points(rp,zp,topo,snow):
    r=topo.loc[0,'param1']*rp**topo.loc[0,'param2']
    sn=topo.loc[1,'param1']*zp**topo.loc[1,'param2']
    if sn>1:
        sn=1
    if sn<=0.35:
        cr=snow.loc[1,'param1']*r**snow.loc[1,'param2']
    else:
        cr=snow.loc[15,'param1']*r + snow.loc[15,'param2']
    StimObj=st.StimSteady()
    [Ks,E,Rc]=StimObj.stim_range_dim(r,cr)
    return r,sn,cr,Ks,E


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

repo_location='/Users/aforte/Documents/GitHub/snowmelt_orography/'
topo=pd.read_csv(repo_location+'stimpy/topo_snow_mean_relationships.csv')
snow=pd.read_csv(repo_location+'stimpy/snowmelt_meanR_cR_relationships.csv')


# Establish color vector
col_vec=colors.Normalize(vmin=0,vmax=4)

# Establish runoff variability relationships
r_vec=np.logspace(-1,1,500)
lowsnow_cr=snow.loc[1,'param1']*r_vec**snow.loc[1,'param2']
highsnow_cr=snow.loc[15,'param1']*r_vec + snow.loc[15,'param2']

# Establish topography relationships
rlf=np.linspace(0,2300,500)
mxz=np.linspace(0,5000,500)
r=topo.loc[0,'param1']*rlf**topo.loc[0,'param2']
sn=topo.loc[1,'param1']*mxz**topo.loc[1,'param2']
sn[sn>1]=1

# Establish hypothetical points
rlf_p=np.linspace(250,2000,5)
mxz_p=np.linspace(1000,4500,5)
r_p=np.zeros(5)
sn_p=np.zeros(5)
cr_p=np.zeros(5)
ksn_p=[]
e_p=[]

for i in range(5):
    [r_p[i],sn_p[i],cr_p[i],ksn,e]=hypo_points(rlf_p[i],mxz_p[i],topo,snow)
    ksn_p.append(ksn)
    e_p.append(e)


f1=plt.figure(figsize=(8,3))
f1.set_dpi(250)

ax2=plt.subplot(1,3,1)
ax2.plot(rlf,r,c='k',label='Local Relief to Runoff')
for i in range(5):
    ax2.scatter(rlf_p[i],r_p[i],color=cm.roma(col_vec(i)),s=40,zorder=2,edgecolor='k')
    
plt.annotate('', xy=(2800, 7), xytext=(1200, 1), 
            arrowprops=dict(facecolor='k', shrink=0.),
            )
plt.text(1600,1.6,'Growing Topography',rotation=58.5)    
ax2.set_xlabel('Topography',labelpad=-10)
ax2.set_ylabel('Runoff',labelpad=-15)
ax2.set_xticks([0,5000])
ax2.set_xticklabels(['Low','High'])
ax2.set_yticks([0,10])
ax2.set_yticklabels(['Low','High'])
ax2.set_xlim((-100,5200))
ax2a=ax2.twinx()
ax2a.plot(mxz,sn,c='k',linestyle='--',label='Maximum Elevation to Snowmelt Fraction')
for i in range(5):
    ax2a.scatter(mxz_p[i],sn_p[i],color=cm.roma(col_vec(i)),s=40,zorder=2,edgecolor='k')
    
ax2a.plot([3000,6000],[0.35,0.35],c='k',linestyle=':') 
ax2a.text(3200,0.1,'Transition to\nHigh Snowmelt\nContribution',horizontalalignment='left')   
ax2a.set_ylabel('Snowmelt Fraction',labelpad=-5)
h1,l1=ax2.get_legend_handles_labels()
h2,l2=ax2a.get_legend_handles_labels()
ax2a.legend(h1+h2,l1+l2,bbox_to_anchor= (-0.1,-0.5),loc='lower left')
ax2a.set_yticks([0,1])
ax2a.set_yticklabels(['0','1'])
ax2.text(0.01, 0.99, 'A',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax2.transAxes,
        fontsize=12,fontweight='extra bold')

ax1=plt.subplot(1,3,2)
plt.plot(r_vec,lowsnow_cr,c='k',label='Low Snowmelt Contribution')
plt.plot(r_vec,highsnow_cr,c='k',linestyle='--',label='High Snowmelt Contribution')
for i in range(5):
    plt.scatter(r_p[i],cr_p[i],color=cm.roma(col_vec(i)),s=40,zorder=2,edgecolor='k')
plt.xlabel('Mean Runoff',labelpad=-10)
plt.ylabel('Shape Parameter',labelpad=-15)
plt.legend(bbox_to_anchor= (-0.1,-0.5),loc='lower left')
ax1.set_xticks([0,10])
ax1.set_xticklabels(['Low','High'])
ax1.set_yticks([0.2,2.8])
ax1.set_yticklabels(['High','Low'])
ax1.text(0.01, 0.99, 'B',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax1.transAxes,
        fontsize=12,fontweight='extra bold')


ax3=plt.subplot(1,3,3)
for i in range(5):
    ax3.plot(e_p[i]/1000,ksn_p[i],c=cm.roma(col_vec(i)),linewidth=2)
ax3.set_xlabel('Erosion Rate',labelpad=-10)
ax3.set_ylabel(r'$k_{sn}$',labelpad=-15)
plt.xlim((0,10))
plt.ylim((0,550))
ax3.set_xticks([0.25,9.75])
ax3.set_xticklabels(['Low','High'])
ax3.set_yticks([10,540])
ax3.set_yticklabels(['Low','High'])
ax3.text(0.01, 0.99, 'C',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax3.transAxes,
        fontsize=12,fontweight='extra bold')

plt.tight_layout(pad=0, w_pad=-2, h_pad=0.1)
plt.rcdefaults()





