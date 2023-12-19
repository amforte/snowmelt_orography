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

def odr_fit(x,y):
    # Filter 0 values
    lx=np.log10(x[(x>0) & (y>0)])
    ly=np.log10(y[(x>0) & (y>0)])
    linmod=odr.Model(linear)
    
    fdlog=odr.Data(lx,ly)
    odrlog=odr.ODR(fdlog,linmod,beta0=[0.1,10])
    outlog=odrlog.run()
    logexp=outlog.beta[0]
    logcoeff=10**outlog.beta[1]

    return logcoeff,logexp


# Define colors
gc_col='black'
alps_col='royalblue'
bc_col='orange'

# Load final ts outupts
df=pd.read_csv('model_final_stats.csv')
dft=pd.read_csv('model_final_theta.csv')

group_list=['gcu','gcl','au','al','bcu','bcl']
mar_list=['o','s','o','s','o','s']
len_list=[50,50,50,50,50,50]
col_list=[gc_col,gc_col,alps_col,alps_col,bc_col,bc_col]
label=['GC Unlinked','GC Linked','Alps Unlinked','Alps Linked','BC Unlinked','BC Linked']
sty_list=['-','--','-','--','-','--']


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

E=np.array([250,500,1000,4000])
kl=np.array([100,200,300,400,600])

f1=plt.figure(1,figsize=(8,4.5))
ax1=plt.subplot(1,3,1)
ax1.set_xlim((100,10000))
ax1.set_ylim((100,700))
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlabel('Erosion Rate [m/Myr]')
ax1.set_ylabel(r'$k_{snQP}$ [m]')
ax1.set_xticks(E,E)
ax1.set_yticks(kl,kl)

ax2=plt.subplot(1,3,2)
ax2.set_xlim((100,10000))
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_ylim((100,700))
ax2.set_xlabel('Erosion Rate [m/Myr]')
ax2.set_ylabel(r'$k_{snQR}$ [m]')
ax2.set_xticks(E,E)
ax2.set_yticks(kl,kl)

ax3=plt.subplot(1,3,3)
ax3.set_xlim((100,10000))
ax3.set_xscale('log')
ax3.set_ylim((1,1.4))
ax3.set_xlabel('Erosion Rate [m/Myr]')
ax3.set_ylabel(r'$k_{snQP}$ / $k_{snQR}$')
ax3.set_xticks(E,E)


E_vec=np.logspace(-1,4,100)

for i in range(len(group_list)):
    idx=df['Group']==group_list[i]
    
    ksn=df.loc[idx,'mn_ksn'].to_numpy()
    E=df.loc[idx,'mn_E'].to_numpy()/1e6
    ksnqp=df.loc[idx,'mn_ksnqp'].to_numpy()
    ksnqr=df.loc[idx,'mn_ksnqr'].to_numpy()
    
    # ksns=df.loc[idx,'std_ksn'].to_numpy()
    # Es=df.loc[idx,'std_E'].to_numpy()
    # ksnqps=df.loc[idx,'std_ksnqp'].to_numpy()
    
    ksns=df.loc[idx,'stde_ksn'].to_numpy()
    Es=df.loc[idx,'stde_E'].to_numpy()
    ksnqps=df.loc[idx,'stde_ksnqp'].to_numpy()
    ksnqrs=df.loc[idx,'stde_ksnqr'].to_numpy()
    
    
    # Fit
    [K,n]=odr_fit(ksn,E)
    [Klp,nlp]=odr_fit(ksnqp,E)
    [Klr,nlr]=odr_fit(ksnqr,E)    
    [C,phi]=odr_fit(E*1e6,ksn)
    [Clp,philp]=odr_fit(E*1e6,ksnqp)
    [Clr,philr]=odr_fit(E*1e6,ksnqr)    
    
    if mar_list[i]=='o':
        ax1.scatter(E*1e6,ksnqp,c=col_list[i],marker=mar_list[i],zorder=2,label=label[i]+'; n = '+str(np.round(nlp,1)))
        ax2.scatter(E*1e6,ksnqr,c=col_list[i],marker=mar_list[i],zorder=2,label=label[i]+'; n = '+str(np.round(nlr,1)))
        ax3.scatter(E*1e6,ksnqp/ksnqr,c=col_list[i],marker=mar_list[i],zorder=2)
        
    else:
        ax1.scatter(E*1e6,ksnqp,c='w',edgecolor=col_list[i],marker=mar_list[i],zorder=2,label=label[i]+'; n = '+str(np.round(nlp,1)))
        ax2.scatter(E*1e6,ksnqr,c='w',edgecolor=col_list[i],marker=mar_list[i],zorder=2,label=label[i]+'; n = '+str(np.round(nlr,1)))
        ax3.scatter(E*1e6,ksnqp/ksnqr,c='w',edgecolor=col_list[i],marker=mar_list[i],zorder=2)
    
    ax1.plot(E_vec,Clp*E_vec**philp,c=col_list[i],linestyle=sty_list[i],zorder=1)
    ax2.plot(E_vec,Clr*E_vec**philr,c=col_list[i],linestyle=sty_list[i],zorder=1)
    
    ax1.errorbar(E*1e6,ksnqp,ksnqps,Es,ecolor=col_list[i],linestyle='',zorder=0,elinewidth=0.5)
    ax2.errorbar(E*1e6,ksnqr,ksnqrs,Es,ecolor=col_list[i],linestyle='',zorder=0,elinewidth=0.5)
    
ax1.legend(loc='upper right',bbox_to_anchor=[0.9,-0.25])
ax2.legend(loc='upper right',bbox_to_anchor=[0.9,-0.25])

ax1.text(0.01, 0.99, 'A',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax1.transAxes,
        fontsize=12,fontweight='extra bold')

ax2.text(0.01, 0.99, 'B',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax2.transAxes,
        fontsize=12,fontweight='extra bold')

ax3.text(0.01, 0.99, 'C',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax3.transAxes,
        fontsize=12,fontweight='extra bold')

plt.tight_layout()
plt.rcdefaults()

f1.savefig('P2_SUP_figureS4.pdf',dpi="figure")