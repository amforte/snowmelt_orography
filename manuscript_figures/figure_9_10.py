#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 13:30:40 2023

@author: aforte
"""

import sys
sys.path.insert(0,'/Users/aforte/Documents/GitHub/snowmelt_orography')
import stimpy as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
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

group_list=['gcu','gcl','gcu_10l','gcl_10l','au','al','bcu','bcl']
mar_list=['o','s','o','s','o','s','o','s']
len_list=[50,50,10,10,50,50,50,50]
col_list=[gc_col,gc_col,gc_col,gc_col,alps_col,alps_col,bc_col,bc_col]
label=['GC Unlinked','GC Linked','GC Unlinked 10 km','GC Linked 10 km',
       'Alps Unlinked','Alps Linked','BC Unlinked','BC Linked']
sty_list=['-','--',':','-.','-','--','-','--']


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

f1=plt.figure(figsize=(8,6.5))
f1.set_dpi(250)

gs1=gridspec.GridSpec(3,2)

ax1=f1.add_subplot(gs1[0:2,0])
ax1.set_xlabel('Erosion Rate [m/Myr]')
ax1.set_ylabel(r'$k_{sn}$ [m]')
ax1.set_xlim((-200,10200))

ax2=f1.add_subplot(gs1[0:2,1])
ax2.set_xlabel('Erosion Rate [m/Myr]')
ax2.set_ylabel(r'$k_{snQ}$ [m]')
ax2.set_xlim((-200,10200))

ax3=f1.add_subplot(gs1[2,0])
ax3.set_xlabel('Erosion Rate [m/Myr]')
ax3.set_ylabel(r'$\theta$ : $\chi$ - z optimization')
ax3.set_ylim((0.25,0.7))
ax3.axhline(0.5,c='k',linestyle=':')
ax3.set_xlim((-200,10200))

ax4=f1.add_subplot(gs1[2,1])
ax4.set_xlabel('Erosion Rate [m/Myr]')
ax4.set_ylabel(r'$\theta$ : $\chi_{Q}$  - z optimization')
ax4.set_ylim((0.25,0.7))
ax4.axhline(0.5,c='k',linestyle=':')
ax4.set_xlim((-200,10200))

f2=plt.figure(figsize=(8,5))
gs2=gridspec.GridSpec(2,2)

ax5=f2.add_subplot(gs2[0,0])
ax5.set_xlabel(r'1/n : $k_{sn}$')
ax5.set_ylabel('Shape Parameter')
ax5.set_xlim((0.01,0.45))

ax6=f2.add_subplot(gs2[0,1])
ax6.set_xlabel(r'1/n : $k_{snQ}$')
ax6.set_ylabel('Shape Parameter')
ax6.set_xlim((0.01,0.45))

ax7=f2.add_subplot(gs2[1,0])
ax7.set_xlabel(r'$K$')
ax7.set_ylabel('Runoff [mm/day]')
ax7.set_xscale('log')

ax8=f2.add_subplot(gs2[1,1])
ax8.set_xlabel(r'$K_{lp}$')
ax8.set_ylabel('Runoff [mm/day]')
ax8.set_xscale('log')



E_vec=np.logspace(-1,4,100)

for i in range(len(group_list)):
    idx=df['Group']==group_list[i]
    
    ksn=df.loc[idx,'mn_ksn'].to_numpy()
    E=df.loc[idx,'mn_E'].to_numpy()/1e6
    ksnqp=df.loc[idx,'mn_ksnqp'].to_numpy()
    
    ksns=df.loc[idx,'std_ksn'].to_numpy()
    Es=df.loc[idx,'std_E'].to_numpy()
    ksnqps=df.loc[idx,'std_ksnqp'].to_numpy()
    
    theta=dft.loc[idx,'theta_chi'].to_numpy()
    thetap=dft.loc[idx,'theta_chi_p'].to_numpy()

    # Fit
    [K,n]=odr_fit(ksn,E)
    [Klp,nlp]=odr_fit(ksnqp,E)
    [C,phi]=odr_fit(E*1e6,ksn)
    [Clp,philp]=odr_fit(E*1e6,ksnqp)
    
    cr=df.loc[idx,'cr'].to_numpy()
    mn_R=df.loc[idx,'mn_R'].to_numpy()
    mn_P=df.loc[idx,'mn_P'].to_numpy()
    
    if len_list[i]==50:
        ax1.scatter(E*1e6,ksn,c=col_list[i],marker=mar_list[i],zorder=2,label=label[i]+'; n = '+str(np.round(n,1)))
        ax2.scatter(E*1e6,ksnqp,c=col_list[i],marker=mar_list[i],zorder=2,label=label[i]+'; n = '+str(np.round(nlp,1)))
        if mar_list[i]=='o':
            ax5.scatter(1/n,np.mean(cr),c=col_list[i],marker=mar_list[i],zorder=1,label=label[i])
        else:
            ax5.scatter(1/n,np.mean(cr),c=col_list[i],marker=mar_list[i],zorder=1)
            
        if mar_list[i]=='s':
            ax6.scatter(1/nlp,np.mean(cr),c=col_list[i],marker=mar_list[i],zorder=1,label=label[i])
        else:
            ax6.scatter(1/nlp,np.mean(cr),c=col_list[i],marker=mar_list[i],zorder=1)            
        ax7.scatter(K,np.mean(mn_R),c=col_list[i],marker=mar_list[i],zorder=1)
        ax8.scatter(Klp,np.mean(mn_R),c=col_list[i],marker=mar_list[i],zorder=1)
        ax3.scatter(E*1e6,theta,c=col_list[i],marker=mar_list[i],zorder=1)
        ax4.scatter(E*1e6,thetap,c=col_list[i],marker=mar_list[i],zorder=1)
    else:
        ax1.scatter(E*1e6,ksn,c='w',edgecolor=col_list[i],marker=mar_list[i],zorder=2,label=label[i]+'; n = '+str(np.round(n,1)))
        ax2.scatter(E*1e6,ksnqp,c='w',edgecolor=col_list[i],marker=mar_list[i],zorder=2,label=label[i]+'; n = '+str(np.round(nlp,1)))
        if mar_list[i]=='o':
            ax5.scatter(1/n,np.mean(cr),c='w',edgecolor=col_list[i],marker=mar_list[i],zorder=1,label=label[i])
        else:
            ax5.scatter(1/n,np.mean(cr),c='w',edgecolor=col_list[i],marker=mar_list[i],zorder=1) 
        if mar_list[i]=='s':
            ax6.scatter(1/nlp,np.mean(cr),c='w',edgecolor=col_list[i],marker=mar_list[i],zorder=1,label=label[i])
        else:
            ax6.scatter(1/nlp,np.mean(cr),c='w',edgecolor=col_list[i],marker=mar_list[i],zorder=1)            
        ax7.scatter(K,np.mean(mn_R),c='w',edgecolor=col_list[i],marker=mar_list[i],zorder=1)
        ax8.scatter(Klp,np.mean(mn_R),c='w',edgecolor=col_list[i],marker=mar_list[i],zorder=1)
        ax3.scatter(E*1e6,theta,c='w',edgecolor=col_list[i],marker=mar_list[i],zorder=1)
        ax4.scatter(E*1e6,thetap,c='w',edgecolor=col_list[i],marker=mar_list[i],zorder=1)
    
    ax1.plot(E_vec,C*E_vec**phi,c=col_list[i],linestyle=sty_list[i],zorder=1)
    ax2.plot(E_vec,Clp*E_vec**philp,c=col_list[i],linestyle=sty_list[i],zorder=1)
    
    ax1.errorbar(E*1e6,ksn,ksns,Es,ecolor=col_list[i],linestyle='',zorder=0,elinewidth=0.5)
    ax2.errorbar(E*1e6,ksnqp,ksnqps,Es,ecolor=col_list[i],linestyle='',zorder=0,elinewidth=0.5)
    ax5.errorbar(1/n,np.mean(cr),np.std(cr),ecolor=col_list[i],linestyle='',zorder=0,elinewidth=0.5)
    ax6.errorbar(1/nlp,np.mean(cr),np.std(cr),ecolor=col_list[i],linestyle='',zorder=0,elinewidth=0.5)
    ax7.errorbar(K,np.mean(mn_R),np.std(mn_R),ecolor=col_list[i],linestyle='',zorder=0,elinewidth=0.5)
    ax8.errorbar(Klp,np.mean(mn_R),np.std(mn_R),ecolor=col_list[i],linestyle='',zorder=0,elinewidth=0.5)
    
    
ax1.legend(loc='best')
ax2.legend(loc='best')
ax5.legend(loc='best')
ax6.legend(loc='best')
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
ax4.text(0.01, 0.99, 'D',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax4.transAxes,
        fontsize=12,fontweight='extra bold')

ax5.text(0.01, 0.99, 'A',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax5.transAxes,
        fontsize=12,fontweight='extra bold')
ax6.text(0.01, 0.99, 'B',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax6.transAxes,
        fontsize=12,fontweight='extra bold')
ax7.text(0.95, 0.99, 'C',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax7.transAxes,
        fontsize=12,fontweight='extra bold')
ax8.text(0.95, 0.99, 'D',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax8.transAxes,
        fontsize=12,fontweight='extra bold')


f1.tight_layout()
f2.tight_layout()
plt.rcdefaults()

f1.savefig('figure_9.pdf',dpi="figure")
f2.savefig('figure_10.pdf',dpi="figure")
