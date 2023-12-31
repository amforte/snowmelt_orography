#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 26 16:17:17 2023

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


    return logcoeff,logexp,outlog.sd_beta[0]

def phi_vs_E(Roi,croi):
    [ksn_st,E_st,rc_st]=stdy.stim_range_dim(np.mean(Roi),np.mean(croi),space_type='log',min_ksn=50,max_ksn=40000,num_points=100)
    logK=np.log10(ksn_st)
    logE=np.log10(E_st.ravel())
    phi=np.diff(logK)/np.diff(logE)
    return phi,E_st[0:-1]
    
    

# Define colors
gc_col='black'
alps_col='royalblue'
bc_col='orange'

# Load final ts outupts
df=pd.read_csv('model_final_stats.csv')
dft=pd.read_csv('model_final_theta.csv')
dfr=pd.read_csv('ksn_e_fit.csv')

repo_location='/Users/aforte/Documents/GitHub/snowmelt-tectonics'

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

stdy=st.StimSteady()
psi=stdy.Psi_c*(1e6*365.25*24*60*60)

## Regime Plot
f1=plt.figure(1,figsize=(8,6),dpi=300)
ax1=plt.subplot(2,1,1)
ax1.set_xscale('log')
ax1.set_xlim((0.001,100))
ax1.set_ylim((0,1))
ax1.axvline(0.1,c='gray')
ax1.axvline(10,c='gray')
ax1.set_xlabel(r'E/$\Psi$')
ax1.set_ylabel(r'$\phi$')
ax1.set_yticks([0,0.2,0.4,0.6,0.8,1.0])
ax1t=ax1.twinx()
ax1t.set_yticks([0,0.2,0.4,0.6,0.8,1.0])
ax1t.set_yticklabels([r'$\infty$','5','2.5','1.7','1.25','1'])
ax1t.get_yticklabels()[0].set_fontsize(15)
ax1t.set_ylabel('n')

ax1a=plt.subplot(2,1,2)
ax1a.set_xscale('log')
ax1a.set_xlim((0.001,0.1))
ax1a.set_ylim((0,0.4))
ax1a.set_xlabel(r'E/$\Psi$')
ax1a.set_ylabel(r'$\phi$')
ax1a.set_yticks([0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4])
ax1at=ax1a.twinx()
ax1at.set_ylim((0,0.4))
ax1at.set_yticks([0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4])
ax1at.set_yticklabels([r'$\infty$','20','10','6.7','5','4','3.3','2.9','2.5'])
ax1at.get_yticklabels()[0].set_fontsize(15)
ax1at.set_ylabel('n')


f2=plt.figure(2,figsize=(8,8),dpi=300,layout='constrained')
ax2=plt.subplot(2,2,1)
ax2.set_xlabel('spatial-STIM E [m/Myr]')
ax2.set_ylabel('STIM E [m/Myr]')
ax2.plot(np.linspace(200,10000),np.linspace(200,10000),c='r',linestyle=':',label='1:1 Line',zorder=0)
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_ylim((10,30000))
ax2.set_xlim((200,10000))
ax2.text(0.01, 0.99, 'A',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax2.transAxes,
        fontsize=12,fontweight='extra bold')

ax3=plt.subplot(2,2,2)
ax3.set_xlabel('spatial-STIM E [m/Myr]')
ax3.set_ylabel('STIM E / spatial-STIM E')
ax3.axhline(1,c='r',linestyle=':',label='1:1 Line',zorder=0)
ax3.set_xscale('log')
ax3.set_xlim((200,10000))
ax3.set_ylim((0.4,1.6))
ax3.text(0.01, 0.99, 'B',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax3.transAxes,
        fontsize=12,fontweight='extra bold')

ax4=plt.subplot(2,2,3)
ax4.set_xlabel('Erosion Rate [m/Myr]')
ax4.set_ylabel(r'$k_{sn}$')
ax4.set_yscale('log')
ax4.set_xscale('log')
ax4.set_xlim((200,10000))
ax4.set_ylim((100,600))
ax4.set_yticks([100,200,300,400,500,600])
ax4.set_yticklabels([100,200,300,400,500,600])
ax4.text(0.01, 0.99, 'C',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax4.transAxes,
        fontsize=12,fontweight='extra bold')

ax5=plt.subplot(2,2,4)
ax5.set_xlabel('spatial-STIM n')
ax5.set_ylabel('STIM n')
ax5.set_xlim((2,20))
ax5.set_ylim((2,20))
ax5.plot(np.linspace(2,20,20),np.linspace(2,20,20),c='r',linestyle=':')
ax5.text(0.01, 0.99, 'D',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax5.transAxes,
        fontsize=12,fontweight='extra bold')

f3=plt.figure(3,figsize=(8,4),dpi=300,layout='constrained')
ax6=plt.subplot(1,2,1)
ax6.set_xlabel('spatial-STIM E [m/Myr]')
ax6.set_ylabel('STIM E / spatial-STIM E')
ax6.axhline(1,c='r',linestyle=':',label='1:1 Line',zorder=0)
ax6.set_xscale('log')
ax6.set_xlim((200,10000))
ax6.set_ylim((0.4,1.6))
ax6.text(0.01, 0.99, 'A',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax6.transAxes,
        fontsize=12,fontweight='extra bold')

ax7=plt.subplot(1,2,2)
ax7.set_xlabel('spatial-STIM n')
ax7.set_ylabel('STIM n')
ax7.set_xlim((2,20))
ax7.set_ylim((2,20))
ax7.plot(np.linspace(2,20,20),np.linspace(2,20,20),c='r',linestyle=':')
ax7.text(0.01, 0.99, 'B',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax7.transAxes,
        fontsize=12,fontweight='extra bold')

for i in range(len(group_list)):
    idx=df['Group']==group_list[i]
    
    ksn=df.loc[idx,'mn_ksn'].to_numpy()
    ksnu=df.loc[idx,'std_ksn'].to_numpy()
    E=df.loc[idx,'mn_E'].to_numpy()
    Eu=df.loc[idx,'std_E'].to_numpy()
    
    cr=df.loc[idx,'cr'].to_numpy()
    mn_R=df.loc[idx,'mn_R'].to_numpy()
    
    stimE=np.zeros(E.shape)
    stimRc=np.zeros(E.shape)
    stimER=np.zeros((2,len(E)))
    for j in range(len(E)):
        [stimE[j],stimRc[j]]=stdy.stim_integrate_one_dim(ksn[j], mn_R[j], cr[j])
        [stimER[0,j],_]=stdy.stim_integrate_one_dim(ksn[j]-ksnu[j],mn_R[j],cr[j])
        [stimER[1,j],_]=stdy.stim_integrate_one_dim(ksn[j]+ksnu[j],mn_R[j],cr[j])
    stimEu=np.copy(stimER)
    stimEu[0,:]=stimE-stimER[0,:]
    stimEu[1,:]=stimER[1,:]-stimE

    
    Ki=np.zeros(E.shape)
    ni=np.zeros(E.shape)
    nistd=np.zeros(E.shape)
    for j in range(len(E)):
        [stimKSrng,stimErng,_]=stdy.stim_range_dim(mn_R[j],cr[j],min_ksn=100,max_ksn=600,space_type='lin')
        
        [Ki[j],ni[j],nistd[j]]=odr_fit(stimKSrng,stimErng.ravel())
        if j==0:
            ax4.plot(stimErng,stimKSrng,c=col_list[i],linestyle=sty_list[i],alpha=0.5,zorder=0,label=label[i])            
        else:           
            ax4.plot(stimErng,stimKSrng,c=col_list[i],linestyle=sty_list[i],alpha=0.5,zorder=0)
    
    [C,phi,phistd]=odr_fit(E,ksn)
    [K,n,nstd]=odr_fit(ksn,E)
    

    if len_list[i]==50:
        ax1.scatter(E/psi,np.ones(E.shape)*phi,c=col_list[i],marker=mar_list[i],zorder=2,label=label[i])
        ax1a.scatter(E/psi,np.ones(E.shape)*phi,c=col_list[i],marker=mar_list[i],zorder=2)
        
        ax2.scatter(E,stimE,c=col_list[i],marker=mar_list[i],label=label[i],zorder=2)
        ax2.errorbar(E,stimE,yerr=stimEu,xerr=Eu,ecolor=col_list[i],elinewidth=0.5,linestyle='',zorder=0)
        ax3.scatter(E,stimE/E,c=col_list[i],marker=mar_list[i],label=label[i],zorder=2)
        
        ax6.scatter(E,stimE/E,c=col_list[i],marker=mar_list[i],label=label[i],zorder=2)

        ax4.scatter(E,ksn,c=col_list[i],marker=mar_list[i],zorder=2)
        ax4.errorbar(E,ksn,yerr=ksnu,xerr=Eu,ecolor=col_list[i],elinewidth=0.5,linestyle='',zorder=1)
        
        ax5.scatter(np.ones(E.shape)*(n),ni,c=col_list[i],marker=mar_list[i],zorder=2)
        ax5.errorbar(np.ones(E.shape)*(n),ni,yerr=np.ones(E.shape)*(nistd),xerr=nistd,ecolor=col_list[i],elinewidth=0.5,linestyle='',zorder=1)
 
        ax7.scatter(np.ones(E.shape)*(n),ni,c=col_list[i],marker=mar_list[i],zorder=2,label=label[i])
        ax7.errorbar(np.ones(E.shape)*(n),ni,yerr=np.ones(E.shape)*(nistd),xerr=nistd,ecolor=col_list[i],elinewidth=0.5,linestyle='',zorder=1)

    else:
        ax1.scatter(E/psi,np.ones(E.shape)*phi,c='w',edgecolor=col_list[i],marker=mar_list[i],zorder=2,label=label[i])
        ax1a.scatter(E/psi,np.ones(E.shape)*phi,c='w',edgecolor=col_list[i],marker=mar_list[i],zorder=2)
        
        ax2.scatter(E,stimE,c='w',edgecolor=col_list[i],marker=mar_list[i],label=label[i],zorder=2)
        ax2.errorbar(E,stimE,yerr=stimEu,xerr=Eu,ecolor=col_list[i],elinewidth=0.5,linestyle='',zorder=0)
        ax3.scatter(E,stimE/E,c='w',edgecolor=col_list[i],marker=mar_list[i],label=label[i],zorder=2)
        ax6.scatter(E,stimE/E,c='w',edgecolor=col_list[i],marker=mar_list[i],label=label[i],zorder=2)
        
        ax4.scatter(E,ksn,c='w',edgecolor=col_list[i],marker=mar_list[i],zorder=2)
        ax4.errorbar(E,ksn,yerr=ksnu,xerr=Eu,ecolor=col_list[i],elinewidth=0.5,linestyle='',zorder=1)
        
        ax5.scatter(np.ones(E.shape)*(n),ni,c='w',edgecolor=col_list[i],marker=mar_list[i],zorder=2)
        ax5.errorbar(np.ones(E.shape)*(n),ni,yerr=np.ones(E.shape)*(nistd),xerr=nistd,ecolor=col_list[i],elinewidth=0.5,linestyle='',zorder=1)

        ax7.scatter(np.ones(E.shape)*(n),ni,c='w',edgecolor=col_list[i],marker=mar_list[i],zorder=2,label=label[i])
        ax7.errorbar(np.ones(E.shape)*(n),ni,yerr=np.ones(E.shape)*(nistd),xerr=nistd,ecolor=col_list[i],elinewidth=0.5,linestyle='',zorder=1)
    
    ax2.plot(E,stimE,c=col_list[i],linestyle=sty_list[i],zorder=0)
    ax3.plot(E,stimE/E,c=col_list[i],linestyle=sty_list[i],zorder=0)
    ax6.plot(E,stimE/E,c=col_list[i],linestyle=sty_list[i],zorder=0)
        
ax2.legend(loc='lower right')  
ax1.legend(loc='upper right',bbox_to_anchor=(1.6,1.1))  
ax4.legend(loc='lower right',bbox_to_anchor=(1.2,0.025))
ax7.legend(loc='lower right')

dfR=pd.read_csv(repo_location+'/stimpy/snowmelt_meanR_cR_relationships.csv')
mnR_vec=np.array([0.25,0.5,1,2.5,5,7.5,10])
lw=np.array([0.3,0.4,0.6,0.8,1,1.5,2])

for i in range(len(mnR_vec)):
    a=dfR.loc[0,'param1']
    b=dfR.loc[0,'param2']
    crR=a*mnR_vec[i]**b
    [phiR,ER]=phi_vs_E(mnR_vec[i],crR)
    ax1.plot(ER/psi,phiR,c='k',linewidth=lw[i],zorder=0)
    ax1a.plot(ER/psi,phiR,c='k',linewidth=lw[i],zorder=0,label=r'$\bar{R}$='+str(mnR_vec[i])+'; $c_{R}$='+str(np.round(crR,2)))
    
for i in range(len(mnR_vec)):
    c=dfR.loc[15,'param1']
    d=dfR.loc[15,'param2']
    crR=a*mnR_vec[i]**b
    crS=c*mnR_vec[i]+d
    [phiS,ES]=phi_vs_E(mnR_vec[i],crS)
    ax1.plot(ES/psi,phiS,c='gray',linewidth=lw[i],linestyle='--',zorder=0)
    ax1a.plot(ES/psi,phiS,c='gray',linewidth=lw[i],linestyle='--',zorder=0,label=r'$\bar{R}$='+str(mnR_vec[i])+'; $c_{R}$='+str(np.round(crS,2)))

ax1a.legend(loc='upper right',bbox_to_anchor=(1.55,1.1))
    
f1.tight_layout()

f2.savefig('P2_SUP_figureS6.pdf',dpi='figure')
f3.savefig('P2_figure11.pdf',dpi='figure')

plt.rcdefaults()