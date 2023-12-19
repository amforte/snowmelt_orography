#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 14 11:31:40 2023

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

# master_location='/Volumes/Choruh/Data/snowmelt_project/'
master_data_location='/Users/aforte/Documents/Python/snowmelt/'
master_location=master_data_location+'model_outputs_v2/'

# Load final ts outupts
df=pd.read_csv('model_final_stats.csv')
dfs=pd.read_csv('model_final_SF.csv')

group_list=['gcu','gcl','gcu_R','gcl_R','bcu','bcl','bcu_R','bcl_R']
mar_list=['o','s','o','s','o','s','o','s'] 
len_list=[50,50,10,10,50,50,10,10]
col_list=[gc_col,gc_col,gc_col,gc_col,bc_col,bc_col,bc_col,bc_col] #,'gray','gray']
label=['GC Unl.','GC L.','GC Unl. - Rain Only','GC L. - Rain Only',
       'BC Unl.','BC L.','BC Unl. - Rain Only','BC L. - Rain Only'] #,'GC Unlinked Area','GC Linked Area']

sty_list=['-','--',':','-.','-','--',':','-.'] 


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

f1=plt.figure(figsize=(8.5,10))
f1.set_dpi(250)

gs1=gridspec.GridSpec(5,2,figure=f1)

ax1=f1.add_subplot(gs1[0:3,0])
# ax1.set_xlabel('Erosion Rate [m/Myr]')
ax1.set_ylabel(r'$k_{sn}$ [m]')
ax1.set_xlim((-200,10200))
ax1.set_ylim((-50,575))


ax3=f1.add_subplot(gs1[3:4,0])
ax3.set_xlabel('Erosion Rate [m/Myr]')
ax3.set_ylabel('Max Snow Fraction')
ax3.set_ylim((-0.1,1.1))
ax3.axhline(0.35,c='k',linestyle=':')
ax3.set_xlim((-200,10200))


E_vec=np.logspace(-1,4,100)


for i in range(len(group_list)):
    idx=df['Group']==group_list[i]
    
    ksn=df.loc[idx,'mn_ksn'].to_numpy()
    E=df.loc[idx,'mn_E'].to_numpy()/1e6

    ksns=df.loc[idx,'stde_ksn'].to_numpy()
    Es=df.loc[idx,'stde_E'].to_numpy()
    
    # Fit
    [K,n]=odr_fit(ksn,E)
    [C,phi]=odr_fit(E*1e6,ksn)

    # Snowmelt
    idxs=dfs['Group']==group_list[i]
    mx_snp=dfs.loc[idxs,'max_snow_fraction']
    
    if len_list[i]==50:
        ax1.scatter(E*1e6,ksn,c=col_list[i],marker=mar_list[i],zorder=2,label=label[i]+'; n = '+str(np.round(n,1)))
        
        ax3.scatter(E*1e6,mx_snp,c=col_list[i],marker=mar_list[i])
        

    else:
        ax1.scatter(E*1e6,ksn,c='w',edgecolor=col_list[i],marker=mar_list[i],zorder=2,label=label[i]+'; n = '+str(np.round(n,1)))
        
                   
    ax1.plot(E_vec,C*E_vec**phi,c=col_list[i],linestyle=sty_list[i],zorder=1)
    ax1.errorbar(E*1e6,ksn,ksns,Es,ecolor=col_list[i],linestyle='',zorder=0,elinewidth=0.5)
    
    

ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlim((200,10000))
ax1.set_ylim((100,500))
ax1.set_xticks(ticks=np.array([250,500,1000,2000,4000,8000]),labels=['250','500','1000','2000','4000','8000'])
ax1.set_yticks(ticks=np.array([100,200,300,400,500]),labels=['100','200','300','400','500'])

ax3.set_xscale('log')
ax3.set_xlim((200,10000))
ax3.set_xticks(ticks=np.array([250,500,1000,2000,4000,8000]),labels=['250','500','1000','2000','4000','8000'])


ax1.legend(loc='best')


ax1.text(0.01, 0.99, 'A',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax1.transAxes,
        fontsize=12,fontweight='extra bold')

ax3.text(0.01, 0.99, 'B',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax3.transAxes,
        fontsize=12,fontweight='extra bold')

# Elevation
ax2a=f1.add_subplot(gs1[0,1])
ax2a.set_ylabel('Elev. Relative to BL [km]')
ax2a.set_ylim((-0.1,3))
ax2a.set_xlim((-1,51))
ax2a.text(0.01, 0.99, 'C',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax2a.transAxes,
        fontsize=12,fontweight='extra bold')
ax2a.set_title('1 mm/yr Models')

# ksn
ax2b=f1.add_subplot(gs1[1,1])
ax2b.set_ylabel(r'$k_{sn}$ [m]')
ax2b.set_xlim((-1,51))
ax2b.set_ylim((100,550))
ax2b.text(0.01, 0.99, 'D',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax2b.transAxes,
        fontsize=12,fontweight='extra bold')

# Snow Fraction
ax2c=f1.add_subplot(gs1[2,1])
ax2c.set_ylabel('Snow Fraction')
ax2c.set_ylim((0,1.1))
ax2c.set_xlim((-1,51))
ax2c.axhline(0.35,c='gray',linestyle=':',zorder=0)
ax2c.text(0.01, 0.99, 'E',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax2c.transAxes,
        fontsize=12,fontweight='extra bold')

# Mean Runoff
ax2d=f1.add_subplot(gs1[3,1])
ax2d.set_ylabel(r'$\bar{R}$ [mm/day]')
ax2d.set_xlim((-1,51))
ax2d.set_ylim((0,10))
ax2d.text(0.01, 0.99, 'F',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax2d.transAxes,
        fontsize=12,fontweight='extra bold')

# Variability
ax2e=f1.add_subplot(gs1[4,1])
ax2e.set_ylabel(r'Variability $c_{R}$')
ax2e.set_xlabel('Stream Distance [km]')
ax2e.set_xlim((-1,51))
ax2e.set_ylim((0.45,1.6))
ax2e.text(0.01, 0.99, 'G',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax2e.transAxes,
        fontsize=12,fontweight='extra bold')

# Set number of timesteps to average across
num_avg=40

mlu=['gc1u','gc1u_RO','bc1u','bc1u_RO','gc1u_1kmBL','gc1u_2kmBL','bc1u_1kmBL','bc1u_2kmBL']
plu=['ts_','ts_','ts_','ts_','ts_','ts_','ts_','ts_']
dlu=['snw','snw','snw','snw','snw','snw','snw','snw']

clu=[gc_col,gc_col,bc_col,bc_col,'darkgray','darkgray','gold','gold']
lsu=['-',':','-',':','-','--','-','--']
llu=['GC Unlinked','GC Unlinked Rain Only','BC Unlinked','BC Unlinked Rain Only',
     'GC Unlinked 1 km BL','GC Unlinked 2 km BL','BC Unlinked 1 km BL','BC Unlinked 2 km BL']
rol=[True,False,True,False,True,True,True,True]
bll=[0,0,0,0,1000,2000,1000,2000]

mcu=st.ModelComparison(master_location,mlu,plu,dlu)
du=mcu.output_result_dict(np.tile([np.inf],4),last_n_ts=num_avg)

sObj=st.Stream(50000,25,dx=100,bin_size=2000)

for i in range(len(mlu)):
    duI=du[i]

    
    # Extract shared x
    x=duI['x']/1000
    xc=duI['x_center']/1000
    
    duIz=np.mean((duI['zout']-bll[i])/1000,axis=0)
    if rol[i]:
        duIsn=np.mean(duI['snp'],axis=0)
    duImr=np.mean(duI['mrout'],axis=0)   
    duIcr=np.mean(duI['crout'],axis=0)


    duIksn=np.mean(np.diff(duI['zout'],axis=1)/np.diff(duI['chi']),axis=0)
    duIksn=np.concatenate(([duIksn[0]],duIksn),axis=0)
    duIksn=np.bincount(sObj.ix,duIksn,sObj.num_bins)[1:sObj.num_bins+1]/np.bincount(sObj.ix,None,sObj.num_bins)[1:sObj.num_bins+1]
    
    print(llu[i]+' : mean ksn = '+str(np.round(np.mean(duIksn),2)))
    
    ax2a.plot(x,duIz,linestyle=lsu[i],c=clu[i],label=llu[i])    
    ax2b.plot(xc,duIksn,linestyle=lsu[i],c=clu[i],label=llu[i])
    
    if rol[i]:
        ax2c.plot(xc,duIsn,linestyle=lsu[i],c=clu[i],label=llu[i])
    
    ax2d.plot(xc,duImr,linestyle=lsu[i],c=clu[i],label=llu[i])
    ax2e.plot(xc,duIcr,linestyle=lsu[i],c=clu[i],label=llu[i])
    

ax2e.legend(loc='lower left',ncol=2,bbox_to_anchor=[-1.4,0.1])
plt.subplots_adjust(hspace=0.25,wspace=0.25)

# f1.tight_layout()
plt.rcdefaults()

f1.savefig('P2_figure12.pdf',dpi="figure")