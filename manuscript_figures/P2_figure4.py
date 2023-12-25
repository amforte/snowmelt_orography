#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 15 08:23:29 2023

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

# Greater Caucasus
gc_col='black'

# Alps
alps_col='royalblue'

# British Columbia
bc_col='orange'

# master_location='/Volumes/Choruh/Data/snowmelt_project/'
master_data_location='/Users/aforte/Documents/Python/snowmelt/'
master_location=master_data_location+'model_outputs_v2/'

# Set uplift rates
u_vec=np.array([0.25,0.5,1,2,4,8])

# Set number of timesteps to average across
num_avg=40

# Generate a stream object for binning ksn
sObj=st.Stream(50000,25,dx=100,bin_size=2000)

#Greater Caucasus
gc_mlu=['gc025u','gc05u','gc1u','gc2u','gc4u','gc8u']
gc_plu=['ts_','ts_','ts_','ts_','ts_','ts_']
gc_dlu=['gcu','gcu','gcu','gcu','gcu','gcu']

gc_mll=['gc025l','gc05l','gc1l','gc2l','gc4l','gc8l']
gc_pll=['ts_','ts_','ts_','ts_','ts_','ts_']
gc_dll=['gcl','gcl','gcl','gcl','gcl','gcl']

gc_mcu=st.ModelComparison(master_location,gc_mlu,gc_plu,gc_dlu)
gc_du=gc_mcu.output_result_dict(np.tile([np.inf],6),last_n_ts=num_avg)

gc_mcl=st.ModelComparison(master_location,gc_mll,gc_pll,gc_dll)
gc_dl=gc_mcl.output_result_dict(np.tile([np.inf],6),last_n_ts=num_avg)

# Alps
a_mlu=['a025u','a05u','a1u','a2u','a4u','a8u']
a_plu=['ts_','ts_','ts_','ts_','ts_','ts_']
a_dlu=['au','au','au','au','au','au']

a_mll=['a025l','a05l','a1l','a2l','a4l','a8l']
a_pll=['ts_','ts_','ts_','ts_','ts_','ts_']
a_dll=['al','al','al','al','al','al']

a_mcu=st.ModelComparison(master_location,a_mlu,a_plu,a_dlu)
a_du=a_mcu.output_result_dict(np.tile([np.inf],6),last_n_ts=num_avg)

a_mcl=st.ModelComparison(master_location,a_mll,a_pll,a_dll)
a_dl=a_mcl.output_result_dict(np.tile([np.inf],6),last_n_ts=num_avg)

# British Columbia
bc_mlu=['bc025u','bc05u','bc1u','bc2u','bc4u','bc8u']
bc_plu=['ts_','ts_','ts_','ts_','ts_','ts_']
bc_dlu=['bcu','bcu','bcu','bcu','bcu','bcu']

bc_mll=['bc025l','bc05l','bc1l','bc2l','bc4l','bc8l']
bc_pll=['ts_','ts_','ts_','ts_','ts_','ts_']
bc_dll=['bcl','bcl','bcl','bcl','bcl','bcl']
    
bc_mcu=st.ModelComparison(master_location,bc_mlu,bc_plu,bc_dlu)
bc_du=bc_mcu.output_result_dict(np.tile([np.inf],6),last_n_ts=num_avg)

bc_mcl=st.ModelComparison(master_location,bc_mll,bc_pll,bc_dll)
bc_dl=bc_mcl.output_result_dict(np.tile([np.inf],6),last_n_ts=num_avg)

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

f2=plt.figure(2,figsize=(8,4),layout='constrained')
ax16=f2.add_subplot(1,3,1)
ax16.set_xlabel(r'$\chi$')
ax16.set_ylabel(r'$\Delta$ Z (spatial-STIM - point-STIM) [km]')
ax16.set_ylim((-0.7,0.9))
ax16.set_title('Greater Caucasus')

ax17=f2.add_subplot(1,3,2)
ax17.set_xlabel(r'$\chi$')
# ax17.set_ylabel(r'$\Delta$ Elevation [km]')
ax17.set_ylim((-0.7,0.9))
ax17.set_title('Alps')

ax18=f2.add_subplot(1,3,3)
ax18.set_xlabel(r'$\chi$')
# ax18.set_ylabel(r'$\Delta$ Elevation [km]')
ax18.set_ylim((-0.7,0.9))
ax18.set_title('British Columbia')



f1=plt.figure(1,figsize=(8,10),dpi=300,layout='constrained')
gs=gridspec.GridSpec(5,7,figure=f1)
# Elevation
ax1=f1.add_subplot(gs[0,0:2])
ax1.set_ylabel('Elevation [km]')
ax1.set_ylim((-0.1,3))
ax1.set_xlim((-1,51))
ax1.set_title('Greater Caucasus')

ax2=f1.add_subplot(gs[0,2:4])
ax2.set_title('Alps')
ax2.set_ylim((-0.1,3))
ax2.set_xlim((-1,51))

ax3=f1.add_subplot(gs[0,4:6])
ax3.set_title('British Columbia')
ax3.set_ylim((-0.1,3))
ax3.set_xlim((-1,51))

ax3a=f1.add_subplot(gs[0,6])
ax3a.set_ylim((-0.1,3))
ax3a.set_xlim((0.7,3.3))
ax3a.set_xticks([1,2,3])
ax3a.set_xticklabels(['GC','Alps','BC'])
ax3a.set_title('Means')

# ksn
ax4=f1.add_subplot(gs[1,0:2])
ax4.set_ylabel(r'$k_{sn}$ [m]')
ax4.set_xlim((-1,51))
ax4.set_ylim((100,550))

ax5=f1.add_subplot(gs[1,2:4])
ax5.set_xlim((-1,51))
ax5.set_ylim((100,550))

ax6=f1.add_subplot(gs[1,4:6])
ax6.set_xlim((-1,51))
ax6.set_ylim((100,550))

ax6a=f1.add_subplot(gs[1,6])
ax6a.set_ylim((100,550))
ax6a.set_xlim((0.7,3.3))
ax6a.set_xticks([1,2,3])
ax6a.set_xticklabels(['GC','Alps','BC'])


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
ax9.set_xlim((-1,51))
ax9.axhline(0.35,c='gray',linestyle=':',zorder=0)

ax9a=f1.add_subplot(gs[2,6])
ax9a.set_ylim((0,1))
ax9a.set_xlim((0.7,3.3))
ax9a.set_xticks([1,2,3])
ax9a.set_xticklabels(['GC','Alps','BC'])
ax9a.axhline(0.35,c='gray',linestyle=':',zorder=0)

# Mean Runoff
ax10=f1.add_subplot(gs[3,0:2])
ax10.set_ylabel(r'$\bar{R}$ [mm/day]')
ax10.set_xlim((-1,51))
ax10.set_ylim((0,10))

ax11=f1.add_subplot(gs[3,2:4])
ax11.set_xlim((-1,51))
ax11.set_ylim((0,10))

ax12=f1.add_subplot(gs[3,4:6])
ax12.set_xlim((-1,51))
ax12.set_ylim((0,10))

ax12a=f1.add_subplot(gs[3,6])
ax12a.set_ylim((0,10))
ax12a.set_xlim((0.7,3.3))
ax12a.set_xticks([1,2,3])
ax12a.set_xticklabels(['GC','Alps','BC'])

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
ax15.set_xlim((-1,51))
ax15.set_ylim((0.45,1.4))

ax15a=f1.add_subplot(gs[4,6])
ax15a.set_ylim((0.45,1.4))
ax15a.set_xlim((0.7,3.3))
ax15a.set_xticks([1,2,3])
ax15a.set_xticklabels(['GC','Alps','BC'])





col_vec=colors.Normalize(vmin=0.25,vmax=8)

for i in range(len(u_vec)):
    gc_duI=gc_du[i]
    gc_dlI=gc_dl[i]
    a_duI=a_du[i]
    a_dlI=a_dl[i]
    bc_duI=bc_du[i]
    bc_dlI=bc_dl[i]
    
    # Extract shared x
    x=gc_duI['x']/1000
    xc=gc_duI['x_center']/1000
    chi=gc_duI['chi']
    
    # GC
    gc_duIz=np.mean(gc_duI['zout']/1000,axis=0)
    gc_dlIz=np.mean(gc_dlI['zout']/1000,axis=0)
    gc_duIsn=np.mean(gc_duI['snp'],axis=0)
    gc_dlIsn=np.mean(gc_dlI['snp'],axis=0)  
    gc_duImr=np.mean(gc_duI['mrout'],axis=0)
    gc_dlImr=np.mean(gc_dlI['mrout'],axis=0)    
    gc_duIcr=np.mean(gc_duI['crout'],axis=0)
    gc_dlIcr=np.mean(gc_dlI['crout'],axis=0)

    gc_duIksn=np.mean(np.diff(gc_duI['zout'],axis=1)/np.diff(gc_duI['chi']),axis=0)
    gc_duIksn=np.concatenate(([gc_duIksn[0]],gc_duIksn),axis=0)
    gc_duIksn=np.bincount(sObj.ix,gc_duIksn,sObj.num_bins)[1:sObj.num_bins+1]/np.bincount(sObj.ix,None,sObj.num_bins)[1:sObj.num_bins+1]
    
    gc_dlIksn=np.mean(np.diff(gc_dlI['zout'],axis=1)/np.diff(gc_dlI['chi']),axis=0)
    gc_dlIksn=np.concatenate(([gc_dlIksn[0]],gc_dlIksn),axis=0)
    gc_dlIksn=np.bincount(sObj.ix,gc_dlIksn,sObj.num_bins)[1:sObj.num_bins+1]/np.bincount(sObj.ix,None,sObj.num_bins)[1:sObj.num_bins+1]
    
    
    # Compare predictions of point STIM for profile shape
    stdy=st.StimSteady()
    [ksl,el,_]=stdy.stim_range_dim(np.mean(gc_dlImr),np.mean(gc_dlIcr),min_ksn=25,space_type='lin')
    ixl=np.argmin(np.abs(el-u_vec[i]*1000))
    [ksu,eu,_]=stdy.stim_range_dim(np.mean(gc_duImr),np.mean(gc_duIcr),min_ksn=25,space_type='lin')
    ixu=np.argmin(np.abs(eu-u_vec[i]*1000)) 
    
    ax16.plot(chi,gc_dlIz-chi*ksl[ixl]/1000,linestyle='--',c=cm.acton_r(col_vec(u_vec[i])))
    ax16.plot(chi,gc_duIz-chi*ksu[ixl]/1000,linestyle='-',c=cm.acton_r(col_vec(u_vec[i])))    
    
    

    if i==5:
        ax1.plot(x,gc_duIz,linestyle='-',c=cm.acton_r(col_vec(u_vec[i])),label='Unlinked')
        ax1.plot(x,gc_dlIz,linestyle='--',c=cm.acton_r(col_vec(u_vec[i])),label='Linked')
        ax1.legend(loc='upper left')
        ax1.text(-0.22, 0.99, 'A',
                horizontalalignment='left',
                verticalalignment='top',
                transform=ax1.transAxes,
                fontsize=12,fontweight='extra bold')
    else:
        ax1.plot(x,gc_duIz,linestyle='-',c=cm.acton_r(col_vec(u_vec[i])))
        ax1.plot(x,gc_dlIz,linestyle='--',c=cm.acton_r(col_vec(u_vec[i])))
        
            
    
    ax4.plot(xc,gc_duIksn,linestyle='-',c=cm.acton_r(col_vec(u_vec[i])))
    ax4.plot(xc,gc_dlIksn,linestyle='--',c=cm.acton_r(col_vec(u_vec[i])))
    if i==0:
        ax4.text(-0.22, 0.99, 'B',
                horizontalalignment='left',
                verticalalignment='top',
                transform=ax4.transAxes,
                fontsize=12,fontweight='extra bold')
    

    ax7.plot(xc,gc_duIsn,linestyle='-',c=cm.acton_r(col_vec(u_vec[i])))
    ax7.plot(xc,gc_dlIsn,linestyle='--',c=cm.acton_r(col_vec(u_vec[i])))
    if i==0:
        ax7.text(-0.22, 0.99, 'C',
                horizontalalignment='left',
                verticalalignment='top',
                transform=ax7.transAxes,
                fontsize=12,fontweight='extra bold')
    
    ax10.plot(xc,gc_duImr,linestyle='-',c=cm.acton_r(col_vec(u_vec[i])))
    ax10.plot(xc,gc_dlImr,linestyle='--',c=cm.acton_r(col_vec(u_vec[i]))) 
    if i==0:
        ax10.text(-0.22, 0.99, 'D',
                horizontalalignment='left',
                verticalalignment='top',
                transform=ax10.transAxes,
                fontsize=12,fontweight='extra bold')    
    
    ax13.plot(xc,gc_duIcr,linestyle='-',c=cm.acton_r(col_vec(u_vec[i])))
    ax13.plot(xc,gc_dlIcr,linestyle='--',c=cm.acton_r(col_vec(u_vec[i])))
    if i==0:
        ax13.text(-0.22, 0.99, 'E',
                horizontalalignment='left',
                verticalalignment='top',
                transform=ax13.transAxes,
                fontsize=12,fontweight='extra bold')
    
     # Alps
    a_duIz=np.mean(a_duI['zout']/1000,axis=0)
    a_dlIz=np.mean(a_dlI['zout']/1000,axis=0)
    a_duIsn=np.mean(a_duI['snp'],axis=0)
    a_dlIsn=np.mean(a_dlI['snp'],axis=0)  
    a_duImr=np.mean(a_duI['mrout'],axis=0)
    a_dlImr=np.mean(a_dlI['mrout'],axis=0)    
    a_duIcr=np.mean(a_duI['crout'],axis=0)
    a_dlIcr=np.mean(a_dlI['crout'],axis=0)

    a_duIksn=np.mean(np.diff(a_duI['zout'],axis=1)/np.diff(a_duI['chi']),axis=0)
    a_duIksn=np.concatenate(([a_duIksn[0]],a_duIksn),axis=0)
    a_duIksn=np.bincount(sObj.ix,a_duIksn,sObj.num_bins)[1:sObj.num_bins+1]/np.bincount(sObj.ix,None,sObj.num_bins)[1:sObj.num_bins+1]
    
    a_dlIksn=np.mean(np.diff(a_dlI['zout'],axis=1)/np.diff(a_dlI['chi']),axis=0)
    a_dlIksn=np.concatenate(([a_dlIksn[0]],a_dlIksn),axis=0)
    a_dlIksn=np.bincount(sObj.ix,a_dlIksn,sObj.num_bins)[1:sObj.num_bins+1]/np.bincount(sObj.ix,None,sObj.num_bins)[1:sObj.num_bins+1]
    
    ax2.plot(x,a_duIz,linestyle='-',c=cm.acton_r(col_vec(u_vec[i])))
    ax2.plot(x,a_dlIz,linestyle='--',c=cm.acton_r(col_vec(u_vec[i])))  
    
    ax5.plot(xc,a_duIksn,linestyle='-',c=cm.acton_r(col_vec(u_vec[i])))
    ax5.plot(xc,a_dlIksn,linestyle='--',c=cm.acton_r(col_vec(u_vec[i])))

    ax8.plot(xc,a_duIsn,linestyle='-',c=cm.acton_r(col_vec(u_vec[i])))
    ax8.plot(xc,a_dlIsn,linestyle='--',c=cm.acton_r(col_vec(u_vec[i])))
    
    ax11.plot(xc,a_duImr,linestyle='-',c=cm.acton_r(col_vec(u_vec[i])))
    ax11.plot(xc,a_dlImr,linestyle='--',c=cm.acton_r(col_vec(u_vec[i])))     
    
    ax14.plot(xc,a_duIcr,linestyle='-',c=cm.acton_r(col_vec(u_vec[i])))
    ax14.plot(xc,a_dlIcr,linestyle='--',c=cm.acton_r(col_vec(u_vec[i])))  
    
    
    # Compare predictions of point STIM for profile shape
    stdy=st.StimSteady()
    [ksl,el,_]=stdy.stim_range_dim(np.mean(a_dlImr),np.mean(a_dlIcr),min_ksn=25,space_type='lin')
    ixl=np.argmin(np.abs(el-u_vec[i]*1000))
    [ksu,eu,_]=stdy.stim_range_dim(np.mean(a_duImr),np.mean(a_duIcr),min_ksn=25,space_type='lin')
    ixu=np.argmin(np.abs(eu-u_vec[i]*1000)) 
    
    ax17.plot(chi,a_dlIz-chi*ksl[ixl]/1000,linestyle='--',c=cm.acton_r(col_vec(u_vec[i])))
    ax17.plot(chi,a_duIz-chi*ksu[ixl]/1000,linestyle='-',c=cm.acton_r(col_vec(u_vec[i])))

    # BC
    bc_duIz=np.mean(bc_duI['zout']/1000,axis=0)
    bc_dlIz=np.mean(bc_dlI['zout']/1000,axis=0)
    bc_duIsn=np.mean(bc_duI['snp'],axis=0)
    bc_dlIsn=np.mean(bc_dlI['snp'],axis=0)  
    bc_duImr=np.mean(bc_duI['mrout'],axis=0)
    bc_dlImr=np.mean(bc_dlI['mrout'],axis=0)    
    bc_duIcr=np.mean(bc_duI['crout'],axis=0)
    bc_dlIcr=np.mean(bc_dlI['crout'],axis=0)

    bc_duIksn=np.mean(np.diff(bc_duI['zout'],axis=1)/np.diff(bc_duI['chi']),axis=0)
    bc_duIksn=np.concatenate(([bc_duIksn[0]],bc_duIksn),axis=0)
    bc_duIksn=np.bincount(sObj.ix,bc_duIksn,sObj.num_bins)[1:sObj.num_bins+1]/np.bincount(sObj.ix,None,sObj.num_bins)[1:sObj.num_bins+1]
    
    bc_dlIksn=np.mean(np.diff(bc_dlI['zout'],axis=1)/np.diff(bc_dlI['chi']),axis=0)
    bc_dlIksn=np.concatenate(([bc_dlIksn[0]],bc_dlIksn),axis=0)
    bc_dlIksn=np.bincount(sObj.ix,bc_dlIksn,sObj.num_bins)[1:sObj.num_bins+1]/np.bincount(sObj.ix,None,sObj.num_bins)[1:sObj.num_bins+1]
    
    ax3.plot(x,bc_duIz,linestyle='-',c=cm.acton_r(col_vec(u_vec[i])))
    ax3.plot(x,bc_dlIz,linestyle='--',c=cm.acton_r(col_vec(u_vec[i])))  
    
    ax6.plot(xc,bc_duIksn,linestyle='-',c=cm.acton_r(col_vec(u_vec[i])))
    ax6.plot(xc,bc_dlIksn,linestyle='--',c=cm.acton_r(col_vec(u_vec[i])))

    ax9.plot(xc,bc_duIsn,linestyle='-',c=cm.acton_r(col_vec(u_vec[i])))
    ax9.plot(xc,bc_dlIsn,linestyle='--',c=cm.acton_r(col_vec(u_vec[i])))
    
    ax12.plot(xc,bc_duImr,linestyle='-',c=cm.acton_r(col_vec(u_vec[i])))
    ax12.plot(xc,bc_dlImr,linestyle='--',c=cm.acton_r(col_vec(u_vec[i])))     
    
    ax15.plot(xc,bc_duIcr,linestyle='-',c=cm.acton_r(col_vec(u_vec[i])))
    ax15.plot(xc,bc_dlIcr,linestyle='--',c=cm.acton_r(col_vec(u_vec[i])))  
    
    # Compare predictions of point STIM for profile shape
    stdy=st.StimSteady()
    [ksl,el,_]=stdy.stim_range_dim(np.mean(bc_dlImr),np.mean(bc_dlIcr),min_ksn=25,space_type='lin')
    ixl=np.argmin(np.abs(el-u_vec[i]*1000))
    [ksu,eu,_]=stdy.stim_range_dim(np.mean(bc_duImr),np.mean(bc_duIcr),min_ksn=25,space_type='lin')
    ixu=np.argmin(np.abs(eu-u_vec[i]*1000)) 
    
    ax18.plot(chi,bc_dlIz-chi*ksl[ixl]/1000,linestyle='--',c=cm.acton_r(col_vec(u_vec[i])))
    ax18.plot(chi,bc_duIz-chi*ksu[ixl]/1000,linestyle='-',c=cm.acton_r(col_vec(u_vec[i])))
    
    if i==0:
        norm=colors.Normalize(vmin=0.25,vmax=8)
        cbar1=plt.colorbar(cmm.ScalarMappable(norm=norm,cmap=cm.acton_r),ax=ax13,orientation='horizontal',shrink=0.8)
        cbar1.ax.set_xlabel('Uplift Rate [m/Myr]')
        
        cbar2=plt.colorbar(cmm.ScalarMappable(norm=norm,cmap=cm.acton_r),ax=ax14,orientation='horizontal',shrink=0.8)
        cbar2.ax.set_xlabel('Uplift Rate [m/Myr]')
        
        cbar3=plt.colorbar(cmm.ScalarMappable(norm=norm,cmap=cm.acton_r),ax=ax15,orientation='horizontal',shrink=0.8)
        cbar3.ax.set_xlabel('Uplift Rate [m/Myr]')
        
        cbar4=plt.colorbar(cmm.ScalarMappable(norm=norm,cmap=cm.acton_r),ax=ax16,orientation='horizontal',shrink=0.8)
        cbar4.ax.set_xlabel('Uplift Rate [m/Myr]')
        
        cbar5=plt.colorbar(cmm.ScalarMappable(norm=norm,cmap=cm.acton_r),ax=ax17,orientation='horizontal',shrink=0.8)
        cbar5.ax.set_xlabel('Uplift Rate [m/Myr]')
        
        cbar6=plt.colorbar(cmm.ScalarMappable(norm=norm,cmap=cm.acton_r),ax=ax18,orientation='horizontal',shrink=0.8)
        cbar6.ax.set_xlabel('Uplift Rate [m/Myr]')        
        
        
    if i==5:    
        ax3a.scatter(0.9,np.mean(gc_duIz),color=cm.acton_r(col_vec(u_vec[i])),marker='o',s=20,label='Unlinked')
        ax3a.scatter(1.1,np.mean(gc_dlIz),edgecolor=cm.acton_r(col_vec(u_vec[i])),facecolor='w',marker='s',s=20,label='Linked')
        ax3a.legend(loc='upper center')
    else:
        ax3a.scatter(0.9,np.mean(gc_duIz),color=cm.acton_r(col_vec(u_vec[i])),marker='o',s=20)
        ax3a.scatter(1.1,np.mean(gc_dlIz),edgecolor=cm.acton_r(col_vec(u_vec[i])),facecolor='w',marker='s',s=20)
    ax3a.scatter(1.9,np.mean(a_duIz),color=cm.acton_r(col_vec(u_vec[i])),marker='o',s=20)
    ax3a.scatter(2.1,np.mean(a_dlIz),edgecolor=cm.acton_r(col_vec(u_vec[i])),facecolor='w',marker='s',s=20)
    ax3a.scatter(2.9,np.mean(bc_duIz),color=cm.acton_r(col_vec(u_vec[i])),marker='o',s=20)
    ax3a.scatter(3.1,np.mean(bc_dlIz),edgecolor=cm.acton_r(col_vec(u_vec[i])),facecolor='w',marker='s',s=20)    

    ax6a.scatter(0.9,np.mean(gc_duIksn),color=cm.acton_r(col_vec(u_vec[i])),marker='o',s=20)
    ax6a.scatter(1.1,np.mean(gc_dlIksn),edgecolor=cm.acton_r(col_vec(u_vec[i])),facecolor='w',marker='s',s=20)
    ax6a.scatter(1.9,np.mean(a_duIksn),color=cm.acton_r(col_vec(u_vec[i])),marker='o',s=20)
    ax6a.scatter(2.1,np.mean(a_dlIksn),edgecolor=cm.acton_r(col_vec(u_vec[i])),facecolor='w',marker='s',s=20)
    ax6a.scatter(2.9,np.mean(bc_duIksn),color=cm.acton_r(col_vec(u_vec[i])),marker='o',s=20)
    ax6a.scatter(3.1,np.mean(bc_dlIksn),edgecolor=cm.acton_r(col_vec(u_vec[i])),facecolor='w',marker='s',s=20)  

    ax9a.scatter(0.9,np.mean(gc_duIsn),color=cm.acton_r(col_vec(u_vec[i])),marker='o',s=20)
    ax9a.scatter(1.1,np.mean(gc_dlIsn),edgecolor=cm.acton_r(col_vec(u_vec[i])),facecolor='w',marker='s',s=20)
    ax9a.scatter(1.9,np.mean(a_duIsn),color=cm.acton_r(col_vec(u_vec[i])),marker='o',s=20)
    ax9a.scatter(2.1,np.mean(a_dlIsn),edgecolor=cm.acton_r(col_vec(u_vec[i])),facecolor='w',marker='s',s=20)
    ax9a.scatter(2.9,np.mean(bc_duIsn),color=cm.acton_r(col_vec(u_vec[i])),marker='o',s=20)
    ax9a.scatter(3.1,np.mean(bc_dlIsn),edgecolor=cm.acton_r(col_vec(u_vec[i])),facecolor='w',marker='s',s=20)

    ax12a.scatter(0.9,np.mean(gc_duImr),color=cm.acton_r(col_vec(u_vec[i])),marker='o',s=20)
    ax12a.scatter(1.1,np.mean(gc_dlImr),edgecolor=cm.acton_r(col_vec(u_vec[i])),facecolor='w',marker='s',s=20)
    ax12a.scatter(1.9,np.mean(a_duImr),color=cm.acton_r(col_vec(u_vec[i])),marker='o',s=20)
    ax12a.scatter(2.1,np.mean(a_dlImr),edgecolor=cm.acton_r(col_vec(u_vec[i])),facecolor='w',marker='s',s=20)
    ax12a.scatter(2.9,np.mean(bc_duImr),color=cm.acton_r(col_vec(u_vec[i])),marker='o',s=20)
    ax12a.scatter(3.1,np.mean(bc_dlImr),edgecolor=cm.acton_r(col_vec(u_vec[i])),facecolor='w',marker='s',s=20)

    ax15a.scatter(0.9,np.mean(gc_duIcr),color=cm.acton_r(col_vec(u_vec[i])),marker='o',s=20)
    ax15a.scatter(1.1,np.mean(gc_dlIcr),edgecolor=cm.acton_r(col_vec(u_vec[i])),facecolor='w',marker='s',s=20)
    ax15a.scatter(1.9,np.mean(a_duIcr),color=cm.acton_r(col_vec(u_vec[i])),marker='o',s=20)
    ax15a.scatter(2.1,np.mean(a_dlIcr),edgecolor=cm.acton_r(col_vec(u_vec[i])),facecolor='w',marker='s',s=20)
    ax15a.scatter(2.9,np.mean(bc_duIcr),color=cm.acton_r(col_vec(u_vec[i])),marker='o',s=20)
    ax15a.scatter(3.1,np.mean(bc_dlIcr),edgecolor=cm.acton_r(col_vec(u_vec[i])),facecolor='w',marker='s',s=20)
    
plt.rcdefaults()
f1.savefig('P2_figure4.pdf')


