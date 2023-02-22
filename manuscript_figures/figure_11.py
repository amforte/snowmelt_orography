#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  8 13:42:31 2023

@author: aforte
"""

import sys
sys.path.insert(0,'/Users/aforte/Documents/GitHub/snowmelt_orography')
import stimpy as st
import os
import matplotlib.ticker
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from cmcrameri import cm
from matplotlib import colors
from scipy.stats import linregress

def response_time(K,n,m,ka,h,xc,L,Ui,fU):
    # Determine beta
    hmn=(h*m)/n 
    if hmn==1:
        bta = (ka**(-m/n)) * np.log(L/xc)
    else:
        term1=ka**(-m/n)
        term2=(1-hmn)**-1
        term3=L**(1-hmn)
        term4=xc**(1-hmn)
        bta = term1 * term2 * (term3 - term4)
    # Calculate response time
    term1=K**(-1/n)
    term2=Ui**(1/n-1)
    term3=(fU**(1/n)-1)
    term4=(fU-1)**-1
    TU = bta * term1 * term2 * term3 * term4
    return TU


def response_time2(K,n,m,ka,h,xc,L,Uf,ksn0):
    Ui=ksn0**n * K
    fU=Uf/Ui
    # Determine beta
    hmn=(h*m)/n 
    if hmn==1:
        bta = (ka**(-m/n)) * np.log(L/xc)
    else:
        term1=ka**(-m/n)
        term2=(1-hmn)**-1
        term3=L**(1-hmn)
        term4=xc**(1-hmn)
        bta = term1 * term2 * (term3 - term4)
    # Calculate response time
    term1=K**(-1/n)
    term2=Ui**(1/n-1)
    term3=(fU**(1/n)-1)
    term4=(fU-1)**-1
    TU = bta * term1 * term2 * term3 * term4
    return TU

master_location='/Volumes/Choruh/Data/snowmelt_project/model_outputs/'

model_list1U=['gc025u','gc05u','gc1u','gc2u','gc4u','gc8u']
model_list1L=['gc025l','gc05l','gc1l','gc2l','gc4l','gc8l']

model_list2U=['a025u','a05u','a1u','a2u','a4u','a8u']
model_list2L=['a025l','a05l','a1l','a2l','a4l','a8l']

model_list3U=['bc025u','bc05u','bc1u','bc2u','bc4u','bc8u']
model_list3L=['bc025l','bc05l','bc1l','bc2l','bc4l','bc8l']


# Set thresholds for determining empirical steady-state
mx_thresh=0.1
mn_thresh=0.005
# Set colors
gc_col='black'
alps_col='royalblue'
bc_col='orange'

# Set contsants for analytical steady-state
# Generate stream  instance
sObj=st.Stream(50000,25,dx=100,bin_size=2000)
x=sObj.length-sObj.x
A=sObj.A - sObj.Ac
[slp,yint,r,p,se]=linregress(np.log10(x[0:-1]),np.log10(A[0:-1]))
h=slp
ka=10**yint
xc=(sObj.Ac/ka)**(1/h)
L=50000 + xc # Accounts for both the length of the channel and the hypothetical hillslope
# Read in fits
fit_df=pd.read_csv('ksn_e_fit.csv')
# Convert C and phi to K and n
K=fit_df['K'].to_numpy()
n=fit_df['n'].to_numpy()
m=n*0.5


ksn0=25

# Greater Caucasus
ts1U=[]
dmxz1U=[]
dmnz1U=[]
mxdz1U=[]
u1U=[]
ss_mxz1U=[]
ss_mnz1U=[]
ss_A1U=[]
idx=fit_df['Group']=='GC Unlinked'
K=fit_df.loc[idx,'K'].to_numpy()
n=fit_df.loc[idx,'n'].to_numpy()
for i in range(len(model_list1U)):
    mObj=st.Stim1D(os.path.join(master_location,model_list1U[i]))
    [ts,u,dmxz,dmnz,mxdz]=mObj.calculate_ss()
    mxix=np.nonzero(dmxz<mx_thresh)[0][:1][0]
    mnix=np.nonzero(dmnz<mn_thresh)[0][:1][0]
    ss_mxz1U.append(ts[mxix])
    ss_mnz1U.append(ts[mnix])
    ts1U.append(ts)
    dmxz1U.append(dmxz)
    dmnz1U.append(dmnz)
    mxdz1U.append(mxdz)
    u1U.append(u)
    TU=response_time2(K,n,n/2,ka,h,xc,L,u,ksn0)/1e6
    ss_A1U.append(TU)
    
ts1L=[]
dmxz1L=[]
dmnz1L=[]
mxdz1L=[]
u1L=[]
ss_mxz1L=[]
ss_mnz1L=[]
ss_A1L=[]
idx=fit_df['Group']=='GC Linked'
K=fit_df.loc[idx,'K'].to_numpy()
n=fit_df.loc[idx,'n'].to_numpy()
for i in range(len(model_list1L)):
    mObj=st.Stim1D(os.path.join(master_location,model_list1L[i]))
    [ts,u,dmxz,dmnz,mxdz]=mObj.calculate_ss()
    mxix=np.nonzero(dmxz<mx_thresh)[0][:1][0]
    mnix=np.nonzero(dmnz<mn_thresh)[0][:1][0]
    ss_mxz1L.append(ts[mxix])
    ss_mnz1L.append(ts[mnix])
    ts1L.append(ts)
    dmxz1L.append(dmxz)
    dmnz1L.append(dmnz)
    mxdz1L.append(mxdz)
    u1L.append(u)
    TU=response_time2(K,n,n/2,ka,h,xc,L,u,ksn0)/1e6
    ss_A1L.append(TU)
    
# Alps
ts2U=[]
dmxz2U=[]
dmnz2U=[]
mxdz2U=[]
u2U=[]
ss_mxz2U=[]
ss_mnz2U=[]
ss_A2U=[]
idx=fit_df['Group']=='Alps Unlinked'
K=fit_df.loc[idx,'K'].to_numpy()
n=fit_df.loc[idx,'n'].to_numpy()
for i in range(len(model_list2U)):
    mObj=st.Stim1D(os.path.join(master_location,model_list2U[i]))
    [ts,u,dmxz,dmnz,mxdz]=mObj.calculate_ss()
    mxix=np.nonzero(dmxz<mx_thresh)[0][:1][0]
    mnix=np.nonzero(dmnz<mn_thresh)[0][:1][0]
    ss_mxz2U.append(ts[mxix])
    ss_mnz2U.append(ts[mnix])
    ts2U.append(ts)
    dmxz2U.append(dmxz)
    dmnz2U.append(dmnz)
    mxdz2U.append(mxdz)
    u2U.append(u)
    TU=response_time2(K,n,n/2,ka,h,xc,L,u,ksn0)/1e6
    ss_A2U.append(TU)
    
ts2L=[]
dmxz2L=[]
dmnz2L=[]
mxdz2L=[]
u2L=[]
ss_mxz2L=[]
ss_mnz2L=[]
ss_A2L=[]
idx=fit_df['Group']=='Alps Linked'
K=fit_df.loc[idx,'K'].to_numpy()
n=fit_df.loc[idx,'n'].to_numpy()
for i in range(len(model_list2L)):
    mObj=st.Stim1D(os.path.join(master_location,model_list2L[i]))
    [ts,u,dmxz,dmnz,mxdz]=mObj.calculate_ss()
    mxix=np.nonzero(dmxz<mx_thresh)[0][:1][0]
    mnix=np.nonzero(dmnz<mn_thresh)[0][:1][0]
    ss_mxz2L.append(ts[mxix])
    ss_mnz2L.append(ts[mnix])
    ts2L.append(ts)
    dmxz2L.append(dmxz)
    dmnz2L.append(dmnz)
    mxdz2L.append(mxdz)
    u2L.append(u)
    TU=response_time2(K,n,n/2,ka,h,xc,L,u,ksn0)/1e6
    ss_A2L.append(TU)

# BC
ts3U=[]
dmxz3U=[]
dmnz3U=[]
mxdz3U=[]
u3U=[]
ss_mxz3U=[]
ss_mnz3U=[]
ss_A3U=[]
idx=fit_df['Group']=='BC Unlinked'
K=fit_df.loc[idx,'K'].to_numpy()
n=fit_df.loc[idx,'n'].to_numpy()
for i in range(len(model_list3U)):
    mObj=st.Stim1D(os.path.join(master_location,model_list3U[i]))
    [ts,u,dmxz,dmnz,mxdz]=mObj.calculate_ss()
    mxix=np.nonzero(dmxz<mx_thresh)[0][:1][0]
    mnix=np.nonzero(dmnz<mn_thresh)[0][:1][0]
    ss_mxz3U.append(ts[mxix])
    ss_mnz3U.append(ts[mnix])
    ts3U.append(ts)
    dmxz3U.append(dmxz)
    dmnz3U.append(dmnz)
    mxdz3U.append(mxdz)
    u3U.append(u)
    TU=response_time2(K,n,n/2,ka,h,xc,L,u,ksn0)/1e6
    ss_A3U.append(TU)
    
ts3L=[]
dmxz3L=[]
dmnz3L=[]
mxdz3L=[]
u3L=[]
ss_mxz3L=[]
ss_mnz3L=[]
ss_A3L=[]
idx=fit_df['Group']=='BC Linked'
K=fit_df.loc[idx,'K'].to_numpy()
n=fit_df.loc[idx,'n'].to_numpy()
for i in range(len(model_list3L)):
    mObj=st.Stim1D(os.path.join(master_location,model_list3L[i]))
    [ts,u,dmxz,dmnz,mxdz]=mObj.calculate_ss()
    mxix=np.nonzero(dmxz<mx_thresh)[0][:1][0]
    mnix=np.nonzero(dmnz<mn_thresh)[0][:1][0]
    ss_mxz3L.append(ts[mxix])
    ss_mnz3L.append(ts[mnix])
    ts3L.append(ts)
    dmxz3L.append(dmxz)
    dmnz3L.append(dmnz)
    mxdz3L.append(mxdz)
    u3L.append(u)
    TU=response_time2(K,n,n/2,ka,h,xc,L,u,ksn0)/1e6
    ss_A3L.append(TU)

# Generate supplemental figure
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

time_ticks=np.arange(0,12,2)

# Start Plot    
f1=plt.figure(figsize=(8,3))
f1.set_dpi(250)
ax1=plt.subplot(1,3,1)
ax1.set_xlabel('Response Time [Myr]')
ax1.set_ylabel('Uplift Rate [mm/yr]')
ax1.set_yscale('log')
ax1.set_yticks([0.25,0.5,1,2,4,8])
ax1.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax1.set_xticks(time_ticks)
ax1.set_xlim((0,11))
plt.title(r'$\Delta$ Maximum Elevation')

ax2=plt.subplot(1,3,2)
ax2.set_xlabel('Response Time [Myr]')
ax2.set_ylabel('Uplift Rate [mm/yr]')
ax2.set_yscale('log')
ax2.set_yticks([0.25,0.5,1,2,4,8])
ax2.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax2.set_xticks(time_ticks)
ax2.set_xlim((0,11))
plt.title(r'$\Delta$ Mean Elevation')

ax3=plt.subplot(1,3,3)
ax3.set_xlabel(r'$\Delta$ Max [Myr]')
ax3.set_ylabel(r'$\Delta$ Mean [Myr]')
plt.title('Response Time')
ax3.plot([0,11],[0,11],c='k',linestyle=':')
ax3.set_xticks(time_ticks)
ax3.set_xlim((0,11))
ax3.set_yticks(time_ticks)
ax3.set_ylim((0,11))


for i in range(len(model_list1L)):
    # Max Z
    ax1.scatter(ss_mxz1L[i],u1L[i]*1000,c='w',marker='s',edgecolor=gc_col)
    ax1.scatter(ss_mxz2L[i],u2L[i]*1000,c='w',marker='s',edgecolor=alps_col)
    ax1.scatter(ss_mxz3L[i],u3L[i]*1000,c='w',marker='s',edgecolor=bc_col)
    # Mean Z
    ax2.scatter(ss_mnz1L[i],u1L[i]*1000,c='w',marker='s',edgecolor=gc_col)
    ax2.scatter(ss_mnz2L[i],u2L[i]*1000,c='w',marker='s',edgecolor=alps_col)
    ax2.scatter(ss_mnz3L[i],u3L[i]*1000,c='w',marker='s',edgecolor=bc_col)
    # Compare
    ax3.scatter(ss_mxz1L[i],ss_mnz1L[i],c='w',marker='s',edgecolor=gc_col)
    ax3.scatter(ss_mxz2L[i],ss_mnz2L[i],c='w',marker='s',edgecolor=alps_col)
    ax3.scatter(ss_mxz3L[i],ss_mnz3L[i],c='w',marker='s',edgecolor=bc_col)

for i in range(len(model_list1U)):
    # Max Z
    ax1.scatter(ss_mxz1U[i],u1U[i]*1000,c=gc_col,marker='o')
    ax1.scatter(ss_mxz2U[i],u2U[i]*1000,c=alps_col,marker='o')
    ax1.scatter(ss_mxz3U[i],u3U[i]*1000,c=bc_col,marker='o')
    # Mean Z
    ax2.scatter(ss_mnz1U[i],u1U[i]*1000,c=gc_col,marker='o')
    ax2.scatter(ss_mnz2U[i],u2U[i]*1000,c=alps_col,marker='o')
    ax2.scatter(ss_mnz3U[i],u3U[i]*1000,c=bc_col,marker='o')
    # Compare
    ax3.scatter(ss_mxz1U[i],ss_mnz1U[i],c=gc_col,marker='o')
    ax3.scatter(ss_mxz2U[i],ss_mnz2U[i],c=alps_col,marker='o')
    ax3.scatter(ss_mxz3U[i],ss_mnz3U[i],c=bc_col,marker='o')

ax1.text(0.92, 0.99, 'A',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax1.transAxes,
        fontsize=12,fontweight='extra bold')

ax2.text(0.92, 0.99, 'B',
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


f2=plt.figure(figsize=(8,3))
f2.set_dpi(250)
ax1=plt.subplot(1,3,1)
ax1.set_xlabel('Response Time [Myr]')
ax1.set_ylabel('Uplift Rate [mm/yr]')
ax1.set_yscale('log')
ax1.set_yticks([0.25,0.5,1,2,4,8])
ax1.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax1.set_xticks(time_ticks)
ax1.set_xlim((0,11))
plt.title(r'$\Delta$ Maximum Elevation')

ax2=plt.subplot(1,3,2)
ax2.set_xlabel('Response Time [Myr]')
ax2.set_yscale('log')
ax2.set_yticks([0.25,0.5,1,2,4,8])
ax2.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax2.set_xticks(time_ticks)
ax2.set_xlim((0,11))
ax2.set_ylabel('Uplift Rate [mm/yr]')
plt.title('Analytical Solution')

ax3=plt.subplot(1,3,3)
ax3.set_xlabel(r'$\Delta$ Max [Myr]')
ax3.set_ylabel('Analtyical Solution [Myr]')
plt.title('Response Time')
ax3.plot([0,11],[0,11],c='k',linestyle=':')
ax3.set_xticks(time_ticks)
ax3.set_xlim((0,11))
ax3.set_yticks(time_ticks)
ax3.set_ylim((0,11))

for i in range(len(model_list1L)):
    if i==0:
        # Max Z
        ax1.scatter(ss_mxz1L[i],u1L[i]*1000,c='w',marker='s',edgecolor=gc_col,label='GC Linked')
        ax1.scatter(ss_mxz2L[i],u2L[i]*1000,c='w',marker='s',edgecolor=alps_col,label='Alps Linked')
        ax1.scatter(ss_mxz3L[i],u3L[i]*1000,c='w',marker='s',edgecolor=bc_col,label='BC Linked')
        # Analytical
        ax2.scatter(ss_A1L[i],u1L[i]*1000,c='w',marker='s',edgecolor=gc_col)
        ax2.scatter(ss_A2L[i],u2L[i]*1000,c='w',marker='s',edgecolor=alps_col)
        ax2.scatter(ss_A3L[i],u3L[i]*1000,c='w',marker='s',edgecolor=bc_col)
        # Compare
        ax3.scatter(ss_mxz1L[i],ss_A1L[i],c='w',marker='s',edgecolor=gc_col)
        ax3.scatter(ss_mxz2L[i],ss_A2L[i],c='w',marker='s',edgecolor=alps_col)
        ax3.scatter(ss_mxz3L[i],ss_A3L[i],c='w',marker='s',edgecolor=bc_col)
        
    else:
        # Max Z
        ax1.scatter(ss_mxz1L[i],u1L[i]*1000,c='w',marker='s',edgecolor=gc_col)
        ax1.scatter(ss_mxz2L[i],u2L[i]*1000,c='w',marker='s',edgecolor=alps_col)
        ax1.scatter(ss_mxz3L[i],u3L[i]*1000,c='w',marker='s',edgecolor=bc_col)
        # Analytical
        ax2.scatter(ss_A1L[i],u1L[i]*1000,c='w',marker='s',edgecolor=gc_col)
        ax2.scatter(ss_A2L[i],u2L[i]*1000,c='w',marker='s',edgecolor=alps_col)
        ax2.scatter(ss_A3L[i],u3L[i]*1000,c='w',marker='s',edgecolor=bc_col)
        # Compare
        ax3.scatter(ss_mxz1L[i],ss_A1L[i],c='w',marker='s',edgecolor=gc_col)
        ax3.scatter(ss_mxz2L[i],ss_A2L[i],c='w',marker='s',edgecolor=alps_col)
        ax3.scatter(ss_mxz3L[i],ss_A3L[i],c='w',marker='s',edgecolor=bc_col)

for i in range(len(model_list1U)):
    if i==0:
        # Max Z
        ax1.scatter(ss_mxz1U[i],u1U[i]*1000,c=gc_col,marker='o')
        ax1.scatter(ss_mxz2U[i],u2U[i]*1000,c=alps_col,marker='o')
        ax1.scatter(ss_mxz3U[i],u3U[i]*1000,c=bc_col,marker='o')
        # Analytical
        ax2.scatter(ss_A1U[i],u1U[i]*1000,c=gc_col,marker='o',label='GC Unlinked')
        ax2.scatter(ss_A2U[i],u2U[i]*1000,c=alps_col,marker='o',label='Alps Unlinked')
        ax2.scatter(ss_A3U[i],u3U[i]*1000,c=bc_col,marker='o',label='BC Unlinked')
        # Compare
        ax3.scatter(ss_mxz1U[i],ss_A1U[i],c=gc_col,marker='o')
        ax3.scatter(ss_mxz2U[i],ss_A2U[i],c=alps_col,marker='o')
        ax3.scatter(ss_mxz3U[i],ss_A3U[i],c=bc_col,marker='o')        
    else:
        # Max Z
        ax1.scatter(ss_mxz1U[i],u1U[i]*1000,c=gc_col,marker='o')
        ax1.scatter(ss_mxz2U[i],u2U[i]*1000,c=alps_col,marker='o')
        ax1.scatter(ss_mxz3U[i],u3U[i]*1000,c=bc_col,marker='o')
        # Analytical
        ax2.scatter(ss_A1U[i],u1U[i]*1000,c=gc_col,marker='o')
        ax2.scatter(ss_A2U[i],u2U[i]*1000,c=alps_col,marker='o')
        ax2.scatter(ss_A3U[i],u3U[i]*1000,c=bc_col,marker='o')
        # Compare
        ax3.scatter(ss_mxz1U[i],ss_A1U[i],c=gc_col,marker='o')
        ax3.scatter(ss_mxz2U[i],ss_A2U[i],c=alps_col,marker='o')
        ax3.scatter(ss_mxz3U[i],ss_A3U[i],c=bc_col,marker='o')

ax1.text(0.92, 0.99, 'A',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax1.transAxes,
        fontsize=12,fontweight='extra bold')
ax1.legend(loc='center right')

ax2.text(0.92, 0.99, 'B',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax2.transAxes,
        fontsize=12,fontweight='extra bold')
ax2.legend(loc='center right')

ax3.text(0.01, 0.99, 'C',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax3.transAxes,
        fontsize=12,fontweight='extra bold') 

plt.tight_layout()

f3=plt.figure(figsize=(8,3))
f3.set_dpi(250)
ax1=plt.subplot(1,3,1)
ax1.set_xlabel('Response Time [Myr]')
ax1.set_ylabel('Uplift Rate [mm/yr]')
ax1.set_yscale('log')
ax1.set_yticks([0.25,0.5,1,2,4,8])
ax1.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax1.set_xticks(time_ticks)
ax1.set_xlim((0,11))
plt.title(r'$\Delta$ Mean Elevation')

ax2=plt.subplot(1,3,2)
ax2.set_xlabel('Response Time [Myr]')
ax2.set_yscale('log')
ax2.set_yticks([0.25,0.5,1,2,4,8])
ax2.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax2.set_xticks(time_ticks)
ax2.set_xlim((0,11))
ax2.set_ylabel('Uplift Rate [mm/yr]')
plt.title('Analytical Solution')

ax3=plt.subplot(1,3,3)
ax3.set_xlabel(r'$\Delta$ Mean [Myr]')
ax3.set_ylabel('Analtyical Solution [Myr]')
plt.title('Response Time')
ax3.plot([0,11],[0,11],c='k',linestyle=':')
ax3.set_xticks(time_ticks)
ax3.set_xlim((0,11))
ax3.set_yticks(time_ticks)
ax3.set_ylim((0,11))

for i in range(len(model_list1L)):
    if i==0:
        # Max Z
        ax1.scatter(ss_mnz1L[i],u1L[i]*1000,c='w',marker='s',edgecolor=gc_col,label='GC-Linked')
        ax1.scatter(ss_mnz2L[i],u2L[i]*1000,c='w',marker='s',edgecolor=alps_col,label='Alps-Linked')
        ax1.scatter(ss_mnz3L[i],u3L[i]*1000,c='w',marker='s',edgecolor=bc_col,label='BC-Linked')
        # Analytical
        ax2.scatter(ss_A1L[i],u1L[i]*1000,c='w',marker='s',edgecolor=gc_col)
        ax2.scatter(ss_A2L[i],u2L[i]*1000,c='w',marker='s',edgecolor=alps_col)
        ax2.scatter(ss_A3L[i],u3L[i]*1000,c='w',marker='s',edgecolor=bc_col)
        # Compare
        ax3.scatter(ss_mnz1L[i],ss_A1L[i],c='w',marker='s',edgecolor=gc_col)
        ax3.scatter(ss_mnz2L[i],ss_A2L[i],c='w',marker='s',edgecolor=alps_col)
        ax3.scatter(ss_mnz3L[i],ss_A3L[i],c='w',marker='s',edgecolor=bc_col)
        
    else:
        # Max Z
        ax1.scatter(ss_mnz1L[i],u1L[i]*1000,c='w',marker='s',edgecolor=gc_col)
        ax1.scatter(ss_mnz2L[i],u2L[i]*1000,c='w',marker='s',edgecolor=alps_col)
        ax1.scatter(ss_mnz3L[i],u3L[i]*1000,c='w',marker='s',edgecolor=bc_col)
        # Analytical
        ax2.scatter(ss_A1L[i],u1L[i]*1000,c='w',marker='s',edgecolor=gc_col)
        ax2.scatter(ss_A2L[i],u2L[i]*1000,c='w',marker='s',edgecolor=alps_col)
        ax2.scatter(ss_A3L[i],u3L[i]*1000,c='w',marker='s',edgecolor=bc_col)
        # Compare
        ax3.scatter(ss_mnz1L[i],ss_A1L[i],c='w',marker='s',edgecolor=gc_col)
        ax3.scatter(ss_mnz2L[i],ss_A2L[i],c='w',marker='s',edgecolor=alps_col)
        ax3.scatter(ss_mnz3L[i],ss_A3L[i],c='w',marker='s',edgecolor=bc_col)

for i in range(len(model_list1U)):
    if i==0:
        # Max Z
        ax1.scatter(ss_mnz1U[i],u1U[i]*1000,c=gc_col,marker='o')
        ax1.scatter(ss_mnz2U[i],u2U[i]*1000,c=alps_col,marker='o')
        ax1.scatter(ss_mnz3U[i],u3U[i]*1000,c=bc_col,marker='o')
        # Analytical
        ax2.scatter(ss_A1U[i],u1U[i]*1000,c=gc_col,marker='o',label='GC-Unlinked')
        ax2.scatter(ss_A2U[i],u2U[i]*1000,c=alps_col,marker='o',label='Alps-Unlinked')
        ax2.scatter(ss_A3U[i],u3U[i]*1000,c=bc_col,marker='o',label='BC-Unlinked')
        # Compare
        ax3.scatter(ss_mnz1U[i],ss_A1U[i],c=gc_col,marker='o')
        ax3.scatter(ss_mnz2U[i],ss_A2U[i],c=alps_col,marker='o')
        ax3.scatter(ss_mnz3U[i],ss_A3U[i],c=bc_col,marker='o')        
    else:
        # Max Z
        ax1.scatter(ss_mnz1U[i],u1U[i]*1000,c=gc_col,marker='o')
        ax1.scatter(ss_mnz2U[i],u2U[i]*1000,c=alps_col,marker='o')
        ax1.scatter(ss_mnz3U[i],u3U[i]*1000,c=bc_col,marker='o')
        # Analytical
        ax2.scatter(ss_A1U[i],u1U[i]*1000,c=gc_col,marker='o')
        ax2.scatter(ss_A2U[i],u2U[i]*1000,c=alps_col,marker='o')
        ax2.scatter(ss_A3U[i],u3U[i]*1000,c=bc_col,marker='o')
        # Compare
        ax3.scatter(ss_mnz1U[i],ss_A1U[i],c=gc_col,marker='o')
        ax3.scatter(ss_mnz2U[i],ss_A2U[i],c=alps_col,marker='o')
        ax3.scatter(ss_mnz3U[i],ss_A3U[i],c=bc_col,marker='o')

ax1.text(0.92, 0.99, 'A',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax1.transAxes,
        fontsize=12,fontweight='extra bold')
ax1.legend(loc='center right')

ax2.text(0.92, 0.99, 'B',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax2.transAxes,
        fontsize=12,fontweight='extra bold')
ax2.legend(loc='center right')

ax3.text(0.01, 0.99, 'C',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax3.transAxes,
        fontsize=12,fontweight='extra bold') 
    
plt.tight_layout()
plt.rcdefaults()

f2.savefig('figure_11.pdf',dpi="figure")

