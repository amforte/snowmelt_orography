#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  6 10:57:10 2023

@author: aforte
"""

import sys
sys.path.insert(0,'/Users/aforte/Documents/GitHub/snowmelt-tectonics')
import stimpy as st
import numpy as np
import matplotlib.pyplot as plt
from cmcrameri import cm
from matplotlib import colors
from matplotlib import cm as cmm

# master_location='/Volumes/Choruh/Data/snowmelt_project/model_outputs_v2/'
master_location='/Users/aforte/Documents/Python/snowmelt/model_outputs_v2'


clip_ts=[np.inf,np.inf,np.inf,np.inf,np.inf,np.inf]
bn=np.linspace(0,2500,50)
odds=np.arange(1,13,2)
evens=np.arange(2,13,2)

f1,ax1=plt.subplots(ncols=2,nrows=6,figsize=(6,10),layout='tight')
f2,ax2=plt.subplots(ncols=2,nrows=6,figsize=(6,10),layout='tight')
f3,ax3=plt.subplots(ncols=2,nrows=6,figsize=(6,10),layout='tight')

U=['025','05','1','2','4','8']
UT=['250 m/Myr','500 m/Myr','1000 m/Myr','2000 m/Myr','4000 m/My','8000 m/Myr']
for i in range(len(U)):
    # Build model lists
    model_list=['gc'+U[i]+'u','a'+U[i]+'u','bc'+U[i]+'u','gc'+U[i]+'l','a'+U[i]+'l','bc'+U[i]+'l']
    prefix_list=['ts_','ts_','ts_','ts_','ts_','ts_']
    descript_list=['GC'+U[i]+'U','A'+U[i]+'U','BC'+U[i]+'U','GC'+U[i]+'L','A'+U[i]+'L','BC'+U[i]+'L']
    mc=st.ModelComparison(master_location,model_list,prefix_list,descript_list)
    d=mc.output_result_dict(clip_ts)
    
    
    
    [c,b]=np.histogram(d[0]['rlf_Z'].ravel(),bn)
    ax1[i,0].hist(b[:-1],b,weights=c/np.sum(c),color='k',alpha=0.5,label='GC')
    [c,b]=np.histogram(d[1]['rlf_Z'].ravel(),bn)
    ax1[i,0].hist(b[:-1],b,weights=c/np.sum(c),color='royalblue',alpha=0.5,label='Alps')
    [c,b]=np.histogram(d[2]['rlf_Z'].ravel(),bn)
    ax1[i,0].hist(b[:-1],b,weights=c/np.sum(c),color='orange',alpha=0.5,label='BC')
    ax1[i,0].set_title(UT[i]+'; Unlinked')
    ax1[i,0].set_ylabel('Normalized Count')
    ax1[i,0].set_yscale('log')
    ax1[i,0].set_ylim((10**-4,1))
    if i==5:
        ax1[i,0].set_xlabel('Local Relief')
    
    [c,b]=np.histogram(d[3]['rlf_Z'].ravel(),bn)
    ax1[i,1].hist(b[:-1],b,weights=c/np.sum(c),color='k',alpha=0.5,label='GC')
    [c,b]=np.histogram(d[4]['rlf_Z'].ravel(),bn)
    ax1[i,1].hist(b[:-1],b,weights=c/np.sum(c),color='royalblue',alpha=0.5,label='Alps')
    [c,b]=np.histogram(d[5]['rlf_Z'].ravel(),bn)
    ax1[i,1].hist(b[:-1],b,weights=c/np.sum(c),color='orange',alpha=0.5,label='BC')
    ax1[i,1].set_title(UT[i]+'; Linked')
    ax1[i,1].legend(loc='best')
    ax1[i,1].set_yscale('log')
    ax1[i,1].set_ylim((10**-4,1))
    if i==5:
        ax1[i,1].set_xlabel('Local Relief')
    
    # Along profile relief change
    [c,b]=np.histogram(np.abs(np.diff(d[0]['rlf_Z'],axis=1)).ravel(),bn)
    ax2[i,0].hist(b[:-1],b,weights=c/np.sum(c),color='k',alpha=0.5,label='GC')
    [c,b]=np.histogram(np.abs(np.diff(d[1]['rlf_Z'],axis=1)).ravel(),bn)
    ax2[i,0].hist(b[:-1],b,weights=c/np.sum(c),color='royalblue',alpha=0.5,label='Alps')
    [c,b]=np.histogram(np.abs(np.diff(d[2]['rlf_Z'],axis=1)).ravel(),bn)
    ax2[i,0].hist(b[:-1],b,weights=c/np.sum(c),color='orange',alpha=0.5,label='BC')    
    ax2[i,0].set_title(UT[i]+'; Unlinked')
    ax2[i,0].set_ylabel('Normalized Count')
    ax2[i,0].set_yscale('log')
    ax2[i,0].set_ylim((10**-4,1))
    if i==5:
        ax2[i,0].set_xlabel(r'$\Delta$ Local Relief in X')
    
    [c,b]=np.histogram(np.abs(np.diff(d[3]['rlf_Z'],axis=1)).ravel(),bn)
    ax2[i,1].hist(b[:-1],b,weights=c/np.sum(c),color='k',alpha=0.5,label='GC')
    [c,b]=np.histogram(np.abs(np.diff(d[4]['rlf_Z'],axis=1)).ravel(),bn)
    ax2[i,1].hist(b[:-1],b,weights=c/np.sum(c),color='royalblue',alpha=0.5,label='Alps')
    [c,b]=np.histogram(np.abs(np.diff(d[5]['rlf_Z'],axis=1)).ravel(),bn)
    ax2[i,1].hist(b[:-1],b,weights=c/np.sum(c),color='orange',alpha=0.5,label='BC')
    ax2[i,1].set_title(UT[i]+'; Linked')
    ax2[i,1].set_yscale('log')
    ax2[i,1].set_ylim((10**-4,1))
    # ax2[i,1].legend(loc='best')
    if i==5:
        ax2[i,1].set_xlabel(r'$\Delta$ Local Relief in X')        

    # Through time relief change
    [c,b]=np.histogram(np.abs(np.diff(d[0]['rlf_Z'],axis=0)).ravel(),bn)
    ax3[i,0].hist(b[:-1],b,weights=c/np.sum(c),color='k',alpha=0.5,label='GC')
    [c,b]=np.histogram(np.abs(np.diff(d[1]['rlf_Z'],axis=0)).ravel(),bn)
    ax3[i,0].hist(b[:-1],b,weights=c/np.sum(c),color='royalblue',alpha=0.5,label='Alps')
    [c,b]=np.histogram(np.abs(np.diff(d[2]['rlf_Z'],axis=0)).ravel(),bn)
    ax3[i,0].hist(b[:-1],b,weights=c/np.sum(c),color='orange',alpha=0.5,label='BC')    
    ax3[i,0].set_title(UT[i]+'; Unlinked')
    ax3[i,0].set_ylabel('Normalized Count')
    ax3[i,0].set_yscale('log')
    ax3[i,0].set_ylim((10**-4,1))
    if i==5:
        ax3[i,0].set_xlabel(r'$\Delta$ Local Relief Through Time')
    
    [c,b]=np.histogram(np.abs(np.diff(d[3]['rlf_Z'],axis=0)).ravel(),bn)
    ax3[i,1].hist(b[:-1],b,weights=c/np.sum(c),color='k',alpha=0.5,label='GC')
    [c,b]=np.histogram(np.abs(np.diff(d[4]['rlf_Z'],axis=0)).ravel(),bn)
    ax3[i,1].hist(b[:-1],b,weights=c/np.sum(c),color='royalblue',alpha=0.5,label='Alps')
    [c,b]=np.histogram(np.abs(np.diff(d[5]['rlf_Z'],axis=0)).ravel(),bn)
    ax3[i,1].hist(b[:-1],b,weights=c/np.sum(c),color='orange',alpha=0.5,label='BC')
    ax3[i,1].set_title(UT[i]+'; Linked')
    ax3[i,1].set_yscale('log')
    ax3[i,1].set_ylim((10**-4,1))
    # ax2[i,1].legend(loc='best')
    if i==5:
        ax3[i,1].set_xlabel(r'$\Delta$ Local Relief Through Time')         
        
    
plt.tight_layout()
f1.savefig('REVIEW_LocalRelief.pdf')
f2.savefig('REVEIW_LocalRelief_DeltaX.pdf')
f3.savefig('REVIEW_LocalRelief_DeltaY.pdf')       
    
    



