#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script written by Adam M. Forte
aforte8@lsu.edu

This script reads the outputs from 'generate_ksn_relief_values.m' and generates
relationships between ksn and relief for use within the stimpy model.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

gc=pd.read_csv('gc_ksn_rlf.csv')
alps=pd.read_csv('alps_ksn_rlf.csv')
bc=pd.read_csv('bc_ksn_rlf.csv')

gb=pd.read_csv('/Users/aforte/Documents/Research/Grants/NSF_2023_Snowmelt/GrtBsn/grt_bsn_sheds.csv')


f1=plt.figure(figsize=(8,10))

plt.subplot(2,1,1)
plt.scatter(gc['mean_ksn'],gc['mean_rlf2500'],s=10,c='k',label='Greater Caucasus')
plt.scatter(alps['mean_ksn'],alps['mean_rlf2500'],s=10,c='royalblue',label='Alps')
plt.scatter(bc['mean_ksn'],bc['mean_rlf2500'],s=10,c='orange',label='BC')
plt.scatter(gb['mean_ksn'],gb['mean_rlf2500'],s=10,c='green',label='Great Basin')
plt.ylim((0,2500))
plt.xlabel(r'Mean $k_{sn}$ [m]')
plt.ylabel('Mean Local 2.5 km Relief [m]')
plt.legend(loc='best')
plt.title('All Basins')

chi_R2=0.90
gc_idx=gc['chi_R_squared']>=chi_R2
alps_idx=alps['chi_R_squared']>=chi_R2
bc_idx=bc['chi_R_squared']>=chi_R2
gb_idx=gb['chi_R_squared']>=chi_R2

[gc_slp,gc_int,gc_r,_,_]=linregress(gc.loc[gc_idx,'mean_ksn'],gc.loc[gc_idx,'mean_rlf2500'])
[alps_slp,alps_int,alps_r,_,_]=linregress(alps.loc[alps_idx,'mean_ksn'],alps.loc[alps_idx,'mean_rlf2500'])
[bc_slp,bc_int,bc_r,_,_]=linregress(bc.loc[bc_idx,'mean_ksn'],bc.loc[bc_idx,'mean_rlf2500'])
[gb_slp,gb_int,gb_r,_,_]=linregress(gb.loc[gb_idx,'mean_ksn'],gb.loc[gb_idx,'mean_rlf2500'])

gc_r2=str(np.round(gc_r**2,2))
alps_r2=str(np.round(alps_r**2,2))
bc_r2=str(np.round(bc_r**2,2))
gb_r2=str(np.round(gb_r**2,2))

# gc_ksn=gc.loc[gc_idx,'mean_ksn'].to_numpy()
# gc_rlf=gc.loc[gc_idx,'mean_rlf2500'].to_numpy()
# alps_ksn=alps.loc[alps_idx,'mean_ksn'].to_numpy()
# alps_rlf=alps.loc[alps_idx,'mean_rlf2500'].to_numpy()
# bc_ksn=bc.loc[bc_idx,'mean_ksn'].to_numpy()
# bc_rlf=bc.loc[bc_idx,'mean_rlf2500'].to_numpy()

# ksn=np.concatenate((gc_ksn,alps_ksn,bc_ksn),axis=0)
# rlf=np.concatenate((gc_rlf,alps_rlf,bc_rlf),axis=0)

# [slp,inte,rval,_,_]=linregress(ksn,rlf)
# r2=str(np.round(rval**2,2))

ksn_vec=np.linspace(0,800,100)

plt.subplot(2,1,2)
plt.scatter(gc.loc[gc_idx,'mean_ksn'],gc.loc[gc_idx,'mean_rlf2500'],s=10,c='k',label='Greater Caucasus')
plt.scatter(alps.loc[alps_idx,'mean_ksn'],alps.loc[alps_idx,'mean_rlf2500'],s=10,c='royalblue',label='Alps')
plt.scatter(bc.loc[bc_idx,'mean_ksn'],bc.loc[bc_idx,'mean_rlf2500'],s=10,c='orange',label='BC')
plt.scatter(gb.loc[gb_idx,'mean_ksn'],gb.loc[gb_idx,'mean_rlf2500'],s=10,c='green',label='Great Basin')

plt.plot(ksn_vec,gc_slp*ksn_vec + gc_int,c='k',label=r'$R^{2}$ = '+gc_r2)
plt.plot(ksn_vec,alps_slp*ksn_vec + alps_int,c='royalblue',label=r'$R^{2}$ = '+alps_r2)
plt.plot(ksn_vec,bc_slp*ksn_vec + bc_int,c='orange',label=r'$R^{2}$ = '+bc_r2)
plt.plot(ksn_vec,gb_slp*ksn_vec + gb_int,c='green',label=r'$R^{2}$ = '+gb_r2)

# plt.plot(ksn_vec,slp*ksn_vec + inte,c='g',linestyle='--',linewidth=2,label=r'All - $R^{2}$ = '+r2)

plt.ylim((0,2500))
plt.xlabel(r'Mean $k_{sn}$ [m]')
plt.ylabel('Mean Local 2.5 km Relief [m]')
plt.legend(loc='best')
plt.title(r'Filtered to $\chi$ $R^{2}$ > 0.90')
plt.tight_layout()

f1.savefig('REVIEW_ksn_to_relief.pdf')












