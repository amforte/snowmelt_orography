#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 07:27:04 2022

@author: aforte
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress


repo_location='/Users/aforte/Documents/GitHub/snowmelt_orography/geospatial_codes/'
gc=pd.read_csv(repo_location+'gc_ksn_rlf.csv')
alps=pd.read_csv(repo_location+'alps_ksn_rlf.csv')
bc=pd.read_csv(repo_location+'bc_ksn_rlf.csv')

gc_col='black'
alps_col='royalblue'
bc_col='orange'

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


f1=plt.figure(figsize=(8,8))
f1.set_dpi(250)

ax1=plt.subplot(2,2,1)
plt.scatter(gc['mean_ksn'],gc['mean_rlf2500'],s=2,c=gc_col,label='Greater Caucasus')
plt.scatter(alps['mean_ksn'],alps['mean_rlf2500'],s=2,c=alps_col,label='Alps')
plt.scatter(bc['mean_ksn'],bc['mean_rlf2500'],s=2,c=bc_col,label='BC')
plt.ylim((0,2500))
plt.xlabel(r'Mean $k_{sn}$ [m]')
plt.ylabel('Mean Local 2.5 km Relief [m]')
plt.legend(loc=4)
plt.title('All Basins')
ax1.text(0.01, 0.99, 'A',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax1.transAxes,
        fontsize=12,fontweight='extra bold')

chi_R2=0.90
gc_idx=gc['chi_R_squared']>=chi_R2
alps_idx=alps['chi_R_squared']>=chi_R2
bc_idx=bc['chi_R_squared']>=chi_R2

[gc_slp,gc_int,gc_r,_,_]=linregress(gc.loc[gc_idx,'mean_ksn'],gc.loc[gc_idx,'mean_rlf2500'])
[alps_slp,alps_int,alps_r,_,_]=linregress(alps.loc[alps_idx,'mean_ksn'],alps.loc[alps_idx,'mean_rlf2500'])
[bc_slp,bc_int,bc_r,_,_]=linregress(bc.loc[bc_idx,'mean_ksn'],bc.loc[bc_idx,'mean_rlf2500'])

gc_r2=str(np.round(gc_r**2,2))
alps_r2=str(np.round(alps_r**2,2))
bc_r2=str(np.round(bc_r**2,2))
gc_slp_s=str(np.round(gc_slp,2))
alps_slp_s=str(np.round(alps_slp,2))
bc_slp_s=str(np.round(bc_slp,2))
gc_int_s=str(np.round(gc_int,2))
alps_int_s=str(np.round(alps_int,2))
bc_int_s=str(np.round(bc_int,2))

gc_label=r'LR = '+gc_slp_s+'$k_{sn}$ + '+gc_int_s+'; $R^{2}$ = '+gc_r2
alps_label=r'LR = '+alps_slp_s+'$k_{sn}$ + '+alps_int_s+'; $R^{2}$ = '+alps_r2
bc_label=r'LR = '+bc_slp_s+'$k_{sn}$ + '+bc_int_s+'; $R^{2}$ = '+bc_r2

ax3=plt.subplot(2,2,3)
ksn_vec=np.linspace(0,800,100)
plt.scatter(gc.loc[gc_idx,'mean_ksn'],gc.loc[gc_idx,'mean_rlf2500'],s=2,c=gc_col)
plt.scatter(alps.loc[alps_idx,'mean_ksn'],alps.loc[alps_idx,'mean_rlf2500'],s=2,c=alps_col)
plt.scatter(bc.loc[bc_idx,'mean_ksn'],bc.loc[bc_idx,'mean_rlf2500'],s=2,c=bc_col)
plt.plot(ksn_vec,gc_slp*ksn_vec + gc_int,c=gc_col,label=gc_label)
plt.plot(ksn_vec,alps_slp*ksn_vec + alps_int,c=alps_col,label=alps_label)
plt.plot(ksn_vec,bc_slp*ksn_vec + bc_int,c=bc_col,label=bc_label)
plt.ylim((0,2500))
plt.xlabel(r'Mean $k_{sn}$ [m]')
plt.ylabel('Mean Local 2.5 km Relief [m]')
plt.legend(loc=4)
plt.title(r'Filtered to $\chi$ $R^{2}$ > 0.90')
ax3.text(0.01, 0.99, 'C',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax3.transAxes,
        fontsize=12,fontweight='extra bold')

ax2=plt.subplot(2,2,2)
plt.scatter(gc['mean_gradient'],gc['mean_rlf2500'],s=2,c=gc_col,label='Greater Caucasus')
plt.scatter(alps['mean_gradient'],alps['mean_rlf2500'],s=2,c=alps_col,label='Alps')
plt.scatter(bc['mean_gradient'],bc['mean_rlf2500'],s=2,c=bc_col,label='BC')
plt.ylim((0,2500))
plt.xlabel('Mean Gradient [m/m]')
plt.ylabel('Mean Local 2.5 km Relief [m]')
# plt.legend(loc=4)
plt.title('All Basins')
ax2.text(0.01, 0.99, 'B',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax2.transAxes,
        fontsize=12,fontweight='extra bold')


ax4=plt.subplot(2,2,4)
plt.scatter(gc.loc[gc_idx,'mean_gradient'],gc.loc[gc_idx,'mean_rlf2500'],s=2,c=gc_col,label='Greater Caucasus')
plt.scatter(alps.loc[alps_idx,'mean_gradient'],alps.loc[alps_idx,'mean_rlf2500'],s=2,c=alps_col,label='Alps')
plt.scatter(bc.loc[bc_idx,'mean_gradient'],bc.loc[bc_idx,'mean_rlf2500'],s=2,c=bc_col,label='BC')
plt.ylim((0,2500))
plt.xlabel('Mean Gradient [m/m]')
plt.ylabel('Mean Local 2.5 km Relief [m]')
# plt.legend(loc=4)
plt.title(r'Filtered to $\chi$ $R^{2}$ > 0.90')
ax4.text(0.01, 0.99, 'D',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax4.transAxes,
        fontsize=12,fontweight='extra bold')


plt.tight_layout()
plt.rcdefaults()
f1.savefig('figure_s2.pdf',dpi="figure")







