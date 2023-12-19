#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 09:56:06 2023

@author: aforte
"""

import sys
sys.path.insert(0,'/Users/aforte/Documents/GitHub/snowmelt_orography')
import stimpy as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from cmcrameri import cm

fit_df=pd.read_csv('ksn_e_fit.csv')


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

# Generate supplemental figure
f1=plt.figure(figsize=(8,4))
f1.set_dpi(250)
e_vec=np.logspace(1,4,100)
gc_ero=pd.read_csv('/Users/aforte/Documents/GitHub/Caucasus_Erosion/data_tables/gc_ero_master_table.csv')

mn_ksn=gc_ero['mean_ksn'].to_numpy()
se_ksn=gc_ero['se_ksn'].to_numpy()
mn_E=gc_ero['St_E_rate_m_Myr'].to_numpy()
se_E=gc_ero['St_Ext_Unc'].to_numpy()
da=gc_ero['drainage_area'].to_numpy()

sc1=plt.scatter(mn_E,mn_ksn,s=30,c=np.log10(da),label='Observed Erosion Rates',zorder=1,cmap=cm.lajolla,edgecolors='k')
plt.errorbar(mn_E,mn_ksn,se_ksn,se_E,ecolor='gray',linestyle='',zorder=0,elinewidth=0.5)
plt.plot(e_vec,fit_df.loc[0,'C']*(e_vec)**fit_df.loc[0,'phi'],c='k',label=r'Unlinked 50 km (LogDA = 2.69 $km^{2}$)')
plt.plot(e_vec,fit_df.loc[1,'C']*(e_vec)**fit_df.loc[1,'phi'],c='k',linestyle='--',label=r'Linked 50 km (LogDA = 2.69 $km^{2}$)')
plt.plot(e_vec,fit_df.loc[2,'C']*(e_vec)**fit_df.loc[2,'phi'],c='gray',linestyle='-',label=r'Unlinked 10 km (LogDA = 1.41 $km^{2}$)')
plt.plot(e_vec,fit_df.loc[3,'C']*(e_vec)**fit_df.loc[3,'phi'],c='gray',linestyle='--',label=r'Linked 10 km (LogDA = 1.41 $km^{2}$)')
plt.plot(e_vec,fit_df.loc[10,'C']*(e_vec)**fit_df.loc[10,'phi'],c='k',linestyle=':',label=r'Unlinked Area (LogDA = 2.69 $km^{2}$)')
plt.plot(e_vec,fit_df.loc[11,'C']*(e_vec)**fit_df.loc[11,'phi'],c='gray',linestyle=':',label=r'Linked Area (LogDA = 2.69 $km^{2}$)')
plt.xlabel('Erosion Rate [m/Myr]')
plt.ylabel(r'$k_{sn}$ [m]')
plt.legend(bbox_to_anchor= (1.3,0.99),loc='upper left')
plt.xscale('log')
plt.yscale('log')
ax1=plt.gca()
cbar1=plt.colorbar(sc1,ax=ax1)
cbar1.ax.set_ylabel(r'Log Drainage Area [$km^{2}$]')

plt.tight_layout()
plt.rcdefaults()

f1.savefig('P2_SUP_figureS5.pdf',dpi='figure')