#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 18 08:38:11 2023

@author: aforte
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
from sklearn.ensemble import RandomForestRegressor
import warnings

# Filter warnings
warnings.filterwarnings('ignore',message='Singular matrix in solving dual problem')

# Define location of data
master_location='/Volumes/Choruh/Data/snowmelt_project/'

## Load global
df=pd.read_csv(master_location+'wrr2_derived_data_v3.csv')
df=df.drop(index=df.index[np.isnan(df['mean_z'])])
df=df.reset_index(drop=True)

# Calc percents
perc_base=df['qsb']/df['mean_runoff']

# Calculate indices
grlf=df['max_z']-df['min_z']

# Set cutoffs
percb_cutoff=0.25

# Index dataset and drop
rem_idx=(grlf<=500) | (df['mean_z']<=250) | (perc_base>percb_cutoff) | (df['mean_rlf']<250)
df_s=df.drop(df.index[rem_idx])
df_s=df_s.reset_index(drop=True)
perc_snow=df_s['qsm']/df_s['mean_runoff']

# Extract values of interest
lt=df_s['latitude'].to_numpy()
mz=df_s['max_z'].to_numpy()
mnz=df_s['mean_z'].to_numpy()
mrl=df_s['mean_rlf'].to_numpy()
cr=df_s['r_c1'].to_numpy()
mr=df_s['mean_runoff'].to_numpy()
cp=df_s['p_c'].to_numpy()
mp=df_s['mean_precip'].to_numpy()
ps=perc_snow.to_numpy()
et=df_s['et_mean'].to_numpy()
pet=df_s['pet_mean'].to_numpy()
mt=df_s['mat'].to_numpy()

dfc=pd.DataFrame(data={'Latitude':lt,
                       'Mean Z':mnz,
                       'Max Z':mz,
                       'Mean Rlf':mrl,
                       'Mean Runoff':mr,
                       'Mean Precip':mp,
                       'Precip Variability':cp,
                       'Runoff Variability':cr,
                       'ET':et,
                       'PET':pet,
                       'MAT':mt,
                       '% Snow':ps})

# Set length of records
l=len(mr)

# Random forest regression
X1=np.concatenate((mnz.reshape(l,1),mz.reshape(l,1),mrl.reshape(l,1),mt.reshape(l,1)),axis=1)
F1=['Mean Elevation', 'Max Elevation', 'Mean Relief','MAT']
X2=np.concatenate((mnz.reshape(l,1),mz.reshape(l,1),mrl.reshape(l,1),mp.reshape(l,1),cp.reshape(l,1),mt.reshape(l,1)),axis=1)
F2=['Mean Elevation', 'Max Elevation', 'Mean Relief','Mean Precip','Precip Shape','MAT']
X3=np.concatenate((mnz.reshape(l,1),mz.reshape(l,1),mrl.reshape(l,1),mr.reshape(l,1),mp.reshape(l,1),cp.reshape(l,1),mt.reshape(l,1),ps.reshape(l,1)),axis=1)
F3=['Mean Elevation', 'Max Elevation', 'Mean Relief','Mean Runoff','Mean Precip','Precip Shape','MAT','Snowmelt Fraction']
X4=np.concatenate((mnz.reshape(l,1),mz.reshape(l,1),mrl.reshape(l,1),mr.reshape(l,1),cr.reshape(l,1),mp.reshape(l,1),cp.reshape(l,1),mt.reshape(l,1)),axis=1)
F4=['Mean Elevation', 'Max Elevation', 'Mean Relief','Mean Runoff','Runoff Shape','Mean Precip','Precip Shape','MAT']
X5=np.concatenate((lt.reshape(l,1),mnz.reshape(l,1),mz.reshape(l,1),mrl.reshape(l,1)),axis=1)
F5=['Latitude','Mean Elevation', 'Max Elevation', 'Mean Relief']


regr1=RandomForestRegressor(random_state=0)
regr1.fit(X1,mp)
yr1=regr1.predict(X1)
rr1=regr1.score(X1,mp)

regr2=RandomForestRegressor(random_state=0)
regr2.fit(X2,mr)
yr2=regr2.predict(X2)
rr2=regr2.score(X2,mr)

regr3=RandomForestRegressor(random_state=0)
regr3.fit(X3,cr)
yr3=regr3.predict(X3)
rr3=regr3.score(X3,cr)

regr4=RandomForestRegressor(random_state=0)
regr4.fit(X4,ps)
yr4=regr4.predict(X4)
rr4=regr4.score(X4,ps)

regr5=RandomForestRegressor(random_state=0)
regr5.fit(X5,mt)
yr5=regr5.predict(X5)
rr5=regr5.score(X5,mt)


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

f1=plt.figure(1,figsize=(8,10))
f1.set_dpi(250)

ax1=plt.subplot(5,2,1)
ax1.scatter(mt,yr5,c='darkred',s=5,alpha=0.25)
[slp,intc,rv,_,_]=linregress(mt,yr5)
x=np.linspace(np.min(mt),np.max(mt))
plt.plot(x,slp*x+intc,c='k',linestyle=':',label=r'$R^{2}$ = '+str(np.round(rv**2,3)))
plt.legend(loc='lower right')
ax1.set_xlabel('SURFEX-TRIP MAT [C]')
ax1.set_ylabel('RFR Predicted MAT [C]')
ax1.text(0.01, 0.99, 'A',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax1.transAxes,
        fontsize=12,fontweight='extra bold')

ax2=plt.subplot(5,2,2)
ax2.barh(F5,regr5.feature_importances_,color='darkred')
ax2.set_xlabel('Importances from Random Forest Regression')
ax2.text(0.95, 0.99, 'B',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax2.transAxes,
        fontsize=12,fontweight='extra bold')

ax3=plt.subplot(5,2,3)
ax3.scatter(mp,yr1,c='cornflowerblue',s=5,alpha=0.25)
[slp,intc,rv,_,_]=linregress(mp,yr1)
x=np.linspace(np.min(mp),np.max(mp))
plt.plot(x,slp*x+intc,c='k',linestyle=':',label=r'$R^{2}$ = '+str(np.round(rv**2,3)))
plt.legend(loc='lower right')
ax3.set_xlabel(r'WaterGAP3  $\bar{P}$ [mm/day]')
ax3.set_ylabel(r'RFR Predicted $\bar{P}$ [mm/day]')
ax3.text(0.01, 0.99, 'C',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax3.transAxes,
        fontsize=12,fontweight='extra bold')

ax4=plt.subplot(5,2,4)
ax4.barh(F1,regr1.feature_importances_,color='cornflowerblue')
ax4.set_xlabel('Importances from Random Forest Regression')
ax4.text(0.95, 0.15, 'D',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax4.transAxes,
        fontsize=12,fontweight='extra bold')

ax5=plt.subplot(5,2,5)
ax5.scatter(mr,yr2,c='darkblue',s=5,alpha=0.25)
[slp,intc,rv,_,_]=linregress(mr,yr2)
x=np.linspace(np.min(mr),np.max(mr))
plt.plot(x,slp*x+intc,c='k',linestyle=':',label=r'$R^{2}$ = '+str(np.round(rv**2,3)))
plt.legend(loc='lower right')
ax5.set_xlabel(r'WaterGap3 $\bar{R}$ [mm/day]')
ax5.set_ylabel(r'RFR Predicted $\bar{R}$ [mm/day]')
ax5.text(0.01, 0.99, 'E',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax5.transAxes,
        fontsize=12,fontweight='extra bold')

ax6=plt.subplot(5,2,6)
ax6.barh(F2,regr2.feature_importances_,color='darkblue')
ax6.set_xlabel('Importances from Random Forest Regression')
ax6.text(0.95, 0.99, 'F',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax6.transAxes,
        fontsize=12,fontweight='extra bold')

ax7=plt.subplot(5,2,7)
ax7.scatter(cr,yr3,c='forestgreen',s=5,alpha=0.25)
[slp,intc,rv,_,_]=linregress(cr,yr3)
x=np.linspace(np.min(cr),np.max(cr))
plt.plot(x,slp*x+intc,c='k',linestyle=':',label=r'$R^{2}$ = '+str(np.round(rv**2,3)))
plt.legend(loc='lower right')
ax7.set_xlabel(r'WaterGAP3 $c_{R}$')
ax7.set_ylabel(r'RFR Predicted $c_{R}$')
ax7.text(0.01, 0.99, 'G',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax7.transAxes,
        fontsize=12,fontweight='extra bold')

ax8=plt.subplot(5,2,8)
ax8.barh(F3,regr3.feature_importances_,color='forestgreen')
ax8.set_xlabel('Importances from Random Forest Regression')
ax8.text(0.95, 0.99, 'H',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax8.transAxes,
        fontsize=12,fontweight='extra bold')

ax9=plt.subplot(5,2,9)
ax9.scatter(ps,yr4,c='darkcyan',s=5,alpha=0.25)
[slp,intc,rv,_,_]=linregress(ps,yr4)
x=np.linspace(np.min(ps),np.max(ps))
plt.plot(x,slp*x+intc,c='k',linestyle=':',label=r'$R^{2}$ = '+str(np.round(rv**2,3)))
plt.legend(loc='lower right')
ax9.set_xlabel('WaterGAP3 SF')
ax9.set_ylabel('RFR Predicted SF')
ax9.text(0.01, 0.99, 'I',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax9.transAxes,
        fontsize=12,fontweight='extra bold')

ax10=plt.subplot(5,2,10)
ax10.barh(F4,regr4.feature_importances_,color='darkcyan')
ax10.set_xlabel('Importances from Random Forest Regression')
ax10.text(0.95, 0.15, 'J',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax10.transAxes,
        fontsize=12,fontweight='extra bold')

plt.tight_layout()
plt.rcdefaults()
f1.savefig('figure_x1.png',dpi='figure')      