#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 24 14:34:53 2023

@author: aforte
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

def survive(ts):
    ts_sort=np.sort(ts)
    tsn=len(ts_sort)
    tsrank=np.arange(1,tsn+1,1)
    ts_freq_excd=(tsn+1-tsrank)/tsn
    return ts_sort,ts_freq_excd

# master_location='/Volumes/Choruh/Data/snowmelt_project/'
master_location='/Users/aforte/Documents/python/snowmelt/'
YR=np.arange(1980,2000,1)

P5=[]
R5=[]
P10=[]
R10=[]
P15=[]
R15=[]
P20=[]
R20=[]
P25=[]
R25=[]
P30=[]
R30=[]
P35=[]
R35=[]

for i in range(len(YR)):
    pdf=pd.read_csv(master_location+'wrr2_events/Precip_Threshold_5_'+str(YR[i])+'.csv')
    rdf=pd.read_csv(master_location+'wrr2_events/Runoff_Threshold_5_'+str(YR[i])+'.csv')
    P5.append(pdf)
    R5.append(rdf)
    pdf=pd.read_csv(master_location+'wrr2_events/Precip_Threshold_10_'+str(YR[i])+'.csv')
    rdf=pd.read_csv(master_location+'wrr2_events/Runoff_Threshold_10_'+str(YR[i])+'.csv')
    P10.append(pdf)
    R10.append(rdf)
    pdf=pd.read_csv(master_location+'wrr2_events/Precip_Threshold_15_'+str(YR[i])+'.csv')
    rdf=pd.read_csv(master_location+'wrr2_events/Runoff_Threshold_15_'+str(YR[i])+'.csv')
    P15.append(pdf)
    R15.append(rdf)
    pdf=pd.read_csv(master_location+'wrr2_events/Precip_Threshold_20_'+str(YR[i])+'.csv')
    rdf=pd.read_csv(master_location+'wrr2_events/Runoff_Threshold_20_'+str(YR[i])+'.csv')
    P20.append(pdf)
    R20.append(rdf)
    pdf=pd.read_csv(master_location+'wrr2_events/Precip_Threshold_25_'+str(YR[i])+'.csv')
    rdf=pd.read_csv(master_location+'wrr2_events/Runoff_Threshold_25_'+str(YR[i])+'.csv')
    P25.append(pdf)
    R25.append(rdf)
    pdf=pd.read_csv(master_location+'wrr2_events/Precip_Threshold_30_'+str(YR[i])+'.csv')
    rdf=pd.read_csv(master_location+'wrr2_events/Runoff_Threshold_30_'+str(YR[i])+'.csv')
    P30.append(pdf)
    R30.append(rdf)
    pdf=pd.read_csv(master_location+'wrr2_events/Precip_Threshold_35_'+str(YR[i])+'.csv')
    rdf=pd.read_csv(master_location+'wrr2_events/Runoff_Threshold_35_'+str(YR[i])+'.csv')
    P35.append(pdf)
    R35.append(rdf)
    
PDF5=pd.concat(P5,ignore_index=True)
RDF5=pd.concat(R5,ignore_index=True)
PDF10=pd.concat(P10,ignore_index=True)
RDF10=pd.concat(R10,ignore_index=True)
PDF15=pd.concat(P15,ignore_index=True)
RDF15=pd.concat(R15,ignore_index=True)
PDF20=pd.concat(P20,ignore_index=True)
RDF20=pd.concat(R20,ignore_index=True)
PDF25=pd.concat(P25,ignore_index=True)
RDF25=pd.concat(R25,ignore_index=True)
PDF30=pd.concat(P30,ignore_index=True)
RDF30=pd.concat(R30,ignore_index=True)
PDF35=pd.concat(P35,ignore_index=True)
RDF35=pd.concat(R35,ignore_index=True)

# Filter
pidx=(np.isnan(PDF5['area_pixels'])) & (PDF5['area_pixels']==1) 
ridx=(np.isnan(RDF5['area_pixels'])) & (RDF5['qs_sum']<0) & (RDF5['area_pixels']==1)
PDF5=PDF5.drop(index=PDF5.index[pidx])
RDF5=RDF5.drop(index=RDF5.index[ridx])
pidx=np.isnan(PDF10['area_pixels'])
ridx=(np.isnan(RDF10['area_pixels'])) & (RDF10['qs_sum']<0)
PDF10=PDF10.drop(index=PDF10.index[pidx])
RDF10=RDF10.drop(index=RDF10.index[ridx])
pidx=np.isnan(PDF15['area_pixels'])
ridx=(np.isnan(RDF15['area_pixels'])) & (RDF15['qs_sum']<0)
PDF15=PDF15.drop(index=PDF15.index[pidx])
RDF15=RDF15.drop(index=RDF15.index[ridx])
pidx=np.isnan(PDF20['area_pixels'])
ridx=(np.isnan(RDF20['area_pixels'])) & (RDF20['qs_sum']<0)
PDF20=PDF20.drop(index=PDF20.index[pidx])
RDF20=RDF20.drop(index=RDF20.index[ridx])
pidx=np.isnan(PDF25['area_pixels'])
ridx=(np.isnan(RDF25['area_pixels'])) & (RDF25['qs_sum']<0)
PDF25=PDF25.drop(index=PDF25.index[pidx])
RDF25=RDF25.drop(index=RDF25.index[ridx])
pidx=np.isnan(PDF30['area_pixels'])
ridx=(np.isnan(RDF30['area_pixels'])) & (RDF30['qs_sum']<0)
PDF30=PDF30.drop(index=PDF30.index[pidx])
RDF30=RDF30.drop(index=RDF30.index[ridx])
pidx=np.isnan(PDF35['area_pixels'])
ridx=(np.isnan(RDF35['area_pixels'])) & (RDF35['qs_sum']<0)
PDF35=PDF35.drop(index=PDF35.index[pidx])
RDF35=RDF35.drop(index=RDF35.index[ridx])

# Calculate metrics
r5_sf = RDF5['qsm_sum']/(RDF5['qs_sum'] + RDF5['qsb_sum'] + RDF5['qsm_sum'])
ix=np.argsort(RDF5['elliptical_area_km'])
r5_sf_sort=r5_sf[ix]
r10_sf = RDF10['qsm_sum']/(RDF10['qs_sum'] + RDF10['qsb_sum'] + RDF10['qsm_sum'])
ix=np.argsort(RDF10['elliptical_area_km'])
r10_sf_sort=r10_sf[ix]
r15_sf = RDF15['qsm_sum']/(RDF15['qs_sum'] + RDF15['qsb_sum'] + RDF15['qsm_sum'])
ix=np.argsort(RDF15['elliptical_area_km'])
r15_sf_sort=r15_sf[ix]
r20_sf = RDF20['qsm_sum']/(RDF20['qs_sum'] + RDF20['qsb_sum'] + RDF20['qsm_sum'])
ix=np.argsort(RDF20['elliptical_area_km'])
r20_sf_sort=r20_sf[ix]
r25_sf = RDF25['qsm_sum']/(RDF25['qs_sum'] + RDF25['qsb_sum'] + RDF25['qsm_sum'])
ix=np.argsort(RDF25['elliptical_area_km'])
r25_sf_sort=r25_sf[ix]
r30_sf = RDF30['qsm_sum']/(RDF30['qs_sum'] + RDF30['qsb_sum'] + RDF30['qsm_sum'])
ix=np.argsort(RDF30['elliptical_area_km'])
r30_sf_sort=r30_sf[ix]
r35_sf = RDF35['qsm_sum']/(RDF35['qs_sum'] + RDF35['qsb_sum'] + RDF35['qsm_sum'])
ix=np.argsort(RDF35['elliptical_area_km'])
r35_sf_sort=r35_sf[ix]

p_sort5,p_freq5=survive(PDF5['elliptical_area_km'])
r_sort5,r_freq5=survive(RDF5['elliptical_area_km'])
p_sort10,p_freq10=survive(PDF10['elliptical_area_km'])
r_sort10,r_freq10=survive(RDF10['elliptical_area_km'])
p_sort15,p_freq15=survive(PDF15['elliptical_area_km'])
r_sort15,r_freq15=survive(RDF15['elliptical_area_km'])
p_sort20,p_freq20=survive(PDF20['elliptical_area_km'])
r_sort20,r_freq20=survive(RDF20['elliptical_area_km'])
p_sort25,p_freq25=survive(PDF25['elliptical_area_km'])
r_sort25,r_freq25=survive(RDF25['elliptical_area_km'])
p_sort30,p_freq30=survive(PDF30['elliptical_area_km'])
r_sort30,r_freq30=survive(RDF30['elliptical_area_km'])
p_sort35,p_freq35=survive(PDF35['elliptical_area_km'])
r_sort35,r_freq35=survive(RDF35['elliptical_area_km'])


f_bin=np.logspace(-7,0,50)

f1=plt.figure(1,figsize=(8,18))
f1.set_dpi(250)
gs = GridSpec(7,4,figure=f1)

ax1=f1.add_subplot(gs[0,0:3])
# plt.scatter(p_sort5,p_freq5,10,'lightblue',label='Precipitation',marker='o')
# plt.scatter(r_sort5[r5_sf_sort>0.35],r_freq5[r5_sf_sort>0.35],10,'cornflowerblue',label='Runoff - Snow Dominated',marker='o')
# plt.scatter(r_sort5[r5_sf_sort<=0.35],r_freq5[r5_sf_sort<=0.35],10,'darkblue',label='Runoff - Rain Dominated',marker='o')
plt.plot(p_sort5,p_freq5,color='lightblue',label='Precipitation')
plt.plot(r_sort5[r5_sf_sort>0.35],r_freq5[r5_sf_sort>0.35],color='cornflowerblue',label='Runoff - Snow Dominated')
plt.plot(r_sort5[r5_sf_sort<=0.35],r_freq5[r5_sf_sort<=0.35],color='darkblue',label='Runoff - Precipitation Dominated')
plt.xscale('log')
plt.yscale('log')
# plt.xlabel(r'Elliptical Event Area [km$^{2}$]')
plt.ylabel('Exccedance Frequency')
plt.legend(loc='best')
plt.title('5 mm/day Threshold')
plt.xlim((70,10**7.1))
plt.ylim((10**-7,2))

ax2=f1.add_subplot(gs[0,3])
fix=np.digitize(r_freq5,f_bin)
for i in range(len(f_bin)-1):
    xvals=r5_sf_sort[fix==i+1]
    yvals=r_freq5[fix==i+1]
    if len(xvals)>0:
        frac=len(xvals[xvals>0.35])/len(xvals)
        plt.scatter(frac*100,np.mean(yvals),c='k',s=5)
plt.ylim((10**-7,2))
plt.yscale('log')
# plt.xlabel('% Snowmelt Dominated')

ax3=f1.add_subplot(gs[1,0:3])
# plt.scatter(p_sort10,p_freq10,10,'lightblue',label='Precipitation',marker='o')
# plt.scatter(r_sort10[r10_sf_sort>0.35],r_freq10[r10_sf_sort>0.35],10,'cornflowerblue',label='Runoff - Snow Dominated',marker='o')
# plt.scatter(r_sort10[r10_sf_sort<=0.35],r_freq10[r10_sf_sort<=0.35],10,'darkblue',label='Runoff - Rain Dominated',marker='o')
plt.plot(p_sort10,p_freq10,color='lightblue',label='Precipitation')
plt.plot(r_sort10[r10_sf_sort>0.35],r_freq10[r10_sf_sort>0.35],color='cornflowerblue',label='Runoff - Snow Dominated')
plt.plot(r_sort10[r10_sf_sort<=0.35],r_freq10[r10_sf_sort<=0.35],color='darkblue',label='Runoff - Rain Dominated')
plt.xscale('log')
plt.yscale('log')
# plt.xlabel(r'Elliptical Event Area [km$^{2}$]')
plt.ylabel('Exccedance Frequency')
plt.legend(loc='best')
plt.title('10 mm/day Threshold')
plt.xlim((10**1.5,10**7.1))
plt.ylim((10**-7,2))

ax4=f1.add_subplot(gs[1,3])
fix=np.digitize(r_freq10,f_bin)
for i in range(len(f_bin)-1):
    xvals=r10_sf_sort[fix==i+1]
    yvals=r_freq10[fix==i+1]
    if len(xvals)>0:
        frac=len(xvals[xvals>0.35])/len(xvals)
        plt.scatter(frac*100,np.mean(yvals),c='k',s=5)
plt.ylim((10**-7,2))
plt.yscale('log')
# plt.xlabel('% Snowmelt Dominated')


ax5=f1.add_subplot(gs[2,0:3])
# plt.scatter(p_sort15,p_freq15,10,'lightblue',label='Precipitation',marker='o')
# plt.scatter(r_sort15[r15_sf_sort>0.35],r_freq15[r15_sf_sort>0.35],10,'cornflowerblue',label='Runoff - Snow Dominated',marker='o')
# plt.scatter(r_sort15[r15_sf_sort<=0.35],r_freq15[r15_sf_sort<=0.35],10,'darkblue',label='Runoff - Rain Dominated',marker='o')
plt.plot(p_sort15,p_freq15,color='lightblue',label='Precipitation')
plt.plot(r_sort15[r15_sf_sort>0.35],r_freq15[r15_sf_sort>0.35],color='cornflowerblue',label='Runoff - Snow Dominated')
plt.plot(r_sort15[r15_sf_sort<=0.35],r_freq15[r15_sf_sort<=0.35],color='darkblue',label='Runoff - Rain Dominated')
plt.xscale('log')
plt.yscale('log')
# plt.xlabel(r'Elliptical Event Area [km$^{2}$]')
plt.ylabel('Exccedance Frequency')
plt.legend(loc='best')
plt.title('15 mm/day Threshold')
plt.xlim((10**1.5,10**7.1))
plt.ylim((10**-7,2))

ax6=f1.add_subplot(gs[2,3])
fix=np.digitize(r_freq15,f_bin)
for i in range(len(f_bin)-1):
    xvals=r15_sf_sort[fix==i+1]
    yvals=r_freq15[fix==i+1]
    if len(xvals)>0:
        frac=len(xvals[xvals>0.35])/len(xvals)
        plt.scatter(frac*100,np.mean(yvals),c='k',s=5)
plt.ylim((10**-7,2))
plt.yscale('log')
# plt.xlabel('% Snowmelt Dominated')


ax7=f1.add_subplot(gs[3,0:3])
# plt.scatter(p_sort20,p_freq20,10,'lightblue',label='Precipitation',marker='o')
# plt.scatter(r_sort20[r20_sf_sort>0.35],r_freq20[r20_sf_sort>0.35],10,'cornflowerblue',label='Runoff - Snow Dominated',marker='o')
# plt.scatter(r_sort20[r20_sf_sort<=0.35],r_freq20[r20_sf_sort<=0.35],10,'darkblue',label='Runoff - Rain Dominated',marker='o')
plt.plot(p_sort20,p_freq20,color='lightblue',label='Precipitation')
plt.plot(r_sort20[r20_sf_sort>0.35],r_freq20[r20_sf_sort>0.35],color='cornflowerblue',label='Runoff - Snow Dominated')
plt.plot(r_sort20[r20_sf_sort<=0.35],r_freq20[r20_sf_sort<=0.35],color='darkblue',label='Runoff - Rain Dominated')
plt.xscale('log')
plt.yscale('log')
# plt.xlabel(r'Elliptical Event Area [km$^{2}$]')
plt.ylabel('Exccedance Frequency')
plt.legend(loc='best')
plt.title('20 mm/day Threshold')
plt.xlim((10**1.5,10**7.1))
plt.ylim((10**-7,2))

ax8=f1.add_subplot(gs[3,3])
fix=np.digitize(r_freq20,f_bin)
for i in range(len(f_bin)-1):
    xvals=r20_sf_sort[fix==i+1]
    yvals=r_freq20[fix==i+1]
    if len(xvals)>0:
        frac=len(xvals[xvals>0.35])/len(xvals)
        plt.scatter(frac*100,np.mean(yvals),c='k',s=5)
plt.ylim((10**-7,2))
plt.yscale('log')
# plt.xlabel('% Snowmelt Dominated')


ax9=f1.add_subplot(gs[4,0:3])
# plt.scatter(p_sort25,p_freq25,10,'lightblue',label='Precipitation',marker='o')
# plt.scatter(r_sort25[r25_sf_sort>0.35],r_freq25[r25_sf_sort>0.35],10,'cornflowerblue',label='Runoff - Snow Dominated',marker='o')
# plt.scatter(r_sort25[r25_sf_sort<=0.35],r_freq25[r25_sf_sort<=0.35],10,'darkblue',label='Runoff - Rain Dominated',marker='o')
plt.plot(p_sort25,p_freq25,color='lightblue',label='Precipitation')
plt.plot(r_sort25[r25_sf_sort>0.35],r_freq25[r25_sf_sort>0.35],color='cornflowerblue',label='Runoff - Snow Dominated')
plt.plot(r_sort25[r25_sf_sort<=0.35],r_freq25[r25_sf_sort<=0.35],color='darkblue',label='Runoff - Rain Dominated')
plt.xscale('log')
plt.yscale('log')
# plt.xlabel(r'Elliptical Event Area [km$^{2}$]')
plt.ylabel('Exccedance Frequency')
plt.legend(loc='best')
plt.title('25 mm/day Threshold')
plt.xlim((10**1.5,10**7.1))
plt.ylim((10**-7,2))

ax10=f1.add_subplot(gs[4,3])
fix=np.digitize(r_freq25,f_bin)
for i in range(len(f_bin)-1):
    xvals=r25_sf_sort[fix==i+1]
    yvals=r_freq25[fix==i+1]
    if len(xvals)>0:
        frac=len(xvals[xvals>0.35])/len(xvals)
        plt.scatter(frac*100,np.mean(yvals),c='k',s=5)
plt.ylim((10**-7,2))
plt.yscale('log')
# plt.xlabel('% Snowmelt Dominated')

ax11=f1.add_subplot(gs[5,0:3])
# plt.scatter(p_sort30,p_freq30,10,'lightblue',label='Precipitation',marker='o')
# plt.scatter(r_sort30[r30_sf_sort>0.35],r_freq30[r30_sf_sort>0.35],10,'cornflowerblue',label='Runoff - Snow Dominated',marker='o')
# plt.scatter(r_sort30[r30_sf_sort<=0.35],r_freq30[r30_sf_sort<=0.35],10,'darkblue',label='Runoff - Rain Dominated',marker='o')
plt.plot(p_sort30,p_freq30,color='lightblue',label='Precipitation')
plt.plot(r_sort30[r30_sf_sort>0.35],r_freq30[r30_sf_sort>0.35],color='cornflowerblue',label='Runoff - Snow Dominated')
plt.plot(r_sort30[r30_sf_sort<=0.35],r_freq30[r30_sf_sort<=0.35],color='darkblue',label='Runoff - Rain Dominated')
plt.xscale('log')
plt.yscale('log')
# plt.xlabel(r'Elliptical Event Area [km$^{2}$]')
plt.ylabel('Exccedance Frequency')
plt.legend(loc='best')
plt.title('30 mm/day Threshold')
plt.xlim((10**1.5,10**7.1))
plt.ylim((10**-7,2))

ax12=f1.add_subplot(gs[5,3])
fix=np.digitize(r_freq30,f_bin)
for i in range(len(f_bin)-1):
    xvals=r30_sf_sort[fix==i+1]
    yvals=r_freq30[fix==i+1]
    if len(xvals)>0:
        frac=len(xvals[xvals>0.35])/len(xvals)
        plt.scatter(frac*100,np.mean(yvals),c='k',s=5)
plt.ylim((10**-7,2))
plt.yscale('log')
# plt.xlabel('% Snowmelt Dominated')

ax13=f1.add_subplot(gs[6,0:3])
# plt.scatter(p_sort35,p_freq30,10,'lightblue',label='Precipitation',marker='o')
# plt.scatter(r_sort35[r35_sf_sort>0.35],r_freq35[r35_sf_sort>0.35],10,'cornflowerblue',label='Runoff - Snow Dominated',marker='o')
# plt.scatter(r_sort35[r35_sf_sort<=0.35],r_freq35[r35_sf_sort<=0.35],10,'darkblue',label='Runoff - Rain Dominated',marker='o')
plt.plot(p_sort35,p_freq35,color='lightblue',label='Precipitation')
plt.plot(r_sort35[r35_sf_sort>0.35],r_freq35[r35_sf_sort>0.35],color='cornflowerblue',label='Runoff - Snow Dominated')
plt.plot(r_sort35[r35_sf_sort<=0.35],r_freq35[r35_sf_sort<=0.35],color='darkblue',label='Runoff - Rain Dominated')
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'Elliptical Event Area [km$^{2}$]')
plt.ylabel('Exccedance Frequency')
plt.legend(loc='best')
plt.title('35 mm/day Threshold')
plt.xlim((10**1.5,10**7.1))
plt.ylim((10**-7,2))

ax14=f1.add_subplot(gs[6,3])
fix=np.digitize(r_freq35,f_bin)
for i in range(len(f_bin)-1):
    xvals=r35_sf_sort[fix==i+1]
    yvals=r_freq35[fix==i+1]
    if len(xvals)>0:
        frac=len(xvals[xvals>0.35])/len(xvals)
        plt.scatter(frac*100,np.mean(yvals),c='k',s=5)
plt.ylim((10**-7,2))
plt.yscale('log')
plt.xlabel('% Snowmelt Dominated')

plt.tight_layout()

###############################

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

f2=plt.figure(2,figsize=(7,9))
f2.set_dpi(250)
gs = GridSpec(4,4,figure=f1)

ax1=f2.add_subplot(gs[0,0:3])
# plt.scatter(p_sort5,p_freq5,10,'lightblue',label='Precipitation',marker='o')
# plt.scatter(r_sort5[r5_sf_sort>0.35],r_freq5[r5_sf_sort>0.35],10,'cornflowerblue',label='Runoff - Snow Dominated',marker='o')
# plt.scatter(r_sort5[r5_sf_sort<=0.35],r_freq5[r5_sf_sort<=0.35],10,'darkblue',label='Runoff - Rain Dominated',marker='o')
plt.plot(p_sort5,p_freq5,color='lightblue',label='Precipitation')
plt.plot(r_sort5[r5_sf_sort>0.35],r_freq5[r5_sf_sort>0.35],color='cornflowerblue',label='Runoff - Snow Dominated')
plt.plot(r_sort5[r5_sf_sort<=0.35],r_freq5[r5_sf_sort<=0.35],color='darkblue',label='Runoff - Precipitation Dominated')
plt.xscale('log')
plt.yscale('log')
# plt.xlabel(r'Elliptical Event Area [km$^{2}$]')
plt.ylabel('Exccedance Frequency')
plt.legend(loc='best')
plt.title('5 mm/day Threshold')
plt.xlim((70,10**7.1))
plt.ylim((10**-7,2))
ax1.text(0.95, 0.99, 'A',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax1.transAxes,
        fontsize=12,fontweight='extra bold')

ax2=f2.add_subplot(gs[0,3])
fix=np.digitize(r_freq5,f_bin)
for i in range(len(f_bin)-1):
    xvals=r5_sf_sort[fix==i+1]
    yvals=r_freq5[fix==i+1]
    if len(xvals)>0:
        frac=len(xvals[xvals>0.35])/len(xvals)
        plt.scatter(frac*100,np.mean(yvals),c='k',s=5)
plt.ylim((10**-7,2))
plt.yscale('log')
# plt.xlabel('% Snowmelt Dominated')
ax2.text(0.85, 0.99, 'B',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax2.transAxes,
        fontsize=12,fontweight='extra bold')


ax5=f2.add_subplot(gs[1,0:3])
# plt.scatter(p_sort15,p_freq15,10,'lightblue',label='Precipitation',marker='o')
# plt.scatter(r_sort15[r15_sf_sort>0.35],r_freq15[r15_sf_sort>0.35],10,'cornflowerblue',label='Runoff - Snow Dominated',marker='o')
# plt.scatter(r_sort15[r15_sf_sort<=0.35],r_freq15[r15_sf_sort<=0.35],10,'darkblue',label='Runoff - Rain Dominated',marker='o')
plt.plot(p_sort15,p_freq15,color='lightblue',label='Precipitation')
plt.plot(r_sort15[r15_sf_sort>0.35],r_freq15[r15_sf_sort>0.35],color='cornflowerblue',label='Runoff - Snow Dominated')
plt.plot(r_sort15[r15_sf_sort<=0.35],r_freq15[r15_sf_sort<=0.35],color='darkblue',label='Runoff - Rain Dominated')
plt.xscale('log')
plt.yscale('log')
# plt.xlabel(r'Elliptical Event Area [km$^{2}$]')
plt.ylabel('Exccedance Frequency')
plt.legend(loc='best')
plt.title('15 mm/day Threshold')
plt.xlim((10**1.5,10**7.1))
plt.ylim((10**-7,2))
ax5.text(0.95, 0.99, 'C',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax5.transAxes,
        fontsize=12,fontweight='extra bold')

ax6=f2.add_subplot(gs[1,3])
fix=np.digitize(r_freq15,f_bin)
for i in range(len(f_bin)-1):
    xvals=r15_sf_sort[fix==i+1]
    yvals=r_freq15[fix==i+1]
    if len(xvals)>0:
        frac=len(xvals[xvals>0.35])/len(xvals)
        plt.scatter(frac*100,np.mean(yvals),c='k',s=5)
plt.ylim((10**-7,2))
plt.yscale('log')
# plt.xlabel('% Snowmelt Dominated')
ax6.text(0.85, 0.99, 'D',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax6.transAxes,
        fontsize=12,fontweight='extra bold')


ax9=f2.add_subplot(gs[2,0:3])
# plt.scatter(p_sort25,p_freq25,10,'lightblue',label='Precipitation',marker='o')
# plt.scatter(r_sort25[r25_sf_sort>0.35],r_freq25[r25_sf_sort>0.35],10,'cornflowerblue',label='Runoff - Snow Dominated',marker='o')
# plt.scatter(r_sort25[r25_sf_sort<=0.35],r_freq25[r25_sf_sort<=0.35],10,'darkblue',label='Runoff - Rain Dominated',marker='o')
plt.plot(p_sort25,p_freq25,color='lightblue',label='Precipitation')
plt.plot(r_sort25[r25_sf_sort>0.35],r_freq25[r25_sf_sort>0.35],color='cornflowerblue',label='Runoff - Snow Dominated')
plt.plot(r_sort25[r25_sf_sort<=0.35],r_freq25[r25_sf_sort<=0.35],color='darkblue',label='Runoff - Rain Dominated')
plt.xscale('log')
plt.yscale('log')
# plt.xlabel(r'Elliptical Event Area [km$^{2}$]')
plt.ylabel('Exccedance Frequency')
plt.legend(loc='best')
plt.title('25 mm/day Threshold')
plt.xlim((10**1.5,10**7.1))
plt.ylim((10**-7,2))
ax9.text(0.95, 0.99, 'E',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax9.transAxes,
        fontsize=12,fontweight='extra bold')

ax10=f2.add_subplot(gs[2,3])
fix=np.digitize(r_freq25,f_bin)
for i in range(len(f_bin)-1):
    xvals=r25_sf_sort[fix==i+1]
    yvals=r_freq25[fix==i+1]
    if len(xvals)>0:
        frac=len(xvals[xvals>0.35])/len(xvals)
        plt.scatter(frac*100,np.mean(yvals),c='k',s=5)
plt.ylim((10**-7,2))
plt.yscale('log')
# plt.xlabel('% Snowmelt Dominated')
ax10.text(0.85, 0.99, 'F',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax10.transAxes,
        fontsize=12,fontweight='extra bold')

ax13=f2.add_subplot(gs[3,0:3])
# plt.scatter(p_sort35,p_freq30,10,'lightblue',label='Precipitation',marker='o')
# plt.scatter(r_sort35[r35_sf_sort>0.35],r_freq35[r35_sf_sort>0.35],10,'cornflowerblue',label='Runoff - Snow Dominated',marker='o')
# plt.scatter(r_sort35[r35_sf_sort<=0.35],r_freq35[r35_sf_sort<=0.35],10,'darkblue',label='Runoff - Rain Dominated',marker='o')
plt.plot(p_sort35,p_freq35,color='lightblue',label='Precipitation')
plt.plot(r_sort35[r35_sf_sort>0.35],r_freq35[r35_sf_sort>0.35],color='cornflowerblue',label='Runoff - Snow Dominated')
plt.plot(r_sort35[r35_sf_sort<=0.35],r_freq35[r35_sf_sort<=0.35],color='darkblue',label='Runoff - Rain Dominated')
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'Elliptical Event Area [km$^{2}$]')
plt.ylabel('Exccedance Frequency')
plt.legend(loc='best')
plt.title('35 mm/day Threshold')
plt.xlim((10**1.5,10**7.1))
plt.ylim((10**-7,2))
ax13.text(0.95, 0.99, 'G',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax13.transAxes,
        fontsize=12,fontweight='extra bold')

ax14=f2.add_subplot(gs[3,3])
fix=np.digitize(r_freq35,f_bin)
for i in range(len(f_bin)-1):
    xvals=r35_sf_sort[fix==i+1]
    yvals=r_freq35[fix==i+1]
    if len(xvals)>0:
        frac=len(xvals[xvals>0.35])/len(xvals)
        plt.scatter(frac*100,np.mean(yvals),c='k',s=5)
plt.ylim((10**-7,2))
plt.yscale('log')
plt.xlabel('% Snowmelt Dominated')
ax14.text(0.85, 0.99, 'H',
        horizontalalignment='left',
        verticalalignment='top',
        transform=ax14.transAxes,
        fontsize=12,fontweight='extra bold')

plt.tight_layout()
plt.rcdefaults()

f2.savefig('P1_figure9.pdf',dpi="figure")

