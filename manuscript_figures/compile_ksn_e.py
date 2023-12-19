#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec  3 14:57:31 2023

@author: aforte
"""

import pandas as pd
import numpy as np
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


# Load final ts outupts
df=pd.read_csv('model_final_stats.csv')
dft=pd.read_csv('model_final_theta.csv')

group_list=['gcu','gcl','gcu_10l','gcl_10l','au','al','bcu','bcl','bcu_R','bcl_R','gcu_A','gcl_A','gcu_R','gcl_R']
label_list=['GC Unlinked','GC Linked','GC Unlinked 10 km','GC Linked 10 km',
       'Alps Unlinked','Alps Linked','BC Unlinked','BC Linked','BC Unlinked RainOnly','BC Linked RainOnly','GC Unlinked Area',
       'GC Linked Area','GC Unlinked RainOnly','GC Linked RainOnly']

Klist=[]
Klplist=[]
nlist=[]
nlplist=[]
Clist=[]
Clplist=[]
philist=[]
philplist=[]


for i in range(len(group_list)):
    idx=df['Group']==group_list[i]
    
    ksn=df.loc[idx,'mn_ksn'].to_numpy()
    E=df.loc[idx,'mn_E'].to_numpy()/1e6
    ksnqp=df.loc[idx,'mn_ksnqp'].to_numpy()
    
    ksns=df.loc[idx,'std_ksn'].to_numpy()
    Es=df.loc[idx,'std_E'].to_numpy()
    ksnqps=df.loc[idx,'std_ksnqp'].to_numpy()
    
    theta=dft.loc[idx,'theta_chi'].to_numpy()
    thetap=dft.loc[idx,'theta_chi_p'].to_numpy()

    # Fit
    [K,n]=odr_fit(ksn,E)
    [Klp,nlp]=odr_fit(ksnqp,E)
    [C,phi]=odr_fit(E*1e6,ksn)
    [Clp,philp]=odr_fit(E*1e6,ksnqp)
    
    Klist.append(K)
    Klplist.append(Klp)
    nlist.append(n)
    nlplist.append(nlp)
    Clist.append(C)
    Clplist.append(Clp)
    philist.append(phi)
    philplist.append(philp)
    
df_out=pd.DataFrame(data={'Group':label_list,
                          'C':Clist,
                          'phi':philist,
                          'K':Klist,
                          'n':nlist,
                          'C_lp':Clplist,
                          'phi_lp':philplist,
                          'K_lp':Klplist,
                          'n_lp':nlplist})
df_out.to_csv('ksn_e_fit.csv',index=False)