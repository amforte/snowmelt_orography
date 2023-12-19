#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 08:06:44 2023

@author: aforte
"""

import sys
sys.path.insert(0,'/Users/aforte/Documents/GitHub/snowmelt-tectonics')
import stimpy as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Greater Caucasus
gc_col='black'

# Alps
alps_col='royalblue'

# British Columbia
bc_col='orange'

# master_location='/Volumes/Choruh/Data/snowmelt_project/'
master_data_location='/Users/aforte/Documents/Python/snowmelt/'

master_location=master_data_location+'model_outputs_v2/'

u_vec=np.array([250,500,1000,2000,4000,8000])

gc_mlu=['gc025u','gc05u','gc1u','gc2u','gc4u','gc8u']
gc_plu=['ts_','ts_','ts_','ts_','ts_','ts_']
gc_dlu=['gcu','gcu','gcu','gcu','gcu','gcu']

gc_mll=['gc025l','gc05l','gc1l','gc2l','gc4l','gc8l']
gc_pll=['ts_','ts_','ts_','ts_','ts_','ts_']
gc_dll=['gcl','gcl','gcl','gcl','gcl','gcl']

gc_mlu10=['gc025u_10l','gc05u_10l','gc1u_10l','gc2u_10l','gc4u_10l','gc8u_10l']
gc_plu10=['ts_','ts_','ts_','ts_','ts_','ts_']
gc_dlu10=['gcu_10l','gcu_10l','gcu_10l','gcu_10l','gcu_10l','gcu_10l']

gc_mll10=['gc025l_10l','gc05l_10l','gc1l_10l','gc2l_10l','gc4l_10l','gc8l_10l']
gc_pll10=['ts_','ts_','ts_','ts_','ts_','ts_']
gc_dll10=['gcl_10l','gcl_10l','gcl_10l','gcl_10l','gcl_10l','gcl_10l']

gc_mluA=['gc025u_Area','gc05u_Area','gc1u_Area','gc2u_Area','gc4u_Area','gc8u_Area']
gc_pluA=['ts_','ts_','ts_','ts_','ts_','ts_']
gc_dluA=['gcu_A','gcu_A','gcu_A','gcu_A','gcu_A','gcu_A']

gc_mllA=['gc025l_Area','gc05l_Area','gc1l_Area','gc2l_Area','gc4l_Area','gc8l_Area']
gc_pllA=['ts_','ts_','ts_','ts_','ts_','ts_']
gc_dllA=['gcl_A','gcl_A','gcl_A','gcl_A','gcl_A','gcl_A']

gc_mluR=['gc025u_RO','gc05u_RO','gc1u_RO','gc2u_RO','gc4u_RO','gc8u_RO']
gc_pluR=['ts_','ts_','ts_','ts_','ts_','ts_']
gc_dluR=['gcu_R','gcu_R','gcu_R','gcu_R','gcu_R','gcu_R']

gc_mllR=['gc025l_RO','gc05l_RO','gc1l_RO','gc2l_RO','gc4l_RO','gc8l_RO']
gc_pllR=['ts_','ts_','ts_','ts_','ts_','ts_']
gc_dllR=['gcl_R','gcl_R','gcl_R','gcl_R','gcl_R','gcl_R']


gc_mcu=st.ModelComparison(master_location,gc_mlu,gc_plu,gc_dlu)
stabil1=gc_mcu.stability_values()

gc_mcl=st.ModelComparison(master_location,gc_mll,gc_pll,gc_dll)
stabil2=gc_mcl.stability_values()

gc_mcu10=st.ModelComparison(master_location,gc_mlu10,gc_plu10,gc_dlu10)
stabil3=gc_mcu10.stability_values()

gc_mcl10=st.ModelComparison(master_location,gc_mll10,gc_pll10,gc_dll10)
stabil4=gc_mcl10.stability_values()

gc_mcuA=st.ModelComparison(master_location,gc_mluA,gc_pluA,gc_dluA)
stabil9=gc_mcuA.stability_values()

gc_mclA=st.ModelComparison(master_location,gc_mllA,gc_pllA,gc_dllA)
stabil10=gc_mclA.stability_values()

gc_mcuR=st.ModelComparison(master_location,gc_mluR,gc_pluR,gc_dluR)
stabil11=gc_mcuR.stability_values()

gc_mclR=st.ModelComparison(master_location,gc_mllR,gc_pllR,gc_dllR)
stabil12=gc_mclR.stability_values()


a_mlu=['a025u','a05u','a1u','a2u','a4u','a8u']
a_plu=['ts_','ts_','ts_','ts_','ts_','ts_']
a_dlu=['au','au','au','au','au','au']

a_mll=['a025l','a05l','a1l','a2l','a4l','a8l']
a_pll=['ts_','ts_','ts_','ts_','ts_','ts_']
a_dll=['al','al','al','al','al','al']

    
a_mcu=st.ModelComparison(master_location,a_mlu,a_plu,a_dlu)
stabil5=a_mcu.stability_values()


a_mcl=st.ModelComparison(master_location,a_mll,a_pll,a_dll)
stabil6=a_mcl.stability_values()

bc_mlu=['bc025u','bc05u','bc1u','bc2u','bc4u','bc8u']
bc_plu=['ts_','ts_','ts_','ts_','ts_','ts_']
bc_dlu=['bcu','bcu','bcu','bcu','bcu','bcu']

bc_mll=['bc025l','bc05l','bc1l','bc2l','bc4l','bc8l']
bc_pll=['ts_','ts_','ts_','ts_','ts_','ts_']
bc_dll=['bcl','bcl','bcl','bcl','bcl','bcl']

    
bc_mcu=st.ModelComparison(master_location,bc_mlu,bc_plu,bc_dlu)
stabil7=bc_mcu.stability_values()

bc_mcl=st.ModelComparison(master_location,bc_mll,bc_pll,bc_dll)
stabil8=bc_mcl.stability_values()

model=gc_mlu+gc_mll+gc_mlu10+gc_mll10+a_mlu+a_mll+bc_mlu+bc_mll+gc_mluA+gc_mllA+gc_mluR+gc_mllR
group=gc_dlu+gc_dll+gc_dlu10+gc_dll10+a_dlu+a_dll+bc_dlu+bc_dll+gc_dluA+gc_dllA+gc_dluR+gc_dllR
stabil=np.concatenate((stabil1,stabil2,stabil3,stabil4,stabil5,stabil6,stabil7,stabil8,stabil9,stabil10,stabil11,stabil12),axis=0)
U=np.tile(u_vec,12)


group_list=['gcu','gcl','gcu_10l','gcl_10l','au','al','bcu','bcl','gcu_A','gcl_A','gcu_R','gcl_R']
mar_list=['o','s','o','s','o','s','o','s','^','v','<','>']
len_list=[50,50,10,10,50,50,50,50,50,50]
col_list=[gc_col,gc_col,gc_col,gc_col,alps_col,alps_col,bc_col,bc_col,'gray','gray','gray','gray']
label=['GC Unlinked','GC Linked','GC Unlinked 10 km','GC Linked 10 km',
       'Alps Unlinked','Alps Linked','BC Unlinked','BC Linked','GC Unlinked Area','GC Linked Area','GC Unlinked Rain','GC Linked Rain']

f1=plt.figure(1,figsize=(5,5))
for i in range(len(group_list)):
    idx=np.array(group)==group_list[i]
    plt.scatter(U[idx],stabil[idx],color=col_list[i],marker=mar_list[i])
    
plt.xscale('log')
plt.xlabel('Uplift Rate [m/Myr]')
plt.ylabel('Max CFL Value During Run')
plt.tight_layout()

f1.savefig('REVIEW_stability.pdf')


model_list=['bc8u','bc8u_HighDef']
prefix_list=['ts_','ts_']
descript_list=['BC 8 mm/yr - Unlinked','BC 8 mm/yr - Unlinked - HighDef']

mc=st.ModelComparison(master_location,model_list,prefix_list,descript_list)
# f1=mc.comp_profile_evol(25,3.25,3.5,35,[4.5e6,2.5e6])
f2=mc.comp_profile_evol2(50,3.25,20,15,[5e5,5e5])
f2.savefig('REVIEW_highdefcomp.pdf')

f3=mc.comp_excd_prob2(0.75,20,[5e5,5e5],['k','gray'])
f3.savefig('REVIEW_highdefcomp2.pdf')
