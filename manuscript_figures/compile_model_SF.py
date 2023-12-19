#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  7 07:43:18 2023

@author: aforte
"""

import sys
sys.path.insert(0,'/Users/aforte/Documents/GitHub/snowmelt_orography')
import stimpy as st
import pandas as pd
import numpy as np



# master_location='/Volumes/Choruh/Data/snowmelt_project/'
master_data_location='/Users/aforte/Documents/Python/snowmelt/'


master_location=master_data_location+'model_outputs_v2/'

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



gc_mcu=st.ModelComparison(master_location,gc_mlu,gc_plu,gc_dlu)
[mn_snp1,std_snp1,min_snp1,max_snp1]=gc_mcu.sf_final_ts()

gc_mcl=st.ModelComparison(master_location,gc_mll,gc_pll,gc_dll)
[mn_snp2,std_snp2,min_snp2,max_snp2]=gc_mcl.sf_final_ts()

gc_mcu10=st.ModelComparison(master_location,gc_mlu10,gc_plu10,gc_dlu10)
[mn_snp3,std_snp3,min_snp3,max_snp3]=gc_mcu10.sf_final_ts()

gc_mcl10=st.ModelComparison(master_location,gc_mll10,gc_pll10,gc_dll10)
[mn_snp4,std_snp4,min_snp4,max_snp4]=gc_mcl10.sf_final_ts()

gc_mcuA=st.ModelComparison(master_location,gc_mluA,gc_pluA,gc_dluA)
[mn_snp9,std_snp9,min_snp9,max_snp9]=gc_mcuA.sf_final_ts()

gc_mclA=st.ModelComparison(master_location,gc_mllA,gc_pllA,gc_dllA)
[mn_snp10,std_snp10,min_snp10,max_snp10]=gc_mclA.sf_final_ts()



a_mlu=['a025u','a05u','a1u','a2u','a4u','a8u']
a_plu=['ts_','ts_','ts_','ts_','ts_','ts_']
a_dlu=['au','au','au','au','au','au']

a_mll=['a025l','a05l','a1l','a2l','a4l','a8l']
a_pll=['ts_','ts_','ts_','ts_','ts_','ts_']
a_dll=['al','al','al','al','al','al']

    
a_mcu=st.ModelComparison(master_location,a_mlu,a_plu,a_dlu)
[mn_snp5,std_snp5,min_snp5,max_snp5]=a_mcu.sf_final_ts()

a_mcl=st.ModelComparison(master_location,a_mll,a_pll,a_dll)
[mn_snp6,std_snp6,min_snp6,max_snp6]=a_mcl.sf_final_ts()

bc_mlu=['bc025u','bc05u','bc1u','bc2u','bc4u','bc8u']
bc_plu=['ts_','ts_','ts_','ts_','ts_','ts_']
bc_dlu=['bcu','bcu','bcu','bcu','bcu','bcu']

bc_mll=['bc025l','bc05l','bc1l','bc2l','bc4l','bc8l']
bc_pll=['ts_','ts_','ts_','ts_','ts_','ts_']
bc_dll=['bcl','bcl','bcl','bcl','bcl','bcl']




bc_mcu=st.ModelComparison(master_location,bc_mlu,bc_plu,bc_dlu)
[mn_snp7,std_snp7,min_snp7,max_snp7]=bc_mcu.sf_final_ts()

bc_mcl=st.ModelComparison(master_location,bc_mll,bc_pll,bc_dll)
[mn_snp8,std_snp8,min_snp8,max_snp8]=bc_mcl.sf_final_ts()


# Package
model=gc_mlu+gc_mll+gc_mlu10+gc_mll10+a_mlu+a_mll+bc_mlu+bc_mll+gc_mluA+gc_mllA
group=gc_dlu+gc_dll+gc_dlu10+gc_dll10+a_dlu+a_dll+bc_dlu+bc_dll+gc_dluA+gc_dllA
mn_snp=np.concatenate((mn_snp1,mn_snp2,mn_snp3,mn_snp4,mn_snp5,mn_snp6,mn_snp7,mn_snp8,mn_snp9,mn_snp10),axis=0)
std_snp=np.concatenate((std_snp1,std_snp2,std_snp3,std_snp4,std_snp5,std_snp6,std_snp7,std_snp8,std_snp9,std_snp10),axis=0)
min_snp=np.concatenate((min_snp1,min_snp2,min_snp3,min_snp4,min_snp5,min_snp6,min_snp7,min_snp8,min_snp9,min_snp10),axis=0)
max_snp=np.concatenate((max_snp1,max_snp2,max_snp3,max_snp4,max_snp5,max_snp6,max_snp7,max_snp8,max_snp9,max_snp10),axis=0)

df=pd.DataFrame(data={'Model':model,
                      'Group':group,
                      'mean_snow_fraction':mn_snp,
                      'std_snow_fraction':std_snp,
                      'min_snow_fraction':min_snp,
                      'max_snow_fraction':max_snp})

df.to_csv('model_final_SF.csv',index=False)







