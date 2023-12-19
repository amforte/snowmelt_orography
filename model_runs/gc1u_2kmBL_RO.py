#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  1 14:56:00 2022

@author: aforte
"""

import sys
sys.path.insert(0,'/Users/aforte/Documents/GitHub/snowmelt_orography')
import stimpy as st

master_location='/Users/aforte/Documents/snowmelt_project/'

## STIM Daily Model ##
output_dir=master_location+'model_outputs/gc1u_2kmBL_RO'
# Generate stream  instance
sObj=st.Stream(50000,25,dx=100,bin_size=2000,zb=2000)
# # Generate runoff instance
rObj=st.GenerateRunoff(sObj,'emp_rain',random_state='unlinked',location='Greater Caucasus',window=1,max_rlf=2500)
# # Generate eroder instance
eObj=st.StimEroder(sObj,1e-3)
# # Generate counter instance
cObj=st.StimCounter(1,3.5e6,5000,100)
# Generate model instance
mObj=st.Stim1D(output_dir)
# Run model
mObj.run_new_model(cObj,sObj,eObj,rObj)
