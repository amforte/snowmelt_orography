#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  1 14:57:40 2022

@author: aforte
"""

import sys
sys.path.insert(0,'/Users/aforte/Documents/GitHub/snowmelt_orography')
import stimpy as st

master_location='/Users/aforte/Documents/snowmelt_project/'
# master_location='/Volumes/Choruh/Data/snowmelt_project/'

## STIM Daily Model ##
output_dir=master_location+'new_model_outputs/gc025l_Area'
# Generate stream  instance
sObj=st.Stream(50000,25,dx=100,bin_size=1.95e7,bin_dim='A')
# Generate runoff instance
rObj=st.GenerateRunoff(sObj,'emp',random_state='linked',location='Greater Caucasus',window=1,max_rlf=2500)
# Generate eroder instance
eObj=st.StimEroder(sObj,0.25e-3)
# Generate counter instance
cObj=st.StimCounter(1,6e6,5000,100)
# Generate model instance
mObj=st.Stim1D(output_dir)
# Run model
mObj.run_new_model(cObj,sObj,eObj,rObj)

