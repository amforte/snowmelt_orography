#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan  7 09:40:08 2023

@author: aforte
"""

import sys
sys.path.insert(0,'/Users/aforte/Documents/GitHub/snowmelt_orography')
import stimpy as st

master_location='/Volumes/Choruh/Data/snowmelt_project/model_outputs/'

model_list=['gc1u','gc1l']
prefix_list=['ts_','ts_']
descript_list=['GC1U - 1 Myr','GC1L - 0.6 Myr']

mc=st.ModelComparison(master_location,model_list,prefix_list,descript_list)
f1=mc.comp_model_setup_plot([1000000,600000],100,500,seed=15,max_z=1250)
f1.savefig('figure_5.pdf',dpi="figure")