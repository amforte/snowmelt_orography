#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan  7 09:40:08 2023

@author: aforte
"""

import sys
sys.path.insert(0,'/Users/aforte/Documents/GitHub/snowmelt_orography')
import stimpy as st

# master_location='/Volumes/Choruh/Data/snowmelt_project/model_outputs/'
master_location='/Users/aforte/Documents/Python/snowmelt/model_outputs_v2'

model_list=['gc1u','gc1l']
prefix_list=['ts_','ts_']
descript_list=['GC1U : GC - 1 mm/yr - Unlinked - 1 Myr','GC1L : GC - 1 mm/yr - Linked  - 0.63 Myr']

mc=st.ModelComparison(master_location,model_list,prefix_list,descript_list)
f1=mc.comp_model_setup_plot([1000000,630000],100,500,seed=15,max_z=1250)
f1.savefig('P2_figure3.pdf',dpi="figure")