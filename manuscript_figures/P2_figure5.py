#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  8 09:26:02 2023

@author: aforte
"""

import sys
sys.path.insert(0,'/Users/aforte/Documents/GitHub/snowmelt_orography')
import stimpy as st
import numpy as np

# master_location='/Volumes/Choruh/Data/snowmelt_project/model_outputs/'
master_location='/Users/aforte/Documents/Python/snowmelt/model_outputs_v2'

model_list=['gc1u','gc1l']
prefix_list=['ts_','ts_']
descript_list=['GC1U : GC - 1 mm/yr - Unlinked','GC1L : GC - 1 mm/yr - Linked']

mc=st.ModelComparison(master_location,model_list,prefix_list,descript_list)
# f1=mc.comp_profile_evol(25,3.25,3.5,35,[4.5e6,2.5e6])
f1=mc.comp_profile_evol2(40,3.25,3.5,10,[np.inf,np.inf])

f1.savefig('P2_figure5.pdf',dpi="figure")