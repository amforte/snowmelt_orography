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
master_location='/Users/aforte/Documents/Python/snowmelt/model_outputs_v2/'

model_list=['gc1u','gc1l']
prefix_list=['ts_','ts_']
descript_list=['GC1U : GC - 1 mm/yr - Unlinked','GC1L : GC - 1 mm/yr - Linked']

mc=st.ModelComparison(master_location,model_list,prefix_list,descript_list)
# mc.comp_excd_prob(0.15,6,[np.inf,2.1e6],['k','gray'])

f1=mc.comp_excd_prob2(0.12,6,[np.inf,np.inf],['k','gray'])
f1.savefig('P2_figure6.pdf',dpi="figure")
