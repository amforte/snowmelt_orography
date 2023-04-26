#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  8 09:26:02 2023

@author: aforte
"""

import sys
sys.path.insert(0,'/Users/aforte/Documents/GitHub/snowmelt_orography')
import stimpy as st

master_location='/Volumes/Choruh/Data/snowmelt_project/model_outputs/'

model_list=['gc1u','gc1l','a1u','a1l','bc1u','bc1l']
prefix_list=['ts_','ts_','ts_','ts_','ts_','ts_']
descript_list=['GC1U','GC1L','A1U','A1L','BC1U','BC1L']
col_list=['black','black','royalblue','royalblue','orange','orange']
style_list=['-','--','-','--','-','--']

mc=st.ModelComparison(master_location,model_list,prefix_list,descript_list)
f1=mc.plot_all2(col_list,style_list,sample_freq=10)

f1.savefig('P2_figure4.pdf',dpi="figure")