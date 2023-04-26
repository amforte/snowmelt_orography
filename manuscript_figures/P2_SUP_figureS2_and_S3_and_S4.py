#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 10 09:55:04 2023

@author: aforte
"""

import sys
sys.path.insert(0,'/Users/aforte/Documents/GitHub/snowmelt_orography')
import stimpy as st

master_location='/Volumes/Choruh/Data/snowmelt_project/model_outputs/'

# Figure S3
model_list=['gc1u','gc1u_5b','gc1u_10b','gc1l']
prefix_list=['ts_','ts_','ts_','ts_']
descript_list=['GC1U','GC1U-5B','GC1U-10B','GC1L']
col_list=['black','olivedrab','limegreen','black']
style_list=['-','--','-','--']

mc=st.ModelComparison(master_location,model_list,prefix_list,descript_list)
f1=mc.plot_all2(col_list,style_list,sample_freq=10)
f1.savefig('P2_SUP_figureS2.pdf',dpi="figure")

# Figure S4
model_list=['gc1u','gc1u_2000r','gc1u_1500r']
prefix_list=['ts_','ts_','ts_']
descript_list=['GC1U','GC1U-2000R','GC1U-1500R']
col_list=['black','tomato','darksalmon']
style_list=['-','--','-']

mc=st.ModelComparison(master_location,model_list,prefix_list,descript_list)
f2=mc.plot_all2(col_list,style_list,sample_freq=10)
f2.savefig('P2_SUP_figureS3.pdf',dpi="figure")

# Figure S5
model_list=['gc1u','gc1u_40l','gc1u_30l','gc1u_20l','gc1u_10l','gc1u_100l']
prefix_list=['ts_','ts_','ts_','ts_','ts_','ts_']
descript_list=['GC1U','GC1U-40L','GC1U-30L','GC1U-20L','GC1U-10L','GC1U-100L']
col_list=['black','aqua','deepskyblue','lightsteelblue','dodgerblue','darkblue']
style_list=['-','--','-','--','-','--']

mc=st.ModelComparison(master_location,model_list,prefix_list,descript_list)
f3=mc.plot_all2(col_list,style_list,sample_freq=10)
f3.savefig('P2_SUP_figureS4.pdf',dpi="figure")