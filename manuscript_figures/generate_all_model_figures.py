#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 15:57:45 2023

@author: aforte
"""

import sys
sys.path.insert(0,'/Users/aforte/Documents/GitHub/snowmelt_orography')
import stimpy as st

master_location='/Volumes/Choruh/Data/snowmelt_project/model_outputs/'
output_location='/Volumes/Choruh/Data/snowmelt_project/model_figures/'

models=['a025u','a05u','a1u','a2u','a4u','a8u',
        'a025l','a05l','a1l','a2l','a4l','a8l',
        'bc025u','bc05u','bc1u','bc2u','bc4u','bc8u',
        'bc025l','bc05l','bc1l','bc2l','bc4l','bc8l',
        'gc025u','gc05u','gc1u','gc2u','gc4u','gc8u',
        'gc025l','gc05l','gc1l','gc2l','gc4l','gc8l',
        'gc025u_10l','gc05u_10l','gc1u_10l','gc2u_10l','gc4u_10l','gc8u_10l',
        'gc025l_10l','gc05l_10l','gc1l_10l','gc2l_10l','gc4l_10l','gc8l_10l',
        'gc1u_10l_1b','gc1u_5b','gc1u_10b','gc1u_1500r','gc1u_2000r']

modelnames=['A025U','A05U','A1U','A2U','A4U','A8U',
            'A025L','A05L','A1L','A2L','A4L','A8L',
            'BC025U','BC05U','BC1U','BC2U','BC4U','BC8U',
            'BC025L','BC05L','BC1L','BC2L','BC4L','BC8L',
            'GC025U','GC05U','GC1U','GC2U','GC4U','GC8U',
            'GC025L','GC05L','GC1L','GC2L','GC4L','GC8L',
            'GC025U_10L','GC05U_10L','GC1U_10L','GC2U_10L','GC4U_10L','GC8U_10L',
            'GC025L_10L','GC05L_10L','GC1L_10L','GC2L_10L','GC4L_10L','GC8L_10L',
            'GC1U_10L_1B','GC1U_5B','GC1U_10B','GC1U_1500R','GC1U_2000R']


for i in range(len(models)):
        print(modelnames[i])
        output_dir=master_location+models[i]

        mObj=st.Stim1D(output_dir)
        f1=mObj.plot_profile_results2(modelnames[i],50)
        f1.savefig(output_location+modelnames[i]+'.png',dpi="figure")