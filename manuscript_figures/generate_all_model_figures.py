#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 15:57:45 2023

@author: aforte
"""

import sys
sys.path.insert(0,'/Users/aforte/Documents/GitHub/snowmelt_orography')
import stimpy as st

# master_location='/Volumes/Choruh/Data/snowmelt_project/model_outputs/'
# output_location='/Volumes/Choruh/Data/snowmelt_project/model_figures/'

master_location='/Users/aforte/Documents/Python/snowmelt/model_outputs_v2/'
output_location='/Users/aforte/Documents/Python/snowmelt/model_figures_v2/'

models=['a025u','a05u','a1u','a2u','a4u','a8u',
        'a025l','a05l','a1l','a2l','a4l','a8l',
        'bc025u','bc05u','bc1u','bc2u','bc4u','bc8u',
        'bc025l','bc05l','bc1l','bc2l','bc4l','bc8l',
        'gc025u','gc05u','gc1u','gc2u','gc4u','gc8u',
        'gc025l','gc05l','gc1l','gc2l','gc4l','gc8l',
        'gc025u_10l','gc05u_10l','gc1u_10l','gc2u_10l','gc4u_10l','gc8u_10l',
        'gc025l_10l','gc05l_10l','gc1l_10l','gc2l_10l','gc4l_10l','gc8l_10l',
        'gc025u_Area','gc05u_Area','gc1u_Area','gc2u_Area','gc4u_Area','gc8u_Area',
        'gc025l_Area','gc05l_Area','gc1l_Area','gc2l_Area','gc4l_Area','gc8l_Area',
        'gc025u_RO','gc05u_RO','gc1u_RO','gc2u_RO','gc4u_RO','gc8u_RO',
        'gc025l_RO','gc05l_RO','gc1l_RO','gc2l_RO','gc4l_RO','gc8l_RO',
        'bc025u_RO','bc05u_RO','bc1u_RO','bc2u_RO','bc4u_RO','bc8u_RO',
        'bc025l_RO','bc05l_RO','bc1l_RO','bc2l_RO','bc4l_RO','bc8l_RO',
        'gc1u_10l_1b','gc1u_5b','gc1u_10b','gc1u_1500r','gc1u_2000r',
        'gc1u_1kmBL','gc1u_2kmBL','bc1u_1kmBL','bc1u_2kmBL','gc_ilr_l','gc_ilr_u']

modelnames=['A025U','A05U','A1U','A2U','A4U','A8U',
            'A025L','A05L','A1L','A2L','A4L','A8L',
            'BC025U','BC05U','BC1U','BC2U','BC4U','BC8U',
            'BC025L','BC05L','BC1L','BC2L','BC4L','BC8L',
            'GC025U','GC05U','GC1U','GC2U','GC4U','GC8U',
            'GC025L','GC05L','GC1L','GC2L','GC4L','GC8L',
            'GC025U_10L','GC05U_10L','GC1U_10L','GC2U_10L','GC4U_10L','GC8U_10L',
            'GC025L_10L','GC05L_10L','GC1L_10L','GC2L_10L','GC4L_10L','GC8L_10L',
            'GC025U_A','GC05U_A','GC1U_A','GC2U_A','GC4U_A','GC8U_A',
            'GC025L_A','GC05L_A','GC1L_A','GC2L_A','GC4L_A','GC8L_A',
            'GC025U_R','GC05U_R','GC1U_R','GC2U_R','GC4U_R','GC8U_R',
            'GC025L_R','GC05L_R','GC1L_R','GC2L_R','GC4L_R','GC8L_R',
            'BC025U_R','BC05U_R','BC1U_R','BC2U_R','BC4U_R','BC8U_R',
            'BC025L_R','BC05L_R','BC1L_R','BC2L_R','BC4L_R','BC8L_R',
            'GC1U_10L_1B','GC1U_5B','GC1U_10B','GC1U_1500R','GC1U_2000R',
            'GC1U_1kmBL','GC1U_2kmBL','BC1U_1kmBL','BC1U_2kmBL','GC_ILR_L','GC_ILR_U']


for i in range(len(models)):
        print(modelnames[i])
        output_dir=master_location+models[i]

        mObj=st.Stim1D(output_dir)
        f1=mObj.plot_profile_results2(modelnames[i],50)
        f1.savefig(output_location+modelnames[i]+'.pdf')