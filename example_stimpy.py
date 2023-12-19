#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Example script for using the stimpy module to run 1D STIM profiles

The stimpy module contains a variety of classes that are used within a single
stimpy run. The classes are:
- Stream - class containing basic stream geometry and methods for updating stream profiles
- GenerateRunoff - class containing runoff and discharge details along profile and methods for 
	calculating and updating details of the runoff as topography evolves
- StimEroder - class containing methods for updating topography through a finite difference solution
	of the stochastic threshold incision equations
- StimCounter - class containing information on timesteps and methods for incrementing timesteps
- Stim1D - model class that incorporates all other classes and methods for plotting results of model
"""

import stimpy as st

# Primary location to store outputs
master_location='/Users/person/model_outputs/'
#  Folder to store individual timesteps
output_dir=master_location+'stim_test'

## Generate stream  instance
# Inputs are:
# 	length of stream in meters (25000)
#	number of bins for calculating runoff and variability (25)
#	various optional inputs (see 'stream.py'), e.g., dx
sObj=st.Stream(25000,25,dx=100)

## Generate runoff instance
# Inputs are:
#	stream instance (sObj)
#	the type of rule for relating topography to runoff and variability ('emp')
#	if rule is 'emp' or 'emp_rain', must specify a valid location for which relationships have been generated
#		(location='Greater Caucasus')
#	control on whether events are linked or unlinked, (random_state='unlinked')
rObj=st.GenerateRunoff(sObj,'emp',random_state='unlinked',location='Greater Caucasus')

## Generate eroder instance
# Inputs are:
#	stream instance (sObj)
#	uplift rate in m/yr (1e-3)
eObj=st.StimEroder(sObj,1e-3)
# Uplift can also spatially vary, but must be same dimensions as x
# u=np.linspace(1e-3,2e-3,len(sObj.x))
# eObj=st.StimEroder(sObj,u)

## Generate counter instance
# Inputs are:
#	timestep in days (1) - using a value other than 1 will produce unexpected results
#	length of run in years (500000)
#	frequency of outputs being saved in years (5000)
#	length of precomputed runoff in years (100) - shorter values will increase computation time
cObj=st.StimCounter(1,500000,5000,100)

## Generate model instance
# Inputs are:
#	location to store outputs (output_dir)
mObj=st.Stim1D(output_dir)
# Various methods within model instance
# Run model
mObj.run_new_model(cObj,sObj,eObj,rObj)
# Restart model
# Inputs are:
#	Time to restart the model (500000)
#	New time to run model to (1000000)
mObj.restart_model(500000,1000000)
# Restart model with different save location and updated timestep (e.g., to run a portion of the
# model at a finer resolution)
# Inputs are:
#	Time to restart model (300000)
#	New time to run model to (310000)
#	Updated save frequency (200)
#	New location to save ouptuts to (master_location+'stim_subsample_test')
mObj.restart_model(300000,310000,200,master_location+'stim_subsample_test')
# Plot profile results
# Inputs are:
#	title for grapth ('STIM Model - Restart')
#	plot every 'n' timesteps (n_step=10)
mObj.plot_profile_results('STIM Model - Restart',n_step=10)
# Run a monte carlo model on the last time step
mc_dict=rObj.monte_carlo_runoff(sObj,100,100)
rObj.plot_monte_carlo(sObj,mc_dict)