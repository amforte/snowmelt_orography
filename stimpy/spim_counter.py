#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 08:37:55 2022

@author: aforte
"""
import numpy as np
import warnings

class SpimCounter:
    def __init__(self,dt,total_time,freq_to_save):
        self.dt = dt
        self.total_time = total_time
        self.freq_to_save = freq_to_save
        self.__check_times()
        self.generate_counts()
        self.i = 0
        print('Counter instance created')
        
        
    def generate_counts(self):
        # Calculate number of timesteps and derivatives
        self.num_steps=np.round((self.total_time)/self.dt,0).astype(int)
        self.num_steps_to_save=np.round((self.total_time/self.freq_to_save)+1,0).astype(int)
        self.steps_to_save=np.linspace(1,(self.total_time)/self.dt,self.num_steps_to_save).astype(int)-1
        self.dif_steps=np.diff(self.steps_to_save)
        
    def __check_times(self):
        if np.mod(self.total_time,self.freq_to_save)!=0:
            self.freq_to_save=self.freq_to_save-np.mod(self.total_time,self.freq_to_save)
            if self.freq_to_save>0:
                print(f'Input frequency to save would produce irregular save intervals, updating to {self.freq_to_save}')
            else:
                raise Exception('Input frequency to save will generate zero saves, please input a frequency to save that can divide the total time without a remainder')
        
        if self.freq_to_save > self.total_time:
            warnings.warn('Frequency to save is less than total time, no timesteps will be saved')
            
    # Update functions
    def update_i(self,new_i):
        self.i = new_i