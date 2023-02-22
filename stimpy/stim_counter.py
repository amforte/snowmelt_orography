#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 08:36:56 2022

@author: aforte
"""
import numpy as np
import warnings

class StimCounter:
    def __init__(self,dt,total_time,freq_to_save,rec_length):
        self.dt = dt
        self.total_time = total_time
        self.freq_to_save = freq_to_save
        self.rec_length = rec_length
        self.__check_times()
        self.generate_counts()
        self.i = 0 # Main loop counter
        self.j = 0 # Random runoff counter
        self.k = 0 # Output counter
        self.gen_rec_i = 0 # Records at which iteration discharge records are generated
        print('Counter instance created')
                
    def generate_counts(self):
        # Calculate number of timesteps and derivatives
        self.num_steps=np.round((self.total_time*365)/self.dt,0).astype(int)
        self.num_steps_to_save=np.round((self.total_time/self.freq_to_save)+1,0).astype(int)
        self.steps_to_save=np.linspace(1,(self.total_time*365)/self.dt,self.num_steps_to_save).astype(int)-1
        self.dif_steps=np.diff(self.steps_to_save)
        self.num_steps_to_gen_record=np.round((self.total_time/self.rec_length)+1,0).astype(int)
        self.steps_to_gen=np.linspace(0,self.total_time*365,self.num_steps_to_gen_record).astype(int)
    
    def __check_times(self):
        if np.mod(self.total_time,self.freq_to_save)!=0:
            self.freq_to_save=self.freq_to_save-np.mod(self.total_time,self.freq_to_save)
            if self.freq_to_save>0:
                print(f'Input frequency to save would produce irregular save intervals, updating to {self.freq_to_save}')
            else:
                raise Exception('Input frequency to save will generate zero saves, please input a frequency to save that can divide the total time without a remainder')
        
        if self.freq_to_save > self.total_time:
            warnings.warn('Frequency to save is less than total time, no timesteps will be saved')
        
        if np.mod(self.freq_to_save,self.rec_length)!=0:
            self.rec_length=self.rec_length-np.mod(self.freq_to_save,self.rec_length)
            if self.freq_to_save>0:
                print(f'Record length is incompatible with frequency to save, setting to {self.rec_length}')
            else:
                raise Exception('Record length will generate zero records, please provide record length that can divide the frequency to save without a remainder')
                
        if self.rec_length>self.freq_to_save:
            self.rec_length=self.freq_to_save
            print(f'Record length must be less than or equal to frequency to save, setting to {self.rec_length}')
            
        if self.dt!=1:
            warnings.warn(f'dt is set {self.dt}, values other than 1 for dt may produce unexpected results')
                
    # Update functions
    def update_i(self,new_i):
        self.i = new_i
        
    def update_j(self,new_j):
        self.j = new_j
    
    def update_k(self,new_k):
        self.k = new_k
    
    def update_gen_rec_i(self,new_gen_rec_i):
        self.gen_rec_i = new_gen_rec_i