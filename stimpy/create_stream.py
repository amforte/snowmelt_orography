#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 08:28:44 2022

@author: aforte
"""
import numpy as np
import scipy.integrate as integrate

class Stream:
    def __init__(self, length, ksn_init, bin_size=5000, bin_dim='x',dx=10, theta=0.5, A0=1, 
                 zb=0, gc=1.5, n=0.54, Ac=1e6, theta0=None):
        self.length = length
        self.ksn_init = ksn_init
        self.bin_size = bin_size
        self.bin_dim = bin_dim
        self.__force_full_x_bins()
        self.dx = dx
        self.theta = theta
        self.theta0 = theta0
        self.A0 = A0
        self.zb = zb
        self.gc = gc
        self.n = n
        self.Ac = Ac
        self.x = np.arange(0,self.length+self.dx,self.dx)
        self.__hack_gravelius()
        # Calculate chi
        self.chi =  integrate.cumtrapz((self.A0/self.A)**self.theta,self.x,initial=0)
        self.__force_full_chi_bins()
        self.__force_full_A_bins()
        # Calculate initial z
        if theta0==None:
            self.z0 = self.zb + (self.chi * self.ksn_init)
        else:
            self.z0 = self.zb + (integrate.cumtrapz((self.A0/self.A)**self.theta0,self.x,initial=0)*self.ksn_init)
        # Duplicate initial z to be updated
        self.z=self.z0
        # Calculate initial slope
        self.slope0=np.concatenate(([0],np.diff(self.z)))/self.dx
        # Duplicate intitial slope to be updated
        self.slope=self.slope0
        # Discretize x
        self.__discretize_x()
        # Calculate initial bin centers
        self.__calc_x_centers()
        self.__calc_z_centers()
        # Calculate differential area along profile
        self.__calculate_diff_area()
        # Calculate Area in each bin
        self.__calc_A_centers()
        print('Stream instance created')
                        
    def __force_full_x_bins(self):
        if self.bin_dim=='x':
            if np.mod(self.length,self.bin_size)==0:
                num_bins=np.round(self.length/self.bin_size,0)
                print('There will be '+str(num_bins)+' runoff bins in the x dimension')
            else:
                self.length=self.length+(self.bin_size-np.mod(self.length,self.bin_size))
                num_bins=np.round(self.length/self.bin_size,0)
                print('The profile length has been updated to '+str(self.length)+', there will be '+str(num_bins)+' runoff bins in the x dimension')
                
            self.num_bins = np.round(self.length/self.bin_size,0).astype(int)

    def __force_full_chi_bins(self):
        if self.bin_dim=='chi':
            if np.mod(np.max(self.chi),self.bin_size)==0:
                num_bins=np.round(np.max(self.chi)/self.bin_size,0)
                print('There will be '+str(num_bins)+' runoff bins in the x dimension')
            else:
                num_bins=np.round(np.max(self.chi)/self.bin_size,0).astype(int)
                bns=np.linspace(0,np.max(self.chi),num_bins+1)
                self.bin_size=bns[1]-bns[0]
                print('The bin length has been updated to '+str(self.bin_size)+' in chi units, there will be '+str(num_bins)+' runoff bins in the x dimension')
            self.num_bins = np.round(np.max(self.chi)/self.bin_size,0).astype(int)            

    def __force_full_A_bins(self):
        if self.bin_dim=='A':
            if np.mod(np.max(self.A),self.bin_size)==0:
                num_bins=np.round(np.max(self.A)/self.bin_size,0)
                print('There will be '+str(num_bins)+' runoff bins in the x dimension')
            else:
                num_bins=np.round(np.max(self.A)/self.bin_size,0).astype(int)
                bns=np.linspace(0,np.max(self.A),num_bins+1)
                self.bin_size=bns[1]-bns[0]
                print('The bin length has been updated to '+str(self.bin_size)+' in meters squared, there will be '+str(num_bins)+' runoff bins in the x dimension')
            self.num_bins = np.round(np.max(self.A)/self.bin_size,0).astype(int)  
        
        
    def __hack_gravelius(self):
        # Calculate area based on Gravelius coeffecient
        c = ((1/2)*self.gc*np.sqrt(np.pi))+((1/4)*np.sqrt((self.gc**2)*np.pi - 4))
        h = 1/self.n
        ka = (c**(-1/self.n))
        self.A = ((ka*((self.length-self.x)/1000)**h)*1e6)+self.Ac

    def __discretize_x(self):
        if self.bin_dim=='x':
            self.x_edges = np.arange(0,self.length,self.bin_size)        
            self.chi_edges = None
            self.A_edges = None
            self.ix = np.digitize(self.x,self.x_edges)
            self.uix = np.unique(self.ix)  
        elif self.bin_dim=='chi':
            self.x_edges = None       
            self.chi_edges = np.arange(0,np.max(self.chi),self.bin_size)
            self.A_edges = None
            self.ix = np.digitize(self.chi,self.chi_edges)
            self.uix = np.unique(self.ix)
        elif self.bin_dim=='A':
            self.x_edges = None     
            self.chi_edges = None
            self.A_edges=np.arange(np.max(self.A),self.Ac,-self.bin_size)
            self.ix = np.digitize(self.A,self.A_edges,right=True)
            self.uix = np.unique(self.ix)
    
    def __calculate_diff_area(self):
        A_diff=np.flip(np.diff(np.flip(self.A)))
        self.A_diff=np.concatenate((A_diff,[self.A[-1]]),0)        
     
    def __calc_A_centers(self):
        self.A_cents=np.bincount(self.ix,self.A_diff,self.num_bins)[1:self.num_bins+1]
    
    def __calc_x_centers(self):
        self.x_cents=np.bincount(self.ix,self.x,self.num_bins)[1:self.num_bins+1]/np.bincount(self.ix,None,self.num_bins)[1:self.num_bins+1]

    def __calc_z_centers(self):
        self.z_cents=np.bincount(self.ix,self.z,self.num_bins)[1:self.num_bins+1]/np.bincount(self.ix,None,self.num_bins)[1:self.num_bins+1]
    
    def update_z(self,f):
        z=np.zeros(self.z.shape)
        z[0]=self.zb
        z[1:len(z)]=self.z[1:len(z)]+f
        self.z = z
        self.slope = np.concatenate(([0],np.diff(self.z)))/self.dx
        self.__calc_z_centers()

    def freeze_z(self):
        z_dict={'z':np.copy(self.z),
             'slope':np.copy(self.slope),
             'z_cents':np.copy(self.z_cents)}
        return z_dict
    
    def replace_z(self,z_dict):
        self.z = z_dict['z']
        self.slope = z_dict['slope']
        self.z_cents = z_dict['z_cents']