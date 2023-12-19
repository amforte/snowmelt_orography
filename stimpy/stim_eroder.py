#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 08:30:10 2022

@author: aforte
"""
import numpy as np

class StimEroder:
    def __init__(self,sobj,uplift,k_e=1e-11,k_w=15,f=0.08313,omega_a=0.50,
                       omega_s=0.25,alpha=2/3,beta=2/3,a=1.5,tau_c=45):
        self.uplift = uplift
        self.k_e = k_e
        self.k_w = k_w
        self.f = f
        self.omega_a = omega_a
        self.omega_s = omega_s
        self.alpha = alpha
        self.beta = beta
        self.a = a
        self.tau_c = tau_c
        self.k_t = 0.5*1000*(9.81**(2/3))*(self.f**(1/3)) # set to 1000 a la Tucker 2004
        self.gmma = self.a*self.alpha*(1-self.omega_s) # gamma exponent
        self.m = self.a*self.alpha*(1-self.omega_a) # m in erosion law
        self.n = self.a*self.beta # n in erosion law
        self.K = self.k_e*(self.k_t**self.a)*(self.k_w**(-self.a*self.alpha)) #erosional efficiency
        self.Psi_c=self.k_e*(self.tau_c**self.a)
        self.cum_E = np.zeros(len(sobj.z0)) # Initiate cumulative erosion value
        self.cum_U = np.zeros(len(sobj.z0)) # Intitiate cumulate uplift value
        self.days_E = np.zeros(len(sobj.z0)) # Initiate days that erosion is exceeded
        self.max_cfl = 0
        print('Eroder instance created')
        
    def incision(self,sobj,q,qbar,dts):
        I=-1*(self.K*(qbar**(self.m-self.gmma))*(q**self.gmma)*(sobj.slope**self.n) - self.Psi_c)*dts
        I[I>0]=0 # Set areas with positive incision to zero as this implies no erosion
        # Calculate stability
        a=-1*(self.K*(qbar**(self.m-self.gmma))*(q**self.gmma)*(sobj.slope**(self.n-1)) - self.Psi_c)
        cfl=(np.max(np.abs(a))*dts)/sobj.dx
        return I,cfl # Inicision in m
        
    def uplift_and_incise(self,sobj,dt,q,qbar):
        # Convert timestep in days to seconds
        dts=dt*(24*60*60) 
        # Convert uplift in m/yr to m/sec
        u=self.uplift/(365.25*24*60*60)
        # Calculate amount of uplift in meters over timestep
        # Apply uplift
        uz=u*dts
        if type(uz)==float:
            self.cum_U = self.cum_U + uz
            sobj.update_z(uz)
        elif type(uz)==np.ndarray:
            self.cum_U[1:] = self.cum_U[1:] + uz[1:]
            sobj.update_z(uz[1:])
        # Calculate incision
        [I,cfl]=StimEroder.incision(self,sobj,q,qbar,dts)
        # Apply incision
        sobj.update_z(I[1:])
        self.days_E[1:]=self.days_E[1:]+ (I[1:]<0).astype(int)
        self.I = I
        self.cum_E[1:] = self.cum_E[1:] + I[1:]
        self.cfl=cfl
        if cfl > self.max_cfl:
            self.max_cfl=cfl
        if cfl>1:
            print('Warning : Stability Criteria Violated')