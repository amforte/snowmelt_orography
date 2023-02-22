#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 08:26:34 2022

@author: aforte
"""
import numpy as np
from scipy.stats import weibull_min
from scipy.special import gamma
import scipy.integrate as integrate

class StimSteady:
    def __init__(self,k_e=1e-11,k_w=15,f=0.08313,omega_a=0.50,
                       omega_s=0.25,alpha=2/3,beta=2/3,a=1.5,tau_c=45):
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
        
    def stim_inst_ero_law_dim(self,Ks,Rbar,R):
        # Assumes input Rbar and R values have already been converted to m/s
        # Multiplier to convert output to erosion in m/Ma
        cnv=(1e6*365.25*24*60*60)
        I=(self.K*(Ks**self.n)*(Rbar**(self.m-self.gmma))*(R**self.gmma))-self.Psi_c
        return I*cnv
    
    def stim_ero_integrand_dim(self,Ks,Rbar,R,cr):
        # Assumes input Rbar and R values have already been converted to m/s
        I=StimSteady.stim_inst_ero_law_dim(self,Ks,Rbar,R)
        # Recast scale to account for units change
        sr_u=Rbar/gamma(1+(1/cr))
        pdf=weibull_min.pdf(R,cr,loc=0,scale=sr_u)
        return I*pdf
    
    def stim_integrate_one_dim(self,Ks,Rbar,cr):
        # Convert R in mm/day to m/s
        Rbar=Rbar/(24*60*60*10*100) 
        # Solve for threshold runoff
        Rc=((self.K*(Rbar**(self.m-self.gmma))*(Ks**self.n))/self.Psi_c)**(-1/self.gmma)
        # Set integration parameters in terms of Qstar
        q_min=0.001
        q_max=10000
        # Redefine intregration paramters in terms of runoff
        r_min=Rbar*q_min
        r_max=Rbar*q_max
        if Rc < r_min:
            Rc = r_min
        elif Rc > r_max:
            Rc = r_max-1
            if Rc<0:
                Rc = r_max - (r_max*1e-1)
        # Generate array over which to integrate
        Rvec=np.logspace(np.log10(Rc),np.log10(r_max),1000)
        Ipdf=StimSteady.stim_ero_integrand_dim(self,Ks,Rbar,Rvec,cr)
        E=integrate.simpson(Ipdf,Rvec)
        return E,Rc*(24*60*60*10*100)
    
    def stim_range_dim(self,Rbar,cr,max_ksn=700,num_points=200,space_type='log'):
        # Initialize variables
        if space_type=='log':
            Ks=np.logspace(0,np.log10(max_ksn),num=num_points)
        elif space_type=='lin':
            Ks=np.linspace(1,max_ksn,num=num_points)        
        E=np.zeros((num_points,1))
        Rc=np.zeros((num_points,1))
        for i in range(len(Ks)):
            [E[i],Rc[i]]=StimSteady.stim_integrate_one_dim(self,Ks[i],Rbar,cr)               
        return Ks,E,Rc
    
    def find_effective_r(self,Ks,Rbar,cr):
        # Convert R in mm/day to m/s
        Rbar=Rbar/(24*60*60*10*100) 
        # Solve for threshold runoff
        Rc=((self.K*(Rbar**(self.m-self.gmma))*(Ks**self.n))/self.Psi_c)**(-1/self.gmma)
        # Set integration parameters in terms of Qstar
        q_min=0.001
        q_max=10000
        # Redefine intregration paramters in terms of runoff
        r_min=Rbar*q_min
        r_max=Rbar*q_max   
        if Rc < r_min:
            Rc = r_min
        elif Rc > r_max:
            Rc = r_max-1
            if Rc<0:
                Rc = r_max - (r_max*1e-1)
        # Generate array over which to integrate
        Rvec=np.logspace(np.log10(Rc),np.log10(r_max),1000)
        Ipdf=StimSteady.stim_ero_integrand_dim(self,Ks,Rbar,Rvec,cr)
        # Geomorphically effective runoff, i.e., max of ero integrand - Wolman & Miller, 1960
        ix=np.argmax(Ipdf)
        Reff=Rvec[ix]*(24*60*60*10*100) # Convert back to mm/day
        return Reff
    
    def find_spim_k(self,uplift,Rbar,cr,n,theta):
        # Calculate m for SPIM values
        m=n*theta
        # Convert uplift in m/yr to m/Myr
        U=uplift*1e6
        # Calculate the ksn and E relationship using STIM
        [Ks,E,_]=StimSteady.stim_range_dim(self,Rbar,cr,space_type='lin')
        # Find Ks closest to the uplift rate
        ix=np.argmin(np.abs(U-E))
        Ksn=Ks[ix]
        # Calculate residual
        res=np.min(np.abs(U-E))
        # Find effective R for this Ksn
        Reff=StimSteady.find_effective_r(self,Ksn,Rbar,cr)
        # Convert Rbar and Reff to m/yr
        RbarC=(Rbar/(10*100))*365.25
        ReffC=(Reff/(10*100))*365.25
        # Solve for "K"
        KR=(Ksn/(uplift**(1/n)))**(-n)
        Kbar=KR/(RbarC**m)
        Keff=KR/(ReffC**m)
        return Kbar,Keff,Reff,Ksn,res # Reff is still in mm/day for routing functions