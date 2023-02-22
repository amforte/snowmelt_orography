#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 08:31:36 2022

@author: aforte
"""
import numpy as np

class SpimEroder:
    def __init__(self,sobj,uplift,n,theta=0.5,r_type='Rbar',algo='explicit',c_max=None):
        self.uplift = uplift
        self.n = n
        self.theta = theta
        self.r_type = r_type
        self.algo = algo
        self.m = self.n * self.theta
        self.cum_E = np.zeros(len(sobj.z0)) # Initiate cumulative erosion value
        self.cum_U = np.zeros(len(sobj.z0)) # Intitiate cumulate uplift value
        self.c_max = c_max
        if (self.c_max==None) & (self.algo=='explicit'):
            self.c_max=1
        elif (self.c_max==None) & (self.algo=='implicit'):
            self.c_max=2
        print('Eroder instance created')        
        
    def incision(self,sobj,K,q,dt):
        # Assumes k, q, and dt are using units of m and years
        I=-1*(K * (q**self.m) * sobj.slope**self.n)*dt
        I[I>0]=0 # Set areas with positive incision to zero as this implies no erosion
        # The above should not be used in a SPIM context
        # Calculate stability
        a=-1*(K * (q**self.m) * sobj.slope**(self.n-1))
        cfl=(np.max(np.abs(a))*dt)/sobj.dx
        return I,cfl
    
    def newtonraphson(self,zt,ztp1d,dx,KQt):
        tempz = zt
        tol = np.inf
        
        while tol >1e-3:
            ztp1 = tempz - (tempz - zt + KQt * ((tempz-ztp1d)/dx)**self.n)/(1+(self.n*KQt/dx) * ((tempz-ztp1d)/dx)**(self.n-1))
            tol=abs(ztp1-tempz)
            tempz = ztp1
        return ztp1
            
    def fastscape_incision(self,sobj,K,q,dt):
        z=np.zeros(sobj.z.shape)
        z[0]=sobj.zb
        for i in range(1,len(sobj.z)):
            KQt=K[i] * q[i]**self.m * dt
            zt=sobj.z[i]
            ztp1d=z[i-1]
            dx=sobj.dx
            if ztp1d < zt:
                ztp1 = SpimEroder.newtonraphson(self,zt,ztp1d,dx,KQt)
            else:
                ztp1 = zt
            z[i]=ztp1    
        return z - sobj.z  # Return as incision to be consistent with explicit algorithm
    
    def fastscape_incision_linear(self,sobj,K,q,dt):
        z=np.zeros(sobj.z.shape)
        z[0]=sobj.zb
        for i in range(1,len(sobj.z)):
            KQtx=(K[i] * q[i]**self.m * dt)/sobj.dx
            zt=sobj.z[i]
            ztp1d=z[i-1]
            z[i]=(zt+ztp1d*KQtx)/(1+KQtx)
        return z -sobj.z # Return as incision to be consistent with explicit algorithm
               
    def sub_timestep(self,sobj,K,q,dt,sub):
        # sub should be an integer
        # Update timestep
        dtn=dt/sub
        cfl_vec=np.zeros(sub)
        if type(self.uplift)==float:
            u_vec=np.zeros((sub,1))
        elif type(self.uplift)==np.ndarray:
            u_vec=np.zeros((sub,sobj.z.shape[0]))   
        i_vec=np.zeros((sub,sobj.z.shape[0]))
        # Freeze z in case it needs to be reverted
        z_dict=sobj.freeze_z()
        fail=False
        for i in range(sub):
           # Calculate amount of uplift per timestep
           uz=self.uplift*dtn
           u_vec[i,:]=uz
           # Apply uplift to profile
           sobj.update_z(uz)
           # Calculate incision
           [I,cfl]=SpimEroder.incision(self,sobj,K,q,dtn)
           sobj.update_z(I[1:len(I)])
           if np.any(np.isinf(sobj.slope)):
               fail=True
               break
           i_vec[i,:]=I
           cfl_vec[i]=cfl
        # Calculate sums
        if fail:
            outcome=False
            sobj.replace_z(z_dict)
        else:
            u_sum=np.sum(u_vec,axis=0)
            i_sum=np.sum(i_vec,axis=0)
            # Determine outcome and act
            if np.any(cfl_vec>self.c_max):
                outcome=False
                # Undo updates
                sobj.replace_z(z_dict)
            else:
                outcome=True
                if len(u_sum)==1:
                    self.cum_U=self.cum_U + u_sum[0]
                else:
                    self.cum_U[1:len(u_sum)] = self.cum_U[1:len(u_sum)] + u_sum[1:len(u_sum)]
                self.cum_E[1:len(i_sum)] = self.cum_E[1:len(i_sum)] + i_sum[1:len(i_sum)]
                self.I=i_sum
        return outcome
            
    def uplift_and_incise(self,sobj,K,q,dt):
        # Record the input K
        self.K = K
        # Freeze z in case it needs to be reverted
        z_dict=sobj.freeze_z()
        # Calculate amount of uplift per timestep
        uz=self.uplift*dt
        # Apply uplift to profile
        if type(uz)==float:
            sobj.update_z(uz)
        else:
            sobj.update_z(uz[1:])
        # Calculate incision
        if self.algo=='explicit':
            [I,cfl]=SpimEroder.incision(self,sobj,K,q,dt)
            if cfl>self.c_max:
                # Revert z
                sobj.replace_z(z_dict)
                # Start adaptive timestep search
                sub=10
                outcome=False
                while outcome==False:
                    if (dt/sub)<(1/365.25):
                        raise Exception('Possible infinite loop in subdivision of timestep')
                    outcome=SpimEroder.sub_timestep(self,sobj,K,q,dt,sub)
                    if outcome==False:
                        sub=sub*10
            else:   
                if type(uz)==float:
                    self.cum_U = self.cum_U + uz
                elif type(uz)==np.ndarray:
                    self.cum_U[1:len(uz)] = self.cum_U[1:len(uz)] + uz[1:len(uz)]                
                sobj.update_z(I[1:len(I)])
                self.I = I
                self.cum_E[1:len(I)] = self.cum_E[1:len(I)] + I[1:len(I)]

        elif self.algo=='implicit':
            if self.n==1:
                I=SpimEroder.fastscape_incision_linear(self,sobj,K,q,dt)
            else:
                I=SpimEroder.fastscape_incision(self,sobj,K,q,dt)
            if type(uz)==float:
                self.cum_U = self.cum_U + uz
            elif type(uz)==np.ndarray:
                self.cum_U[1:len(uz)] = self.cum_U[1:len(uz)] + uz[1:len(uz)]                
            sobj.update_z(I[1:len(I)])
            a=-1*(K * (q**self.m) * sobj.slope**(self.n-1))
            cfl=(np.max(np.abs(a))*dt)/sobj.dx
            if cfl>self.c_max:
                print(f'Warning : stability criteria exceeded : C = {cfl:.4}')
            self.I = I
            self.cum_E[1:len(I)] = self.cum_E[1:len(I)] + I[1:len(I)]