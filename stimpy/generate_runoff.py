#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 08:32:47 2022

@author: aforte
"""
import numpy as np
import pandas as pd
from scipy.stats import weibull_min
from scipy.special import gamma
import matplotlib.pyplot as plt
from cmcrameri import cm
import matplotlib.gridspec as gridspec

import os

class GenerateRunoff:
    def __init__(self, sobj, pattern, random_state='linked',location=None,
                 r_constant=1,r_bottom=0.5,r_top=1,z_bottom=100,z_top=2500,
                 cr_a=0.109,cr_b=0.21,pexp=1,cpos=0.5,r_center=2,window=1,
                 min_runoff=None,kr_slope=None,kr_int=None,max_rlf=None): 
        self.pattern = pattern
        self.random_state = random_state
        self.location = location
        self.r_constant = r_constant
        self.r_bottom = r_bottom
        self.r_top = r_top
        self.z_bottom = z_bottom
        self.z_top = z_top
        self.cr_a = cr_a
        self.cr_b = cr_b
        self.pexp = pexp
        self.cpos = cpos
        self.r_center = r_center
        self.window = window
        self.min_runoff = min_runoff
        self.kr_slope = kr_slope
        self.kr_int = kr_int
        self.max_rlf = max_rlf
        self.__parse_pattern(sobj)
        print('Runoff instance created')
        
    def __parse_pattern(self,sobj):
        if self.pattern=='emp':
            GenerateRunoff.read_emp(self)
        elif self.pattern=='emp_rain':
            GenerateRunoff.read_emp_rain(self)
        GenerateRunoff.update_r(self,sobj)
        if np.mod(self.window,2)==1:
            pass
        else:
            raise Exception("Value provide to 'window' must be odd")
        
    def __moving_average_con(self,a, w):
        # This moving average solution has problematic edge handling, especially on log-transformed
        return np.convolve(a,np.ones(w),'same')/w
    
    def __moving_average(self,a,n):
        # Rolling average with edge handling
        N = len(a)
        cumsum_vec = np.cumsum(np.insert(np.pad(a,(n-1,n-1),'constant'), 0, 0)) 
        d = np.hstack((np.arange(n//2+1,n),np.ones(N-n)*n,np.arange(n,n//2,-1)))  
        return (cumsum_vec[n+n//2:-n//2+1] - cumsum_vec[n//2:-n-n//2]) / d    
        
    def read_emp(self):
        pt = os.path.dirname(os.path.realpath(__file__))
        top = pd.read_csv(os.path.join(pt,'topo_snow_mean_relationships.csv'))
        sno = pd.read_csv(os.path.join(pt,'snowmelt_meanR_cR_relationships.csv'))
        snosr = pd.read_csv(os.path.join(pt,'snowmelt_approx_sR_fit_sR_relationships.csv'))
        kr = pd.read_csv(os.path.join(pt,'ksn_relief_relationships.csv'))
        # Generate indices to isolate relationships of interest
        idx1 = (top['location']==self.location) & (top['relation_between']=='rlf to mR')
        idx2 = (top['location']==self.location) & (top['relation_between']=='maxZ to snowP')
        idx3 = kr['location']==self.location
        # Extract coefficients and exponents for the topography
        rlf_coeff=top.loc[idx1,'param1'].to_numpy()[0]
        rlf_expo=top.loc[idx1,'param2'].to_numpy()[0]
        mxz_coeff=top.loc[idx2,'param1'].to_numpy()[0]
        mxz_expo=top.loc[idx2,'param2'].to_numpy()[0]    
        topo_func_vals=np.array([rlf_coeff,rlf_expo,mxz_coeff,mxz_expo])
        # Convert snow to bins
        snowPL=sno['snowP_left'].to_numpy()
        snowPR=sno['snowP_right'].to_numpy()
        snow_bins=np.concatenate((snowPL,[snowPR[-1]+0.001]),axis=0)
        # Extract relationships
        snow_form=sno['relationship'].to_list()
        snow_1=sno['param1'].to_numpy()
        snow_2=sno['param2'].to_numpy()
        snow_func_vals=np.concatenate((snow_1.reshape(len(snow_1),1),snow_2.reshape(len(snow_1),1)),axis=1)
        snowsr_form=snosr['relationship'].to_list()
        snowsr_1=snosr['param1'].to_numpy()
        snowsr_2=snosr['param2'].to_numpy()
        snowsr_func_vals=np.concatenate((snowsr_1.reshape(len(snowsr_1),1),snowsr_2.reshape(len(snowsr_1),1)),axis=1)
        # Extract ksn relief relationships
        kr_slope=kr.loc[idx3,'param1'].to_numpy()[0]
        kr_inter=kr.loc[idx3,'param2'].to_numpy()[0]
        kr_func_vals=np.array([kr_slope,kr_inter])
        # Package and return
        self.emp_values=[topo_func_vals,snow_bins,snow_form,snow_func_vals,kr_func_vals,snowsr_form,snowsr_func_vals]

    def read_emp_rain(self):
        pt = os.path.dirname(os.path.realpath(__file__))
        top = pd.read_csv(os.path.join(pt,'topo_snow_mean_relationships.csv'))
        kr = pd.read_csv(os.path.join(pt,'ksn_relief_relationships.csv'))
        crsr = pd.read_csv(os.path.join(pt,'rainfall_relationships.csv'))
        # Generate indices to isolate relationships of interest
        idx1 = (top['location']==self.location) & (top['relation_between']=='rlf to mR')
        idx21 = crsr['function']=='rbar_to_cR'
        idx22 = crsr['function']=='sRm_to_sRf'
        idx3 = kr['location']==self.location
        # Extract coefficients and exponents for the topography
        rlf_coeff=top.loc[idx1,'param1'].to_numpy()[0]
        rlf_expo=top.loc[idx1,'param2'].to_numpy()[0]
        topo_func_vals=np.array([rlf_coeff,rlf_expo])
        # Extract coefficients and exponents for cr and sr
        cr_coeff=crsr.loc[idx21,'param1'].to_numpy()[0]
        cr_exp=crsr.loc[idx21,'param2'].to_numpy()[0]
        sr_slp=crsr.loc[idx22,'param1'].to_numpy()[0]
        sr_int=crsr.loc[idx22,'param2'].to_numpy()[0]
        crsr_func_vals=np.array([cr_coeff,cr_exp,sr_slp,sr_int])
        # Extract ksn relief relationships
        kr_slope=kr.loc[idx3,'param1'].to_numpy()[0]
        kr_inter=kr.loc[idx3,'param2'].to_numpy()[0]
        kr_func_vals=np.array([kr_slope,kr_inter])
        # Package and return
        self.emp_rain_values=[topo_func_vals,crsr_func_vals,kr_func_vals]  

    def z_to_runoff_empirical_rain(self,sobj):
        topo_func_vals=self.emp_rain_values[0]
        kr_func_vals=self.emp_rain_values[2]
        crsr_func_vals=self.emp_rain_values[1]
        # Calculate ksn over whole profile
        ksn=np.diff(sobj.z-sobj.zb)/np.diff(sobj.chi)
        ksn=np.concatenate(([ksn[0]],ksn),axis=0)
        # Estimate relief by bins based on the empirical relationship between ksn and relief
        mean_ksn=np.bincount(sobj.ix,ksn,sobj.num_bins)[1:sobj.num_bins+1]/np.bincount(sobj.ix,None,sobj.num_bins)[1:sobj.num_bins+1]
        if (self.kr_slope==None) & (self.kr_int==None):
            rlf_Z=kr_func_vals[0]*mean_ksn + kr_func_vals[1]
        elif (self.kr_slope!=None) & (self.kr_int!=None):
            rlf_Z=self.kr_slope*mean_ksn + self.kr_int
        # Apply threshold relief if provided
        if self.max_rlf!=None:
            rlf_Z[rlf_Z>self.max_rlf]=self.max_rlf
        # Find index of first element in each bin
        fix=np.concatenate(([1],np.diff(sobj.ix))).astype('bool')
        # Use rlf_Z + base elevation of each segment to get max elevation
        max_Z=rlf_Z+sobj.z[fix]        
        # Convert these to mean runoff and snow percentages
        r=topo_func_vals[0]*rlf_Z**topo_func_vals[1]  
        # Perform checks and replacements - this should no longer be necessary
        if self.min_runoff==None:
            r[r<=0]=1e-5 # Set floor to avoid negative runoffs
        else:
            r[r<self.min_runoff]=self.min_runoff  
        # Calculate cr and sr
        cr=crsr_func_vals[0]*r**crsr_func_vals[1]
        srm=r/gamma(1+(1/cr))
        sr=crsr_func_vals[2]*srm + crsr_func_vals[3]
        # At very low values of cR, sR can end up negative because of the linear fit
        sr[sr<0]=srm[sr<0]
        # Update values within instance
        self.r_bin = r
        self.cr_bin =cr
        self.sr_bin = sr
        self.max_Z_bin = max_Z
        self.rlf_Z_bin = rlf_Z        

    def z_to_runoff_empirical(self,sobj):
        # Unpack list
        topo_func_vals=self.emp_values[0]
        snow_bins=self.emp_values[1]
        snow_form=self.emp_values[2]
        snow_func_vals=self.emp_values[3]
        kr_func_vals=self.emp_values[4]
        snowsr_form=self.emp_values[5]
        snowsr_func_vals=self.emp_values[6]
        # Calculate ksn over whole profile
        ksn=np.diff(sobj.z-sobj.zb)/np.diff(sobj.chi)
        ksn=np.concatenate(([ksn[0]],ksn),axis=0)
        # Estimate relief by bins based on the empirical relationship between ksn and relief
        mean_ksn=np.bincount(sobj.ix,ksn,sobj.num_bins)[1:sobj.num_bins+1]/np.bincount(sobj.ix,None,sobj.num_bins)[1:sobj.num_bins+1]
        if (self.kr_slope==None) & (self.kr_int==None):
            rlf_Z=kr_func_vals[0]*mean_ksn + kr_func_vals[1]
        elif (self.kr_slope!=None) & (self.kr_int!=None):
            rlf_Z=self.kr_slope*mean_ksn + self.kr_int
        # Apply threshold relief if provided
        if self.max_rlf!=None:
            rlf_Z[rlf_Z>self.max_rlf]=self.max_rlf
        # Find index of first element in each bin
        fix=np.concatenate(([1],np.diff(sobj.ix))).astype('bool')
        # Use rlf_Z + base elevation of each segment to get max elevation
        max_Z=rlf_Z+sobj.z[fix]        
        # Convert these to mean runoff and snow percentages
        r=topo_func_vals[0]*rlf_Z**topo_func_vals[1]
        snp=topo_func_vals[2]*max_Z**topo_func_vals[3]
        # Perform checks and replacements - this should no longer be necessary
        if self.min_runoff==None:
            r[r<=0]=1e-5 # Set floor to avoid negative runoffs
        else:
            r[r<self.min_runoff]=self.min_runoff
        snp[snp>1]=1 # Set ceiling for snow percent
        snp[snp<0]=0 # Set floor for snow percent
        # Use mean r and snowmelt to predict cr
        cr=np.zeros(len(sobj.uix))
        sr=np.zeros(len(sobj.uix))
        # Find which snow bin each segment belongs to
        nix=np.digitize(snp,snow_bins)-1
        for i in range(len(sobj.uix)):
            # Determine cr from r and relationship
            if snow_form[nix[i]]=='linear':
                cr[i]=snow_func_vals[nix[i],0]*r[i] + snow_func_vals[nix[i],1]
            elif snow_form[nix[i]]=='power':
                cr[i]=snow_func_vals[nix[i],0]*r[i] ** snow_func_vals[nix[i],1]
            # Estimate sr from cr and mean
            srm=r[i]/gamma(1+(1/cr[i]))
            # Apply correction to scale
            sr0=snowsr_func_vals[nix[i],0]*srm + snowsr_func_vals[nix[i],1]
            # At very low values of cR, sR can end up negative because of the linear
            # fit
            if sr0>0:
                sr[i]=sr0
            else:
                sr[i]=srm
        # Update values within instance
        self.r_bin = r
        self.cr_bin =cr
        self.sr_bin = sr
        self.snp_bin = snp
        self.max_Z_bin = max_Z
        self.rlf_Z_bin = rlf_Z

    # def runoff_to_shape(self):
    #     # r- runoff in mm_day
    #     # c = a*mar^b
    #     # default a and b from CONUS Rossi et al., 2016
    #     mar=self.r_bin*365.25
    #     c=self.cr_a*mar**self.cr_b
    #     s=self.r_bin/gamma(1+(1/c))
    #     self.cr_bin=c
    #     self.sr_bin=s  

    # def z_to_runoff_linear(self,sobj):
    #     # Define slope and y intercept based on inputs
    #     m=(self.r_top-self.r_bottom)/(self.z_top-self.z_bottom)
    #     b=self.r_top-(m*self.z_top)
    #     # Generate binned runoffs
    #     r=m*sobj.z_cents + b
    #     if self.min_runoff==None:
    #         r[r<=0]=1e-5 # Set floor to avoid negative runoffs
    #     else:
    #         r[r<self.min_runoff]=self.min_runoff
    #     self.r_bin=r
        
    # def z_to_runoff_power(self,sobj):
    #     # Define constants
    #     a=(self.r_bottom-self.r_top)/(self.z_bottom**self.pexp - self.z_top**self.pexp)
    #     c=self.r_bottom-a*self.z_bottom**self.pexp
    #     # Generate binned runoffs
    #     r=a*sobj.z_cents**self.pexp + c
    #     if self.min_runoff==None:
    #         r[r<=0]=1e-5 # Set floor to avoid negative runoffs
    #     else:
    #         r[r<self.min_runoff]=self.min_runoff
    #     self.r_bin=r
        
    # def z_to_runoff_peaked(self,sobj):        
    #     # Find x and z position of peak
    #     xpos=np.max(sobj.x)*self.cpos
    #     xix=np.argmin(np.abs(sobj.x-xpos))
    #     z_peak=sobj.z[xix]
    #     # Find closest pin of peak
    #     zix=np.argmin(np.abs(sobj.z_cents-z_peak))
    #     zb_peak=sobj.z_cents[zix]
    #     # Determine scale if rc based on elevation of profile
    #     sf=(np.max(sobj.z_cents)/self.z_top)
    #     if sf>1:
    #         sf=1
    #     r_center=self.r_center*sf
    #     # Determine slopes and y-intercepts for each part
    #     m1=(r_center-self.r_bottom)/(zb_peak-self.z_bottom)
    #     b1=r_center-(m1*zb_peak)
    #     m2=(self.r_top-r_center)/(self.z_top-zb_peak)
    #     b2=self.r_top-(m2*self.z_top)
    #     # Fill runoff in piecewise
    #     r=np.zeros(sobj.z_cents.shape)
    #     r[0:zix]=m1*sobj.z_cents[0:zix]+b1
    #     r[zix:len(r)]=m2*sobj.z_cents[zix:len(r)]+b2
    #     if self.min_runoff==None:
    #         r[r<=0]=1e-5 # Set floor to avoid negative runoffs
    #     else:
    #         r[r<self.min_runoff]=self.min_runoff
    #     self.r_bin = r
        
    def constant_runoff(self,sobj):
        r=np.ones(sobj.uix.shape)*self.r_constant
        self.r_bin=r
        
    def route_discharge(self, sobj,r):
        # r in mm/day - convert to m/sec
        rmsec=r/(10*100*24*60*60)
        # Calculate the discharge added per node
        qdelt=rmsec*sobj.A_diff
        q=np.flip(np.cumsum(np.flip(qdelt)))
        return q # returns q in m3/sec
    
    def route_record_discharge(self,sobj,rts,count):
        rts=np.transpose(rts)
        rday=np.repeat(rts,count,axis=1)/(10*100*24*60*60) 
        qdelt=rday*sobj.A_diff
        qday=np.fliplr(np.cumsum(np.fliplr(qdelt),axis=1))
        return qday
        
    def basin_average_runoff(self, sobj,q):
        # q in m3/sec
        rmsec=q[0]/sobj.A[0]
        return rmsec*(10*100*24*60*60) # returns runoff in mm/day
 
    def area_weighted_mean(self,sobj,val):
        # val is value along stream and representative for delta upstream
        return np.sum(sobj.A_diff*val)/np.max(sobj.A)
    
    def area_exp_weighted_mean(self,sobj,val):
        # val is value along stream and representative for delta upstream
        return np.log(np.sum(sobj.A_diff*np.exp(val))/np.max(sobj.A))    
 
    def find_area_weighted_means(self,sobj):
        mr=np.zeros(sobj.ix.shape)
        x_cr=np.zeros(sobj.ix.shape)
        x_sr=np.zeros(sobj.ix.shape)
        for i in range(len(sobj.uix)):
            mr[sobj.ix==sobj.uix[i]]=self.r_bin[i]
            x_cr[sobj.ix==sobj.uix[i]]=self.cr_bin[i]
            x_sr[sobj.ix==sobj.uix[i]]=self.sr_bin[i]
        aw_mr=self.area_weighted_mean(sobj,mr)
        aw_cr=self.area_exp_weighted_mean(sobj,x_cr)
        aw_sr=self.area_weighted_mean(sobj,x_sr)
        aw_z=self.area_weighted_mean(sobj,sobj.z)
        return aw_mr,aw_cr,aw_sr,aw_z

    def route_binned(self,sobj):
        r=np.zeros(sobj.z.shape)
        for i in range(len(sobj.uix)):
            r[sobj.ix==sobj.uix[i]]=self.r_bin[i]
        if self.window>1:
            r=GenerateRunoff.__moving_average(self,r,self.window)
        self.r=r
        self.qbar=GenerateRunoff.route_discharge(self, sobj,r)
    
    def update_r(self,sobj):
        if self.pattern=='emp':
            GenerateRunoff.z_to_runoff_empirical(self,sobj)
        elif self.pattern=='emp_rain':
            GenerateRunoff.z_to_runoff_empirical_rain(self,sobj)
        # elif self.pattern=='linear':
        #     GenerateRunoff.z_to_runoff_linear(self,sobj)
        #     GenerateRunoff.runoff_to_shape(self)
        # elif self.pattern=='power':
        #     GenerateRunoff.z_to_runoff_power(self,sobj)
        #     GenerateRunoff.runoff_to_shape(self)
        # elif self.pattern=='peaked':
        #     GenerateRunoff.z_to_runoff_peaked(self,sobj)
        #     GenerateRunoff.runoff_to_shape(self)
        # elif self.pattern=='constant':
        #     GenerateRunoff.constant_runoff(self,sobj)
        #     GenerateRunoff.runoff_to_shape(self)
        GenerateRunoff.route_binned(self,sobj)
        
    def spim_calc_K_and_Q(self,spobj,sobj,eobj):
        # Determine number bins
        num_bins=sobj.uix.shape[0]
        # Generate empty arrays
        K_vec=np.zeros(sobj.z.shape)
        R_vec=np.zeros(sobj.z.shape)
        # Start loop
        for i in range(num_bins):
            if type(eobj.uplift)==float:
                [Kbar,Keff,Reff,Ksn,res]=spobj.find_spim_k(eobj.uplift,self.r_bin[i],
                                                           self.cr_bin[i],eobj.n,eobj.theta)
            else:
                mn_u=np.mean(eobj.uplift[sobj.ix==sobj.uix[i]])
                [Kbar,Keff,Reff,Ksn,res]=spobj.find_spim_k(mn_u,self.r_bin[i],
                                                           self.cr_bin[i],eobj.n,eobj.theta)                
                
            if eobj.r_type=='Rbar':
                K_vec[sobj.ix==sobj.uix[i]]=Kbar
                R_vec[sobj.ix==sobj.uix[i]]=self.r_bin[i]
            elif eobj.r_type=='Reff':
                K_vec[sobj.ix==sobj.uix[i]]=Keff
                R_vec[sobj.ix==sobj.uix[i]]=Reff
        
        # Blur
        if self.window>1:
            R_vec=GenerateRunoff.__moving_average(self,R_vec,self.window)
            # Perform rolling average on log transformed K
            K_vec=10**GenerateRunoff.__moving_average(self,np.log10(K_vec),self.window)
        # Route and convert to m3/year from m3/sec
        Q_vec=GenerateRunoff.route_discharge(self,sobj,R_vec)*(60*60*24*365.25)
        # Update q_bar (output as m3/sec to be consistent)
        self.qbar=Q_vec/(60*60*24*365.25)
        return K_vec,Q_vec
            
    def generate_random_runoff(self,rec_length,seed):
        # Determine lenght of record in days
        rec_day=np.round(rec_length*365,0).astype(int)
        rts=np.zeros((len(self.r_bin),rec_day))
        for i in range(len(self.r_bin)):
            if self.random_state=='unlinked':
                rts[i,:]=weibull_min.rvs(self.cr_bin[i],loc=0,scale=self.sr_bin[i],size=rec_day,
                                          random_state=seed+i)
            elif self.random_state=='linked':
                rts[i,:]=weibull_min.rvs(self.cr_bin[i],loc=0,scale=self.sr_bin[i],size=rec_day,
                                          random_state=seed)
        return rts
        
    def route_one_day(self,sobj,rts,j):
        rday=np.zeros(sobj.ix.shape)
        for i in range(len(sobj.uix)):
            rday[sobj.ix==sobj.uix[i]]=rts[i,j]
        self.qd=GenerateRunoff.route_discharge(self,sobj,rday)

    def survive(self,ts):
        ts_sort=np.sort(ts)
        tsn=len(ts_sort)
        tsrank=np.arange(1,tsn+1,1)
        ts_freq_excd=(tsn+1-tsrank)/tsn
        return ts_sort,ts_freq_excd
    
    def weibull_tail_fit(self,x,y,thresh):
        ix=np.nonzero(y<thresh)[0][:1][0]
        xtrim=x[ix:]
        ytrim=y[ix:]
        xts=np.log(xtrim)
        yts=np.log(-np.log(ytrim))      
        [lin,r,rnk,sng,V]=np.polyfit(xts,yts,1,full=True)
        c=lin[0]
        s=np.exp(-1*lin[1]/c)
        return c,s 
            
    def monte_carlo_runoff(self,sobj,rec_length,num_trials,flood_recur=2.5,print_freq=10,verbose=True):
        # Calculate area weighted means
        [aw_mean_runoff,aw_cr,aw_sr,aw_z]=self.find_area_weighted_means(sobj)
        # Generate empty arrays
        tail_fit_cr=np.zeros(num_trials)
        tail_fit_sr=np.zeros(num_trials)
        mm_fit_cr=np.zeros(num_trials)
        mm_fit_sr=np.zeros(num_trials)      
        runoff_means=np.zeros(num_trials)
        obs_tail=np.zeros(num_trials)
        tail_tail=np.zeros(num_trials)
        mm_tail=np.zeros(num_trials)
        
        [_,count]=np.unique(sobj.ix,return_counts=True)
        
        # Begin main trial loop
        for i in range(num_trials):
            if (np.mod(i,print_freq)==0) & (verbose):
                print(f'Iteration {i} of {num_trials} total')
            # Generate the record
            rts=self.generate_random_runoff(rec_length,i)
            # Vectorized solution, much faster            
            qday=self.route_record_discharge(sobj,rts,count)
            r_daily=(qday[:,0]/sobj.A[0])*(10*100*24*60*60)
            #########################################################
            ## Old Loop Solution - SLOW
            # num_days=rts.shape[1]
            # # Generate empty array for output
            # r_daily=np.zeros(num_days)
            # for j in range(num_days):
            #     runoffs=rts[:,j]
            #     rday=np.zeros(sobj.ix.shape)
            #     for k in range(len(sobj.uix)):
            #         rday[sobj.ix==sobj.uix[k]]=runoffs[k]
            #     qdaily=self.route_discharge(sobj,rday)
            #     # max_q_daily[j]=qdaily[0]
            #     r_daily[j]=self.basin_average_runoff(sobj,qdaily)
            #     # r_daily=(max_q_daily/np.max(sobj.A))*(10*100*60*60*24)
            ############################################################
            runoff_means[i]=np.mean(r_daily)
            [r_sort,r_freq]=self.survive(r_daily)
            [tail_fit_cr[i],tail_fit_sr[i]]=self.weibull_tail_fit(r_sort,r_freq,0.01)
            [mm_fit_cr[i],loc,mm_fit_sr[i]]=weibull_min.fit(r_daily,floc=0)
            
            # Find magnitude of recurring flood
            targ=(1/(flood_recur*365.25))
            rix=np.argmin(np.abs(r_freq-targ))
            obs_tail[i]=r_sort[rix]
            tail_tail[i]=weibull_min.isf(targ,tail_fit_cr[i],loc=0,scale=tail_fit_sr[i])
            mm_tail[i]=weibull_min.isf(targ,mm_fit_cr[i],loc=0,scale=mm_fit_sr[i])
        
        # Package as a dictionary for output
        mc_dict={'aw_mean_runoff':aw_mean_runoff,
                  'aw_cr':aw_cr,
                  'aw_sr':aw_sr,
                  'aw_z':aw_z,
                  'tail_fit_cr':tail_fit_cr,
                  'tail_fit_sr':tail_fit_sr,
                  'mm_fit_cr':mm_fit_cr,
                  'mm_fit_sr':mm_fit_sr,
                  'mean_runoffs':runoff_means,
                  'obs_flood':obs_tail,
                  'tail_fit_flood':tail_tail,
                  'mm_fit_flood':mm_tail,
                  'flood_recur':flood_recur}
        return mc_dict
    
    def plot_monte_carlo(self,sobj,mc_dict):
        f1=plt.figure(figsize=(25,12))
        gsm=gridspec.GridSpec(1,2)
        gs0=gsm[0].subgridspec(4,3,wspace=0.4,hspace=0.3)
        gs1=gsm[1].subgridspec(4,4,wspace=0.4,hspace=0.3)
        
        ax1=f1.add_subplot(gs0[0,:])
        sc1=plt.scatter((sobj.x)/1000,sobj.A/1e6,c=self.r,cmap=cm.vik_r)
        cbar1=plt.colorbar(sc1,ax=ax1)
        cbar1.ax.set_ylabel('Runoff [mm/day]')    
        plt.xlabel(r'Distance from Mouth [$km$]')
        plt.ylabel(r'Drainage Area [$km^{2}$]')
        
        ax2=f1.add_subplot(gs0[1,:])
        sc2=plt.scatter((sobj.x)/1000,sobj.z,c=self.r,cmap=cm.vik_r)
        cbar2=plt.colorbar(sc2,ax=ax2)
        cbar2.ax.set_ylabel('Runoff [mm/day]') 
        plt.xlabel(r'Distance from Mouth [$km$]')
        plt.ylabel(r'Elevation [$m$]')
        
        ax3=f1.add_subplot(gs0[2,:])
        sc3=plt.scatter((sobj.x)/1000,self.qbar,c=self.r,cmap=cm.vik_r)
        cbar3=plt.colorbar(sc3,ax=ax3)
        cbar3.ax.set_ylabel('Runoff [mm/day]')
        plt.xlabel(r'Distance from Mouth [$km$]')
        plt.ylabel(r'Discharge [$m^{3}/s$]')
        
        ax4a=f1.add_subplot(gs0[3,0])
        plt.scatter(self.r_bin,self.cr_bin,c='k')
        ax4a.axvline(mc_dict['aw_mean_runoff'],c='k',linestyle=':')
        ax4a.axhline(mc_dict['aw_cr'],c='k',linestyle=':')
        plt.xlabel('Mean Runoff [mm/day]')
        plt.ylabel('Shape')
        
        ax4b=f1.add_subplot(gs0[3,1])
        plt.scatter(self.r_bin,sobj.z_cents,c='k')
        ax4b.axvline(mc_dict['aw_mean_runoff'],c='k',linestyle=':')
        ax4b.axhline(mc_dict['aw_z'],c='k',linestyle=':')
        plt.xlabel('Mean Runoff [mm/day]')
        plt.ylabel('Elevation [m]')
        
        ax4c=f1.add_subplot(gs0[3,2])
        plt.scatter(self.cr_bin,sobj.z_cents,c='k')
        ax4c.axvline(mc_dict['aw_cr'],c='k',linestyle=':')
        ax4c.axhline(mc_dict['aw_z'],c='k',linestyle=':')
        plt.xlabel('Shape')
        plt.ylabel('Elevation [m]')
        
        # Set limits
        xmin=np.min(np.array([np.min(mc_dict['tail_fit_cr']),np.min(self.cr_bin),np.min(mc_dict['mm_fit_cr'])]))
        xmax=np.max(np.array([np.max(mc_dict['tail_fit_cr']),np.max(self.cr_bin),np.max(mc_dict['mm_fit_cr'])]))
        ymin=np.min(np.array([np.min(mc_dict['tail_fit_sr']),np.min(self.sr_bin),np.min(mc_dict['mm_fit_cr'])]))
        ymax=np.max(np.array([np.max(mc_dict['tail_fit_sr']),np.max(self.sr_bin),np.max(mc_dict['mm_fit_cr'])]))
        
        ax5=f1.add_subplot(gs1[1:,3])
        plt.hist(mc_dict['tail_fit_sr'],100,color='gray',zorder=1,label='Tail Fits of Trials',alpha=0.5,orientation='horizontal')
        plt.hist(mc_dict['mm_fit_sr'],100,color='blue',zorder=0,label='MM Fit of Trials',alpha=0.5,orientation='horizontal')
        plt.ylim((ymin-ymin*0.1,ymax+ymax*0.1))
        plt.xlabel('Count')
        
        ax6=f1.add_subplot(gs1[0,:-1])
        plt.hist(mc_dict['tail_fit_cr'],100,color='gray',zorder=0,label='Tail Fits of Trials',alpha=0.5)
        plt.hist(mc_dict['mm_fit_cr'],100,color='blue',zorder=0,label='MM Fit of Trials',alpha=0.5)
        plt.xlim(xmin-xmin*0.1,xmax+xmax*0.1)
        plt.ylabel('Count')
        
        ax7=f1.add_subplot(gs1[1:,0:-1])
        plt.scatter(self.cr_bin,self.sr_bin,c=self.r_bin,cmap=cm.vik_r,s=25,zorder=1,marker='s',label='Elevation Bands')
        plt.scatter(mc_dict['aw_cr'],mc_dict['aw_sr'],c='k',s=60,zorder=1,marker='s',label='Area Weighted Mean')
        plt.scatter(np.median(mc_dict['tail_fit_cr']),np.median(mc_dict['tail_fit_sr']),c='k',s=50,zorder=2,label='Median of Trials')
        plt.plot([np.median(mc_dict['tail_fit_cr']),np.median(mc_dict['tail_fit_cr'])],
                  [np.percentile(mc_dict['tail_fit_sr'],25),np.percentile(mc_dict['tail_fit_sr'],75)],c='k',linewidth=0.5,zorder=1)
        plt.plot([np.percentile(mc_dict['tail_fit_cr'],25),np.percentile(mc_dict['tail_fit_cr'],75)],
                  [np.median(mc_dict['tail_fit_sr']),np.median(mc_dict['tail_fit_sr'])],c='k',linewidth=0.5,zorder=1,label='Interquartile Range')
        plt.scatter(mc_dict['tail_fit_cr'],mc_dict['tail_fit_sr'],c='gray',s=5,alpha=0.5,zorder=0,label='All Trials Tail Fit')
        plt.scatter(mc_dict['mm_fit_cr'],mc_dict['mm_fit_sr'],c='blue',s=5,alpha=0.5,zorder=0,label='All Trials MM Fit')
        plt.xlim(xmin-xmin*0.1,xmax+xmax*0.1)
        plt.ylim((ymin-ymin*0.1,ymax+ymax*0.1))
        plt.legend(loc='best')
        plt.xlabel('Shape')
        plt.ylabel('Scale')
        
        ax8=f1.add_subplot(gs1[0,-1])
        all_vals=np.concatenate((mc_dict['obs_flood'],mc_dict['tail_fit_flood'],mc_dict['mm_fit_flood']))
        plt.plot([np.min(all_vals),np.max(all_vals)],[np.min(all_vals),np.max(all_vals)],c='k',linewidth=1,linestyle=':',zorder=0)
        plt.scatter(mc_dict['obs_flood'],mc_dict['tail_fit_flood'],c='gray',s=5,alpha=0.5,zorder=0,label='All Trials Tail Fit')
        plt.scatter(mc_dict['obs_flood'],mc_dict['mm_fit_flood'],c='blue',s=5,alpha=0.5,zorder=0,label='All Trials MM Fit')
        plt.xlabel(f"Observed {mc_dict['flood_recur']} yr Flood [mm/day]")
        plt.ylabel(f"Fit {mc_dict['flood_recur']} yr Flood [mm/day]")
        
        
