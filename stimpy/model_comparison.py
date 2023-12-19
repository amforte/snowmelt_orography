#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  5 13:14:04 2022

@author: aforte
"""

import os
import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import weibull_min
from scipy.stats import linregress
from scipy.optimize import minimize_scalar
import scipy.integrate as integrate
from . import StimCounter
from . import SpimCounter
from . import Stim1D
from . import Spim1D
from cmcrameri import cm
from matplotlib import colors
from matplotlib import cm as cmm
import matplotlib.gridspec as gridspec
import glob
import re
import string
from scipy import odr

def linear(B,x):
    return B[0]*x + B[1]

def odr_fit(x,y):
    # Filter 0 values
    lx=np.log10(x[(x>0) & (y>0)])
    ly=np.log10(y[(x>0) & (y>0)])
    linmod=odr.Model(linear)
    
    fdlog=odr.Data(lx,ly)
    odrlog=odr.ODR(fdlog,linmod,beta0=[0.1,10])
    outlog=odrlog.run()
    logexp=outlog.beta[0]
    logcoeff=10**outlog.beta[1]
    return logcoeff,logexp

def theta_min(theta,*args):
    A,x,z=args
    c=integrate.cumtrapz((1/A)**theta,x,initial=0)
    res=linregress(c,z)
    return np.abs(1-res.rvalue**2)

class ModelComparison:
    def __init__(self,parent_directory,model_list,prefix_list,descript_list,
                 lower_ts=0,upper_ts=np.inf):
        self.parent_directory = parent_directory
        self.model_list = model_list
        self.prefix_list = prefix_list
        self.descript_list = descript_list
        self.num_models = len(model_list)
        self.lower_ts = lower_ts
        self.upper_ts = upper_ts

    
    def determine_type(self,yr_list):
        type_list=[]
        for i in range(len(yr_list)):
            yr=int(yr_list[i])
            fname=os.path.join(self.parent_directory,self.model_list[i],self.prefix_list[i]+str(yr)+'.pkl')
            with open(fname,'rb') as f:
                cObj=pickle.load(f)
            if isinstance(cObj, StimCounter):
                type_list.append('STIM')
            elif isinstance(cObj, SpimCounter):
                type_list.append('SPIM')
        return type_list
    
    def recover_stim_state(self,fname):
        with open(fname,'rb') as f:
            cObj=pickle.load(f)
            sObj=pickle.load(f)
            eObj=pickle.load(f)
            rObj=pickle.load(f)
        return cObj,sObj,eObj,rObj
    
    def recover_spim_state(self,fname):
        with open(fname,'rb') as f:
            cObj=pickle.load(f)
            spObj=pickle.load(f)
            sObj=pickle.load(f)
            eObj=pickle.load(f)
            rObj=pickle.load(f)
        return cObj,spObj,sObj,eObj,rObj

    def output_result_dict(self,clip_ts,last_n_ts=None):
        # Build ts list
        ts_list=[]
        ts0_list=[]
        ts1_list=[]
        for i in range(self.num_models):
            ts=self.sorted_time_steps(self.model_list[i],self.prefix_list[i])
            ts_list.append(ts)
            ts0_list.append(ts[0])
            ts1_list.append(ts[-1])
        # Determine type
        type_list=self.determine_type(ts0_list)
        
        d_list=[]
        for i in range(self.num_models):
            if type_list[i]=='STIM':
                mObj=Stim1D(os.path.join(self.parent_directory,self.model_list[i]))
            elif type_list[i]=='SPIM':
                mObj=Spim1D(os.path.join(self.parent_directory,self.model_list[i]))
                
            if last_n_ts==None:    
                d=mObj.parse_results(0,clip_ts[i],False,1e-1)
            else:
                tsOI=ts_list[i]
                d=mObj.parse_results(tsOI[-last_n_ts],tsOI[-1],False,1e-1)
            d_list.append(d)
        return d_list    
    
    def comp_model_setup_plot(self,yr_list,rec_length,num_trials,n_skip=1,stream='long profile',seed=1,max_z=None,max_ksn=550):
        letters=list(string.ascii_uppercase)
        # Set figure details
        SMALL_SIZE = 8
        MEDIUM_SIZE = 10
        BIGGER_SIZE = 12
        plt.rc('font', size=SMALL_SIZE,family='Futura')          # controls default text sizes
        plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
        plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
        plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
        plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
        plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
        plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
        # Determine types
        type_list=self.determine_type(yr_list)
        # Initiate figure
        f1=plt.figure(figsize=(8,9))
        f1.set_dpi(250)
        gs=gridspec.GridSpec(3,self.num_models)
        # Start plot
        letter_order=np.arange(0,self.num_models*3,1).reshape((3,self.num_models))
        for i in range(self.num_models):
            yr=int(yr_list[i])
            fname=os.path.join(self.parent_directory,self.model_list[i],self.prefix_list[i]+str(yr)+'.pkl')
            if type_list[i]=='STIM':
                [cObj,sObj,eObj,rObj]=self.recover_stim_state(fname)
            elif type_list[i]=='SPIM':
                [cObj,spObj,sObj,eObj,rObj]=self.recover_spim_state(fname)
            
            ksn=np.diff(sObj.z-sObj.zb)/np.diff(sObj.chi)
            ksn=np.concatenate(([ksn[0]],ksn),axis=0)
            # Estimate relief by bins based on the empirical relationship between ksn and relief
            mean_ksn=np.bincount(sObj.ix,ksn,sObj.num_bins)[1:sObj.num_bins+1]/np.bincount(sObj.ix,None,sObj.num_bins)[1:sObj.num_bins+1]
            
            col_vec=colors.Normalize(vmin=0,vmax=len(sObj.uix)-1)
            
            if stream=='long profile':
                ax1=f1.add_subplot(gs[0,i])
                ax1.set_xlabel('Distance [km]')
                if i==0:
                    ax1.set_ylabel('Elevation [km]')
                ax1.set_title(self.descript_list[i])
                ax1a=ax1.twinx()
                if i==self.num_models-1:
                    ax1a.set_ylabel(r'Mean k$_{sn}$ of Bin')
                ax1a.set_ylim(0,max_ksn)
                
            elif stream=='chi':
                ax1=f1.add_subplot(gs[0,i])
                ax1.set_xlabel(r'$\chi$')
                if i==0:
                    ax1.set_ylabel('Elevation [km]')
                ax1.set_title(self.descript_list[i])
            if max_z!=None:
                ax1.set_ylim((-10/1000,max_z/1000))
            
            ax2=f1.add_subplot(gs[1,i])
            ax2.set_xlabel('Mean Runoff [mm/day]')
            if i==0:
                ax2.set_ylabel('Shape Parameter')
            ax2.set_xlim((-0.25,12.1))
            ax2.set_ylim((0,2.5))
            ax2.set_facecolor('none')
            ax2a=ax2.twinx()
            if i==self.num_models-1:
                ax2a.set_ylabel('Count of Mean Runoff Bin Values')
            ax2a.set_ylim((0,18))
            ax2a.set_zorder(-1)
            ax2a.set_facecolor('none')
            # ax2b=ax2.twiny()
            # ax2b.set_xlabel('Count of Shape Bin Values')
            # ax2b.set_xlim((0,16))
            # ax2b.set_zorder(-2)
            
            # ax2.set_xscale('log')
            
            ax3=f1.add_subplot(gs[2,i])
            ax3.set_xlabel('Runoff [mm/day]')
            if i==0:
                ax3.set_ylabel('Exceedance Frequency')
            ax3.set_yscale('log')
            ax3.set_xscale('log')
            r_vec=np.logspace(-2,2,100)
            ax3.set_ylim((1/(rec_length*365.25),1))
            ax3.set_xlim((10**-2,10**2))
            
            if type_list[i]=='STIM':
                mc_dict=rObj.monte_carlo_runoff(sObj,cObj.rec_length,num_trials,verbose=True)
                print('Using original record length for STIM models')
            elif type_list[i]=='SPIM':
                mc_dict=rObj.monte_carlo_runoff(sObj,rec_length,num_trials,verbose=True)
            
            [_,count]=np.unique(sObj.ix,return_counts=True)
            rts=rObj.generate_random_runoff(1,seed)
            # qday=rObj.route_record_discharge(sObj,rts,count)
            # r_daily=(qday[:,0]/sObj.A[0])*(10*100*24*60*60)
            # mn_day=r_daily[0]
            
            ax2a.hist(rObj.r_bin,np.linspace(0,12,13),zorder=0,color='lightblue',edgecolor='k',alpha=0.5)
            # ax2b.hist(rObj.cr_bin,np.linspace(0,2.5,13),color='lightblue',edgecolor='k',orientation='horizontal',alpha=0.5)

            ax2.scatter(mc_dict['aw_mean_runoff'],mc_dict['aw_cr'],c='w',s=40,zorder=2,marker='s',label='Area Weighted Mean',edgecolor='k')
            ax2.scatter(np.median(mc_dict['mean_runoffs']),np.median(mc_dict['tail_fit_cr']),c='k',s=40,zorder=3,label='Median of Trials')
            
            ax2.plot([np.median(mc_dict['mean_runoffs']),np.median(mc_dict['mean_runoffs'])],
                      [np.percentile(mc_dict['tail_fit_cr'],25),np.percentile(mc_dict['tail_fit_cr'],75)],c='k',linewidth=0.5,zorder=2)
            ax2.plot([np.percentile(mc_dict['mean_runoffs'],25),np.percentile(mc_dict['mean_runoffs'],75)],
                      [np.median(mc_dict['tail_fit_cr']),np.median(mc_dict['tail_fit_cr'])],c='k',linewidth=0.5,zorder=2,label='Interquartile Range')
            ax2.scatter(mc_dict['mean_runoffs'],mc_dict['tail_fit_cr'],c='gray',s=1,alpha=0.5,zorder=1,label='All Trials Tail Fit')
            
            
            
            ax3.plot(r_vec,weibull_min.sf(r_vec,mc_dict['aw_cr'],loc=0,scale=mc_dict['aw_sr']),
                      c='k',linewidth=2,label='Area Weighted Mean',zorder=1,linestyle='--')
            ax3.plot(r_vec,weibull_min.sf(r_vec,np.median(mc_dict['tail_fit_cr']),loc=0,scale=np.median(mc_dict['tail_fit_sr'])),
                    c='k',linewidth=2,label='Median of All Trials',zorder=1)
            # ax3.scatter(mn_day,weibull_min.sf(mn_day,mc_dict['aw_cr'],loc=0,scale=mc_dict['aw_sr']),
            #             marker='s',c='w',zorder=2,edgecolor='k',s=40)
            # ax3.scatter(mn_day,weibull_min.sf(mn_day,np.median(mc_dict['tail_fit_cr']),loc=0,scale=np.median(mc_dict['tail_fit_sr'])),
            #             marker='s',c='k',zorder=2,edgecolor='k',s=40)
            
            if i>0:
                ax3.legend(loc='lower left')

            for j in range(len(sObj.uix)):
                x=sObj.x[sObj.ix==sObj.uix[j]]
                z=sObj.z[sObj.ix==sObj.uix[j]]
                chi=sObj.chi[sObj.ix==sObj.uix[j]]
                if stream=='long profile':
                    # Long Profile
                    ax1.plot(x/1000,z/1000,c=cm.roma(col_vec(j)),linewidth=2,zorder=0)
                    if j==0:
                        ax1a.scatter(np.mean(x/1000),mean_ksn[j],color=cm.roma(col_vec(j)),s=20,marker='o',edgecolor='k',label=r'Mean k$_{sn}$')
                    else:
                        ax1a.scatter(np.mean(x/1000),mean_ksn[j],color=cm.roma(col_vec(j)),s=20,marker='o',edgecolor='k')
                    if j==0:
                        ax1.scatter(x[0]/1000,z[0]/1000,c='k',s=5,zorder=1,marker='s')
                        ax1.scatter(x[-1]/1000,z[-1]/1000,c='k',s=5,zorder=1,marker='s')
                    else:
                        ax1.scatter(x[-1]/1000,z[-1]/1000,c='k',s=5,zorder=1,marker='s')
                elif stream=='chi':
                    # Chi Elevation
                    ax1.plot(chi,z/1000,c=cm.roma(col_vec(j)),linewidth=2,zorder=0)
                    if j==0:
                        ax1.scatter(chi[0],z[0]/1000,c='k',s=5,zorder=1,marker='s')
                        ax1.scatter(chi[-1],z[-1]/1000,c='k',s=5,zorder=1,marker='s')
                    else:
                        ax1.scatter(chi[-1],z[-1]/1000,c='k',s=5,zorder=1,marker='s')
                  # Mean Runoff - Variability  
                ax2.scatter(rObj.r_bin[j],rObj.cr_bin[j],color=cm.roma(col_vec(j)),edgecolor='k',marker='s',s=20,zorder=0)
                # Exceed freq
            for j in range(0,len(sObj.uix),n_skip):
                efreq=weibull_min.sf(r_vec,rObj.cr_bin[j],loc=0,scale=rObj.sr_bin[j])
                ax3.plot(r_vec,efreq,c=cm.roma(col_vec(j)),zorder=0,linewidth=0.5)
                ax3.scatter(rts[j,0],weibull_min.sf(rts[j,0],rObj.cr_bin[j],loc=0,scale=rObj.sr_bin[j]),
                            color=cm.roma(col_vec(j)),zorder=2,edgecolor='k',marker='s',s=20)
            if i>0:
                ax2.legend(loc='upper right')
                ax1a.legend(loc='upper right')
            
            ax1.text(0.01, 0.99, letters[letter_order[0,i]],
                    horizontalalignment='left',
                    verticalalignment='top',
                    transform=ax1.transAxes,
                    fontsize=12,fontweight='extra bold')
            ax2.text(0.01, 0.99, letters[letter_order[1,i]],
                    horizontalalignment='left',
                    verticalalignment='top',
                    transform=ax2.transAxes,
                    fontsize=12,fontweight='extra bold')
            ax3.text(0.96, 0.99, letters[letter_order[2,i]],
                    horizontalalignment='left',
                    verticalalignment='top',
                    transform=ax3.transAxes,
                    fontsize=12,fontweight='extra bold')
        
        plt.tight_layout()
        plt.rcdefaults()
        return f1
        
    def model_setup_plot(self,yr_list,rec_length,num_trials):
        type_list=self.determine_type(yr_list)

        for i in range(self.num_models):
            yr=int(yr_list[i])
            fname=os.path.join(self.parent_directory,self.model_list[i],self.prefix_list[i]+str(yr)+'.pkl')
            if type_list[i]=='STIM':
                [cObj,sObj,eObj,rObj]=self.recover_stim_state(fname)
            elif type_list[i]=='SPIM':
                [cObj,spObj,sObj,eObj,rObj]=self.recover_spim_state(fname)
            
            
            f1=plt.figure(figsize=(25,12))
            gs=gridspec.GridSpec(3,2)
            col_vec=colors.Normalize(vmin=0,vmax=len(sObj.uix)-1)
            
            ax1=f1.add_subplot(gs[0,0])
            ax1.set_xlabel('Distance [m]')
            ax1.set_ylabel('Elevation [m]')
            ax1.set_title(self.descript_list[i])
            
            ax2=f1.add_subplot(gs[1,0])
            ax2.set_xlabel(r'$\chi$')
            ax2.set_ylabel('Elevation [m]')
            
            ax3=f1.add_subplot(gs[2,0])
            ax3.set_xlabel('Mean Runoff [mm/day]')
            ax3.set_ylabel('Variability')
            ax3.set_xscale('log')
            # ax3.set_yscale('log')
            
            ax4=f1.add_subplot(gs[:,1])
            ax4.set_xlabel('Runoff [mm/day]')
            ax4.set_ylabel('Exceedance Frequency')
            ax4.set_yscale('log')
            ax4.set_xscale('log')
            r_vec=np.logspace(-2,2,100)
            ax4.set_ylim((1/(100*365.25),1))
            ax4.set_xlim((10**-2,10**2))
            
            [_,count]=np.unique(sObj.ix,return_counts=True)
            rts=rObj.generate_random_runoff(1,1)
            qday=rObj.route_record_discharge(sObj,rts,count)
            r_daily=(qday[:,0]/sObj.A[0])*(10*100*24*60*60)
            mn_day=r_daily[0]
        
            if type_list[i]=='STIM':
                mc_dict=rObj.monte_carlo_runoff(sObj,cObj.rec_length,num_trials,verbose=True)
                print('Using original record length for STIM models')
            elif type_list[i]=='SPIM':
                mc_dict=rObj.monte_carlo_runoff(sObj,rec_length,num_trials,verbose=True)
            
            ax3.scatter(mc_dict['aw_mean_runoff'],mc_dict['aw_cr'],c='k',s=60,zorder=1,marker='s',label='Area Weighted Mean')
            ax3.scatter(np.median(mc_dict['mean_runoffs']),np.median(mc_dict['tail_fit_cr']),c='k',s=50,zorder=2,label='Median of Trials')
            
            ax3.plot([np.median(mc_dict['mean_runoffs']),np.median(mc_dict['mean_runoffs'])],
                      [np.percentile(mc_dict['tail_fit_cr'],25),np.percentile(mc_dict['tail_fit_cr'],75)],c='k',linewidth=0.5,zorder=1)
            ax3.plot([np.percentile(mc_dict['mean_runoffs'],25),np.percentile(mc_dict['mean_runoffs'],75)],
                      [np.median(mc_dict['tail_fit_cr']),np.median(mc_dict['tail_fit_cr'])],c='k',linewidth=0.5,zorder=1,label='Interquartile Range')
            ax3.scatter(mc_dict['mean_runoffs'],mc_dict['tail_fit_cr'],c='gray',s=5,alpha=0.5,zorder=0,label='All Trials Tail Fit')
            
            ax4.plot(r_vec,weibull_min.sf(r_vec,mc_dict['aw_cr'],loc=0,scale=mc_dict['aw_sr']),
                      c='k',linewidth=4,label='Area Weighted Mean',zorder=1,linestyle='--')
            ax4.plot(r_vec,weibull_min.sf(r_vec,np.median(mc_dict['tail_fit_cr']),loc=0,scale=np.median(mc_dict['tail_fit_sr'])),
                    c='k',linewidth=4,label='Median of All Trials',zorder=1)
            ax4.scatter(mn_day,weibull_min.sf(mn_day,mc_dict['aw_cr'],loc=0,scale=mc_dict['aw_sr']),
                        marker='s',c='w',zorder=2,edgecolor='k',s=100)
            ax4.scatter(mn_day,weibull_min.sf(mn_day,np.median(mc_dict['tail_fit_cr']),loc=0,scale=np.median(mc_dict['tail_fit_sr'])),
                        marker='s',c='k',zorder=2,edgecolor='k',s=100)
            
            
            ax4.legend(loc='best')

            for i in range(len(sObj.uix)):
                x=sObj.x[sObj.ix==sObj.uix[i]]
                z=sObj.z[sObj.ix==sObj.uix[i]]
                chi=sObj.chi[sObj.ix==sObj.uix[i]]
                # Long Profile
                ax1.plot(x,z,c=cm.roma(col_vec(i)),linewidth=2,zorder=0)
                if i==0:
                    ax1.scatter(x[0],z[0],c='k',s=20,zorder=1,marker='s')
                    ax1.scatter(x[-1],z[-1],c='k',s=20,zorder=1,marker='s')
                else:
                    ax1.scatter(x[-1],z[-1],c='k',s=20,zorder=1,marker='s')
                # Chi Elevation
                ax2.plot(chi,z,c=cm.roma(col_vec(i)),linewidth=2,zorder=0)
                if i==0:
                    ax2.scatter(chi[0],z[0],c='k',s=20,zorder=1,marker='s')
                    ax2.scatter(chi[-1],z[-1],c='k',s=20,zorder=1,marker='s')
                else:
                    ax2.scatter(chi[-1],z[-1],c='k',s=20,zorder=1,marker='s')
                  # Mean Runoff - Variability  
                ax3.scatter(rObj.r_bin[i],rObj.cr_bin[i],color=cm.roma(col_vec(i)),s=20,zorder=0)
                # Exceed freq
                efreq=weibull_min.sf(r_vec,rObj.cr_bin[i],loc=0,scale=rObj.sr_bin[i])
                ax4.plot(r_vec,efreq,c=cm.roma(col_vec(i)),zorder=0)
                ax4.scatter(rts[i,0],weibull_min.sf(rts[i,0],rObj.cr_bin[i],loc=0,scale=rObj.sr_bin[i]),
                            color=cm.roma(col_vec(i)),zorder=2,edgecolor='k',marker='s')
            ax3.legend(loc='best')
            
    
    def comp_excd_prob(self,max_ep,max_e,clip_ts,col_list):
        # Build ts list
        ts_list=[]
        ts0_list=[]
        ts1_list=[]
        for i in range(self.num_models):
            ts=self.sorted_time_steps(self.model_list[i],self.prefix_list[i])
            ts_list.append(ts)
            ts0_list.append(ts[0])
            ts1_list.append(ts[-1])
        # Determine type
        type_list=self.determine_type(ts0_list)
        max_ts=np.max(np.array(ts1_list))
        
        # Generate figure details
        letters=list(string.ascii_uppercase)
        # Set figure details
        SMALL_SIZE = 8
        MEDIUM_SIZE = 10
        BIGGER_SIZE = 12
        plt.rc('font', size=SMALL_SIZE,family='Futura')          # controls default text sizes
        plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
        plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
        plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
        plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
        plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
        plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
        # Initiate figure
        f1=plt.figure(figsize=(8,8))
        f1.set_dpi(250)
        gs=gridspec.GridSpec(6,self.num_models)
        num_panels=3
        letter_order=np.arange(0,self.num_models*num_panels,1).reshape((num_panels,self.num_models))


        ax2=f1.add_subplot(gs[4:,0])
        ax2.set_xlabel('Model Time [Myr]')
        ax2.set_ylabel('Frequency of E. Threshold Exceedance')
        
        ax3=f1.add_subplot(gs[4:,1])
        ax3.set_xlabel('Avgerage Erosion Rate [mm/yr]')
        ax3.set_ylabel('Frequency of E. Threshold Exceedance')
     
        
        for i in range(self.num_models):
            if type_list[i]=='STIM':
                mObj=Stim1D(os.path.join(self.parent_directory,self.model_list[i]))
            elif type_list[i]=='SPIM':
                mObj=Spim1D(os.path.join(self.parent_directory,self.model_list[i]))
            d=mObj.parse_results(0,clip_ts[i],False,1e-1) 
            
            x=d['x']
            ts=d['ts']
            fts=d['freq_to_save']
            eCout=d['eCout']
            dout=d['dout']
            
            # Calculate average erosion rate
            avgE=np.diff(eCout,axis=0)
            avgE=np.concatenate((eCout[0,:].reshape((1,len(eCout[0,:]))),avgE),axis=0)
            avgE=((avgE*(-1))/fts)*(10*100)
            
            # Calculate percentage of days where erosion threshold is exceeded
            excdP=np.diff(dout,axis=0)
            excdP=np.concatenate((dout[0,:].reshape((1,len(dout[0,:]))),excdP),axis=0)
            excdP=excdP/(fts*365)
            
            cumP=dout[-1,:]/(ts[-1]*365)
            
            ax1a=f1.add_subplot(gs[3,i])
            ax1a.plot(x[1:]/1000,cumP[1:],c=col_list[i])
            ax1a.set_xlabel('Stream Distance [km]')
            ax1a.set_ylabel('Cum. Ex.')
            ax1a.set_ylim((0,max_ep))
            ax1a.set_xlim((0,np.max(x)/1000))

            
            ax1=f1.add_subplot(gs[0:3,i])
            ax1.set_title(self.descript_list[i])
            im1=ax1.imshow(np.flipud(excdP),extent=[x[0]/1000,x[-1]/1000,ts[0]/1e6,ts[-1]/1e6],
                       cmap=cm.hawaii_r,norm=colors.Normalize(vmin=0,vmax=np.max(excdP.ravel())),aspect='auto')
            cbar1=f1.colorbar(im1,ax=ax1,orientation='horizontal')
            ax1.set_xlabel('Stream Distance [km]')
            ax1.set_ylabel('Model Time [Myr]')
            cbar1.ax.set_xlabel('Frequency of Erosion Threshold Exceedance') 

            
            ax2.plot(ts/1e6,np.mean(excdP,axis=1),label=self.descript_list[i],
                     c=col_list[i],linewidth=1.5)
            ax2.plot(ts/1e6,np.min(excdP,axis=1),
                     linestyle='--',c=col_list[i],linewidth=0.5)
            ax2.plot(ts/1e6,np.max(excdP,axis=1),
                     linestyle='--',c=col_list[i],linewidth=0.5)
            
            ax3.scatter(avgE.ravel(),excdP.ravel(),c=col_list[i],s=1)
            
            
            ax1.text(0.01, 0.99, letters[letter_order[0,i]],
                    horizontalalignment='left',
                    verticalalignment='top',
                    transform=ax1.transAxes,
                    fontsize=12,fontweight='extra bold')
            
            ax1a.text(0.95, 0.99, letters[letter_order[1,i]],
                    horizontalalignment='left',
                    verticalalignment='top',
                    transform=ax1a.transAxes,
                    fontsize=12,fontweight='extra bold')            
                            
        ax2.legend(loc='best')
        ax2.text(0.01, 0.99, letters[letter_order[2,0]],
                horizontalalignment='left',
                verticalalignment='top',
                transform=ax2.transAxes,
                fontsize=12,fontweight='extra bold')
        ax3.text(0.01, 0.99, letters[letter_order[2,1]],
                horizontalalignment='left',
                verticalalignment='top',
                transform=ax3.transAxes,
                fontsize=12,fontweight='extra bold')

        
        f1.tight_layout()
        plt.rcdefaults()
        return f1

    def comp_excd_prob2(self,max_ep,max_e,clip_ts,col_list):
        # Build ts list
        ts_list=[]
        ts0_list=[]
        ts1_list=[]
        for i in range(self.num_models):
            ts=self.sorted_time_steps(self.model_list[i],self.prefix_list[i])
            ts_list.append(ts)
            ts0_list.append(ts[0])
            ts1_list.append(ts[-1])
        # Determine type
        type_list=self.determine_type(ts0_list)
        max_ts=np.max(np.array(ts1_list))
        
        # Generate figure details
        letters=list(string.ascii_uppercase)
        # Set figure details
        SMALL_SIZE = 8
        MEDIUM_SIZE = 10
        BIGGER_SIZE = 12
        plt.rc('font', size=SMALL_SIZE,family='Futura')          # controls default text sizes
        plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
        plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
        plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
        plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
        plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
        plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
        # Initiate figure
        f1=plt.figure(figsize=(8,8),layout='tight')
        f1.set_dpi(250)
        gs=gridspec.GridSpec(6,self.num_models)
        num_panels=4
        letter_order=np.arange(0,self.num_models*num_panels,1).reshape((num_panels,self.num_models))
    

     
        
        for i in range(self.num_models):
            if type_list[i]=='STIM':
                mObj=Stim1D(os.path.join(self.parent_directory,self.model_list[i]))
            elif type_list[i]=='SPIM':
                mObj=Spim1D(os.path.join(self.parent_directory,self.model_list[i]))
            d=mObj.parse_results(0,clip_ts[i],False,1e-1) 
            
            x=d['x']
            ts=d['ts']
            fts=d['freq_to_save']
            eCout=d['eCout']
            dout=d['dout']
            
            # Calculate average erosion rate
            avgE=np.diff(eCout,axis=0)
            avgE=np.concatenate((eCout[0,:].reshape((1,len(eCout[0,:]))),avgE),axis=0)
            avgE=((avgE*(-1))/fts)*(10*100)
            
            # Calculate percentage of days where erosion threshold is exceeded
            excdP=np.diff(dout,axis=0)
            excdP=np.concatenate((dout[0,:].reshape((1,len(dout[0,:]))),excdP),axis=0)
            excdP=excdP/(fts*365)
            
            cumP=dout[-1,:]/(ts[-1]*365)
            
            ax1a=f1.add_subplot(gs[3,i])
            ax1a.plot(x[1:]/1000,cumP[1:],c=col_list[i])
            ax1a.set_xlabel('Stream Distance [km]')
            ax1a.set_ylabel('Cum. Ex.')
            ax1a.set_ylim((0,max_ep))
            ax1a.set_xlim((0,np.max(x)/1000))

            ax23=f1.add_subplot(gs[4:,i])
            ax23.spines['top'].set_color('none')
            ax23.spines['bottom'].set_color('none')
            ax23.spines['left'].set_color('none')
            ax23.spines['right'].set_color('none')
            ax23.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
            ax23.set_ylabel('Frequency of E. Threshold Exceedance')            
            
            ax2=f1.add_subplot(gs[4,i])
            ax2.set_xlabel('Model Time [Myr]')
            
            ax3=f1.add_subplot(gs[5,i])
            ax3.set_xlabel('Avgerage Erosion Rate [mm/yr]')
    
            ax1=f1.add_subplot(gs[0:3,i])
            ax1.set_title(self.descript_list[i])
            im1=ax1.imshow(np.flipud(excdP),extent=[x[0]/1000,x[-1]/1000,ts[0]/1e6,ts[-1]/1e6],
                       cmap=cm.hawaii_r,norm=colors.Normalize(vmin=0,vmax=max_ep),aspect='auto')
            cbar1=f1.colorbar(im1,ax=ax1,orientation='horizontal')
            ax1.set_xlabel('Stream Distance [km]')
            ax1.set_ylabel('Model Time [Myr]')
            cbar1.ax.set_xlabel('Frequency of Erosion Threshold Exceedance') 
    
            
            ax2.plot(ts/1e6,np.mean(excdP,axis=1),label=self.descript_list[i],
                     c=col_list[i],linewidth=1.5)
            ax2.plot(ts/1e6,np.min(excdP,axis=1),
                     linestyle='--',c=col_list[i],linewidth=0.5)
            ax2.plot(ts/1e6,np.max(excdP,axis=1),
                     linestyle='--',c=col_list[i],linewidth=0.5)
            ax2.set_ylim((-0.01,max_ep*1.5))
            

            xb=np.linspace(0,max_e,50)
            yb=np.linspace(0,max_ep,50)
            im2=ax3.hist2d(avgE.ravel(),excdP.ravel(),[xb,yb],norm=colors.LogNorm(vmin=1,vmax=10000),cmap=cm.grayC)
            cbar2=plt.colorbar(im2[3],ax=ax3)
            cbar2.ax.set_ylabel('Density')
            
            ax1.text(0.01, 0.99, letters[letter_order[0,i]],
                    horizontalalignment='left',
                    verticalalignment='top',
                    transform=ax1.transAxes,
                    fontsize=12,fontweight='extra bold')
            
            ax1a.text(0.95, 0.99, letters[letter_order[1,i]],
                    horizontalalignment='left',
                    verticalalignment='top',
                    transform=ax1a.transAxes,
                    fontsize=12,fontweight='extra bold')
            ax2.text(0.01, 0.99, letters[letter_order[2,i]],
                    horizontalalignment='left',
                    verticalalignment='top',
                    transform=ax2.transAxes,
                    fontsize=12,fontweight='extra bold')
            ax3.text(0.01, 0.99, letters[letter_order[3,i]],
                    horizontalalignment='left',
                    verticalalignment='top',
                    transform=ax3.transAxes,
                    fontsize=12,fontweight='extra bold')            
                            
        # ax2.legend(loc='best')

        # f1.tight_layout()
        plt.rcdefaults()
        return f1
            
    def comp_profile_evol(self,num_ts,max_z,max_e,max_q,clip_ts):
        # Build ts list
        ts_list=[]
        ts0_list=[]
        ts1_list=[]
        for i in range(self.num_models):
            ts=self.sorted_time_steps(self.model_list[i],self.prefix_list[i])
            ts_list.append(ts)
            ts0_list.append(ts[0])
            ts1_list.append(ts[-1])
        # Determine type
        type_list=self.determine_type(ts0_list)
        max_ts=np.max(np.array(ts1_list))
        # Account for zero time
        num_ts=num_ts-1
        
        # Generate figure details
        letters=list(string.ascii_uppercase)
        # Set figure details
        SMALL_SIZE = 8
        MEDIUM_SIZE = 10
        BIGGER_SIZE = 12
        plt.rc('font', size=SMALL_SIZE,family='Futura')          # controls default text sizes
        plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
        plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
        plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
        plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
        plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
        plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
        # Initiate figure
        f1=plt.figure(figsize=(8,9))
        f1.set_dpi(250)
        num_panels=4
        gs=gridspec.GridSpec(num_panels,self.num_models)
        letter_order=np.arange(0,self.num_models*num_panels,1).reshape((num_panels,self.num_models))
        col_vec=colors.Normalize(vmin=0,vmax=max_ts)
        
        for i in range(self.num_models):
            if type_list[i]=='STIM':
                mObj=Stim1D(os.path.join(self.parent_directory,self.model_list[i]))
            elif type_list[i]=='SPIM':
                mObj=Spim1D(os.path.join(self.parent_directory,self.model_list[i]))
            d=mObj.parse_results(0,clip_ts[i],False,1e-1)
            
            ts=d['ts']
            x=d['x']
            x_center=d['x_center']
            z0=d['z0']
            dt=d['dt']
            fts=d['freq_to_save']
            A=d['A']
            slp0=d['slp0']
            slp=d['sout']
            zout=d['zout']
            zcout=d['zcout']
            chi=d['chi']
            eCout=d['eCout']
            eout=d['eout']*(-1)*(10*100)*(dt*365)
            qout=d['qout']
            mrout=d['mrout']
            crout=d['crout']
            rlf_Z=d['rlf_Z']
            
            avgE=np.diff(eCout,axis=0)
            avgE=np.concatenate((eCout[0,:].reshape((1,len(eCout[0,:]))),avgE),axis=0)
            avgE=((avgE*(-1))/fts)*(10*100)
            
            # Determine step
            step=int(len(ts)/num_ts)
            
            # Chi-Elevation
            ax1=f1.add_subplot(gs[0,i])
            ax1.set_xlabel(r'$\chi$')
            ax1.set_title(self.descript_list[i])
            if i==0:
                ax1.set_ylabel('Elevation [km]')
            ax1.plot(chi,z0/1000,c='k',linestyle=':')
            for j in range(0,len(ts),step):
                ax1.plot(chi,zout[j,:]/1000,c=cm.batlowK_r(col_vec(ts[j])))
            norm=colors.Normalize(vmin=0,vmax=np.max(max_ts)/1e6)
            cbar1=plt.colorbar(cmm.ScalarMappable(norm=norm,cmap=cm.batlowK_r),ax=ax1)
            if i>0:
                cbar1.ax.set_ylabel('Time [Myrs]')
            ax1.set_ylim((-0.01,max_z))
            
            # Stream profile
            ax2=f1.add_subplot(gs[1,i])
            ax2.set_xlabel('Stream Distance [km]')
            if i==0:
                ax2.set_ylabel('Elevation [km]')
            ax2.plot(x/1000,z0/1000,c='k',linestyle=':')
            for j in range(0,len(ts),step):
                ax2.plot(x/1000,zout[j,:]/1000,c=cm.batlowK_r(col_vec(ts[j])))
            ax2.set_ylim((-0.01,max_z))
            
            # Erosion rate
            ax3=f1.add_subplot(gs[2,i])
            ax3.set_xlabel('Stream Distance [km]')
            if i==0:
                ax3.set_ylabel('Erosion Rate [mm/yr]')
            for j in range(0,len(ts),step):
                ax3.plot(x[1:]/1000,avgE[j,1:],c=cm.batlowK_r(col_vec(ts[j])))
            ax3.set_ylim((-0.1,max_e))
            
            # Discharge    
            ax4=f1.add_subplot(gs[3,i])
            ax4.set_xlabel('Stream Distance [km]')
            if i==0:
                ax4.set_ylabel('Mean Discharge [$m^{3}/s$]') 
            for j in range(0,len(ts),step):
                ax4.plot(x/1000,qout[j,:],c=cm.batlowK_r(col_vec(ts[j])))
            ax4.set_ylim((-0.1,max_q))
                
            ax1.text(0.01, 0.99, letters[letter_order[0,i]],
                    horizontalalignment='left',
                    verticalalignment='top',
                    transform=ax1.transAxes,
                    fontsize=12,fontweight='extra bold')
            ax2.text(0.01, 0.99, letters[letter_order[1,i]],
                    horizontalalignment='left',
                    verticalalignment='top',
                    transform=ax2.transAxes,
                    fontsize=12,fontweight='extra bold')
            ax3.text(0.01, 0.99, letters[letter_order[2,i]],
                    horizontalalignment='left',
                    verticalalignment='top',
                    transform=ax3.transAxes,
                    fontsize=12,fontweight='extra bold')
            ax4.text(0.95, 0.99, letters[letter_order[3,i]],
                    horizontalalignment='left',
                    verticalalignment='top',
                    transform=ax4.transAxes,
                    fontsize=12,fontweight='extra bold')               
            
        plt.tight_layout()
        plt.rcdefaults()
        return f1

    def comp_profile_evol2(self,num_ts,max_z,max_e,max_q,clip_ts):
        # Build ts list
        ts_list=[]
        ts0_list=[]
        ts1_list=[]
        for i in range(self.num_models):
            ts=self.sorted_time_steps(self.model_list[i],self.prefix_list[i])
            ts_list.append(ts)
            ts0_list.append(ts[0])
            ts1_list.append(ts[-1])
        # Determine type
        type_list=self.determine_type(ts0_list)
        max_ts=np.max(clip_ts)
        if max_ts==np.inf:
            max_ts=np.max(np.array(ts1_list))
        
        # Account for zero time
        num_ts=num_ts-1
        
        # Generate figure details
        letters=list(string.ascii_uppercase)
        # Set figure details
        SMALL_SIZE = 8
        MEDIUM_SIZE = 10
        BIGGER_SIZE = 12
        plt.rc('font', size=SMALL_SIZE,family='Futura')          # controls default text sizes
        plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
        plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
        plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
        plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
        plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
        plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
        # Initiate figure
        f1=plt.figure(figsize=(8,9))
        f1.set_dpi(250)
        num_panels=4
        gs=gridspec.GridSpec(num_panels,self.num_models)
        letter_order=np.arange(0,self.num_models*num_panels,1).reshape((num_panels,self.num_models))
        col_vec=colors.Normalize(vmin=0,vmax=max_ts)
        
        for i in range(self.num_models):
            if type_list[i]=='STIM':
                mObj=Stim1D(os.path.join(self.parent_directory,self.model_list[i]))
            elif type_list[i]=='SPIM':
                mObj=Spim1D(os.path.join(self.parent_directory,self.model_list[i]))
            d=mObj.parse_results(0,clip_ts[i],False,1e-1)
            
            ts=d['ts']
            x=d['x']
            x_center=d['x_center']
            z0=d['z0']
            dt=d['dt']
            fts=d['freq_to_save']
            A=d['A']
            slp0=d['slp0']
            slp=d['sout']
            zout=d['zout']
            zcout=d['zcout']
            chi=d['chi']
            eCout=d['eCout']
            eout=d['eout']*(-1)*(10*100)*(dt*365)
            qout=d['qout']
            mrout=d['mrout']
            crout=d['crout']
            rlf_Z=d['rlf_Z']
            
            avgE=np.diff(eCout,axis=0)
            avgE=np.concatenate((eCout[0,:].reshape((1,len(eCout[0,:]))),avgE),axis=0)
            avgE=((avgE*(-1))/fts)*(10*100)
            
            # Determine step
            step=int(len(ts)/num_ts)
            
            # Chi-Elevation
            ax1=f1.add_subplot(gs[0,i])
            ax1.set_xlabel(r'$\chi$')
            ax1.set_title(self.descript_list[i])
            if i==0:
                ax1.set_ylabel('Elevation [km]')
            ax1.plot(chi,z0/1000,c='k',linestyle=':')
            for j in range(0,len(ts),step):
                ax1.plot(chi,zout[j,:]/1000,c=cm.batlowK_r(col_vec(ts[j])))
            ax1.set_ylim((-0.01,max_z))
            
            # Stream profile
            ax2=f1.add_subplot(gs[1,i])
            ax2.set_xlabel('Stream Distance [km]')
            if i==0:
                ax2.set_ylabel('Elevation [km]')
            ax2.plot(x/1000,z0/1000,c='k',linestyle=':')
            for j in range(0,len(ts),step):
                ax2.plot(x/1000,zout[j,:]/1000,c=cm.batlowK_r(col_vec(ts[j])))
            norm=colors.Normalize(vmin=0,vmax=np.max(max_ts)/1e6)
            cbar1=plt.colorbar(cmm.ScalarMappable(norm=norm,cmap=cm.batlowK_r),ax=ax2)
            if i>0:
                cbar1.ax.set_ylabel('Model Time [Myrs]')
            ax2.set_ylim((-0.01,max_z))
            ax2.set_xlim((0,np.max(x/1000)))
            
            # Erosion rate
            ax3=f1.add_subplot(gs[2,i])
            ax3.set_xlabel('Stream Distance [km]')
            if i==0:
                ax3.set_ylabel('Model Time [Myrs]')
            im1=ax3.imshow(np.flipud(avgE),extent=[x[0]/1000,x[-1]/1000,ts[0]/1e6,ts[-1]/1e6],
                       cmap=cm.tokyo_r,norm=colors.Normalize(vmin=0,vmax=max_e),aspect='auto')
            cbar3=plt.colorbar(im1,ax=ax3)
            if i>0:
                cbar3.ax.set_ylabel('Erosion Rate [mm/yr]')
            
           
            # # Discharge    
            # ax4=f1.add_subplot(gs[3,i])
            # ax4.set_xlabel('Stream Distance [km]')
            # if i==0:
            #     ax4.set_ylabel('Model Time [Myrs]')
            # im2=ax4.imshow(np.flipud(qout),extent=[x[0]/1000,x[-1]/1000,ts[0]/1e6,ts[-1]/1e6],
            #            cmap=cm.davos_r,norm=colors.Normalize(vmin=0,vmax=max_q),aspect='auto')
            # cbar4=plt.colorbar(im2,ax=ax4)
            # if i>0:
            #     cbar4.ax.set_ylabel('Discharge [$m^{3}/s$]')

            # Discharge    
            ax4=f1.add_subplot(gs[3,i])
            ax4.set_xlabel('Stream Distance [km]')
            if i==0:
                ax4.set_ylabel('Model Time [Myrs]')
            im2=ax4.imshow(np.flipud(mrout),extent=[x_center[0]/1000,x_center[-1]/1000,ts[0]/1e6,ts[-1]/1e6],
                       cmap=cm.davos_r,norm=colors.Normalize(vmin=0,vmax=max_q),aspect='auto')
            cbar4=plt.colorbar(im2,ax=ax4)
            if i>0:
                cbar4.ax.set_ylabel('Mean Runoff [mm/day]')
            

                
            ax1.text(0.01, 0.99, letters[letter_order[0,i]],
                    horizontalalignment='left',
                    verticalalignment='top',
                    transform=ax1.transAxes,
                    fontsize=12,fontweight='extra bold')
            ax2.text(0.01, 0.99, letters[letter_order[1,i]],
                    horizontalalignment='left',
                    verticalalignment='top',
                    transform=ax2.transAxes,
                    fontsize=12,fontweight='extra bold')
            ax3.text(0.01, 0.99, letters[letter_order[2,i]],
                    horizontalalignment='left',
                    verticalalignment='top',
                    transform=ax3.transAxes,
                    fontsize=12,fontweight='extra bold')
            ax4.text(0.01, 0.99, letters[letter_order[3,i]],
                    horizontalalignment='left',
                    verticalalignment='top',
                    transform=ax4.transAxes,
                    fontsize=12,fontweight='extra bold')               
            
        plt.tight_layout()
        plt.rcdefaults()
        return f1
        
        
    def time_to_ss(self,grp_num,line_style):
        # Build ts list
        yr_list=[]
        for i in range(self.num_models):
            ts=self.sorted_time_steps(self.model_list[i],self.prefix_list[i])
            yr_list.append(ts[-1])
        type_list=self.determine_type(yr_list)
        
        # Set figure details
        SMALL_SIZE = 8
        MEDIUM_SIZE = 10
        BIGGER_SIZE = 12
        plt.rc('font', size=SMALL_SIZE,family='Futura')          # controls default text sizes
        plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
        plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
        plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
        plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
        plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
        plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
        # Initiate figure
        f1=plt.figure(figsize=(11,8))
        f1.set_dpi(250)
        # Set color normalization
        col_vec=colors.Normalize(vmin=0.25,vmax=8)
        
        gs=gridspec.GridSpec(3,3)
        
        ax1a=f1.add_subplot(gs[0,0])
        ax1a.set_ylabel(r'$\Delta$ Maximum Elevation [m]')
        ax1a.set_yscale('log')

        ax1b=f1.add_subplot(gs[0,1])
        ax1b.set_yscale('log')

        ax1c=f1.add_subplot(gs[0,2])
        ax1c.set_yscale('log') 

        ax2a=f1.add_subplot(gs[1,0])
        ax2a.set_ylabel(r'$\Delta$ Mean Elevation [m]')
        ax2a.set_yscale('log')

        ax2b=f1.add_subplot(gs[1,1])
        ax2b.set_yscale('log')

        ax2c=f1.add_subplot(gs[1,2])
        ax2c.set_yscale('log')

        ax3a=f1.add_subplot(gs[2,0])
        ax3a.set_ylabel(r'Max $\Delta$ Elevation [m]')
        ax3a.set_xlabel('Model Time [Myr]')
        ax3a.set_yscale('log')

        ax3b=f1.add_subplot(gs[2,1])
        ax3b.set_xlabel('Model Time [Myr]')
        ax3b.set_yscale('log')

        ax3c=f1.add_subplot(gs[2,2])
        ax3c.set_xlabel('Model Time [Myr]')
        ax3c.set_yscale('log')                 
        
        for i in range(self.num_models):
            if type_list[i]=='STIM':
                mObj=Stim1D(os.path.join(self.parent_directory,self.model_list[i]))
            elif type_list[i]=='SPIM':
                mObj=Spim1D(os.path.join(self.parent_directory,self.model_list[i]))
            d=mObj.parse_results(0,np.inf,False,1e-1)
            # Extract from dict
            ts=d['ts']/1e6
            zout=d['zout']
            z0=d['z0']
            U=d['uplift']*1000
            # Setup containers
            delta_maxZ=np.zeros(len(ts))
            delta_mnZ=np.zeros(len(ts))
            max_deltaZ=np.zeros(len(ts))
            # Calcuate SS
            for j in range(len(ts)):
                if j==0:
                    zf=zout[j,:]
                    zi=z0
                else:
                    zf=zout[j,:]
                    zi=zout[j-1,:]
                
                delta_maxZ[j]=np.abs(np.max(zf)-np.max(zi))
                delta_mnZ[j]=np.abs(np.mean(zf)-np.mean(zi))
                max_deltaZ[j]=np.max(np.abs(zf-zi))
            if grp_num[i]==0:
                ax1a.plot(ts,delta_maxZ,color=cm.vik_r(col_vec(U)),linestyle=line_style[i])
                ax2a.plot(ts,delta_mnZ,color=cm.vik_r(col_vec(U)),linestyle=line_style[i])
                ax3a.plot(ts,max_deltaZ,color=cm.vik_r(col_vec(U)),linestyle=line_style[i])
            elif grp_num[i]==1:
                ax1b.plot(ts,delta_maxZ,color=cm.vik_r(col_vec(U)),linestyle=line_style[i])
                ax2b.plot(ts,delta_mnZ,color=cm.vik_r(col_vec(U)),linestyle=line_style[i])
                ax3b.plot(ts,max_deltaZ,color=cm.vik_r(col_vec(U)),linestyle=line_style[i])
            elif grp_num[i]==2:
                ax1c.plot(ts,delta_maxZ,color=cm.vik_r(col_vec(U)),linestyle=line_style[i])
                ax2c.plot(ts,delta_mnZ,color=cm.vik_r(col_vec(U)),linestyle=line_style[i])
                ax3c.plot(ts,max_deltaZ,color=cm.vik_r(col_vec(U)),linestyle=line_style[i])
                
        plt.tight_layout()
        plt.rcdefaults()
    
    def stability_values(self):
        # Build ts list
        yr_list=[]
        for i in range(self.num_models):
            ts=self.sorted_time_steps(self.model_list[i],self.prefix_list[i])
            yr_list.append(ts[-1])
              
        type_list=self.determine_type(yr_list)

        stabil=np.zeros(self.num_models)

        for i in range(self.num_models):
            yr=int(yr_list[i])
            ts=self.sorted_time_steps(self.model_list[i],self.prefix_list[i])
            ix=np.argmin(np.abs(ts-yr))
            
            fname1=os.path.join(self.parent_directory,self.model_list[i],self.prefix_list[i]+str(ts[ix])+'.pkl')                
            if type_list[i]=='STIM':
                [cObj1,sObj1,eObj1,rObj1]=self.recover_stim_state(fname1)
            elif type_list[i]=='SPIM':
                [cObj1,spObj1,sObj1,eObj1,rObj1]=self.recover_spim_state(fname1)

            stabil[i]=eObj1.max_cfl

        return stabil


    def ksn_final_ts(self,rp_slp,rp_yint):
        # Build ts list
        yr_list=[]
        for i in range(self.num_models):
            ts=self.sorted_time_steps(self.model_list[i],self.prefix_list[i])
            yr_list.append(ts[-1])
              
        type_list=self.determine_type(yr_list)
        
        mn_ksn=np.zeros(self.num_models)
        std_ksn=np.zeros(self.num_models)
        stde_ksn=np.zeros(self.num_models)
        mn_ksnqr=np.zeros(self.num_models)
        std_ksnqr=np.zeros(self.num_models)
        stde_ksnqr=np.zeros(self.num_models)
        mn_ksnqp=np.zeros(self.num_models)
        std_ksnqp=np.zeros(self.num_models)
        stde_ksnqp=np.zeros(self.num_models)
        mn_E=np.zeros(self.num_models)
        std_E=np.zeros(self.num_models)
        stde_E=np.zeros(self.num_models)
        p25_E=np.zeros(self.num_models)
        p75_E=np.zeros(self.num_models)
        
        for i in range(self.num_models):
            yr=int(yr_list[i])
            ts=self.sorted_time_steps(self.model_list[i],self.prefix_list[i])
            ix=np.argmin(np.abs(ts-yr))
            
            fname1=os.path.join(self.parent_directory,self.model_list[i],self.prefix_list[i]+str(ts[ix])+'.pkl')
            fname0=os.path.join(self.parent_directory,self.model_list[i],self.prefix_list[i]+str(ts[ix-1])+'.pkl')
            if type_list[i]=='STIM':
                [cObj0,sObj0,eObj0,rObj0]=self.recover_stim_state(fname0)
                [cObj1,sObj1,eObj1,rObj1]=self.recover_stim_state(fname1)
            elif type_list[i]=='SPIM':
                [cObj0,spObj0,sObj0,eObj0,rObj0]=self.recover_spim_state(fname0)
                [cObj1,spObj1,sObj1,eObj1,rObj1]=self.recover_spim_state(fname1)
                
            # Calculate ksn-qr along profile
            AR=rObj1.qbar * (60*60*24*365.25)
            chir=integrate.cumtrapz((1/AR)**sObj1.theta,sObj1.x,initial=0)
            ksnqr=np.diff(sObj1.z)/np.diff(chir)
            
            # Calculate ksn-qp along profile
            # Convert binned runoff to precip equivalent            
            p_bin=rp_slp[i]*rObj1.r_bin + rp_yint[i]
            # Propagate binned precip
            p=np.zeros(sObj1.z.shape)
            for j in range(len(sObj1.uix)):
                p[sObj1.ix==sObj1.uix[j]]=p_bin[j]
            # Route binned precip
            pq=rObj1.route_discharge(sObj1,p)
            AP=pq  * (60*60*24*365.25)
            # Calcuate chip and ksn-qp
            chip=integrate.cumtrapz((1/AP)**sObj1.theta,sObj1.x,initial=0)
            ksnqp=np.diff(sObj1.z)/np.diff(chip)
            
            # Calculate ksn along profile
            ksn=np.diff(sObj1.z)/np.diff(sObj1.chi)
            # Calculate erosion rate along profile         
            fts=cObj0.freq_to_save
            # Erate=(((eObj1.cum_E - eObj0.cum_E)*-1)/fts)*(10*100)*1000
        
            # Modified erosion calculation
            ts_slice=ts[-10:]
            er=np.zeros((len(ts_slice),len(sObj1.x)))
            for j in range(len(ts_slice)):
                fname0=os.path.join(self.parent_directory,self.model_list[i],self.prefix_list[i]+str(ts_slice[j])+'.pkl')
                [cObj0,sObj0,eObj0,rObj0]=self.recover_stim_state(fname0)
                er[j,:]=eObj0.cum_E
                
            Erate=np.diff(er,axis=0)
            Erate=((Erate*(-1))/fts)*(10*100)*1000
            
            # Calculate averages and std
            mn_ksn[i]=np.mean(ksn)
            std_ksn[i]=np.std(ksn)
            stde_ksn[i]=np.std(ksn)/np.sqrt(ksn.size)
            mn_ksnqr[i]=np.mean(ksnqr)
            std_ksnqr[i]=np.std(ksnqr)
            stde_ksnqr[i]=np.std(ksnqr)/np.sqrt(ksnqr.size)
            mn_ksnqp[i]=np.mean(ksnqp)
            std_ksnqp[i]=np.std(ksnqp)
            stde_ksnqp[i]=np.std(ksnqp)/np.sqrt(ksnqp.size)
            mn_E[i]=np.mean(Erate)
            std_E[i]=np.std(Erate)
            stde_E[i]=np.std(Erate)/np.sqrt(Erate.size)
            p25_E[i]=np.percentile(Erate,25)
            p75_E[i]=np.percentile(Erate,75)


        return mn_ksn,std_ksn,stde_ksn,mn_ksnqr,std_ksnqr,stde_ksnqr,mn_ksnqp,std_ksnqp,stde_ksnqp,mn_E,std_E,stde_E,p25_E,p75_E
    
    def rpv_final_ts(self,rp_slp,rp_yint):
        
        # Build ts list
        yr_list=[]
        for i in range(self.num_models):
            ts=self.sorted_time_steps(self.model_list[i],self.prefix_list[i])
            yr_list.append(ts[-1])
              
        type_list=self.determine_type(yr_list)
        
        mnR=np.zeros(self.num_models)
        mnP=np.zeros(self.num_models)
        cr=np.zeros(self.num_models)
        sr=np.zeros(self.num_models)
        
        for i in range(self.num_models):
            yr=int(yr_list[i])
            ts=self.sorted_time_steps(self.model_list[i],self.prefix_list[i])
            ix=np.argmin(np.abs(ts-yr))
            
            fname1=os.path.join(self.parent_directory,self.model_list[i],self.prefix_list[i]+str(ts[ix])+'.pkl')
            if type_list[i]=='STIM':
                [cObj1,sObj1,eObj1,rObj1]=self.recover_stim_state(fname1)
            elif type_list[i]=='SPIM':
                [cObj1,spObj1,sObj1,eObj1,rObj1]=self.recover_spim_state(fname1)
                
            Q_m3_day=rObj1.qbar * (60*60*24)
            mnR[i]=(Q_m3_day[0]/sObj1.A[0])*(100*10)
            mnP[i]=rp_slp[i]*mnR[i] + rp_yint[i]
            
            print('Performing monte-carlo for '+str(i+1)+' of '+str(self.num_models))
            mcd=rObj1.monte_carlo_runoff(sObj1,500,1,verbose=False)
            cr[i]=mcd['tail_fit_cr']
            sr[i]=mcd['tail_fit_cr']
            
        return mnR,mnP,cr,sr
    
    def theta_final_ts(self,rp_slp,rp_yint):
        # Build ts list
        yr_list=[]
        for i in range(self.num_models):
            ts=self.sorted_time_steps(self.model_list[i],self.prefix_list[i])
            yr_list.append(ts[-1])
              
        type_list=self.determine_type(yr_list)
        
        theta_chi=np.zeros(self.num_models)
        theta_slp=np.zeros(self.num_models)
        theta_chir=np.zeros(self.num_models)
        theta_chip=np.zeros(self.num_models)
        
        for i in range(self.num_models):
            yr=int(yr_list[i])
            ts=self.sorted_time_steps(self.model_list[i],self.prefix_list[i])
            ix=np.argmin(np.abs(ts-yr))
            
            fname1=os.path.join(self.parent_directory,self.model_list[i],self.prefix_list[i]+str(ts[ix])+'.pkl')
            if type_list[i]=='STIM':
                [cObj1,sObj1,eObj1,rObj1]=self.recover_stim_state(fname1)
            elif type_list[i]=='SPIM':
                [cObj1,spObj1,sObj1,eObj1,rObj1]=self.recover_spim_state(fname1)
                
            # sObj1.A, sObj1.x, sObj1.z, sObj1.slope
            # Find theta through optimization
            args=(sObj1.A,sObj1.x,sObj1.z)
            res=minimize_scalar(theta_min,args=args,bounds=(0.1,1.5),method='bounded')
            theta_chi[i]=res.x
            
            # Find theta through slope area fit
            idx=(sObj1.A==0) | (sObj1.slope==0) 
            lres=linregress(np.log10(sObj1.A[~idx]),np.log10(sObj1.slope[~idx]))
            theta_slp[i]=-1*lres.slope
            
            # Calculate runoff weighted area for use in optimization
            AR=rObj1.qbar * (60*60*24*365.25)
            # Calculate precip weigthed area for use in optimization
            # Convert binned runoff to precip equivalent            
            p_bin=rp_slp[i]*rObj1.r_bin + rp_yint[i]
            # Propagate binned precip
            p=np.zeros(sObj1.z.shape)
            for j in range(len(sObj1.uix)):
                p[sObj1.ix==sObj1.uix[j]]=p_bin[j]
            # Route binned precip
            pq=rObj1.route_discharge(sObj1,p)
            AP=pq  * (60*60*24*365.25)
            
            # Perform optimizations
            args=(AR,sObj1.x,sObj1.z)
            res=minimize_scalar(theta_min,args=args,bounds=(0.1,1.5),method='bounded')
            theta_chir[i]=res.x  
            
            args=(AP,sObj1.x,sObj1.z)
            res=minimize_scalar(theta_min,args=args,bounds=(0.1,1.5),method='bounded')
            theta_chip[i]=res.x 

        return theta_chi,theta_slp,theta_chir,theta_chip

    def sf_final_ts(self):
        # Build ts list
        yr_list=[]
        for i in range(self.num_models):
            ts=self.sorted_time_steps(self.model_list[i],self.prefix_list[i])
            yr_list.append(ts[-1])
              
        type_list=self.determine_type(yr_list)
        
        mn_snp=np.zeros(self.num_models)
        min_snp=np.zeros(self.num_models)
        max_snp=np.zeros(self.num_models)
        std_snp=np.zeros(self.num_models)
        
        for i in range(self.num_models):
            yr=int(yr_list[i])
            ts=self.sorted_time_steps(self.model_list[i],self.prefix_list[i])
            ix=np.argmin(np.abs(ts-yr))
            
            fname1=os.path.join(self.parent_directory,self.model_list[i],self.prefix_list[i]+str(ts[ix])+'.pkl')
            if type_list[i]=='STIM':
                [cObj1,sObj1,eObj1,rObj1]=self.recover_stim_state(fname1)
            elif type_list[i]=='SPIM':
                [cObj1,spObj1,sObj1,eObj1,rObj1]=self.recover_spim_state(fname1)
                
            try:
                snp=rObj1.snp_bin
                mn_snp[i]=np.mean(snp)
                min_snp[i]=np.mean(snp)
                max_snp[i]=np.max(snp)
                std_snp[i]=np.std(snp)
            except:
                print('No snow fraction data')
                mn_snp[i]=np.nan
                min_snp[i]=np.nan
                max_snp[i]=np.nan
                std_snp[i]=np.nan
                
        return mn_snp,std_snp,min_snp,max_snp
            

            
    def comp_final_ts(self,group_list,group_col,group_shape,group_filled,group_line,grp,rp_slp,rp_yint):
        # Build ts list
        yr_list=[]
        for i in range(self.num_models):
            ts=self.sorted_time_steps(self.model_list[i],self.prefix_list[i])
            yr_list.append(ts[-1])
              
        type_list=self.determine_type(yr_list)
        
        mn_ksn=np.zeros(self.num_models)
        std_ksn=np.zeros(self.num_models)
        mn_ksnqr=np.zeros(self.num_models)
        std_ksnqr=np.zeros(self.num_models)
        mn_ksnqp=np.zeros(self.num_models)
        std_ksnqp=np.zeros(self.num_models)
        mn_E=np.zeros(self.num_models)
        std_E=np.zeros(self.num_models)
        
        for i in range(self.num_models):
            yr=int(yr_list[i])
            ts=self.sorted_time_steps(self.model_list[i],self.prefix_list[i])
            ix=np.argmin(np.abs(ts-yr))
            
            fname1=os.path.join(self.parent_directory,self.model_list[i],self.prefix_list[i]+str(ts[ix])+'.pkl')
            fname0=os.path.join(self.parent_directory,self.model_list[i],self.prefix_list[i]+str(ts[ix-1])+'.pkl')
            if type_list[i]=='STIM':
                [cObj0,sObj0,eObj0,rObj0]=self.recover_stim_state(fname0)
                [cObj1,sObj1,eObj1,rObj1]=self.recover_stim_state(fname1)
            elif type_list[i]=='SPIM':
                [cObj0,spObj0,sObj0,eObj0,rObj0]=self.recover_spim_state(fname0)
                [cObj1,spObj1,sObj1,eObj1,rObj1]=self.recover_spim_state(fname1)
                
            # Calculate ksn-qr along profile
            AR=rObj1.qbar * (60*60*24*365.25)
            chir=integrate.cumtrapz((1/AR)**sObj1.theta,sObj1.x,initial=0)
            ksnqr=np.diff(sObj1.z)/np.diff(chir)
            
            # Calculate ksn-qp along profile
            # Convert binned runoff to precip equivalent            
            p_bin=rp_slp[i]*rObj1.r_bin + rp_yint[i]
            # Propagate binned precip
            p=np.zeros(sObj1.z.shape)
            for j in range(len(sObj1.uix)):
                p[sObj1.ix==sObj1.uix[j]]=p_bin[j]
            # Route binned precip
            pq=rObj1.route_discharge(sObj1,p)
            AP=pq  * (60*60*24*365.25)
            # Calcuate chip and ksn-qp
            chip=integrate.cumtrapz((1/AP)**sObj1.theta,sObj1.x,initial=0)
            ksnqp=np.diff(sObj1.z)/np.diff(chip)
            
            # Calculate ksn along profile
            ksn=np.diff(sObj1.z)/np.diff(sObj1.chi)
            # Calculate erosion rate along profile         
            fts=cObj0.freq_to_save
            Erate=(((eObj1.cum_E - eObj0.cum_E)*-1)/fts)*(10*100)*1000
            
            # Calculate averages and std
            mn_ksn[i]=np.mean(ksn)
            std_ksn[i]=np.std(ksn)
            mn_ksnqr[i]=np.mean(ksnqr)
            std_ksnqr[i]=np.std(ksnqr)
            mn_ksnqp[i]=np.mean(ksnqp)
            std_ksnqp[i]=np.std(ksnqp)
            mn_E[i]=np.mean(Erate)
            std_E[i]=np.std(Erate)
            
        # Determine groups
        grp=np.array(grp)
        group_list=np.array(group_list)
        num_grp=len(grp)
        
        # Set figure details
        SMALL_SIZE = 8
        MEDIUM_SIZE = 10
        BIGGER_SIZE = 12
        plt.rc('font', size=SMALL_SIZE,family='Futura')          # controls default text sizes
        plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
        plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
        plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
        plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
        plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
        plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
        # Initiate figure
        f1=plt.figure(figsize=(8,4))
        f1.set_dpi(250)
        ax1=f1.add_subplot(1,2,1)
        ax1.set_xlabel('Erosion Rate [m/Myr]')
        ax1.set_ylabel(r'$k_{sn}$ [m]')

        ax2=f1.add_subplot(1,2,2)
        ax2.set_xlabel('Erosion Rate [m/Myr]')
        ax2.set_ylabel(r'$k_{snQ}$ [m]')        
        
        
        f2=plt.figure(figsize=(8,6))
        f2.set_dpi(250)

        ax3=f2.add_subplot(1,3,1)
        ax3.set_xlabel('Erosion Rate [m/Myr]')
        ax3.set_ylabel(r'$k_{snQP}$ [m]')
        
        ax4=f2.add_subplot(1,3,2)
        ax4.set_xlabel('Erosion Rate [m/Myr]')
        ax4.set_ylabel(r'$k_{snQR}$ [m]')
        
        ax5=f2.add_subplot(1,3,3)
        ax5.set_xlabel('Erosion Rate [m/Myr]')
        ax5.set_ylabel(r'$k_{snQP}$ / $k_{snQR}$')        
        
        e_vec=np.logspace(2,4,100)
        
        c_out=np.zeros(num_grp)
        phi_out=np.zeros(num_grp)
        K_out=np.zeros(num_grp)
        n_out=np.zeros(num_grp)
        
        for i in range(num_grp):
            idx=group_list==grp[i]
            
            # Convert e to m/yr
            ef=mn_E[idx]/1e6
            
            [K,n]=odr_fit(mn_ksn[idx],ef)
            K_out[i]=K
            n_out[i]=n
            
            [c,phi]=odr_fit(ef,mn_ksn[idx])
            c_out[i]=c
            phi_out[i]=phi
            if group_filled[i]=='filled':
                ax1.scatter(mn_E[idx],mn_ksn[idx],color=group_col[i],marker=group_shape[i],s=50,zorder=2,label=grp[i]+'; n = '+str(np.round(1/phi,1)) )
            elif group_filled[i]=='open':
                ax1.scatter(mn_E[idx],mn_ksn[idx],color='w',edgecolor=group_col[i],marker=group_shape[i],s=50,zorder=2,label=grp[i]+'; n = '+str(np.round(1/phi,1)) )
            ax1.errorbar(mn_E[idx],mn_ksn[idx],std_ksn[idx],std_E[idx],linestyle='',ecolor=group_col[i],zorder=1,elinewidth=0.5)
            ax1.plot(e_vec,c*(e_vec/1e6)**phi,c=group_col[i],linestyle=group_line[i],zorder=0)

            [cp,phip]=odr_fit(ef,mn_ksnqp[idx])
            if group_filled[i]=='filled':
                ax2.scatter(mn_E[idx],mn_ksnqp[idx],color=group_col[i],marker=group_shape[i],s=50,zorder=2,label=grp[i]+'; n = '+str(np.round(1/phip,1)) )
            elif group_filled[i]=='open':
                ax2.scatter(mn_E[idx],mn_ksnqp[idx],color='w',edgecolor=group_col[i],marker=group_shape[i],s=50,zorder=2,label=grp[i]+'; n = '+str(np.round(1/phip,1)) )
            ax2.errorbar(mn_E[idx],mn_ksnqp[idx],std_ksnqp[idx],std_E[idx],linestyle='',ecolor=group_col[i],zorder=1,elinewidth=0.5)
            ax2.plot(e_vec,cp*(e_vec/1e6)**phip,c=group_col[i],linestyle=group_line[i],zorder=0) 
            
            ##
            [cr,phir]=odr_fit(ef,mn_ksnqr[idx])
            if group_filled[i]=='filled':
                ax4.scatter(mn_E[idx],mn_ksnqr[idx],color=group_col[i],marker=group_shape[i],s=50,zorder=2,label=grp[i]+'; n = '+str(np.round(1/phir,1)) )
            elif group_filled[i]=='open':
                ax4.scatter(mn_E[idx],mn_ksnqr[idx],color='w',edgecolor=group_col[i],marker=group_shape[i],s=50,zorder=2,label=grp[i]+'; n = '+str(np.round(1/phir,1)) )
            ax4.errorbar(mn_E[idx],mn_ksnqr[idx],std_ksnqr[idx],std_E[idx],linestyle='',ecolor=group_col[i],zorder=1,elinewidth=0.5)
            ax4.plot(e_vec,cr*(e_vec/1e6)**phir,c=group_col[i],linestyle=group_line[i],zorder=0) 
            
            [cp,phip]=odr_fit(ef,mn_ksnqp[idx])
            if group_filled[i]=='filled':
                ax3.scatter(mn_E[idx],mn_ksnqp[idx],color=group_col[i],marker=group_shape[i],s=50,zorder=2,label=grp[i]+'; n = '+str(np.round(1/phip,1)) )
            elif group_filled[i]=='open':
                ax3.scatter(mn_E[idx],mn_ksnqp[idx],color='w',edgecolor=group_col[i],marker=group_shape[i],s=50,zorder=2,label=grp[i]+'; n = '+str(np.round(1/phip,1)) )                
            ax3.errorbar(mn_E[idx],mn_ksnqp[idx],std_ksnqp[idx],std_E[idx],linestyle='',ecolor=group_col[i],zorder=1,elinewidth=0.5)
            ax3.plot(e_vec,cp*(e_vec/1e6)**phip,c=group_col[i],linestyle=group_line[i],zorder=0)
            
            if group_filled[i]=='filled':
                ax5.scatter(mn_E[idx],mn_ksnqp[idx]/mn_ksnqr[idx],color=group_col[i],marker=group_shape[i],s=50,zorder=2)
            elif group_filled[i]=='open':
                ax5.scatter(mn_E[idx],mn_ksnqp[idx]/mn_ksnqr[idx],color='w',edgecolor=group_col[i],marker=group_shape[i],s=50,zorder=2)
            
            
        # ax1.legend(bbox_to_anchor= (-0.1,-0.75),loc='lower left')
        ax1.legend(loc='best')
        ax1.text(0.01, 0.99, 'A',
                horizontalalignment='left',
                verticalalignment='top',
                transform=ax1.transAxes,
                fontsize=12,fontweight='extra bold')   
        # ax2.legend(bbox_to_anchor= (-0.1,-0.75),loc='lower left')
        ax2.legend(loc='best')
        ax2.text(0.01, 0.99, 'B',
                horizontalalignment='left',
                verticalalignment='top',
                transform=ax2.transAxes,
                fontsize=12,fontweight='extra bold')
        
        ax3.legend(bbox_to_anchor= (-0.1,-0.75),loc='lower left')
        ax3.text(0.01, 0.99, 'A',
                horizontalalignment='left',
                verticalalignment='top',
                transform=ax3.transAxes,
                fontsize=12,fontweight='extra bold')
        ax4.legend(bbox_to_anchor= (-0.1,-0.75),loc='lower left')
        ax4.text(0.01, 0.99, 'B',
                horizontalalignment='left',
                verticalalignment='top',
                transform=ax4.transAxes,
                fontsize=12,fontweight='extra bold')
        ax5.text(0.92, 0.99, 'C',
                horizontalalignment='left',
                verticalalignment='top',
                transform=ax5.transAxes,
                fontsize=12,fontweight='extra bold')
        plt.tight_layout()

        
        df=pd.DataFrame({'Group':grp,
                          'C':c_out,
                          'phi':phi_out,
                          'K':K_out,
                          'n':n_out})
        
        plt.rcdefaults()
        
        return df,f1,f2
        
    
    def plot_individual_ts(self,yr_list):
        type_list=self.determine_type(yr_list)
        slist1=[]
        elist1=[]
        rlist1=[]
        slist0=[]
        elist0=[]
        rlist0=[]
        for i in range(self.num_models):
            yr=int(yr_list[i])
            ts=self.sorted_time_steps(self.model_list[i],self.prefix_list[i])
            ix=np.argmin(np.abs(ts-yr))
            
            fname1=os.path.join(self.parent_directory,self.model_list[i],self.prefix_list[i]+str(ts[ix])+'.pkl')
            fname0=os.path.join(self.parent_directory,self.model_list[i],self.prefix_list[i]+str(ts[ix-1])+'.pkl')
            if type_list[i]=='STIM':
                [cObj0,sObj0,eObj0,rObj0]=self.recover_stim_state(fname0)
                [cObj1,sObj1,eObj1,rObj1]=self.recover_stim_state(fname1)
            elif type_list[i]=='SPIM':
                [cObj0,spObj0,sObj0,eObj0,rObj0]=self.recover_spim_state(fname0)
                [cObj1,spObj1,sObj1,eObj1,rObj1]=self.recover_spim_state(fname1)
            slist0.append(sObj0)
            elist0.append(eObj0)
            rlist0.append(rObj0)            
            slist1.append(sObj1)
            elist1.append(eObj1)
            rlist1.append(rObj1)
            
        # Begin plotting
        col_vec=colors.Normalize(vmin=0,vmax=self.num_models-1)
        
        f1=plt.figure(figsize=(25,20))
        ax1=plt.subplot(4,2,1)
        ax1.set_xlabel('Stream Distance [km]')
        ax1.set_ylabel('Elevation [km]')      
        
        ax2=plt.subplot(4,2,2)
        ax2.set_xlabel('Stream Distance [km]')
        ax2.set_ylabel('Mean Discharge [$m^{3}/s$]')
    
        ax3=plt.subplot(4,2,3)
        ax3.set_xlabel('Stream Distance [km]')
        ax3.set_ylabel('$k_{sn}$ [m]') 
        
        ax4=plt.subplot(4,2,4)
        ax4.set_xlabel(r'$\chi$')
        ax4.set_ylabel('Elevation [km]')
        
        ax5=plt.subplot(4,2,5)
        ax5.set_xlabel('Stream Distance [km]')
        ax5.set_ylabel('Erosion Rate [mm/yr]')
        
        ax6=plt.subplot(4,2,6)
        ax6.set_xscale('log')
        ax6.set_yscale('log')
        ax6.set_xlabel(r'Drainage Area [$m^{2}$]')
        ax6.set_ylabel('Slope [m/m]')
        
        ax7=plt.subplot(4,2,7)
        ax7.set_xlabel('Mean Runoff [mm/day]')
        ax7.set_ylabel('Variablity')

        ax8=plt.subplot(4,2,8)
        ax8.set_xlabel('Erosion Rate [mm/yr]')
        ax8.set_ylabel(r'$k_{sn}$ [m]')

        for i in range(self.num_models):
            sObj1=slist1[i]
            eObj1=elist1[i]
            rObj1=rlist1[i]
            sObj0=slist0[i]
            eObj0=elist0[i]
            rObj0=rlist0[i]
            ax1.plot(sObj1.x/1000,sObj1.z/1000,c=cm.roma(col_vec(i)))
            ax2.plot(sObj1.x/1000,rObj1.qbar,c=cm.roma(col_vec(i)))
            # Calculate ksn
            ksn=np.diff(sObj1.z)/np.diff(sObj1.chi)
            ax3.plot(sObj1.x[1:len(sObj1.x)]/1000,ksn,c=cm.roma(col_vec(i)))
            ax4.plot(sObj1.chi,sObj1.z/1000,c=cm.vik_r(col_vec(i)))
            # Calculate erosion rate         
            fts=cObj0.freq_to_save
            Erate=(((eObj1.cum_E - eObj0.cum_E)*-1)/fts)*(10*100)
            ax5.plot(sObj1.x,Erate,c=cm.roma(col_vec(i)))
            ax6.plot(sObj1.A[1:len(sObj1.A)],sObj1.slope[1:len(sObj1.A)],c=cm.roma(col_vec(i)),label=self.descript_list[i])
            ax7.scatter(rObj1.r_bin,rObj1.cr_bin,color=cm.roma(col_vec(i)))
            # Calculate averages
            U=np.mean(eObj1.uplift)*(10*100)
            E=np.mean(Erate)
            mnksn=np.mean(ksn)
            if i==0:
                ax8.scatter(E,mnksn,color=cm.roma(col_vec(i)),marker='o',label='Mean Erosion Rate')
                ax8.scatter(U,mnksn,color=cm.roma(col_vec(i)),marker='s',label='Mean Uplift Rate')
            else:
                ax8.scatter(E,mnksn,color=cm.roma(col_vec(i)),marker='o')
                ax8.scatter(U,mnksn,color=cm.roma(col_vec(i)),marker='s')                
        ax6.legend(loc='best')
        ax8.legend(loc='best')

    def sorted_time_steps(self,model,prefix):
        # Generate unsorted file list
        p_files=glob.glob(os.path.join(self.parent_directory,model,prefix+'*.pkl'))
        # Extract numbers
        p_nums=np.zeros(len(p_files))
        regex=re.compile(r'\d+')
        for i in range(len(p_files)):
            p_nums[i]=[int(j) for j in regex.findall(p_files[i])][-1]
        # Convert file list to array to ease indexing
        p_files=np.array([p_files])
        p_files=p_files.reshape((p_files.shape[1]))
        # Trim timesteps if necessary
        p_idx=(p_nums>=self.lower_ts) & (p_nums<=self.upper_ts)
        p_nums=p_nums[p_idx]
        # Created sorted index
        pix=np.argsort(p_nums)
        # Apply index to list and generate list of model times
        # p_files_sorted=[p_files[i] for i in pix]
        ts=p_nums[pix].astype(int)
        return ts
    
    def build_stim_dict(self,j,ts,sample_freq):
        fname=os.path.join(self.parent_directory,self.model_list[j],self.prefix_list[j]+str(ts[0])+'.pkl')
        # Load static values
        [cObj0,sObj0,eObj0,rObj0]=self.recover_stim_state(fname)
        x=sObj0.x
        z0=sObj0.z0
        chi=sObj0.chi
        A=sObj0.A
        slp0=sObj0.slope0
        x_center=sObj0.x_cents
        runoff_pattern=rObj0.pattern
        U=np.mean(eObj0.uplift)
        
        # Generate empty array for storage
        zout=np.zeros((len(ts),len(z0)))
        eout=np.zeros((len(ts),len(z0)))
        eCout=np.zeros((len(ts),len(z0)))
        qout=np.zeros((len(ts),len(z0)))
        sout=np.zeros((len(ts),len(z0)))
        dout=np.zeros((len(ts),len(z0)))
        zcout=np.zeros((len(ts),len(x_center)))
        mrout=np.zeros((len(ts),len(x_center)))
        crout=np.zeros((len(ts),len(x_center)))
        srout=np.zeros((len(ts),len(x_center)))
        cr_wavg=np.zeros(len(ts))
        cr_fit=np.zeros(len(ts))

        if runoff_pattern=='emp':
            snp=np.zeros((len(ts),len(x_center)))
            max_Z=np.zeros((len(ts),len(x_center)))
            rlf_Z=np.zeros((len(ts),len(x_center)))
            
        if sample_freq>0:
            print('Starting calculation for '+self.descript_list[j])

        # Read in each
        for i in range(len(ts)):
            fname=os.path.join(self.parent_directory,self.model_list[j],self.prefix_list[j]+str(ts[i])+'.pkl')
            [_,sObjOI,eObjOI,rObjOI]=self.recover_stim_state(fname)
            zout[i,:]=sObjOI.z
            sout[i,:]=sObjOI.slope
            eout[i,:]=eObjOI.I
            eCout[i,:]=eObjOI.cum_E
            dout[i,:]=eObjOI.days_E
            qout[i,:]=rObjOI.qbar
            zcout[i,:]=sObjOI.z_cents
            mrout[i,:]=rObjOI.r_bin
            crout[i,:]=rObjOI.cr_bin
            srout[i,:]=rObjOI.sr_bin
            if sample_freq==0:
                cr_wavg[i]=np.mean(rObjOI.cr_bin) 
                cr_fit[i]=np.nan
            else:
                if np.mod(i,sample_freq)==0:
                    print(f'Routing runoff for timestep {i} of {len(ts)}')
                    mcd=rObjOI.monte_carlo_runoff(sObjOI,500,1,verbose=False)
                    cr_wavg[i]=mcd['aw_cr']
                    cr_fit[i]=mcd['tail_fit_cr']
                else:
                    cr_wavg[i]=np.nan
                    cr_fit[i]=np.nan
            
            if runoff_pattern=='emp':
                snp[i,:]=rObjOI.snp_bin
                max_Z[i,:]=rObjOI.max_Z_bin
                rlf_Z[i,:]=rObjOI.rlf_Z_bin
                
        # Package in dictionary for output
        if runoff_pattern=='emp':
            out_dict={'runoff_pattern':runoff_pattern,
                      'x':x,
                      'x_center':x_center,
                      'dt':cObj0.dt,
                      'freq_to_save':cObj0.freq_to_save,
                      'A':A,
                      'chi':chi,
                      'z0':z0,
                      'slp0':slp0,
                      'ts':ts,
                      'U':U,
                      'zout':zout,
                      'sout':sout,
                      'eout':eout,
                      'eCout':eCout,
                      'dout':dout,
                      'qout':qout,
                      'zcout':zcout,
                      'mrout':mrout,
                      'crout':crout,
                      'awcr':cr_wavg,
                      'fitcr':cr_fit,
                      'srout':srout,
                      'snp':snp,
                      'max_Z':max_Z,
                      'rlf_Z':rlf_Z}            
        else:
            out_dict={'runoff_pattern':runoff_pattern,
                      'x':x,
                      'x_center':x_center,
                      'dt':cObj0.dt,
                      'freq_to_save':cObj0.freq_to_save,
                      'A':A,
                      'chi':chi,
                      'z0':z0,
                      'slp0':slp0,
                      'ts':ts,
                      'U':U,
                      'zout':zout,
                      'sout':sout,
                      'eout':eout,
                      'eCout':eCout,
                      'dout':dout,
                      'qout':qout,
                      'zcout':zcout,
                      'mrout':mrout,
                      'crout':crout,
                      'awcr':cr_wavg,
                      'fitcr':cr_fit,
                      'srout':srout}
        return out_dict
    
    def build_spim_dict(self,j,ts,sample_freq):
        fname=os.path.join(self.parent_directory,self.model_list[j],self.prefix_list[j]+str(ts[0])+'.pkl')
        # Load static values
        [cObj0,_,sObj0,eObj0,rObj0]=self.recover_spim_state(fname)
        x=sObj0.x
        z0=sObj0.z0
        chi=sObj0.chi
        A=sObj0.A
        slp0=sObj0.slope0
        x_center=sObj0.x_cents
        runoff_pattern=rObj0.pattern
        U=np.mean(eObj0.uplift)
        
        # Generate empty array for storage
        zout=np.zeros((len(ts),len(z0)))
        eout=np.zeros((len(ts),len(z0)))
        eCout=np.zeros((len(ts),len(z0)))
        qout=np.zeros((len(ts),len(z0)))
        sout=np.zeros((len(ts),len(z0)))
        kout=np.zeros((len(ts),len(z0)))
        zcout=np.zeros((len(ts),len(x_center)))
        mrout=np.zeros((len(ts),len(x_center)))
        crout=np.zeros((len(ts),len(x_center)))
        srout=np.zeros((len(ts),len(x_center)))
        cr_wavg=np.zeros(len(ts))
        cr_fit=np.zeros(len(ts))
        
        if runoff_pattern=='emp':
            snp=np.zeros((len(ts),len(x_center)))
            max_Z=np.zeros((len(ts),len(x_center)))
            rlf_Z=np.zeros((len(ts),len(x_center)))
            
        if sample_freq>0:
            print('Starting calculation for '+self.descript_list[j])

        # Read in each
        for i in range(len(ts)):
            fname=os.path.join(self.parent_directory,self.model_list[j],self.prefix_list[j]+str(ts[i])+'.pkl')
            [_,_,sObjOI,eObjOI,rObjOI]=self.recover_spim_state(fname)
            zout[i,:]=sObjOI.z
            sout[i,:]=sObjOI.slope
            eout[i,:]=eObjOI.I
            eCout[i,:]=eObjOI.cum_E
            qout[i,:]=rObjOI.qbar
            kout[i,:]=eObjOI.K
            zcout[i,:]=sObjOI.z_cents
            mrout[i,:]=rObjOI.r_bin
            crout[i,:]=rObjOI.cr_bin
            srout[i,:]=rObjOI.sr_bin
            if sample_freq==0:
                cr_wavg[i]=np.mean(rObjOI.cr_bin) 
                cr_fit[i]=np.nan
            else:
                if np.mod(i,sample_freq)==0:
                    print(f'Routing runoff for timestep {i} of {len(ts)}')
                    mcd=rObjOI.monte_carlo_runoff(sObjOI,500,1,verbose=False)
                    cr_wavg[i]=mcd['aw_cr']
                    cr_fit[i]=mcd['tail_fit_cr']
                else:
                    cr_wavg[i]=np.nan
                    cr_fit[i]=np.nan

            if runoff_pattern=='emp':
                snp[i,:]=rObjOI.snp_bin
                max_Z[i,:]=rObjOI.max_Z_bin
                rlf_Z[i,:]=rObjOI.rlf_Z_bin
        
        # Package in dictionary for output
        if runoff_pattern=='emp':
            out_dict={'runoff_pattern':runoff_pattern,
                      'dt':cObj0.dt,
                      'x':x,
                      'x_center':x_center,
                      'A':A,
                      'chi':chi,
                      'z0':z0,
                      'slp0':slp0,
                      'U':U,
                      'ts':ts,
                      'freq_to_save':cObj0.freq_to_save,
                      'zout':zout,
                      'sout':sout,
                      'eout':eout,
                      'eCout':eCout,
                      'qout':qout,
                      'kout':kout,
                      'zcout':zcout,
                      'mrout':mrout,
                      'crout':crout,
                      'awcr':cr_wavg,
                      'fitcr':cr_fit,
                      'srout':srout,
                      'snp':snp,
                      'max_Z':max_Z,
                      'rlf_Z':rlf_Z}            
        else:
            out_dict={'runoff_pattern':runoff_pattern,
                      'dt':cObj0.dt,
                      'x':x,
                      'x_center':x_center,
                      'A':A,
                      'chi':chi,
                      'z0':z0,
                      'slp0':slp0,
                      'U':U,
                      'ts':ts,
                      'freq_to_save':cObj0.freq_to_save,
                      'zout':zout,
                      'sout':sout,
                      'eout':eout,
                      'eCout':eCout,
                      'qout':qout,
                      'kout':kout,
                      'zcout':zcout,
                      'mrout':mrout,
                      'crout':crout,
                      'awcr':cr_wavg,
                      'fitcr':cr_fit,
                      'srout':srout}
        return out_dict

    def parse_multiple(self,sample_freq):
        d_list=[]
        for i in range(self.num_models):
            ts=self.sorted_time_steps(self.model_list[i],self.prefix_list[i])
            mod_type=self.determine_type([ts[0]])
            if mod_type[0]=='STIM':
                out_dict=self.build_stim_dict(i,ts,sample_freq)
                d_list.append(out_dict)
            elif mod_type[0]=='SPIM':
                out_dict=self.build_spim_dict(i,ts,sample_freq)
                d_list.append(out_dict)
        return d_list
    
    def derive_values(self,d):
        ksn=np.diff(d['zout'],axis=1)/np.diff(d['chi'])
        mean_ksn=np.mean(ksn,axis=1)
        mean_z=np.mean(d['zout'],axis=1)
        max_Q=d['qout'][:,0]
        
        mean_cr=d['awcr']
        fit_cr=d['fitcr']
        ts_cr=d['ts']
        idx=np.isnan(mean_cr)
        mean_cr=mean_cr[~idx]
        fit_cr=fit_cr[~idx]
        ts_cr=ts_cr[~idx]
        
        eCout=d['eCout']
        fts=d['freq_to_save']
        avgE=np.diff(eCout,axis=0)
        avgE=np.concatenate((eCout[0,:].reshape((1,len(eCout[0,:]))),avgE),axis=0)
        avgE=((avgE*(-1))/fts)*(10*100)
        mean_E=np.mean(avgE,axis=1)
        
        if d['runoff_pattern']=='emp':
            max_snp=np.max(d['snp'],axis=1)
            return mean_ksn,mean_z,mean_E,max_Q,mean_cr,fit_cr,ts_cr,max_snp
        else:
            return mean_ksn,mean_z,mean_E,max_Q,mean_cr,fit_cr,ts_cr
        
    
    def plot_all(self,sample_freq=1):
        d_list=self.parse_multiple(sample_freq)
        col_vec=colors.Normalize(vmin=0,vmax=self.num_models-1)
        
        plt.figure(figsize=(15,15))
        ax1=plt.subplot(3,2,1)
        ax1.set_xlabel('Time [Myr]')
        ax1.set_ylabel(r'Mean $k_{sn}$ [m]')

        ax2=plt.subplot(3,2,2)
        ax2.set_xlabel('Time [Myr]')
        ax2.set_ylabel('Mean Elevation [m]')        
 
        ax3=plt.subplot(3,2,3)
        ax3.set_xlabel('Time [Myr]')
        ax3.set_ylabel('Mean Erosion Rate [mm/yr]')
        ax3.axhline(d_list[0]['U']*(10*100),c='k',linestyle=':',label='Uplift Rate')
        
        ax4=plt.subplot(3,2,4)
        ax4.set_xlabel('Time [Myr]')
        ax4.set_ylabel(r'Discharge at Outlet [$m^{3}$/s]')
        
        ax5=plt.subplot(3,2,5)
        ax5.set_xlabel('Time [Myr]')
        ax5.set_ylabel('Mean Variability')
 
        for i in range(self.num_models):
            if d_list[i]['runoff_pattern']=='emp':
                [mean_ksn,mean_z,mean_E,max_Q,mean_cr,fit_cr,ts_cr,max_snp]=self.derive_values(d_list[i])
                ax6=plt.subplot(3,2,6)
                ax6.set_xlabel('Time [Myr]')
                ax6.set_ylabel('Max Snow Fraction')
                ax6.plot(d_list[i]['ts']/1e6,max_snp,c=cm.roma(col_vec(i)))
            else:
                [mean_ksn,mean_z,mean_E,max_Q,mean_cr,fit_cr,ts_cr]=self.derive_values(d_list[i])                
            ax1.plot(d_list[i]['ts']/1e6,mean_ksn,c=cm.roma(col_vec(i)))
            ax2.plot(d_list[i]['ts']/1e6,mean_z,c=cm.roma(col_vec(i)))
            ax3.plot(d_list[i]['ts']/1e6,mean_E,c=cm.roma(col_vec(i)),label=self.descript_list[i])
            ax4.plot(d_list[i]['ts']/1e6,max_Q,c=cm.roma(col_vec(i)))
            
            ax5.plot(ts_cr/1e6,mean_cr,c=cm.roma(col_vec(i)))
            ax5.plot(ts_cr/1e6,fit_cr,c=cm.roma(col_vec(i)),linestyle='--')

        ax3.legend(loc='best')
        
    def plot_all2(self,col_list,style_list,sample_freq=1):
        d_list=self.parse_multiple(sample_freq)
        col_vec=colors.Normalize(vmin=0,vmax=self.num_models-1)
        
        SMALL_SIZE = 8
        MEDIUM_SIZE = 10
        BIGGER_SIZE = 12
        plt.rc('font', size=SMALL_SIZE,family='Futura')          # controls default text sizes
        plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
        plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
        plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
        plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
        plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
        plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
        # Initiate figure
        f1=plt.figure(figsize=(8,8))
        f1.set_dpi(250)
        
        ax1=plt.subplot(3,2,1)
        # ax1.set_xlabel('Time [Myr]')
        ax1.set_ylabel(r'Mean $k_{sn}$ [m]')

        ax2=plt.subplot(3,2,2)
        # ax2.set_xlabel('Time [Myr]')
        ax2.set_ylabel('Mean Elevation [m]')        
 
        ax3=plt.subplot(3,2,3)
        # ax3.set_xlabel('Time [Myr]')
        ax3.set_ylabel('Mean Erosion Rate [mm/yr]')
        ax3.axhline(d_list[0]['U']*(10*100),c='k',linestyle=':',label='Uplift Rate')
        
        ax4=plt.subplot(3,2,4)
        # ax4.set_xlabel('Time [Myr]')
        ax4.set_ylabel(r'Discharge at Outlet [$m^{3}$/s]')
        
        ax5=plt.subplot(3,2,5)
        ax5.set_xlabel('Time [Myr]')
        ax5.set_ylabel('Shape Parameter')
 
        for i in range(self.num_models):
            if d_list[i]['runoff_pattern']=='emp':
                [mean_ksn,mean_z,mean_E,max_Q,mean_cr,fit_cr,ts_cr,max_snp]=self.derive_values(d_list[i])
                ax6=plt.subplot(3,2,6)
                ax6.set_xlabel('Time [Myr]')
                ax6.set_ylabel('Max Snow Fraction')
                ax6.plot(d_list[i]['ts']/1e6,max_snp,c=col_list[i],linestyle=style_list[i],linewidth=0.75)
            else:
                [mean_ksn,mean_z,mean_E,max_Q,mean_cr,fit_cr,ts_cr]=self.derive_values(d_list[i])                
            ax1.plot(d_list[i]['ts']/1e6,mean_ksn,c=col_list[i],linestyle=style_list[i],linewidth=0.75)
            ax2.plot(d_list[i]['ts']/1e6,mean_z,c=col_list[i],linestyle=style_list[i],linewidth=0.75)
            ax3.plot(d_list[i]['ts']/1e6,mean_E,c=col_list[i],linestyle=style_list[i],linewidth=0.75,label=self.descript_list[i])
            ax4.plot(d_list[i]['ts']/1e6,max_Q,c=col_list[i],linestyle=style_list[i],linewidth=0.75)
            
            if i==0:
                ax5.plot(ts_cr/1e6,mean_cr,c=col_list[i],linestyle=style_list[i],linewidth=0.5,label='Area Weighted')
                ax5.plot(ts_cr/1e6,fit_cr,c=col_list[i],linestyle=style_list[i],linewidth=2,label='Routed')
            else:
                ax5.plot(ts_cr/1e6,mean_cr,c=col_list[i],linestyle=style_list[i],linewidth=0.5)
                ax5.plot(ts_cr/1e6,fit_cr,c=col_list[i],linestyle=style_list[i],linewidth=2)
        ax5.legend(loc='best')
        ax3.legend(loc='best')
        
        ax1.text(0.01, 0.99, 'A',
                horizontalalignment='left',
                verticalalignment='top',
                transform=ax1.transAxes,
                fontsize=12,fontweight='extra bold')
        ax2.text(0.01, 0.99, 'B',
                horizontalalignment='left',
                verticalalignment='top',
                transform=ax2.transAxes,
                fontsize=12,fontweight='extra bold')
        ax3.text(0.01, 0.99, 'C',
                horizontalalignment='left',
                verticalalignment='top',
                transform=ax3.transAxes,
                fontsize=12,fontweight='extra bold')
        ax4.text(0.01, 0.99, 'D',
                horizontalalignment='left',
                verticalalignment='top',
                transform=ax4.transAxes,
                fontsize=12,fontweight='extra bold') 
        ax5.text(0.01, 0.99, 'E',
                horizontalalignment='left',
                verticalalignment='top',
                transform=ax5.transAxes,
                fontsize=12,fontweight='extra bold') 
        ax6.text(0.01, 0.99, 'F',
                horizontalalignment='left',
                verticalalignment='top',
                transform=ax6.transAxes,
                fontsize=12,fontweight='extra bold')               
        
        plt.tight_layout()
        plt.rcdefaults()
        return f1

        
            


                     

    
            
            
    
        

        
    