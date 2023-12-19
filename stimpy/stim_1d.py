#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 08:42:43 2022

@author: aforte
"""
import numpy as np
import os
import pickle
import time
import glob
import re
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from cmcrameri import cm
from matplotlib import colors
from matplotlib import cm as cmm

class Stim1D:
    def __init__(self,output_dir,file_prefix='ts_'):
        self.output_dir = output_dir
        self.file_prefix = file_prefix
        self.__gen_folder()
        
    def __gen_folder(self):
        # Creates folder if it doesn't exist
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
            
    def save_state(self,cObj,sObj,eObj,rObj,yr):
        fname=os.path.join(self.output_dir,self.file_prefix+str(yr)+'.pkl')
        with open(fname,'wb') as f:
            pickle.dump(cObj,f)
            pickle.dump(sObj,f)
            pickle.dump(eObj,f)
            pickle.dump(rObj,f)
            
    def recover_state(self,yr):
        fname=os.path.join(self.output_dir,self.file_prefix+str(yr)+'.pkl')
        with open(fname,'rb') as f:
            cObj=pickle.load(f)
            sObj=pickle.load(f)
            eObj=pickle.load(f)
            rObj=pickle.load(f)
        return cObj,sObj,eObj,rObj

    def __parse_time(self,cObj,sObj,tic,toc):
        yr=np.round(((cObj.i+1)*cObj.dt)/365,0).astype(int)
        tic_toc=(toc-tic)/(1e9*60)
        # Time_per step elpased
        if cObj.k==0:
            time_per_step=tic_toc/cObj.dif_steps[0]
        else: 
            time_per_step=tic_toc/cObj.dif_steps[1]
        # Build status string
        # Choose relevant time period
        stat_time=time_per_step*(cObj.num_steps-1-cObj.i)
        if stat_time<=60:
            time_string=str(np.round(stat_time,2))+ ' minutes remaining'
        elif (stat_time>60) & (stat_time<=1440):
            time_string=str(np.round(stat_time/60,2))+ ' hours remaining'
        elif stat_time>1440:
            time_string=str(np.round(stat_time/(60*24),2))+ ' days remaining'

        stat_str=(str(yr)+ ' model years completed : ~' +
                  time_string + 
                  ' : Current Max Elevation - ' +
                  str(np.round(np.max(sObj.z),2)))
        return stat_str,yr 

    def run_new_model(self,cObj,sObj,eObj,rObj,restart=False,restart_from=0,restart_i=0):
        # Start main model loop
        ##
        [_,counts]=np.unique(sObj.ix,return_counts=True)
        ##
        while cObj.i < cObj.num_steps:
            # Generate random runoff            
            if np.any(cObj.i==cObj.steps_to_gen):
                rts=rObj.generate_random_runoff(cObj.rec_length,cObj.i)
                rtsl=rts.shape[1]
                ##
                qdays=rObj.route_record_discharge(sObj,rts,counts)
                ##
                cObj.update_gen_rec_i(cObj.i)
            cObj.update_j(0) # Set inner loop counter to 0 at start of loop
            while (cObj.j < rtsl) & (cObj.i < cObj.num_steps):
                # Start timer
                if (cObj.i==0):
                    # Start loop timer
                    tic = time.time_ns()
                elif (restart) & (cObj.i==restart_i):
                    # Start loop timer
                    tic=time.time_ns() 
                # Route one day of random runoff
                # This is an inefficient way to do this, vectorized
                # solution results in performance boost
                # rObj.route_one_day(sObj,rts,cObj.j)
                ##
                rObj.qd=qdays[cObj.j,:]
                ##
                # Uplift and incise
                eObj.uplift_and_incise(sObj,cObj.dt,rObj.qd,rObj.qbar)
                # Output if necessary
                if (np.any(cObj.i==cObj.steps_to_save)) & (cObj.i>0) & (cObj.i>restart_i):
                    # Stop timer
                    toc=time.time_ns()
                    [stat_str,yr]=self.__parse_time(cObj,sObj,tic,toc)
                    self.save_state(cObj,sObj,eObj,rObj,yr)
                    cObj.update_k(cObj.k+1)
                    print(stat_str)
                    # Restart timer
                    tic=time.time_ns()
                # Update runoff based on updated elevations
                rObj.update_r(sObj)
                # Increment counters
                cObj.update_i(cObj.i+1)
                cObj.update_j(cObj.j+1)
                
    def restart_model(self,restart_from=0,new_total_time=0,new_freq_to_save=None,new_output_dir=None):
        # Load in prior model
        [cObjR,sObjR,eObjR,rObjR]=self.recover_state(restart_from)
        # Perform check and reset time
        if new_total_time <= 0:
            raise Exception('Must set new_total_time to value greater than zero')
        elif new_total_time <= restart_from:
            raise Exception('Must set new_total_time parameter to a value greater than the restart_from time')
        else:
            cObjR.total_time=new_total_time

        # Check for any necessary updates
        if new_freq_to_save!=None:
            if new_freq_to_save<cObjR.rec_length:
                raise Exception('Updated frequency to save would be less than the record length')
            else:
                cObjR.freq_to_save=new_freq_to_save
        if new_output_dir!=None:
            self.output_dir=new_output_dir
            self.__gen_folder()

        # Regenerate counting variables based on new time
        cObjR.generate_counts()
        # Generate restart_i
        cObjR.update_i(cObjR.i+1)
        cObjR.update_j(cObjR.j+1)
        restart_i=cObjR.i
        # Call run_new_model and pass updated values
        self.run_new_model(cObjR,sObjR,eObjR,rObjR,restart=True,
                            restart_from=restart_from,restart_i=restart_i)
        
    def parse_results(self,lower_ts,upper_ts,trim_to_ss,ss_value):
        # Generate unsorted file list
        p_files=glob.glob(os.path.join(self.output_dir,self.file_prefix+'*.pkl'))
        # Extract numbers
        p_nums=np.zeros(len(p_files))
        regex=re.compile(r'\d+')
        for i in range(len(p_files)):
            p_nums[i]=[int(j) for j in regex.findall(p_files[i])][-1]
        # Convert file list to array to ease indexing
        p_files=np.array([p_files])
        p_files=p_files.reshape((p_files.shape[1]))
        # Trim timesteps if necessary
        p_idx=(p_nums>=lower_ts) & (p_nums<=upper_ts)
        p_nums=p_nums[p_idx]
        p_files=p_files[p_idx]
        # Created sorted index
        pix=np.argsort(p_nums)
        # Apply index to list and generate list of model times
        # p_files_sorted=[p_files[i] for i in pix]
        ts=p_nums[pix]
        # Load static values
        [cObj0,sObj0,eObj0,rObj0]=self.recover_state(ts[0].astype(int))
        x=sObj0.x
        z0=sObj0.z0
        chi=sObj0.chi
        A=sObj0.A
        slp0=sObj0.slope0
        x_center=sObj0.x_cents
        runoff_pattern=rObj0.pattern
        
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

        if runoff_pattern=='emp':
            snp=np.zeros((len(ts),len(x_center)))
            max_Z=np.zeros((len(ts),len(x_center)))
            rlf_Z=np.zeros((len(ts),len(x_center)))
        elif runoff_pattern=='emp_rain':
            max_Z=np.zeros((len(ts),len(x_center)))
            rlf_Z=np.zeros((len(ts),len(x_center)))

        # Read in each
        for i in range(len(ts)):
            [_,sObjOI,eObjOI,rObjOI]=self.recover_state(ts[i].astype(int))
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
            if runoff_pattern=='emp':
                snp[i,:]=rObjOI.snp_bin
                max_Z[i,:]=rObjOI.max_Z_bin
                rlf_Z[i,:]=rObjOI.rlf_Z_bin
            elif runoff_pattern=='emp_rain':
                max_Z[i,:]=rObjOI.max_Z_bin
                rlf_Z[i,:]=rObjOI.rlf_Z_bin
                
              
        if trim_to_ss:
            # Calculate steady-state status
            ss=np.zeros(len(ts))
            for i in range(len(ts)):
                if i==0:
                    ss[i]=np.abs(np.max(zout[i,:])-np.max(z0))
                else:
                    ss[i]=np.abs(np.max(zout[i,:])-np.max(zout[i-1,:]))
            # Find index and trim all variables       
            ss_ix=np.nonzero(ss<ss_value)[0][:1]
            if len(ss_ix)>0:
                ss_ix=ss_ix[0]+1
            else:
                ss_ix=len(ts)
            ts=ts[0:ss_ix]
            zout=zout[0:ss_ix,:]
            eout=eout[0:ss_ix,:]
            eCout=eCout[0:ss_ix,:]
            dout=dout[0:ss_ix,:]
            qout=qout[0:ss_ix,:]
            zcout=zcout[0:ss_ix,:]
            mrout=mrout[0:ss_ix,:]
            crout=crout[0:ss_ix,:]
            srout=srout[0:ss_ix,:]
            if runoff_pattern=='emp':
                snp=snp[0:ss_ix,:]
                max_Z=max_Z[0:ss_ix,:]
                rlf_Z=rlf_Z[0:ss_ix,:]
            elif runoff_pattern=='emp_rain':
               max_Z=max_Z[0:ss_ix,:]
               rlf_Z=rlf_Z[0:ss_ix,:]
        
        # Package in dictionary for output
        if runoff_pattern=='emp':
            out_dict={'runoff_pattern':runoff_pattern,
                      'uplift':eObj0.uplift,
                      'x':x,
                      'x_center':x_center,
                      'dt':cObj0.dt,
                      'freq_to_save':cObj0.freq_to_save,
                      'A':A,
                      'chi':chi,
                      'z0':z0,
                      'slp0':slp0,
                      'ts':ts,
                      'zout':zout,
                      'sout':sout,
                      'eout':eout,
                      'eCout':eCout,
                      'dout':dout,
                      'qout':qout,
                      'zcout':zcout,
                      'mrout':mrout,
                      'crout':crout,
                      'srout':srout,
                      'snp':snp,
                      'max_Z':max_Z,
                      'rlf_Z':rlf_Z} 
        elif runoff_pattern=='emp_rain':
            out_dict={'runoff_pattern':runoff_pattern,
                      'uplift':eObj0.uplift,
                      'x':x,
                      'x_center':x_center,
                      'dt':cObj0.dt,
                      'freq_to_save':cObj0.freq_to_save,
                      'A':A,
                      'chi':chi,
                      'z0':z0,
                      'slp0':slp0,
                      'ts':ts,
                      'zout':zout,
                      'sout':sout,
                      'eout':eout,
                      'eCout':eCout,
                      'dout':dout,
                      'qout':qout,
                      'zcout':zcout,
                      'mrout':mrout,
                      'crout':crout,
                      'srout':srout,
                      'max_Z':max_Z,
                      'rlf_Z':rlf_Z}
        else:
            out_dict={'runoff_pattern':runoff_pattern,
                      'uplift':eObj0.uplift,
                      'x':x,
                      'x_center':x_center,
                      'dt':cObj0.dt,
                      'freq_to_save':cObj0.freq_to_save,
                      'A':A,
                      'chi':chi,
                      'z0':z0,
                      'slp0':slp0,
                      'ts':ts,
                      'zout':zout,
                      'sout':sout,
                      'eout':eout,
                      'eCout':eCout,
                      'dout':dout,
                      'qout':qout,
                      'zcout':zcout,
                      'mrout':mrout,
                      'crout':crout,
                      'srout':srout}
        return out_dict
    
    def calculate_ss(self):
        d=self.parse_results(0,np.inf,False,1e-11)
        ts=d['ts']/1e6
        zout=d['zout']
        z0=d['z0']
        U=d['uplift']
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
        return ts,U,delta_maxZ,delta_mnZ,max_deltaZ
    
    def plot_profile_results(self,title,n_step=1,lower_ts=0,upper_ts=np.inf,
                             trim_to_ss=False,ss_value=1e-1):
        
        d=self.parse_results(lower_ts,upper_ts,trim_to_ss,ss_value)
        
        # Unpack
        runoff_pattern=d['runoff_pattern']
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
        dout=d['dout']
        eout=d['eout']*(-1)*(10*100)*(dt*365)
        qout=d['qout']
        mrout=d['mrout']
        crout=d['crout']
        rlf_Z=d['rlf_Z']
        
        # Calculate average erosion rate
        avgE=np.diff(eCout,axis=0)
        avgE=np.concatenate((eCout[0,:].reshape((1,len(eCout[0,:]))),avgE),axis=0)
        avgE=((avgE*(-1))/fts)*(10*100)
        
        # Calculate percentage of days where erosion threshold is exceeded
        excdP=np.diff(dout,axis=0)
        excdP=np.concatenate((dout[0,:].reshape((1,len(dout[0,:]))),excdP),axis=0)
        excdP=excdP/(fts*365)
        
        # Begin plotting
        col_vec=colors.Normalize(vmin=0,vmax=len(ts)-1)
        
        f1=plt.figure(figsize=(25,15))
        ax1a=plt.subplot(3,2,1)
        plt.title(title)
        ax1a.plot(x/1000,z0/1000,c='k',linestyle=':')
        for i in range(0,len(ts),n_step):
            ax1a.plot(x/1000,zout[i,:]/1000,c=cm.vik_r(col_vec(i)))
        ax1a.set_xlabel('Stream Distance [km]')
        ax1a.set_ylabel('Elevation [km]')      
        
        
        plt.subplot(3,2,3)
        plt.plot(x[1:len(x)]/1000,np.diff(z0)/np.diff(chi),c='k',linestyle=':')
        for i in range(0,len(ts),n_step):
            plt.plot(x[1:len(x)]/1000,np.diff(zout[i,:])/np.diff(chi),c=cm.vik_r(col_vec(i)))
        plt.xlabel('Stream Distance [km]')
        plt.ylabel('$k_{sn}$ [m]')        
        
        plt.subplot(3,2,5)
        for i in range(0,len(ts),n_step):
            plt.plot(x/1000,eCout[i,:],c=cm.vik_r(col_vec(i)))
        plt.xlabel('Stream Distance [km]')
        plt.ylabel('Cumulative Erosion [m]')        

        ax1b=plt.subplot(3,2,2)
        for i in range(0,len(ts),n_step):
            ax1b.plot(x/1000,qout[i,:],c=cm.vik_r(col_vec(i)))
        ax1b.set_xlabel('Stream Distance [km]')
        ax1b.set_ylabel('Mean Discharge [$m^{3}/s$]')    

        plt.subplot(3,2,4)
        plt.plot(chi,z0/1000,c='k',linestyle=':')
        for i in range(0,len(ts),n_step):
            plt.plot(chi,zout[i,:]/1000,c=cm.vik_r(col_vec(i)))
        plt.xlabel(r'$\chi$')
        plt.ylabel('Elevation [km]')
        
        ax3=plt.subplot(3,2,6)
        plt.plot(A[1:len(A)],slp0[1:len(A)],c='k',linestyle=':')
        for i in range(0,len(ts),n_step):
            plt.plot(A[1:len(A)],slp[i,1:len(A)],c=cm.vik_r(col_vec(i)))
        plt.xlabel(r'Drainage Area [$m^{2}$]')
        plt.ylabel('Slope [m/m]')
        plt.xscale('log')
        plt.yscale('log')
        norm=colors.Normalize(vmin=np.min(ts)/1e6,vmax=np.max(ts)/1e6)
        cbar1=plt.colorbar(cmm.ScalarMappable(norm=norm,cmap=cm.vik_r),ax=ax3)
        cbar1.ax.set_ylabel('Time [Myrs]')
   
        f2=plt.figure(figsize=(20,15))
        ax1=plt.subplot(3,2,1)
        im1=plt.imshow(np.flipud(mrout),extent=[x_center[0]/1000,x_center[-1]/1000,ts[0]/1e6,ts[-1]/1e6],
                   cmap=cm.batlow_r,norm=colors.Normalize(vmin=np.min(mrout.ravel()),vmax=np.max(mrout.ravel())),aspect='auto')
        cbar1=plt.colorbar(im1,ax=ax1)
        plt.xlabel('Stream Distance [km]')
        plt.ylabel('Model Time [Myr]')
        cbar1.ax.set_ylabel('Mean Runoff [mm/day]')
        plt.title(title)
                
        ax2=plt.subplot(3,2,2)
        im2=plt.imshow(np.flipud(crout),extent=[x_center[0]/1000,x_center[-1]/1000,ts[0]/1e6,ts[-1]/1e6],
                   cmap=cm.lapaz,norm=colors.Normalize(vmin=np.min(crout.ravel()),vmax=np.max(crout.ravel())),aspect='auto')
        cbar2=plt.colorbar(im2,ax=ax2)
        plt.xlabel('Stream Distance [km]')
        plt.ylabel('Model Time [Myr]')
        cbar2.ax.set_ylabel('Variability') 
       
        ax3=plt.subplot(3,2,3)
        im3=plt.imshow(np.flipud(excdP),extent=[x[0]/1000,x[-1]/1000,ts[0]/1e6,ts[-1]/1e6],
                   cmap=cm.hawaii_r,norm=colors.Normalize(vmin=0,vmax=np.max(excdP.ravel())),aspect='auto')
        cbar3=plt.colorbar(im3,ax=ax3)
        plt.xlabel('Stream Distance [km]')
        plt.ylabel('Model Time [Myr]')
        cbar3.ax.set_ylabel('Frequency of Erosion Threshold Exceedance') 

        
        ax4=plt.subplot(3,2,4)
        im4=plt.imshow(np.flipud(avgE),extent=[x[0]/1000,x[-1]/1000,ts[0]/1e6,ts[-1]/1e6],
                    cmap=cm.hawaii_r,norm=colors.Normalize(vmin=np.min(avgE.ravel()),vmax=np.max(avgE.ravel())),aspect='auto')
        cbar4=plt.colorbar(im4,ax=ax4)
        plt.xlabel('Stream Distance [km]')
        plt.ylabel('Model Time [Myr]')
        cbar4.ax.set_ylabel('Average Erosion Rate [mm/yr]')          
        
        if runoff_pattern=='emp':  
            
            ax5=plt.subplot(3,2,5)
            im5=plt.imshow(np.flipud(rlf_Z),extent=[x_center[0]/1000,x_center[-1]/1000,ts[0]/1e6,ts[-1]/1e6],
                       cmap=cm.bamako,norm=colors.Normalize(vmin=np.min(rlf_Z.ravel()),vmax=np.max(rlf_Z.ravel())),aspect='auto')
            cbar5=plt.colorbar(im5,ax=ax5)
            plt.xlabel('Stream Distance [km]')
            plt.ylabel('Model Time [Myr]')
            cbar5.ax.set_ylabel('Relief [m]')
        
            snp=d['snp']
            ax6=plt.subplot(3,2,6)
            im6=plt.imshow(np.flipud(snp),extent=[x_center[0]/1000,x_center[-1]/1000,ts[0]/1e6,ts[-1]/1e6],
                       cmap=cm.devon,norm=colors.Normalize(vmin=np.min(snp.ravel()),vmax=np.max(snp.ravel())),aspect='auto')
            cbar6=plt.colorbar(im6,ax=ax6)
            plt.xlabel('Stream Distance [km]')
            plt.ylabel('Model Time [Myr]')
            cbar6.ax.set_ylabel('Percent Runoff as Snow') 
        
         
        f1.canvas.draw()
        f2.canvas.draw()
        
    def plot_profile_results2(self,title,num_ts,lower_ts=0,upper_ts=np.inf,
                             trim_to_ss=False,ss_value=1e-1):
        
        d=self.parse_results(lower_ts,upper_ts,trim_to_ss,ss_value)
        
        # Unpack
        runoff_pattern=d['runoff_pattern']
        ts=d['ts']
        x=d['x']
        x_center=d['x_center']
        z0=d['z0']
        # dt=d['dt']
        fts=d['freq_to_save']
        A=d['A']
        slp0=d['slp0']
        slp=d['sout']
        zout=d['zout']
        # zcout=d['zcout']
        chi=d['chi']
        eCout=d['eCout']
        dout=d['dout']
        # eout=d['eout']*(-1)*(10*100)*(dt*365)
        # qout=d['qout']
        mrout=d['mrout']
        crout=d['crout']
        # rlf_Z=d['rlf_Z']
        if runoff_pattern=='emp':
            snp=d['snp']
        
        # Calculate average erosion rate
        avgE=np.diff(eCout,axis=0)
        avgE=np.concatenate((eCout[0,:].reshape((1,len(eCout[0,:]))),avgE),axis=0)
        avgE=((avgE*(-1))/fts)*(10*100)
        
        
        # Calculate percentage of days where erosion threshold is exceeded
        excdP=np.diff(dout,axis=0)
        excdP=np.concatenate((dout[0,:].reshape((1,len(dout[0,:]))),excdP),axis=0)
        excdP=excdP/(fts*365)
        
        n_step=int(len(ts)/num_ts)
        
        # Begin plotting
        col_vec=colors.Normalize(vmin=0,vmax=len(ts)-1)
        
        SMALL_SIZE = 8
        MEDIUM_SIZE = 10
        BIGGER_SIZE = 12
        plt.rc('font', size=SMALL_SIZE,family='Futura')          # controls default text sizes
        plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
        plt.rc('axes', labelsize=SMALL_SIZE)    # fontsize of the x and y labels
        plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
        plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
        plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
        plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
        # Initiate figure
        f1=plt.figure(figsize=(14,6),layout='tight')
        f1.set_dpi(250)

        if runoff_pattern=='emp_rain':
            gs=gridspec.GridSpec(2,12,figure=f1)
    
            ax3=f1.add_subplot(gs[0,0:4])
            ax3.plot(x/1000,z0/1000,c='k',linestyle=':')
            for i in range(0,len(ts),n_step):
                ax3.plot(x/1000,zout[i,:]/1000,c=cm.batlowK_r(col_vec(i)),linewidth=0.75)
            ax3.set_xlabel('Stream Distance [km]')
            ax3.set_ylabel('Elevation [km]')  
            plt.xlim((0,np.max(x)/1000))        
    
            ax1=f1.add_subplot(gs[0,4:8])
            plt.title(title)
            ax1.plot(chi,z0/1000,c='k',linestyle=':')
            for i in range(0,len(ts),n_step):
                ax1.plot(chi,zout[i,:]/1000,c=cm.batlowK_r(col_vec(i)),linewidth=0.75)
            plt.xlabel(r'$\chi$')
            plt.ylabel('Elevation [km]')
            
            ax2=f1.add_subplot(gs[0,8:12])
            ax2.plot(A[1:len(A)],slp0[1:len(A)],c='k',linestyle=':')
            for i in range(0,len(ts),n_step):
                ax2.plot(A[1:len(A)],slp[i,1:len(A)],c=cm.batlowK_r(col_vec(i)),linewidth=0.75)
            plt.xlabel(r'Drainage Area [$m^{2}$]')
            plt.ylabel('Slope [m/m]')
            plt.xscale('log')
            plt.yscale('log')
            norm=colors.Normalize(vmin=np.min(ts)/1e6,vmax=np.max(ts)/1e6)
            cbar1=plt.colorbar(cmm.ScalarMappable(norm=norm,cmap=cm.batlowK_r),ax=ax2)
            cbar1.ax.set_ylabel('Time [Myrs]')
    
            ax4=f1.add_subplot(gs[1,0:3])
            im4=plt.imshow(np.flipud(avgE),extent=[x[0]/1000,x[-1]/1000,ts[0]/1e6,ts[-1]/1e6],
                        cmap=cm.acton,norm=colors.Normalize(vmin=np.min(avgE.ravel()),vmax=np.max(avgE.ravel())),aspect='auto')
            cbar4=plt.colorbar(im4,ax=ax4)
            plt.xlabel('Stream Distance [km]')
            plt.ylabel('Model Time [Myr]')
            ax4.set_xticks(np.linspace(0,np.max(x)/1000,6))
            cbar4.ax.set_ylabel('Average Erosion Rate [mm/yr]')
            
            ax7=f1.add_subplot(gs[1,3:6])
            im3=plt.imshow(np.flipud(excdP),extent=[x[0]/1000,x[-1]/1000,ts[0]/1e6,ts[-1]/1e6],
                        cmap=cm.hawaii_r,norm=colors.Normalize(vmin=0,vmax=np.max(excdP.ravel())),aspect='auto')
            cbar3=plt.colorbar(im3,ax=ax7)
            plt.xlabel('Stream Distance [km]')
            ax7.set_xticks(np.linspace(0,np.max(x)/1000,6))
            # plt.ylabel('Model Time [Myr]')
            cbar3.ax.set_ylabel('Ero Threshold Excd Freq')  
            
            ax6=f1.add_subplot(gs[1,6:9])
            im2=plt.imshow(np.flipud(crout),extent=[x_center[0]/1000,x_center[-1]/1000,ts[0]/1e6,ts[-1]/1e6],
                        cmap=cm.lapaz,norm=colors.Normalize(vmin=np.min(crout.ravel()),vmax=np.max(crout.ravel())),aspect='auto')
            cbar2=plt.colorbar(im2,ax=ax6)
            plt.xlabel('Stream Distance [km]')
            ax6.set_xticks(np.linspace(0,np.max(x)/1000,6))
            # plt.ylabel('Model Time [Myr]')
            cbar2.ax.set_ylabel('Shape Parameter') 
    
            ax8=f1.add_subplot(gs[1,9:12])
            im1=plt.imshow(np.flipud(mrout),extent=[x_center[0]/1000,x_center[-1]/1000,ts[0]/1e6,ts[-1]/1e6],
                        cmap=cm.vik_r,norm=colors.Normalize(vmin=np.min(mrout.ravel()),vmax=np.max(mrout.ravel())),aspect='auto')
            cbar1=plt.colorbar(im1,ax=ax8)
            plt.xlabel('Stream Distance [km]')
            ax8.set_xticks(np.linspace(0,np.max(x)/1000,6))
            # plt.ylabel('Model Time [Myr]')
            cbar1.ax.set_ylabel('Mean Runoff [mm/day]')
            
        elif runoff_pattern=='emp':
            gs=gridspec.GridSpec(2,15,figure=f1)
    
            ax3=f1.add_subplot(gs[0,0:5])
            ax3.plot(x/1000,z0/1000,c='k',linestyle=':')
            for i in range(0,len(ts),n_step):
                ax3.plot(x/1000,zout[i,:]/1000,c=cm.batlowK_r(col_vec(i)),linewidth=0.75)
            ax3.set_xlabel('Stream Distance [km]')
            ax3.set_ylabel('Elevation [km]')  
            plt.xlim((0,np.max(x)/1000))        
    
            ax1=f1.add_subplot(gs[0,5:10])
            plt.title(title)
            ax1.plot(chi,z0/1000,c='k',linestyle=':')
            for i in range(0,len(ts),n_step):
                ax1.plot(chi,zout[i,:]/1000,c=cm.batlowK_r(col_vec(i)),linewidth=0.75)
            plt.xlabel(r'$\chi$')
            plt.ylabel('Elevation [km]')
            
            ax2=f1.add_subplot(gs[0,10:15])
            ax2.plot(A[1:len(A)],slp0[1:len(A)],c='k',linestyle=':')
            for i in range(0,len(ts),n_step):
                ax2.plot(A[1:len(A)],slp[i,1:len(A)],c=cm.batlowK_r(col_vec(i)),linewidth=0.75)
            plt.xlabel(r'Drainage Area [$m^{2}$]')
            plt.ylabel('Slope [m/m]')
            plt.xscale('log')
            plt.yscale('log')
            norm=colors.Normalize(vmin=np.min(ts)/1e6,vmax=np.max(ts)/1e6)
            cbar1=plt.colorbar(cmm.ScalarMappable(norm=norm,cmap=cm.batlowK_r),ax=ax2)
            cbar1.ax.set_ylabel('Time [Myrs]')
    
            ax4=f1.add_subplot(gs[1,0:3])
            im4=plt.imshow(np.flipud(avgE),extent=[x[0]/1000,x[-1]/1000,ts[0]/1e6,ts[-1]/1e6],
                        cmap=cm.acton,norm=colors.Normalize(vmin=np.min(avgE.ravel()),vmax=np.max(avgE.ravel())),aspect='auto')
            cbar4=plt.colorbar(im4,ax=ax4)
            plt.xlabel('Stream Distance [km]')
            plt.ylabel('Model Time [Myr]')
            ax4.set_xticks(np.linspace(0,np.max(x)/1000,6))
            cbar4.ax.set_ylabel('Average Erosion Rate [mm/yr]')
            
            ax7=f1.add_subplot(gs[1,3:6])
            im3=plt.imshow(np.flipud(excdP),extent=[x[0]/1000,x[-1]/1000,ts[0]/1e6,ts[-1]/1e6],
                        cmap=cm.hawaii_r,norm=colors.Normalize(vmin=0,vmax=np.max(excdP.ravel())),aspect='auto')
            cbar3=plt.colorbar(im3,ax=ax7)
            plt.xlabel('Stream Distance [km]')
            ax7.set_xticks(np.linspace(0,np.max(x)/1000,6))
            # plt.ylabel('Model Time [Myr]')
            cbar3.ax.set_ylabel('Ero Threshold Excd Freq')  
            
            ax6=f1.add_subplot(gs[1,6:9])
            im2=plt.imshow(np.flipud(crout),extent=[x_center[0]/1000,x_center[-1]/1000,ts[0]/1e6,ts[-1]/1e6],
                        cmap=cm.lapaz,norm=colors.Normalize(vmin=np.min(crout.ravel()),vmax=np.max(crout.ravel())),aspect='auto')
            cbar2=plt.colorbar(im2,ax=ax6)
            plt.xlabel('Stream Distance [km]')
            ax6.set_xticks(np.linspace(0,np.max(x)/1000,6))
            # plt.ylabel('Model Time [Myr]')
            cbar2.ax.set_ylabel('Shape Parameter') 
    
            ax8=f1.add_subplot(gs[1,9:12])
            im1=plt.imshow(np.flipud(mrout),extent=[x_center[0]/1000,x_center[-1]/1000,ts[0]/1e6,ts[-1]/1e6],
                        cmap=cm.vik_r,norm=colors.Normalize(vmin=np.min(mrout.ravel()),vmax=np.max(mrout.ravel())),aspect='auto')
            cbar1=plt.colorbar(im1,ax=ax8)
            plt.xlabel('Stream Distance [km]')
            ax8.set_xticks(np.linspace(0,np.max(x)/1000,6))
            # plt.ylabel('Model Time [Myr]')
            cbar1.ax.set_ylabel('Mean Runoff [mm/day]')
            
            ax9=f1.add_subplot(gs[1,12:15])
            im5=plt.imshow(np.flipud(snp),extent=[x_center[0]/1000,x_center[-1]/1000,ts[0]/1e6,ts[-1]/1e6],
                        cmap=cm.imola_r,norm=colors.Normalize(vmin=np.min(snp.ravel()),vmax=np.max(snp.ravel())),aspect='auto')
            cbar5=plt.colorbar(im5,ax=ax9)
            plt.xlabel('Stream Distance [km]')
            ax8.set_xticks(np.linspace(0,np.max(x)/1000,6))
            # plt.ylabel('Model Time [Myr]')
            cbar5.ax.set_ylabel('Snowmelt Fraction')

        # plt.tight_layout()
        plt.rcdefaults()
        return f1
        
