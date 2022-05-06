# -*- coding: utf-8 -*-
"""
Drifter Tracking using velocity field from FVCOM GOM3 model

NOTE: THIS VERSION WAS DERIVED FROM VITALII'S dtr_32_mean.py

2019-06-19 dtr_v32.py
split fields into separate npy files
optimizing reads from disk using
ua=np.load(FNU,mmap_mode='r')
get_uv in line

2021-11-23 dtr_JiM.py
JiM's minor modifications taken from Vitalii's "dtr_v33_c.py"
modified format so all function lines are indented
got rid of unused functions

2022-04-17 using climatology fields with to output single dataframe

2022-05-04 
put all functions in a module called "pt_functions.py"
getting back to actual monthly fields (not climatology)
moved all files from "FATE" folder to "PT" folder
"""

import numpy as np
import matplotlib.pyplot as plt
import os
from datetime import datetime,timedelta
import multiprocessing as mp
import sys
import pandas as pd
from get_fvcom_gom3_grid import get_fvcom_gom3_grid
from pt_functions import *
########################
####################HARDCODES ###################
NCPU=1
NUMDAYS=30 # NUMBER OF DAYS TO TRACK
D='a' # a means vertically averaged
run='201608' # use 'JUNE_MEAN', for example, for climatology fields
YEAR0='2016'# use '2001' for climatology fields
year0=int(YEAR0)
CASE='VERTICALLY_AVE_'+run+'_OVER_'+str(NUMDAYS)+'_DAYS'# text for plot title
FNUV='gom3_uv_'+run+'.npy' # output of get_fvcom_mon_nc.py that has a dictionary
# set up days of a year to start runs
# note that 2001 is a regular year, the year value will be overwritten later 
#t0dates=np.arange(datetime(year0,6,1),datetime(year0,6,30),timedelta(days=31)).astype(datetime) 
t0dates=np.arange(datetime(year0,1,1),datetime(year0,1,30),timedelta(days=31)).astype(datetime) 
####################### HARDCODES END HERE ###################################
############################################################################
    
# FVCOM GOM3 triangular grid
############################################################

#Grid=get_fvcom_gom3_grid('server') # download a fresh copy of Grid
Grid=get_fvcom_gom3_grid('dict') # load from local file
#Grid=np.load('Grid.npy').item()
lon=Grid['lon']
lat=Grid['lat']
lonc=Grid['lonc']
latc=Grid['latc']
bathy=Grid['h']
############################################################
"""
# test values of initial conditions
lond0=np.arange(-69.,-66.,0.25)
latd0=np.arange(41.,43.,0.25)

lond0=np.arange(-70.93,-70.87,0.01)
latd0=np.arange(42.31,42.36,0.01)

#lond0=np.array([-66.])
#latd0=np.array([ 42.])

llond0,llatd0=np.meshgrid(lond0,latd0)
llond0=llond0.flatten()
llatd0=llatd0.flatten()
ND=llond0.size
"""
############################################################################
# delete files if you want new initial conditions
try:
    #llond0=np.load('llond0_BH.npy')# where BH=Boston Harbor
    #llatd0=np.load('llatd0_BH.npy')
    # Mass Disposal is at 42° 25.55’ N, 70° 35.13
    llond0=np.array([-70.625,-70.575,-70.575,-70.625])
    llatd0=np.array([42.375,42.375,42.425,42.425])
    #llond0=llond0[0:1] #JiM subsampling start positions for testing
    #llatd0=llatd0[0:1]
    print('loading init positions')
    ND=llond0.size
    print('ND =',ND)
    print(llond0.min(),llond0.max())
    print(llatd0.min(),llatd0.max())

except:
    # calculate intial conditions for some cases    
    print('ipython run Vitaliis dtr_init_positions.py to make llond0.npy llatd0.npy')
    sys.exit()
    
#############################################################################


for t0date in t0dates: # loop through a list of datetimes to start tracking
        
    mo=int(t0date.month)
    da=int(t0date.day)
    MO=str(mo).zfill(2)
    DA=str(da).zfill(2)
    print('Drifter Launch:',YEAR0,MO,DA,D)
    
    t0=RataDie(year0,mo,da)
    FTS=YEAR0+MO+DA
    
    # specify the length of calculation, typically 62 days
    #t1=t0+62
    #t1=t0+30
    t1=t0+NUMDAYS
    
    # specify time step: standard 1h
    # number of time steps per hour
    #    MH=5 # 1h/10=6min step
    MH=1 # 1h step  only 1h step is implemented 
    
            #    if MH==1:
            #        
            #        if D=='a':
            #            get_uv=get_uv1_from_gom3_uva_YYYYMM_npy_tRD
            #        else:
            #            get_uv=get_uv1_from_gom3_uvs_YYYYMM_npy_tRD
            #        
            #    else:
            #        get_uv=get_uv2 # needed for step size equal to a fraction of an hour 
        
    dtsec=60*60./MH
    tau=dtsec/111111. # deg per (velocityunits*dt)
    # dt in seconds
    # vel units m/s
    # in other words v*tau -> deg 
    
    dtday=1./24./MH
    tt=np.arange(t0,t1+dtday,dtday)# where t1 is defined above as the end of the tracking time
    NT=len(tt)
    lont=np.zeros((NT,ND))
    latt=np.zeros((NT,ND))
    #tempt=np.zeros((NT,ND))
    
    # initial positions
    lont[0,:]=llond0
    latt[0,:]=llatd0
    kd=0;kt=0; print(lont[kt,kd],latt[kt,kd])

    dEpochMJD_RD=(datetime(1858,11,17)-datetime(1,1,1)).days+1.    #Out[5]: 678576
    t0dati=datetime(1858,11,17) # MJD=JD-2400000.5
    
    #time dependent u,v
    kt=0
    tRD=tt[kt]    
    #u1,v1=get_uv(tRD)
    # inline get_uv #############################################
    tMJD=tRD-dEpochMJD_RD
    tdati=t0dati+timedelta(seconds=tMJD*86400.)
    YEAR=str(tdati.year) # current year during the run, may become YEAR0+1
    MO=str(tdati.month).zfill(2)
    YEARMO=YEAR+MO
    
    #if YEARMO != YEARMO0: # load new monthly file
    print('loading uv',YEARMO)
    if run[-4:]=='MEAN':
        # took the following from "dtr_32_mean.py" and files from vitalii's "FVCOM_GOM3_Climate_mean_19782016" folder
        FNU='gom3_ua_mean.npy'
        FNV='gom3_va_mean.npy'
        ua=np.load(FNU,mmap_mode='r')
        va=np.load(FNV,mmap_mode='r')
        u1=ua*1.
        v1=va*1.
        it=0
    else: # case where we have a time-varying monthly field like "gom3_uv_197801.npy" 
        dict=np.load(FNUV,allow_pickle=True).item()
        ua=dict['ua']
        va=dict['va']
        tFVCOM=dict['time']
        it=np.argmin(np.abs(tFVCOM-tMJD))
        u1=ua[it,:].flatten() # do we need flatten?
        v1=va[it,:].flatten() # 
    
    for kt in range(NT-1): # loop though the total number of time steps NT
        # time dependent u,v at current time level (from previous step)    
        u0=u1*1.0;v0=v1*1.0
        YEARMO0=YEAR+MO 
   
        #time dependent u,v at next time level
        tRD=tt[kt+1]
        #datet=datet.append(t0dates[0]+timedelta(hours=kt))
        #u1,v1=get_uv(tRD)
        # inline get_uv #############################################
        tMJD=tRD-dEpochMJD_RD

        tdati=t0dati+timedelta(seconds=tMJD*86400.)
        YEAR=str(tdati.year)
        MO=str(tdati.month).zfill(2)
        YEARMO=YEAR+MO
        
        if YEARMO != YEARMO0: # load new monthly file where not necessary if < 30 days
            print('loading uv',YEARMO)
            FNT='/home/vsheremet/GOM3_DATAMon_uva/gom3_time_'+YEARMO+'.npy'
            FNU='/home/vsheremet/GOM3_DATAMon_uva/gom3_ua_'+YEARMO+'.npy'
            FNV='/home/vsheremet/GOM3_DATAMon_uva/gom3_va_'+YEARMO+'.npy'
            tFVCOM=np.load(FNT,mmap_mode='r')
            ua=np.load(FNU,mmap_mode='r')
            va=np.load(FNV,mmap_mode='r')
        
        if run[-4:]!='MEAN':
            it=np.argmin(np.abs(tFVCOM-tMJD))
            t1dati=t0dati+timedelta(seconds=tFVCOM[it]*86400.)
            u1=ua[it,:].flatten() # do we need flatten?
            v1=va[it,:].flatten() # 
        ui=(u0+u1)*0.5;vi=(v0+v1)*0.5 # velocity at the middle of time step, linear interpolation
    
        for kd in range(ND): # where "ND" is the number of particles
                # for each drifter make one time step using classic 4th order Runge-Kutta method        
                lont[kt+1,kd],latt[kt+1,kd]=RungeKutta4_lonlattime(lont[kt,kd],latt[kt,kd],Grid,u0,v0,ui,vi,u1,v1,tau)
    
    datets=np.arange(t0date,t0date+timedelta(hours=len(tt)),timedelta(hours=1)).astype(datetime)
    if t0date==t0dates[0]:
        df=pd.DataFrame([tt,latt,lont]).T
    else:
        df1=pd.DataFrame([tt,latt,lont]).T
        df=pd.concat([df,df1])
    df['datet']=datets
    df['launch_date']=t0date
    df=df.set_index('datet')
    df=df.drop(0,axis=1) 
    df.to_csv(CASE+'.csv')

    plt.figure()
    plt.axes().set_aspect(1.3)
    for kt in [0,1,2,3]:
            plt.plot(lont[:,kt],latt[:,kt])
    #plt.tricontour(lon,lat,bathy,[0.,100.],colors='k')
    #basemap_region('wv')
    c=pd.read_csv('capecod.txt',header=None,names=['lon','lat']) # this is a quick coastline option
    plt.plot(c.lon,c.lat)
    plt.plot(np.mean(llond0),np.mean(llatd0),'mo',markersize=30,zorder=1)
    plt.text(np.mean(llond0),np.mean(llatd0),'MBDS',color='w',verticalalignment='center',horizontalalignment='center')
    plt.ylim(np.min(c.lat),np.max(latt)+.05)
    plt.xlim(np.min(lont)-.2,np.max(lont)+.1)
    plt.title(CASE,fontsize=10,fontweight='bold')
    plt.savefig(CASE+'.png')
