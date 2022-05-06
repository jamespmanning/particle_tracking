# -*- coding: utf-8 -*-
"""
s01_get_fvcom_mon_nc.py 
    get fvcom data from month-long nc files and save to local disk selected fields as dictionary gom3_uvzt_YYYYMM.npy
    Further used for drifter dtracking and tide analysis

    note that u,v,lonc,latc,uwindstress,vwindstress defined at triangle centers 
              zeta,temp,lon,lat at triangle vertices

@author: Vitalii Sheremet, FATE Project, 2012-2019, vsheremet@whoi.edu

Modifications by JiM in May 2022 to a) reduce the code and b) simplified output 
"""

import numpy as np
from datetime import *
import os
from netCDF4 import Dataset        # netCDF4 version
 

URL0='http://www.smast.umassd.edu:8080/thredds/dodsC/models/fvcom/NECOFS/Archive/'

# see catalog at 
#http://www.smast.umassd.edu:8080/thredds/catalog/models/fvcom/NECOFS/Archive/Seaplan_33_Hindcast_v1/catalog.html

#http://www.smast.umassd.edu:8080/thredds/dodsC/models/fvcom/NECOFS/Archive/Seaplan_33_Hindcast_v1/catalog.html?dataset=models/NECOFS/Archive/Seaplan_33_Hindcast_v1/gom3_197801.nc

# opendap url
#http://www.smast.umassd.edu:8080/thredds/dodsC/models/fvcom/NECOFS/Archive/Seaplan_33_Hindcast_v1/gom3_197801.nc?time[0:1:744],zeta[0:1:744][0:1:48450],u[0:1:744][0:1:44][0:1:90414],v[0:1:744][0:1:44][0:1:90414],ua[0:1:744][0:1:90414],va[0:1:744][0:1:90414],uwind_stress[0:1:744][0:1:90414],vwind_stress[0:1:744][0:1:90414]

# sometimes fails 199407005959
# need to round to 1 hour
#tt=np.round(tRD*24.)/24.
#ti=datetime.fromordinal(int(tt))
#YEAR=str(ti.year)
#MO=str(ti.month).zfill(2)
#DA=str(ti.day).zfill(2)
#hr=(tt-int(tt))*24
#HR=str(int(np.round(hr))).zfill(2)            
#TS=YEAR+MO+DA+HR+'0000'
#FN2=TS

#T1=datetime.strftime(t1,'%Y%m%dT%H:%M') # '%Y%m%dT%H:%M'
YEARMOS=np.array([])
#for yr in range(1978,2017):
for yr in [2016]:
    for mo in range(1,1+12):
#    for mo in range(2,2+1):
        YEARMO=datetime.strftime(datetime(yr,mo,1),'%Y%m') # # '%Y%m%dT%H:%M'
        YEARMOS=np.append(YEARMOS,YEARMO)
YEARMOS=['201608']
for k in range(0,len(YEARMOS)):
    YEARMO=YEARMOS[k]
    print (k,YEARMO) 

    # YEARMO='197801' # test
    #http://www.smast.umassd.edu:8080/thredds/dodsC/models/fvcom/NECOFS/Archive/Seaplan_33_Hindcast_v1/gom3_197801.nc
    URL='http://www.smast.umassd.edu:8080/thredds/dodsC/models/fvcom/NECOFS/Archive/Seaplan_33_Hindcast_v1/gom3_'+YEARMO+'.nc'
    #xxx=ds['xxx']; np.save('gom3.xxx.npy',np.array(xxx))
    ds = Dataset(URL,'r').variables   # netCDF4 version
    time=np.array(ds['time'][:])
    lon=np.array(ds['lon'][:])
    lat=np.array(ds['lat'][:])
    lonc=np.array(ds['lonc'][:])
    latc=np.array(ds['latc'][:])
    #us=np.array(ds['u'][:,0,:]) # surface current
    #vs=np.array(ds['v'][:,0,:])
    #ub=np.array(ds['u'][:,-1,:]) # bottom current
    #vb=np.array(ds['v'][:,-1,:])
    ua=np.array(ds['ua'][:,:]) # depth avg current
    va=np.array(ds['va'][:,:])
    #uwind_stress=np.array(ds['uwind_stress'][:,:]) # wind stress
    #vwind_stress=np.array(ds['vwind_stress'][:,:])
    #zeta=np.array(ds['zeta'][:,:])
    #temps=np.array(ds['temp'][:,0,:]) # surface temperature
    #tempb=np.array(ds['temp'][:,-1,:]) # bottom temperature
    # All relevant data for dynamical analysis
    '''    
    D={'time':time,'lon':lon,'lat':lat,'lonc':lonc,'latc':latc,
       'zeta':zeta,'ua':ua,'va':va,'us':us,'vs':vs,'ub':ub,'vb':vb,
       'uwind_stress':uwind_stress,'vwind_stress':vwind_stress,
       'temps':temps,'tempb':tempb
       }
    '''
    D={'time':time,'lon':lon,'lat':lat,'lonc':lonc,'latc':latc,'ua':ua,'va':va}
    FNOUT='gom3_uv_'+YEARMO+'.npy'
    #FNOUT='gom3_uvwizte_'+YEARMO+'.npy'
    np.save(FNOUT,D)



