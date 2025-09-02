#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  1 06:29:15 2025

@author: user
"""
#!conda install -c conda-forge opendrift #this takes a few minutes
#!conda install -c conda-forge gdal
from osgeo import gdal
from osgeo import ogr
from osgeo import osr

from datetime import datetime,timedelta
from opendrift.models.oceandrift import OceanDrift
from opendrift.readers import reader_ROMS_native

lat=[42.0,42.5]  # default position can be changed as you wish (41.8 for CC Bay)
lon=[-70.,-69.75]

o = OceanDrift(loglevel=50) #suppressed nearly all log to screen
datet=datetime.utcnow()-timedelta(days=3) # these single DOPPIO files appear online a few days prior to now
mth=str(datet.month).zfill(2);day=str(datet.day).zfill(2); # added leading zeros to mth and day
url='https://tds.marine.rutgers.edu/thredds/dodsC/roms/doppio/2017_da/his/files/doppio_his_2024'+mth+day+'_0000_0001.nc' # most recent output
reader_doppio = reader_ROMS_native.Reader(url)
o.add_reader([reader_doppio])
o.seed_elements(lon=lon, lat=lat,time=reader_doppio.start_time)
o.run()
o.plot()
