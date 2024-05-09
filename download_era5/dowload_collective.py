#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 19 17:37:28 2023
https://codes.ecmwf.int/grib/param-db/

@author: doan
"""


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os, sys, glob
from collections import namedtuple
import xarray as xr
import mpl_toolkits.basemap

#dates = pd.date_range('2010-01-01','2020-02-01', freq='MS')
dates = pd.date_range('2012-12-01','2014-01-01', freq='MS')

#==================================================================:::
# This script is to download ERA5 data used for HRLDAS downscaling :::
#==== monthly basic

if 1: # download era5

    print('download_x')
    import cdsapi
    c = cdsapi.Client()


    #dom = 'global'
    #area = '90/-180/-90/180' # NWSE
    
    dom = 'dbsh'
    area = '25/98/17/110' # NWSE
    
    by = '1' # each hour
    odir0 = 'era5_hrldas_new/' +dom
        
    n = 0
    for i, (d, d1) in enumerate( zip( dates[n:-1], dates[n+1:])):
        ed = d1 - pd.Timedelta('1 day')
        print(d)    
        date = '/'.join([ a.strftime('%Y%m%d') for a in pd.date_range(d,ed, freq='D')][:2])
        dateland = d.strftime('%Y%m%d') +'/'+ed.strftime('%Y%m%d')

        if 0:
            
            variables = ["surface_solar_radiation_downwards", 
                         "surface_thermal_radiation_downwards", 
                         "surface_pressure",  
                         "total_precipitation",  
                         'skin_temperature', 
                         'snow_depth',
                         'soil_temperature_level_1', 
                         'soil_temperature_level_2', 
                         'soil_temperature_level_3',  
                         'soil_temperature_level_4', 
                         'volumetric_soil_water_layer_1',  
                         'volumetric_soil_water_layer_2',  
                         'volumetric_soil_water_layer_3',  
                         'volumetric_soil_water_layer_4', 
                         '2m_dewpoint_temperature',
                         '2m_temperature'
                         ]
            
            odir = odir0 + '/land/'  
            if not os.path.exists(odir): os.makedirs(odir) 
            ofile = odir  + '/'+d.strftime('%Y%m.nc')
                
            print(ofile, dateland)
        
            c.retrieve(
            'reanalysis-era5-single-levels',
            {
                'product_type':'reanalysis',
                'format': 'netcdf', #dformat,
                'variable':variables, 
                'date':dateland,
                'area':area,
                'time':'00/to/23/by/'+by,
            },
            ofile )            
            
            
            
        if 1:
            variables_air = ['t', 'u', 'v', 'q']
            param = ['130', '131', '132', '133']
            
            odir = odir0 + '/air/'
            if not os.path.exists(odir): os.makedirs(odir) 
                
            levelist = '136'
            lev = levelist
            ofile = odir  + '/lev'+lev+'_'+d.strftime('%Y%m.nc')
                
            print(ofile, date)
            c.retrieve('reanalysis-era5-complete', { # Requests follow MARS syntax
                                                         # Keywords 'expver' and 'class' can be dropped. They are obsolete
                                                         # since their values are imposed by 'reanalysis-era5-complete'
                        'date'    : date, 
                        'levelist': levelist,               # 1 is top level, 137 the lowest model level in ERA5. Use '/' to separate values.
                        'levtype' : 'ml',
                        #'variable': variables_air,
                        'param'   : param,                 # Full information at https://apps.ecmwf.int/codes/grib/param-db/
                                                            # The native representation for temperature is spherical harmonics
                        'stream'  : 'oper',                  # Denotes ERA5. Ensemble members are selected by 'enda'
                        'time'    : '00/to/23/by/1',         # You can drop :00:00 and use MARS short-hand notation, instead of '00/06/12/18'
                        'type'    : 'an',
                        'grid'    : '0.1/0.1',               # Latitude/longitude. Default: spherical harmonics or reduced Gaussian grid
                        'format'  : 'netcdf',                # Output needs to be regular lat-lon, so only works in combination with 'grid'!
                        'area':area,
                    }, 
                           ofile)

















