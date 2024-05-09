#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 15:10:21 2022

@author: doan
"""
import glob, sys, os
import numpy as np
import pandas as pd
import xarray as xr
import metpy.calc as mpcalc
from metpy.units import units
import  mpl_toolkits.basemap  
import matplotlib.pyplot as plt
#=====






#==============================================================================#
def extract_geo_em(geo_em_file, gxdir):
    '''
    This function is to extract necessary variables from geo_em.d0?.nc file
    Variables include: LAI, VEGFRA, ...
    
    Parameters
    ----------
    geo_em_file : link to geo_em file
    gxdir : link to output directory
    '''
    dg = xr.open_dataset(geo_em_file)
    xlat = dg.XLAT_M.values[0]
    xlon1 = dg.XLONG_M[0]
    xlon = xlon1.where(xlon1>0, xlon1+360).values                  
        
    
    vv1 = {"LAI366": {'units': ' ', 'geoname': 'LAI12M'} ,       
                 'VEGFRA366':{'unites': '', 'geoname': 'GREENFRAC' }}
    
    
    for i, var in enumerate(list(vv1.keys())[:]):
    
        v = vv1[var]['geoname']
        
        lg = xr.Dataset()
        x = dg[v][0].values
        x = np.append( x, [x[0]], axis=0)        
        x = xr.DataArray(x , dims=('t', 'y', 'x'), coords= {'t': ('t', pd.date_range( '2010', periods=13,freq='MS' ) )   })
        x = x.resample(t = '1D').interpolate(kind="cubic")
        if v == 'GREENFRAC': x = xr.where(x<=0 ,0.01, x) * 100
        lg[var] = x #.sel(t = time)
    
        lg.coords['t'] = ( ('t'), range(366) )
        lg.coords['lat'] = (  ('y', 'x'), xlat  )
        lg.coords['lon'] = (  ('y', 'x'), xlon  )
    
    
        lg.to_netcdf(gxdir + '/'+var+'.nc')
    #=====================================
    
    
    #=====================================
    vv2 = {
            # from geo_em file
      "XLAT": {'units': 'degree_north', 'geoname': 'XLAT_M'} , 
      "XLONG": {'units': 'degree_east', 'geoname': 'XLONG_M'} , 
      "HGT": {'units': 'm', 'geoname': 'HGT_M'} , 
      "MAPFAC_MX": {'units': '', 'geoname': 'MAPFAC_MX'} , 
      "MAPFAC_MY": {'units': '', 'geoname': 'MAPFAC_MY'} , 
      "IVGTYP": {'units': '', 'geoname': 'LU_INDEX'} , 
      # edit from geo_em file
      "TMN": {'units': 'K', 'geoname': 'SOILTEMP'} ,  # adjust to elevation
      "SHDMAX": {'units': '%', 'geoname': 'GREENFRAC'} ,  # max(100*GREENFRAC) 
      "SHDMIN": {'units': '%', 'geoname': 'GREENFRAC'} ,  # min(100*GREENFRAC) 
      #"LAI": {'units': 'm^2/m^2', 'geoname': 'LAI12M'} , # LAI12M after interpolated
      "XLAND": {'units': '', 'geoname': 'LU_INDEX'} ,  # if LU_INDEX==iswater or islake,2; else 1.  
      "ISLTYP": {'units': '', 'geoname': 'SOILCTOP'} ,  
      }
    
    iswater = int(dg.attrs['ISWATER'])
    islake = int(dg.attrs['ISLAKE'])
    issoilwater = int(dg.attrs['ISOILWATER'])
    
    for i, var in enumerate(list(vv2.keys())[:]):
        print('\p', var)    
        v = vv2[var]['geoname']

        do = xr.Dataset()
        
        if var in ['XLAT', 'XLONG', 'HGT', "MAPFAC_MX", "MAPFAC_MY", "IVGTYP"]: 
            dat = dg[v].values

        if var == 'TMN':
            dat = dg[v].values  # - 0.0065 * dg['HGT_M'].values
            #data_var = xr.where((dg.LU_INDEX==iswater)|(dg.LU_INDEX==islake), -1.e36, data_var).values

        if var == 'SHDMAX': 
            dat = dg['GREENFRAC'].max(axis=1).values*100
        
        if var == 'SHDMIN': 
            dat = dg['GREENFRAC'].min(axis=1).values*100
            
        # if LU_INDEX==iswater or islake,2; else 1.  
        if var == 'XLAND':
            LU_data = dg[v]
            dat = xr.where((LU_data==iswater)|(LU_data==islake), 2, 1).values
            xland = dat
        
        elif var == 'ISLTYP':
            dominant_index = dg['SOILCTOP'].argmax(dim='soil_cat') + 1
            dominant_value = dg['SOILCTOP'].max(dim='soil_cat')
            dominant_index_corrected = xr.where(dominant_value<0.01, 8, dominant_index)
            dat = xr.where(xland ==2, issoilwater, dominant_index_corrected)
            dat = xr.where((xland !=2)&(dat==14), 8, dat).values                       


        do[var] = ( ( 'y', 'x'), dat[0]  )
        do.coords['lat'] = (  ('y', 'x'), xlat  )
        do.coords['lon'] = (  ('y', 'x'), xlon  )
        do.attrs = dg.attrs
        
        do.to_netcdf( gxdir + '/'+var+'.nc')        
    return 
#==================================================================================





#==============================================================================
def bdcon(era_link, gxdir, sdate, edate, hidir, dom='1'):
    '''
    This function is to generate boundary conditions for hrldas program, based 
    ERA5 data
    
    Parameters
    ----------
    era_link : need to specify link to ERA5
    gxdir : link to geo_em parameters
    sdate : start date (YYYYMMDD HH)
    edate : end date
    hidir : hrldas input directory
    dom : optional, if geo_em.d?.nc is other than d01, need to specify the dom number  

    '''
    
    xlat = xr.open_dataset(gxdir+'/XLAT.nc')['XLAT'].values
    xlon = xr.open_dataset(gxdir+'/XLONG.nc')['XLONG'].values
    
    
    print('****')
    vv3 = { 
        'T2D':    {  'eraname':   't', 'units': 'K'},          
        'Q2D':    {  'eraname':   'q', 'units': 'kg/kg'},  
        'U2D':    {  'eraname':   'u', 'units': 'm/s'}, 
        'V2D':    {  'eraname':   'v', 'units': 'm/s'},
        'PSFC':   {  'eraname':   'sp', 'units': 'Pa'},
        'LWDOWN':  {  'eraname':   'strd', 'units': 'W/m^2'},
        'SWDOWN':  {  'eraname':   'ssrd', 'units': 'W/m^2'}, 
        'RAINRATE':  {  'eraname': 'tp', 'units'  : 'kg/m^2/s'},
        'LAI': {'geoname':'LAI366', 'units':'-'},
        'VEGFRA': {'geoname':'VEGFRA366', 'units':'-'},
        }
    
    
    for time in pd.date_range(sdate, edate, freq= 'H')[:]:
        print(time)
        
        do = xr.Dataset()
        for i, var in enumerate( list(vv3.keys())[:] ):
            print(var)
            if not var in ['LAI', 'VEGFRA']:
                v_era = vv3[var]['eraname']
                unit = vv3[var]['units']
                
                idir = era_link + '/*/' + v_era+ '/'
                lev = '*'
                ifile = idir + v_era + '_lev'+lev+'_'+time.strftime('%Y%m.nc')        
    
                print(glob.glob(ifile))
            
                ds = xr.open_dataset(glob.glob(ifile)[0])        
                d1 = ds[v_era].sel(time=time, method='nearest')[::-1].values
                raw_lat, raw_lon = ds.latitude.values[::-1], ds.longitude.values
                d2 = mpl_toolkits.basemap.interp(d1, raw_lon, raw_lat, xlon, xlat, checkbounds=False, masked=False, order=1)
           
        
                if var in ['LWDOWN', 'SWDOWN']: d2 = d2 / 3600
                if var in ['RAINRATE']: d2 = d2 / 3600 * 1000                
       
            else:
                print('****')
                v_geo = vv3[var]['geoname']
                unit = vv3[var]['units']
                ifile = gxdir + '/' + v_geo+ '.nc'
                ds = xr.open_dataset(ifile)[v_geo]
                d2 = ds.loc[ time.day_of_year - 1].values
      
            do[var ] = (('Time','south_north','west_east'), [d2])
            do[var].attrs['units'] = unit      
            
        ofile = hidir+"/"+time.strftime('%Y%m%d%H')+'.LDASIN_DOMAIN'+dom
        encoding=[{var: {'_FillValue': None}} for var in do.variables]                
        do.to_netcdf( ofile, encoding=encoding[0])
        print(ofile)                
         
    return 
#==============================================================================





#==============================================================================
#==============================================================================
def bdcon_mpi(era_link, gxdir, sdate, edate, hidir, dom='1'):
    '''
    This function is to generate boundary conditions for hrldas program, based 
    ERA5 data
    
    Parameters
    ----------
    era_link : need to specify link to ERA5
    gxdir : link to geo_em parameters
    sdate : start date (YYYYMMDD HH)
    edate : end date
    hidir : hrldas input directory
    dom : optional, if geo_em.d?.nc is other than d01, need to specify the dom number  

    '''
    
    from mpi4py import MPI
        
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()
        
        
    xlat = xr.open_dataset(gxdir+'/XLAT.nc')['XLAT'].values
    xlon = xr.open_dataset(gxdir+'/XLONG.nc')['XLONG'].values
    
    
    print('****')
    vv3 = { 
        'T2D':    {  'eraname':   't', 'units': 'K'},          
        'Q2D':    {  'eraname':   'q', 'units': 'kg/kg'},  
        'U2D':    {  'eraname':   'u', 'units': 'm/s'}, 
        'V2D':    {  'eraname':   'v', 'units': 'm/s'},
        'PSFC':   {  'eraname':   'sp', 'units': 'Pa'},
        'LWDOWN':  {  'eraname':   'strd', 'units': 'W/m^2'},
        'SWDOWN':  {  'eraname':   'ssrd', 'units': 'W/m^2'}, 
        'RAINRATE':  {  'eraname': 'tp', 'units'  : 'kg/m^2/s'},
        'LAI': {'geoname':'LAI366', 'units':'-'},
        'VEGFRA': {'geoname':'VEGFRA366', 'units':'-'},
        }
    
    
    for impi, time in enumerate(pd.date_range(sdate, edate, freq= 'H')[:]):
        
        if np.mod(impi, size) == rank:
            print(rank, size, impi, np.mod(impi, size))
            
            print(time)
            
            do = xr.Dataset()
            for i, var in enumerate( list(vv3.keys())[:] ):
                print(var)
                if not var in ['LAI', 'VEGFRA']:
                    v_era = vv3[var]['eraname']
                    unit = vv3[var]['units']
                    
                    idir = era_link + '/*/' + v_era+ '/'
                    lev = '*'
                    ifile = idir + v_era + '_lev'+lev+'_'+time.strftime('%Y%m.nc')        
        
                    print( glob.glob(ifile) )
                
                    ds = xr.open_dataset(glob.glob(ifile)[0])        
                    d1 = ds[v_era].sel(time=time, method='nearest')[::-1].values
                    raw_lat, raw_lon = ds.latitude.values[::-1], ds.longitude.values
                    d2 = mpl_toolkits.basemap.interp(d1, raw_lon, raw_lat, xlon, xlat, checkbounds=False, masked=False, order=1)
               
            
                    if var in ['LWDOWN', 'SWDOWN']: d2 = d2 / 3600
                    if var in ['RAINRATE']: d2 = d2 / 3600 * 1000                
           
                else:
                    print('****')
                    v_geo = vv3[var]['geoname']
                    unit = vv3[var]['units']
                    ifile = gxdir + '/' + v_geo+ '.nc'
                    ds = xr.open_dataset(ifile)[v_geo]
                    d2 = ds.loc[ time.day_of_year - 1].values
          
                do[var ] = (('Time','south_north','west_east'), [d2])
                do[var].attrs['units'] = unit      
                
            ofile = hidir+"/"+time.strftime('%Y%m%d%H')+'.LDASIN_DOMAIN'+dom
            encoding=[{var: {'_FillValue': None}} for var in do.variables]                
            do.to_netcdf( ofile, encoding=encoding[0])
            print(ofile)                
#==============================================================================





#==============================================================================
def setupfile(era_link,gxdir, stdir, sdate):
    '''
    This function to prepare data to build setup file for hrldas run

    Parameters
    ----------
    era_link : link to ERA5 data
    gxdir : link to geo_em parameters
    stdir : link to directory to save setup variables
    sdate : state date (YYYYMMDD HH)

    '''
    
    xlat = xr.open_dataset(gxdir+'/XLAT.nc')['XLAT'].values
    xlon = xr.open_dataset(gxdir+'/XLONG.nc')['XLONG'].values
            
    time = pd.to_datetime(sdate)
    lev = '*'
    for iv, varlist in enumerate( [ ['skt'], 
                     ['swvl1', 'swvl2', 'swvl3', 'swvl4'], 
                     ['stl1', 'stl2', 'stl3', 'stl4'], 
                     ['sde'] 
                    ]):
        print(iv, varlist)
        dat = []
        vvv = ['TSK','SMOIS', 'TSLB', 'SNOW']
        for var in varlist:
            print(var)
            idir = era_link + '/*/' + var + '/'
            ifile = idir + var + '_lev'+lev+'_'+time.strftime('%Y%m.nc')
            print(ifile)
            ifile_0 = glob.glob(ifile)[0]
            
            raw_data_file = xr.open_dataset( ifile_0 )
            
            raw_lat, raw_lon = raw_data_file.latitude.values[::-1], raw_data_file.longitude.values       
            if var == 'sde': var = 'sd'
            data_var = raw_data_file[var][0][::-1].values #********
            data_var_interpolated = mpl_toolkits.basemap.interp(data_var, 
                                                                raw_lon, raw_lat, 
                                                                xlon, xlat,  
                                                                checkbounds=False, masked=False, order=1)
            dat.append(data_var_interpolated)     
        do = xr.Dataset()
        do[ vvv[iv] ] = xr.DataArray(dat) #( dat  )
        do.to_netcdf(stdir+'/'+vvv[iv]+'.nc')
#==============================================================================                






#==============================================================================
def writesetup(gxdir, stdir, sdate, hidir, dom='1'):
    '''
    This function is to create HRLDAS setup file

    Parameters
    ----------
    gxdir : link to geo_em variables
    stdir : link to setup variables 
    sdate : state date
    hidir : link to hrldas input directory
    dom : Toptional
          The default is '1'.

    '''
    vv4 = {
        "Times": {'units':''} ,
        # from geo_em file
        "XLAT": {'units': 'degree_north', 'geoname': 'XLAT_M'} , 
        "XLONG": {'units': 'degree_east', 'geoname': 'XLONG_M'} , 
        "HGT": {'units': 'm', 'geoname': 'HGT_M'} , 
        "MAPFAC_MX": {'units': '', 'geoname': 'MAPFAC_MX'} , 
        "MAPFAC_MY": {'units': '', 'geoname': 'MAPFAC_MY'} , 
        "IVGTYP": {'units': '', 'geoname': 'LU_INDEX'} , 
        # edit from geo_em file
        "TMN": {'units': 'K', 'geoname': 'SOILTEMP'} ,  # adjust to elevation
        "SHDMAX": {'units': '%'} ,  # max(100*GREENFRAC) 
        "SHDMIN": {'units': '%'} ,  # min(100*GREENFRAC) 
        "LAI": {'units': 'm^2/m^2'} , # LAI12M after interpolated
        "XLAND": {'units': '', 'geoname': 'LU_INDEX'} ,  # if LU_INDEX==iswater or islake,2; else 1.  
        "ISLTYP": {'units': '', 'geoname': 'SOILCTOP'} ,
        # from raw data file
        "TSK": {'units': 'K'} ,  # skin temperature from ERA5 
        "TSLB": {'units': 'K'} ,  # soil layer temp 
        "SMOIS": {'units': 'm^3/m^3'},  # layer volumetric total water content [m3/m3] !!!
        "DZS": {'units': 'm'} ,  # each soil layer depth
        "ZS": {'units': 'm'} ,   # soil layer 
        "SNOW": {'units': 'kg/m^2'} , #snow depth
        #'SNODEP':{'units': 'kg/m^2'} , #snow depth
        # add
        "SEAICE": {'units': ''} ,       # sea ice fraction (=0 for a land point)
        "CANWAT": {'units': 'kg/m^2'} , # set CANWAT = 0
    } 

    
    time = pd.to_datetime(sdate)
    

    do = xr.Dataset()
    for i, var in enumerate(list(vv4.keys())[:]):
        print('****  '+var)   
        unit = vv4[var]['units']
        dims = ('Time', 'south_north', 'west_east')
        dim4 = ('Time', 'soil_layers_stag', 'south_north', 'west_east')
        dim2 = ('Time', 'soil_layers_stag')        
        
        if var == 'Times': 
            dat, dims =  [pd.to_datetime(sdate).strftime('%Y-%m-%d_%H:%M:%S')], ( 'Time' )        
        
        if var in ['XLAT', 'XLONG', 'HGT', "MAPFAC_MX", "MAPFAC_MY", "IVGTYP", 'TMN', 'SHDMAX', 'SHDMIN', 'XLAND', 'ISLTYP']: 
            ifile = gxdir + '/'+var+'.nc' 
            ds = xr.open_dataset(ifile)
            dat = np.expand_dims(ds[var].values, axis=0)
            attrs = ds.attrs
            shape = dat.shape
            
            
        if var == 'LAI':
            ifile = gxdir + '/LAI366.nc' 
            ds = xr.open_dataset(ifile)
            dat = np.expand_dims(ds['LAI366'].sel(t =time.day_of_year).values, axis=0)
            
        ####################
        # from raw data file
        ####################
        if var in ['TSK', 'TSLB', 'SMOIS', 'SNOW']:  
            ifile = stdir + '/'+var+'.nc' 
            x = xr.open_dataset(ifile)[var]
            print(x.shape)
            dat = x.values
            if var == 'SNOW': dat = dat*1000
            if x.shape[0] == 4: dat, dims = [dat], dim4
        
        if var == 'ZS': 
            dat, dims = [ [0.035, 0.175, 0.64, 1.945] ], dim2

        if var == 'DZS': 
            dat, dims = [ [0.07 , 0.21 , 0.72, 1.89 ] ], dim2

        
        elif var == 'SEAICE':
            dat = np.zeros(shape)

        elif var == 'CANWAT': 
            dat = np.zeros(shape)
            
            
            
        print(dims, np.array(dat).shape)        
        do[var] = ( dims, dat )

        do[var].attrs['units'] = unit
    
    
    do = do.fillna({'SNOW': -999})
    do = do.fillna({'SMOIS': -999})
    do.attrs = attrs
    print(attrs)
    
    ofile = hidir+'/HRLDAS_setup_'+pd.to_datetime(sdate).strftime('%Y%m%d')+'00_d'+dom
    do.to_netcdf(ofile)        
    return do
#==============================================================================


    


#==============================================================================
def copyprog(hrldas_prog, rdir):
    '''
    Prepare hrldas program

    Parameters
    ----------
    hrldas_prog : link to compile hrldas program
    rdir : link to run directory
    '''
    os.system('cp '+hrldas_prog+'/run/hrldas.exe '+hrldas_prog+'/run/*.TBL  '+rdir)
    os.system('cp '+hrldas_prog+'/run/examples/GLDAS/namelist.hrldas.GLDAS '+ rdir + '/namelist.hrldas')


def writenamelist(namelist_link, sdate, edate, hidir,hodir1, rdir, update, dom='1'):
    '''
    Modiry namelist

    Parameters
    ----------
    namelist_link : link to standard namelist file
    sdate : state date 
    edate : end date
    hidir : link to hrldas input directory
    hodir1 : link to hrldas output directory
    rdir : link to run directory
    update : updated scheme option (under dictionary type)
    dom : optional
        The default is '1'.
    '''    
    sd, ed = pd.to_datetime(sdate), pd.to_datetime(edate)
    su_file = os.path.abspath(hidir + 'HRLDAS_setup_'+pd.to_datetime(sdate).strftime('%Y%m%d')+'00_d'+dom)
    
    if not 'kday' in list(update.keys()):
        nuday = (ed - sd ).days + 1
    else:
        nuday = update['kday']
    
    
    
    for v1, v2 in zip( ['hrldas_setup_file','indir', 'outdir', 'start_year', 'start_month', 'start_day','start_hour','kday'  ], 
                      ['"'+ su_file +'"', '"'+os.path.abspath(hidir)+'"', '"'+os.path.abspath(hodir1)+'"', sd.strftime('%Y'),  sd.strftime('%m'), sd.strftime('%d'), sd.strftime('%H'), str(nuday)]):
        update[v1] = v2
                    
    ll = open(namelist_link).readlines()
    ll = [l for l in ll if l.strip() != '' if l.strip()[0]!='!' ]
    # namelist: https://github.com/NCAR/hrldas/blob/master/hrldas/run/README.namelist
    
    
    
    
    nn = {
        
        #'hrldas_setup_file':'"'+ su_file +'"',
        #'indir':'"'+hidir+'"',
        #'outdir':'"'+hodir1+'"',
        #'start_year':sd.strftime('%Y'),
        #'start_month':sd.strftime('%m'),
        #'start_day':sd.strftime('%d'),
        #'start_hour':sd.strftime('%H'),
        'start_min':'00',
        #'kday':str(nuday),
        'khour': '0',
        'spinup_loops':'0',
        'forcing_name_pr': '"RAINRATE"',
        #=======================
        'dynamic_veg_option':                 '1',
        'canopy_stomatal_resistance_option':  '1',
        'btr_option':                         '1',
        'surface_runoff_option':              '3',
        'dvic_infiltration_option':           '1',
        'surface_drag_option':                '1',
        'frozen_soil_option':                 '1',
        'supercooled_water_option':           '1',
        'radiative_transfer_option':          '3',
        'snow_albedo_option':                 '1',
        'pcp_partition_option':               '1',
        'tbot_option':                        '2',
        'temp_time_scheme_option':            '1',
        'glacier_option':                     '1',
        'surface_resistance_option':'1',
        'soil_data_option':'1',
        'pedotransfer_option':'1',
        'crop_option':'0',
        'irrigation_option':'0',
        'irrigation_method':'0',
        'tile_drainage_option':'0',
        #=========================
        'sf_urban_physics':'0',
        'use_wudapt_lcz':'0',
        'forcing_timestep':'3600',
        'noah_timestep':'1800',
        'output_timestep':'3600',
        'split_output_count':'1',
        'skip_first_output':'.false.',
        'restart_frequency_hours':'0',
        'nsoil':'4',
        'soil_thick_input(1)':'0.07',
        'soil_thick_input(2)':'0.21',
        'soil_thick_input(3)':'0.72',
        'soil_thick_input(4)':'1.89',
        'zlvl':'30.0',
        }

    for n1 in list(set(update.keys()) or set(nn.keys())): nn[n1] = update[n1]
    
    for i, l in enumerate(ll):
        for n1, n2 in nn.items():
            
            if l.lower().split('=')[0].strip() == n1: 
                #print(n1,n2)
                ll[i] =  n1 + '=' + n2 + '\n'
              
    if 'khour' in list(update.keys()):
        ikday = [i for i, l in enumerate(ll) if l.lower().split('=')[0].strip() == 'kday'][0]
        ll = ll[:ikday] + ['khour = '+update['khour']+'\n'] + ll[ikday:]
    
    ofile = rdir + '/namelist.hrldas'
    print(nn)
    print(ofile)
    open(ofile, 'w' ).write( ''.join(ll))
    return ll
#==============================================================================







#==============================================================================
def extract_hrldas(geo_em_file, hidir, hodir1, exdir1, sdate, edate, varnames = ['TEMP'], dom='1'):
    '''
    Extract hrldas output

    Parameters
    ----------
    geo_em_file : link to gem_em file
    hidir : link to hrldas input file
    hodir1 : link to hrldas output files
    exdir1 : link to extract directory
    sdate : state date
    edate : end date
    varnames : list of variables for now only TEMP and RH
    dom : optional
        DESCRIPTION. The default is '1'.
    '''    
    def filt(x): return x.where(x>=0)
    def filt2(x): return x.where(x<=1)
    
    #su_file = hidir + 'HRLDAS_setup_'+pd.to_datetime(sdate).strftime('%Y%m%d')+'00_d1'          
    #g = xr.open_dataset(su_file)
    urbfrac = xr.open_dataset(geo_em_file).LANDUSEF[0,12].values
    
    from metpy.calc import relative_humidity_from_mixing_ratio
    from metpy.units import units
    
    for date in pd.date_range(sdate,edate, freq='H')[:]:
        
        try:
        
            time = date #- pd.Timedelta('7 hours')
            print(time) 
    
            #if1 = hidir+"/"+time.strftime('%Y%m%d%H')+'.LDASIN_DOMAIN'+dom
            if2 = hodir1+"/"+time.strftime('%Y%m%d%H')+'.LDASOUT_DOMAIN'+dom
            print(if2)
            #di = xr.open_dataset(if1), 
            ds = xr.open_dataset(if2)
            urbscheme = 'TH2_URB2D' in list(ds.variables.keys())
            if urbscheme: print('        urban scheme is used')
            for var in varnames[:]:
                
                do = xr.Dataset()
                if var in ['TEMP', 'RH']:
                    x1 = (ds['T2MV']*ds['FVEG']+ds['T2MB']*(1-ds['FVEG'])) - 273.15
                    
                    if urbscheme:
                        x2 = ds['TH2_URB2D'] - 273.15
                        a1 = urbfrac 
                        #a1 = ds['FRC_URB']
                        isurban = ds.IVGTYP.values == 13
                        x3 = (x2*a1 + x1*(1-a1)) 
                        temp = np.where(isurban,x3, x1)
                        print('        urban scheme is used')
                    else:
                        temp = x1.values
                    
                    if var == 'TEMP': 
                        do[var] =  ( ('Time',  'south_north', 'west_east'), temp )
                    
                    if var == 'RH':
                        m1 = ds['Q2MV']*ds['FVEG']+ds['Q2MB']*(1-ds['FVEG'])
                        m2 = 0
                        #if urbscheme:
                        #    #a1 = ds['FRC_URB']
                        #    mir = (m2*a1 + m1*(1-a1))
                        #else:
                        mir = m1
                        p = ds.FORCPLSM
                        rh = relative_humidity_from_mixing_ratio(p.values * units.Pa, temp * units.degC, mir.values)*100 #.to('percent')
                        do[var] =  ( ('Time',  'south_north', 'west_east'), rh )
                    
                exdir2 = mdir(exdir1 + '/' + var + '/')
                ofile = exdir2 + os.path.basename(if2) + '.nc'
                print('    ', ofile)
                do.to_netcdf(ofile)
        except:
            print('something wring')
#==============================================================================
def mdir(odir): 
    if not os.path.exists(odir): os.makedirs(odir)
    return odir
#==============================================================================












if __name__ == "__main__":
# global path:
    
    era_link = '../download_era5/era5_hrldas/dbsh/'
    hrldas_prog = '/mnt/work/working/HRLDAS/uHRLDAS/hrldas/'
    hrldas_prog = '/Users/doan/working/HRLDAS/20240509/hrldas/hrldas/'

    #hrldas_prog = '/Users/doan/Documents/GitHub/imhen_agri_dbsh/run_hrldas/hrldas_github/hrldas/hrldas/'

    sdate, edate = '2013-01-01', '2013-01-01 05:00:00'
    geo_em_file = 'geo_em.d03.nc'
    
    odir00 = 'dataout_cordex/'
    
    rdir0  = odir00 + '/era5/'
    #schemes = {'set1':{'dynamic_veg_option': '2', 'crop_option':'1'},
    #            'set2':{'dynamic_veg_option': '1', 'crop_option':'1'},
    #            'set3':{'dynamic_veg_option': '4', 'crop_option':'1'},
    #            'set4':{'dynamic_veg_option': '2', 'crop_option':'0'}}
    
    schemes = {}
    for idyn in range(1, 10):
        print(idyn)
        schemes['set'+str(idyn)] = { 'dynamic_veg_option': str(idyn)} 


    gxdir, stdir, hidir, hodir, exdir = [ mdir(rdir0 + a + '/') for a in ['geo_em', 'setup', 'datain', 'dataout', 'extract']]

    if 1:    
        extract_geo_em(geo_em_file, gxdir)
        setupfile(era_link,gxdir, stdir, sdate)
        writesetup(gxdir, stdir, sdate, hidir, dom='3')
        bdcon(era_link, gxdir, sdate, edate, hidir, dom='3')
            
    
        for iss, sset  in enumerate(list(schemes.keys())[:1]):
                
            print(sset)     
            rdir = mdir(hodir+'/'+sset) 
            hodir1 = mdir(hodir+'/'+sset+'/dataout/') 
            exdir1 = mdir(exdir+'/'+sset) 
                        
            update = schemes[sset]
            namelist_link = hrldas_prog+'/run/examples/GLDAS/namelist.hrldas.GLDAS'
            
            update['khour'] = '5'
            update['kday'] = '0'
            
            copyprog(hrldas_prog, rdir )
            ll = writenamelist(namelist_link, sdate, edate, hidir, hodir1, rdir, update, dom='3')       
            os.system('source /home/doan/env_wrf.sh ; cd  '+rdir+'; ./hrldas.exe; cd -')    
            
            varnames = ['TEMP', 'RH'][:]
            extract_hrldas(geo_em_file, hidir, hodir1, exdir1, sdate, edate, varnames, dom='3' )
        
    if 0:
        xx = []
        for iss, sset  in enumerate(list(schemes.keys())[:]):
                
            print(sset)     
            rdir = hodir+'/'+sset 
            hodir1 = hodir+'/'+sset+'/dataout/'
            exdir1 = exdir+'/'+sset 
                                
            f = exdir1 + '/TEMP/2013010105.LDASOUT_DOMAIN3.nc'
            x = xr.open_dataset(f).TEMP[0]
            xx.append(x)
            
            
            
            
    




















