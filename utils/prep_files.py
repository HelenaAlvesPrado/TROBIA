#!/usr/bin/env python

import os
import glob
import numpy as np
import gdal
from netCDF4 import Dataset as dt


mask_fpath = './mask12.npy'
NO_DATA = [-9999.0, -9999.0]
lsmk = np.load(mask_fpath)

fdir = ['../outputs_pft', '../outputs']



def flt_attrs():
    return {'header'  : ['long_name',                 'unit',           'standart_name'],
            'rsds'    : ['short_wav_rad_down',        'W m-2',                   'rsds'], 
            'wind'    : ['wind_velocity',             'm s-1',                   'wind'],
            'ps'      : ['sur_pressure',              'Pa',                        'ps'],
            'tas'     : ['sur_temperature_2m',        'celcius',                  'tas'],  
            'pr'      : ['precipitation',             'Kg m-2 month-1',            'pr'],      
            'wsoil'   : ['soil_water_content-wsoil',  'kg m-2',                  'mrso'],
            'evaptr'  : ['evapotranpiration',         'kg m-2 day-1',              'et'],
            'runom'   : ['total_runoff',              'kg m-2 day-1',            'mrro'],
            'ar'      : ['autothrophic respiration',  'kg m-2 year-1',             'ar'],
            'ph'      : ['photosynthesis',            'kg m-2 year-1',             'ph'],
            'npp'     : ['net primary productivity',  'kg m-2 year-1',            'npp'],
            'lai'     : ['Leaf Area Index',           'm-2 m-2',                  'LAI'],
            'rcm'     : ['stomatal resistence',       's m-1',                    'rcm'],
            'hr'      : ['heterothrophic respiration','kg m-2 year-1',             'hr'],
            'clit'    : ['Litter Carbon',             'Kg m-2',                  'clit'],
            'csoil'   : ['Soil Carbon',               'Kg m-2',                 'csoil'],
            'rm'      : ['maintenance respiration',   'kg m-2 year-1',             'rm'],
            'rms'     : ['sapwood rm',                'kg m-2 year-1',            'rms'],
            'rmf'     : ['fine root rm',              'kg m-2 year-1',            'rmf'],
            'rml'     : ['leaf rm',                   'kg m-2 year-1',            'rml'],
            'rg'      : ['growth respiration',        'kg m-2 year-1',             'rg'],
            'rgs'     : ['sapwood rg',                'kg m-2 year-1',            'rgs'],
            'rgf'     : ['fine root rg',              'kg m-2 year-1',            'rgf'],
            'rgl'     : ['leaf rg',                   'kg m-2 year-1',            'rgl'],
            'cawood'  : ['C in abovewgrownd wood',    'kg m-2',                'cawood'],
            'cfroot'  : ['C in fine roots',           'kg m-2',                'cfroot'],
            'cleaf'   : ['C in leaves',               'kg m-2',                 'cleaf'],
          'cawood_ini': ['init C in abovewgrownd wood','kg m-2',               'cawood'],
          'cfroot_ini': ['init C in fine roots',       'kg m-2',               'cfroot'],
          'cleaf_ini' : ['init C in leaves',           'kg m-2 ',               'cleaf']}

def read_raster(fpath):
    """Returns the raster file in fpath as a masked numpy array
    with shape as raster file
    formal arguments:

    fpath :: string: complete file path

    """
    ds = gdal.Open(fpath).ReadAsArray()
    nlayers = ds.shape[0]
    lsmk_internal = lsmk[0:nlayers,:,:]
    np.place(arr=ds, mask=lsmk_internal, vals=NO_DATA)
    return ds


def write_CAETE_output(nc_filename, arr, var):

    t, la, lo,  = arr.shape
    nlayers = arr.shape[0]
    lsmk_internal = lsmk[0:nlayers,:,:]
    # create netcdf file
    rootgrp = dt(nc_filename, mode='w', format='NETCDF3_CLASSIC')

    #dimensions
    rootgrp.createDimension("time", None)
    rootgrp.createDimension("latitude", la)
    rootgrp.createDimension("longitude", lo)


    #variables
    time      = rootgrp.createVariable(varname="time", datatype=np.float32, dimensions=("time",))
    latitude  = rootgrp.createVariable(varname="latitude", datatype=np.float32,dimensions=("latitude",))
    longitude = rootgrp.createVariable(varname="longitude", datatype=np.float32, dimensions=("longitude",))
    var_      = rootgrp.createVariable(varname = 'annual_cycle_mean_of_'+str(flt_attrs()[var][2]), datatype=np.float32,
                                       dimensions=("time","latitude","longitude",),
                                       fill_value=NO_DATA[0])

    #attributes
    ## rootgrp
    rootgrp.description = flt_attrs()[var][0] + " from CAETE_1981-2010--> annual cycle"
    rootgrp.source = "CAETE model outputs"
    ## time
    time.units = "days since 1850-01-01 00:00:00.0"
    time.calendar = "noleap"
    time.axis='T'
    ## lat
    latitude.units = u"degrees_north"
    latitude.long_name=u"latitude"
    latitude.standart_name =u"latitude"
    latitude.axis = u'Y'
    ## lon
    longitude.units = "degrees_east"
    longitude.long_name = "longitude"
    longitude.standart_name = "longitude"
    longitude.axis = 'X'
    ## var
    var_.long_name=flt_attrs()[var][0]
    var_.units = flt_attrs()[var][1]
    var_.standard_name=flt_attrs()[var][2]
    var_.missing_value=NO_DATA[0]

    ## WRITING DATA
    times_fill = np.array([15.5, 45., 74.5, 105., 135.5, 166.,
                           196.5, 227.5, 258., 288.5, 319., 349.5])
    
    time[:] = times_fill[0:nlayers]
    longitude[:] = np.arange(-179.75, 180, 0.5)
    latitude[:] =  np.arange(-89.75, 90, 0.5)

    var_[:,:,:] = np.fliplr(np.ma.masked_array(arr, lsmk_internal))
    rootgrp.close()


for folder in range(len(fdir)):
    files = glob.glob1(fdir[folder], '*.bin')
    for fl in range(len(files)):
        fpath = fdir[folder] + os.sep + files[fl]
        varname = files[fl].split('.')[0]
        if len(files[fl].split('.'))  == 3:
            caete_name = fdir[folder] + os.sep + varname + '_pft_' + files[fl].split('.')[1] + '_' + 'annual_cycle_mean_CAETE.nc'
        else:
            caete_name = fdir[folder] + os.sep + varname + '_' + '1981_2010_annual_cycle_mean_CAETE.nc'
        if varname in flt_attrs().keys():
            arr = read_raster(fpath)
            write_CAETE_output(caete_name, arr, varname)

