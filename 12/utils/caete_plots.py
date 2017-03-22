#!/usr/env python


import os
from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset as dts




def read_data_bin12(file_in):
    """FILE IN = caete netcdf result file """
    
    # prepare cell centers

    prefix = file_in.split('_')[0]
    if prefix == 'runom':
        prefix ='mrro'
    if prefix == 'wsoil':
        prefix ='mrso'
    if prefix == 'evaptr':
        prefix ='et'
#    if prefix == 'runom':
#        prefix ='mrro'
    data_name = 'annual_cycle_mean_of_' + prefix

    dt = dts(file_in, 'r')
    unit  = dt.variables[data_name].units
    lons = dt.variables['longitude'][:]
    lats =  dt.variables['latitude'][:]
    lons, lats = np.meshgrid(lons,lats)

    #
    dtt = np.fliplr(dt.variables[data_name][:,:,:])
    mean_dta = np.mean(dtt, axis=0,)

    return (mean_dta,(lons,lats), prefix, unit)

def plot_map(raster, meshed_coord, prefix, unit):
    lons,lats = meshed_coord
    # create figure, axes instances.
    fig = plt.figure()
    ax = fig.add_axes([0.05,0.05,0.9,0.9])
    # create Basemap instance.
    # coastlines not used, so resolution set to None to skip
    # continent processing (this speeds things up a bit)
    m = Basemap(projection='kav7',lon_0=0,resolution=None)
    # draw line around map projection limb.
    # color background of map projection region.
    # missing values over land will show up this color.
    m.drawmapboundary(fill_color='0.3')
    # plot sst, then ice with pcolor
    im1 = m.pcolormesh(lons,lats,np.flipud(raster),shading='flat',cmap=plt.cm.jet,latlon=True)
    # draw parallels and meridians, but don't bother labelling them.
    m.drawparallels(np.arange(-90.,99.,30.))
    m.drawmeridians(np.arange(-180.,180.,60.))
    # add colorbar
    cb = m.colorbar(im1,"bottom", size="5%", pad="2%")
    # add a title.
    ax.set_title('caete_annual_mean %s - units = %s' %(prefix, unit))
#    plt.show()

    fh = open(prefix + '.png', mode='w' )
    plt.savefig(fh, dpi=350)
    fh.close()
