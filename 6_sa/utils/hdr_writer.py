#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
# writen by jpdarela 16-10-2016

# script to write header files to .bin files (bsq) with matrix data

#GLOBAL VARIABLES
nx = 120
ny = 160

linesep = '\r\n'

lizt = []

pixel_depht = 32

pixel_type  = 'FLOAT' # || SIGNINT || UNSIGNINT

bnd_layout = 'BSQ' # || BIL || BIP

NO_DATA = -9999.0

FILE_EXT = ['flt','bin','bsq']
# END GLOBAL

def catch_nt(input_file, nx, ny, pixel_depht):

    """Get the number of layers in input_file

    input_file = flat binary filename

    nx = (int) number of columns

    ny = (int) number of rows

    pixel_depth = (int) stride length in bits

    returns nt = number of layers stored in input_file"""

    image_size = (nx * ny * (pixel_depht / 8)) / 1024 # in bytes

    num_lay = (os.path.getsize(input_file)) / 1024 / image_size

    return int(num_lay)


def catch_data(input_file, layers, nx, ny):


    """Loads the input_file as a np.array once you know
    the number of layers in input_file

    input_file = flat binary filename (.bin extension mandatory)

    nx = (int) number of columns

    ny = (int) number of rows

    layers = (int) number of layers in input_file * ease with catch_nt()

    returns np.array shape(layers,nx,ny)"""
    if layers == 1:
        Bcount = nx * ny
        return np.fromfile(input_file, count=Bcount,
                    dtype=np.float32).reshape((ny,nx))
    else:
        Bcount = nx * ny * layers
        return np.fromfile(input_file, count=Bcount,
                    dtype=np.float32).reshape((layers,ny,nx))


def write_header(file_conn, NBANDS, nx=nx, ny=ny, xllcorner=-90.,
                yllcorner=-57.5,byteOrder='LSBFIRST'):

    """ Cria um cabeçalho.hdr nos padrões dos arquivos.flt """

    cellsize = 0.5

    write = ['NCOLS %i%s'%(nx, linesep),
            'NROWS %i%s'%(ny, linesep),
            'NBANDS %i%s'%(NBANDS, linesep),
            'NBITS %i%s'%(pixel_depht, linesep),
            'PIXELTYPE %s%s'%(pixel_type, linesep),
            'LAYOUT %s%s'%(bnd_layout, linesep),
            'XLLCORNER %f%s'%(xllcorner,linesep),
            'YLLCORNER %f%s'%(yllcorner,linesep),
            'CELLSIZE %f%s'%(cellsize,linesep),
            'NODATA_VALUE %f%s'%(NO_DATA,linesep),
            'BYTEORDER %s%s'%(byteOrder,linesep)
            ]
    with open(file_conn, 'w') as fh:
        for line in write:
            fh.write(line)


def main():
    bin_files_path = ['../outputs', '../outputs_pft', '../inputs', '../spinup' ] 
    for j in range(len(bin_files_path)):
        raw_list =[ i for i in os.listdir(bin_files_path[j]) if i.split('.')[-1] in FILE_EXT]
    
        for i in raw_list:
        
            path_in = os.path.join(bin_files_path[j],i)
            if len(i.split('.')) == 3:
                path_out = os.path.join(bin_files_path[j], (i.split('.')[0] + '.' + i.split('.')[1] + str('.img')))
            else:
                path_out = os.path.join(bin_files_path[j], (i.split('.')[0] + str('.img')))
            nlayers = catch_nt(path_in, nx, ny, pixel_depht)
            # write_header(path_out, nlayers, nx, ny )
            if nlayers == 1:
                dt = read_bin_ob(path_in)
            else:
                dt = read_bin_mb(path_in)
            array2raster(path_out, origin, 0.5, -0.5, nlayers, dt)
        
            print(nlayers)
            print(path_in)
            print(path_out, '\n')

if __name__ == '__main__':
    from array2raster import *
#    import gdal
    main()
