import os
import glob

os.chdir('../outputs')
files = os.listdir()
files = glob.glob1(os.getcwd(), '*.bin')

files_to  = [fl.split('.')[0] + '.tif' for fl in files]

z = zip(files, files_to)

for x in z:
    os.system('gdal_translate %s %s' %x) 
