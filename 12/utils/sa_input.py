# cortar o nparray
# salvar en texto 
# converter em bin 
import os 
import numpy as np
import gdal

txts = []
ylow = -57.25
yupp = 22.75
xleft = -89.75
xright = -29.75
llcorner = (xleft, ylow)

#mask0 = gdal.Open('../inputs/lsmk.bin').ReadAsArray()
npp0 = gdal.Open('../inputs/npp.bin').ReadAsArray()
#rsds = gdal.Open('../inputs/rsds.bin').ReadAsArray()
#tas = gdal.Open('../inputs/tas.bin').ReadAsArray()
#pr = gdal.Open('../inputs/pr.bin').ReadAsArray()
#ps = gdal.Open('../inputs/ps.bin').ReadAsArray()
#hurs = gdal.Open('../inputs/hurs.bin').ReadAsArray()
#txt_file = 'la_mask.txt'
#np.savetxt(txt_file, mask0[135:295,180:300], fmt='%.10f')

def crop_array(arr, filename):
    global txts
    with open(filename, mode='a') as fh:
        index = 0
        for array in arr:
            cropped = arr[index,135:295,180:300]
            print(cropped.shape)
            txt_file = 'np_array_calc_avg_py.txt'
            np.savetxt(txt_file, cropped, fmt='%.10f')
            with open(txt_file, newline='\n') as fh1:
                reader = fh1.readlines()
            os.remove(txt_file)
            for line in reader:
                fh.write(line)
            index += 1
    txts.append(filename)
    return None
    
        
#crop_array(pr, 'la_pr.txt')
#crop_array(ps, 'la_ps.txt')
#crop_array(rsds,'la_rsds.txt')
#crop_array(tas,'la_tas.txt')
crop_array(npp0,'la_npp.txt')
#crop_array(hurs,'la_hurs.txt')
                
