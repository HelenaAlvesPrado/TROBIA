#!usr/bin/env python3

import os
import random as rd
import numpy as np
import gdal
import matplotlib.pyplot as plt


def random_npp():
    """CREATE A RANDOM NPP GLOBAL MAP... AS AN np.array - 0.5°resolution """
    mask = np.load('mask3.npy')[0]
    rnpp = np.zeros(shape=(360,720),dtype=np.float32)
    
    for j in range(rnpp.shape[0]):
        for i in range(rnpp.shape[1]):
            if mask[j][i]:
                rnpp[j][i] = -9999.0
            else:
                if j <= 60 or j >= 270:
                    rnpp[j][i] = rd.random()*(rd.random() + rd.randrange(0,1))
                elif j > 60 or j < 270:
                    if j < 140 or j > 210 :
                        rnpp[j][i] = rd.random() * (rd.random() + rd.randrange(0,2))
                    else:
                        rnpp[j][i] = rd.random() * (rd.random() + rd.randrange(0,3))
    return rnpp


# Valores que vão sair de uma distribuição lognormal
#g1 = np.linspace(1.6, 7.1, 10)
#vcmax = np.linspace(3.e-5,25e-5,10)
#jmax = np.linspace(1e-4,3e-4,10)

#tleaf = np.linspace(1,100,50)/12 # years
#twood = np.linspace(0,80,90)
#troot = np.linspace(1,100,50)/12

# estes valores combinados devem somar 100% 

#aleaf = np.linspace(25,90,33)
#aroot = np.linspace(25,90,33)
#awood = np.linspace(0,90,45)
g1 = np.linspace(1.6, 7.1, 10)
vcmax = np.linspace(3.e-5,25e-5,10)
jmax = np.linspace(1e-4,3e-4,10)

tleaf = np.arange(1,100,12)/12 # years
twood = np.arange(1,80,5)
troot = np.arange(1,100,12)/12

aleaf = np.arange(20,81,5)
aroot = np.arange(20,81,5)
awood = np.arange(20,81,5)

pls_list = []
colnames_a = ['aleaf','awood','aroot']


pls_grass1 = [[a/100,0.0,c/100] for a in aleaf for c in aroot if abs(a + 0.0 + c) == 100.]
pls_grass2 = [[c/100,0.0,a/100] for a in aleaf for c in aroot if abs(c + 0.0 + a) == 100.]
pls_grass3 = [[a/100,0.0,c/100] for c in aroot for a in aleaf if abs(a + 0.0 + c) == 100.]
pls_grass4 = [[c/100,0.0,a/100] for c in aroot for a in aleaf if abs(c + 0.0 + a) == 100.]

pls_woody1 = [[a/100,b/100,c/100] for a in aleaf for b in awood for c in aroot if ((a + b + c) == 100.) and (b > 19)]
pls_woody2 = [[a/100,b/100,c/100] for c in aroot for b in awood for a in aleaf if ((c + b + a) == 100.) and (b > 19)]
pls_woody3 = [[a/100,b/100,c/100] for b in awood for c in aroot for a in aleaf if ((c + b + a) == 100.) and (b > 19)]
pls_woody4 = [[a/100,b/100,c/100] for a in aleaf for c in aroot for b in awood if ((c + b + a) == 100.) and (b > 19)]
pls_woody5 = [[a/100,b/100,c/100] for c in aroot for b in awood for a in aleaf if ((c + b + a) == 100.) and (b > 19)]
pls_woody6 = [[a/100,b/100,c/100] for b in awood for a in aleaf for c in aroot if ((c + b + a) == 100.) and (b > 19)]
plsa =(pls_grass1 + pls_grass2 + pls_grass3 + pls_grass4 + pls_woody1 +\
      pls_woody2 +pls_woody3+pls_woody4+pls_woody5+pls_woody6)

plsa_wood =(pls_woody1 + pls_woody2 + pls_woody3 + pls_woody4 + pls_woody5 + pls_woody6)

plsa_grass =(pls_grass1 + pls_grass2 + pls_grass3 + pls_grass4)


# CREATING ALLOCATION COMBINATIONS
for i in range(len(plsa_grass)):
    x = plsa_grass.pop()
    if x in plsa_grass:
        pass
    else:
        plsa_grass.insert(0,x)
        
for i in range(len(plsa_wood)):
    x = plsa_wood.pop()
    if x in plsa_wood:
        pass
    else:
        plsa_wood.insert(0,x)

# CREATING TURNOVER COMBINATIONS
colnames_t = ['tleaf','twood','troot']
turnover_wood = [[a,b,c] for a in tleaf for b in twood for c in troot]
turnover_grass = [[a,0.0,c] for a in tleaf for c in troot]
turnover = turnover_grass + turnover_wood

for i in range(len(turnover_grass)):
    x = turnover_grass.pop()
    if x in turnover_grass:
        pass
    else:
        turnover_grass.insert(0,x)

for i in range(len(turnover_wood)):
    x = turnover_wood.pop()
    if x in turnover_wood:
        pass
    else:
        turnover_wood.insert(0,x)
        
# CREATING PHYSIOLOGICAL COMBINATIONS
colenames_p = ['g1','vcmax','jmax']
phys = [[a,b,c] for a in g1 for b in vcmax for c in jmax]

sec_hand_wood = [a + b  for a in turnover_wood for b in phys]
sec_hand_grass = [a + b  for a in turnover_grass for b in phys]

sec_hand_arr_grass = np.array(sec_hand_grass)
sec_hand_arr_wood = np.array(sec_hand_wood)

# juntando as combinações de turnover + g1 + vcmax etc temos
# mais de 1300000 possíveis combinações

# selecionando randomicamente (usando uma distribuição discreta uniforme) 10000
plss_wood = sec_hand_arr_wood[np.random.random_integers(0,sec_hand_arr_wood.shape[0],10000)][:]
plss_grass = sec_hand_arr_grass[np.random.random_integers(0,sec_hand_arr_grass.shape[0],10000)][:]

pls_list = []

for alloc_pls in plsa_grass:
    for x in range(10):
        plst = alloc_pls + list(sec_hand_arr_grass[np.random.randint(0,10000)][:])
        pls_list.append(plst)
    
for alloc_pls in plsa_wood:
    for x in range(10):
        plst = alloc_pls + list(sec_hand_arr_wood[np.random.randint(0,10000)][:])
        pls_list.append(plst)

out_arr = np.array(pls_list).T

np.savetxt('pls_580.txt', out_arr, fmt='%.12f')



