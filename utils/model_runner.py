#!/usr/bin/env python3
# -*- coding: utf-8 -*-


# simple script to compile and run caete

import os
import time

fh = open('caete_exec_log.txt', mode='w')

curdir = os.getcwd()

os.chdir('../source1')
fh.write('compiling sources from ./source1\n')
os.system('gfortran -g -Wall env5_TR.f wbm4_TR.f productivity1.f  allocation1_pft.f -o env.out')
fh.write('\nrunning executable')
fh.write('\n\n\n\tExecutando caete- inicio em: %s\n' %time.ctime())
print('\n\n\n\tExecutando caete- inicio em', end='--->')
print(time.ctime())
init = time.time()

os.system('./env.out')

print('\n\n\n\tFinalizando caete - em', end='--->')
fh.write('\n\n\n\tFinalizando caete - em: %s' %time.ctime())
print(time.ctime())
end = time.time()
spend_time = end - init
fh.write('\n\nTempo de execução: %.0f minutos e %.2f segundos' %(spend_time // 60, spend_time % 60))
print('\n\nTempo de execução: %.0f minutos e %.2f segundos' %(spend_time // 60, spend_time % 60))
fh.close()
os.chdir(curdir)
os.system('./hdr_writer.py')
os.system('./prep_files.py')

