#!/usr/bin/python

import os

nomes = ['1','2','3','4','5','6','7']

path_to_inputs = '../Inputs/DNA/'
os.system('mkdir -p Results/DNA')


for i in range(1,len(nomes)):

    print '##'+nomes[i-1]

    for j in range(i+1,len(nomes)+1):

        cmd = './mosal '+path_to_inputs+str(i)+'.txt '+path_to_inputs+str(j)+'.txt indels dpp -b=15 --no-traceback'
        cmd += ' ' + '> out.txt; ./check_bounds.py > Results/DNA/'+str(i)+'_'+str(j)+'.txt;'
        print cmd
        os.system(cmd)
