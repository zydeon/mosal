#!/usr/bin/python

import os

nomesRNA=['Banathema','HepatocystisCRNP11' ,'HepatocystisPP1' ,'Pjuxtanucleare', 'Tparva' ,'Pinfectans' ,'Chorrida' ,'Acatenella' ,'Atriacantha']
#nomesDNA=['CalbicansPAPa', 'CalbicansPAPalpha', 'CdubliniensisPAPa' ,'CtropicalisPAPalpha','CtropicalisPAPa' ,'PstipitisPAPalpha', 'ScerevisiaePAP']

#nomesDNA_=['1','2','3']#,'4','5','6','7']
nomesRNA_=['1','2','3','4','5','6','7','8','9']

nomes = nomesRNA_
total_scores = []

MIDs = [10] #,10,15,20]

path_to_inputs = '../../../../../../Inputs/RNA/'

os.system('mkdir -p Results/RNA')

for M in MIDs:

    print '\tMID = '+str(M)

    times = []
    total_scores = []  

    for i in range(1,len(nomes)):

        print '##'+nomesRNA[i-1]

        for j in range(i+1,len(nomes)+1):

            cmd = './prog '+path_to_inputs+str(i)+'.txt '+path_to_inputs+str(j)+'.txt ' + '../../../../../../Tables/BLOSUM85.tbl ' + str(M)
            cmd += ' ' + '> Results/RNA/'+nomesRNA[i-1]+'_'+nomesRNA[j-1]
            print cmd
            os.system(cmd)


    print '\n'