import random

A = ['A','G','C','T']

N = ['300','600','900','1200','1500','1800','2100','2400','2700','3000','3300','3600','3900' ]
seq_per_N = 5

for n in range(len(N)):
	seqs = [ '' for i in range( seq_per_N ) ]
	
	for k in range( seq_per_N ):

		for i in range( int(N[n]) ):
			seqs[k] += A[ random.randint(0,len(A)-1) ]

		file = open( 'Alpha_4/'+N[n]+'/'+N[n]+'_'+str(k+1) , 'w')
		file.write('>'+N[n]+'_'+str(k+1)+'\n')
		file.write(seqs[k])
		file.close()