import math

for m in range(1,14):

	for k in range(1,6):

		f = open(str(m*300)+'/'+str(m*300)+'_'+str(k),'r')

		header = f.readline()
		seq = f.readline()

		f.close()

		f = open(str(m*300)+'/'+str(m*300)+'_'+str(k),'w')
		f.write(header)

		for i in range( len(seq)/80 ) :
			f.write( seq[ i*80 : (i+1)*80 ] + '\n' )

		f.write(seq[ len(seq)/80*80 : ] + '\n')

		f.close()