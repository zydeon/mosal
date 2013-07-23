import math

for k in range(1,8):

	f = open(str(k)+'.txt','r')

	header = f.readline()
	seq = f.readline()

	f.close()

	f = open(str(k)+'.txt','w')
	f.write(header)

	for i in range( len(seq)/80 ) :
		f.write( seq[ i*80 : (i+1)*80 ] + '\n' )

	f.write(seq[ len(seq)/80*80 : ] + '\n')

	f.close()