import copy

for i in range(1,10):

	f = open(str(i)+'.txt','r')
	f_ = open(str(i)+'_.txt','w')
	f_.write(f.readline())

	string = ''

	for line in f:
		string += copy.copy(line[0:-2])
	
	f_.write( string )
	f.close()
	f_.close()