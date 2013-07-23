#!/usr/bin/python
#
fout = open('out.txt', 'r')
fbounds = open('bounds.txt','r')

points = [ [int(x.split()[0]), int(x.split()[1])] for x in fout.readlines()]
bounds = [ [int(y.split()[0]), int(y.split()[1])] for y in fbounds.readlines()]

errors = False

for b in bounds:
	count = points.count(b)
	if count > 1:
		errors = True
		print 'Point (' + str(b[0]) + ','+str(b[1])+') appears ' + count + ' times!'
	elif count == 0:
		errors = True
		print 'Point (' + str(b[0]) + ','+str(b[1])+') does not appear in output'

if errors == False:
	print '\nBounds are correct!'
else:
	print '\nBounds incorrect!'	