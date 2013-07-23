#!/usr/bin/python

import os

files = os.listdir('DP/Results/RNA')

for f in files:
	cmd = 'diff '+'DP/Results/RNA/'+f+ ' DP_prune/Results/RNA/'+f
	print cmd
	os.system(cmd)