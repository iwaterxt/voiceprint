#!/bin/python


import sys

file = open(sys.argv[1],'r')
lines = file.readlines()
r = 0
w = 0
for line in lines:
    line_split = line.split() 
    if len(line_split) == 10:
	if line_split[0] == line_split[5]:
	    r = r + 1	
	else:
	    w = w + 1
file.close()
print "the EER is: ", float(r)/float(r+w)
