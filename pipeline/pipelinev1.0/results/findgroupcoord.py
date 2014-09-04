#!/usr/bin/python

#Import
import argparse
import re
import json
import collections
import copy

#Define arguments
parser = argparse.ArgumentParser()
parser.add_argument('groupedfile', help = 'Location and file for grouped file')
args = parser.parse_args()

firstline = open(args.groupedfile).readline().split('\t')
fcontig, fgroup = firstline[0], firstline[14]
startlist = []
endlist = []
startlist.append(float(firstline[6]))
endlist.append(float(firstline[7]))

#Read in file
for astrline in open(args.groupedfile):
	aastrline = astrline.strip().split('\t')
	scontig, sgroup = aastrline[0], aastrline[14]
	start, end = float(aastrline[6]), float(aastrline[7])
	#print fcontig, fgroup, scontig, sgroup, start, end
	if fcontig == scontig:
		if fgroup == sgroup:
			#print 'adding'
			startlist.append(start)
			endlist.append(end)
		else:
			newstart = min(startlist)
			newend = max(endlist)
			startlist = []
			endlist = []
			startlist.append(start)
			endlist.append(end)
			print fcontig, fgroup, newstart, newend, newend-newstart+1
		fcontig = scontig
		fgroup = sgroup
	else:
		newstart = min(startlist)
		newend = max(endlist)
		print fcontig, fgroup, newstart, newend, newend-newstart+1
		startlist = []
		endlist = []
		startlist.append(start)
		endlist.append(end)
		fcontig = scontig
		fgroup = sgroup
		
			
	
