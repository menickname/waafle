#!/usr/bin/python

'''
This script is to get the coordinates of the synthetic contigs so we can compare them against the contigs from BLAST.
'''

# Import
import sys
import re
# Define arguments

# Read in answerkey
for astrline in open(sys.argv[1]):
	aastrline = astrline.strip().split('\t')
	contigname, genelist = aastrline[0], aastrline[7]
	listgenes = genelist.split(', [')
	for gene in listgenes:
		mygene = gene.replace(']', '').replace('[', '').split(',')
		status, start, end = mygene[0], mygene[5], mygene[6]
		print contigname, status.strip("'"), start.strip(), end.strip()
	
