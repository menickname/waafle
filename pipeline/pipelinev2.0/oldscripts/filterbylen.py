#!/usr/bin/python

'''
This script will filter all BLAST hits by a length filter.
'''

#Import
import argparse
import hgtmodules

#Define arguments
parser = argparse.ArgumentParser()
parser.add_argument('--groupmerged', help = 'Location and file of merged group BLAST output.')
parser.add_argument('--length', type = float, help = 'Groups below this length will be removed.')
args = parser.parse_args()

#Organize BLAST output
dictContigHits = hgtmodules.contigHits(open(args.groupmerged))
ddictContigGroupHits = hgtmodules.contigGroupHits(dictContigHits)	

#For each group in each contig, find the smallest and largest start/end sites among the alignments
for contig in ddictContigGroupHits.iterkeys():
        for group in ddictContigGroupHits[contig].iterkeys():
                startlist = []
                endlist = []
                for i in range(len(ddictContigGroupHits[contig][group])):
                        startlist.append(float(ddictContigGroupHits[contig][group][i][6]))
                        endlist.append(float(ddictContigGroupHits[contig][group][i][7]))
                minstart, maxend = min(startlist), max(endlist) 
        	grouplen = maxend - minstart
		if grouplen < args.length:
			del ddictContigGroupHits[contig][group]

#Renumber the groups before reprinting the output
for contig2 in ddictContigGroupHits.iterkeys():
	totalgroups = len(ddictContigGroupHits[contig2].iterkeys()) #This is the number of new groups there should be.
	startlist, endlist = [], []
	for group2 in ddictContigGroupHits.iterkeys():
		startlist.append(float(ddictContigGroupHits[contig2][group2][6]))
		endlist.append(float(ddictContigGroupHits[contig2][group2][7]))
	minStart = min(startlist)
	maxEnd = max(endlist)
	
	
			for j in range(len(ddictContigGroupHits[contig][group])):
				print '\t'.join(ddictContigGroupHits[contig][group][j])
