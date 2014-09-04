#!/usr/bin/python

'''
This script should create groups based on length of contig. Before the script is run, the BLAST output should be sorted first by contig name, second by length, and third by bit score. 

This script differs from 'aggregatebylen_revised.py' in that it does not merge groups as it runs. Instead, the script aggregatebylen_comboverlaps.py should be run to combine overlapping groups at the end. 
'''

#Import
import argparse
import sys
import re
import collections


#Define arguments
parser = argparse.ArgumentParser()
parser.add_argument('blast_output', help = 'Location and file of sorted BLAST output')
parser.add_argument('overlap', type = float, help = 'Percentage of overlap between two reads')
args = parser.parse_args()


#Create a dictionary which has the contigs as keys, and all the BLAST output as a value (key is in list format)
dictContigHits = {}
for line in open(args.blast_output):
        aItems = line.strip().split('\t')
        dictContigHits.setdefault(aItems[0],[]).append(aItems)


#Define function for calculating overlap
def calcOverlap(iStart1, iEnd1, iStart2, iEnd2):
	listcoord = [float(iStart1), float(iEnd1), float(iStart2), float(iEnd2)]
	listcoord.sort()
	contig1len = iEnd1 - iStart1 + 1
	contig2len = iEnd2 - iStart2 + 1
	divisor = min(contig1len, contig2len)
	overlap = (listcoord[2]-listcoord[1] + 1)/divisor
	return overlap


#Group by length. Anything that overlaps by more than 'overlap' will be grouped together.
groups = {}

for contig, alignList in dictContigHits.iteritems():
	count = 0

	for align in alignList: #For each alignment in this contig, record the start & end sites
		#print align
		
		count = count + 1
		start, end = float(align[6]), float(align[7])
		counter = 0
		overlapcounter = 0
		overlaplist = []		

		if count == 1: #For the first alignment, add it to a new group.
			groupname = 'Group1'
			groups.setdefault(groupname, [])
			groups[groupname].append(align)
		else: #For all other alignments, see if it goes in any of the groups.
			for group in groups.iterkeys():
				gstart = float(groups[group][0][6]) #Record the group start/end sites.
				gend  = float(groups[group][0][7])
				
				if gstart-end > 0 or start - gend > 0: #If the two do not overlap at all, count number of times this happens
					counter = counter + 1

				else: #If there is any overlap, calculate the %overlap.
					overlap = calcOverlap(gstart, gend, start, end)
					
					if overlap > args.overlap: #Add to group if it overlaps with one, then continue to compare groups. 
						groups[group].append(align)
						#overlaplist.append(group) 
						#overlapcounter = overlapcounter + 1
					else: #If they do not overlap by enough, count and continue to the next group. 
						counter = counter + 1
		if counter == len(groups): #If we've bypassed all the groups, and none overlap or overlap by enough, add it as a new group.
			groupnum = len(groups) + 1
			groupname = 'Group' + str(groupnum)
			groups.setdefault(groupname, [])
			groups[groupname].append(align)
	
	#Print groups without sorting groups
	for groupname in groups.iterkeys():
		for i in range(len(groups[groupname])):
			linelist = groups[groupname][i][:]
			linelist.append(groupname)
			#print '\t'.join(str(linelist[j]) for j in range(len(linelist)))
			
	#Find out order of groups by start site
	startlist = []
	for groupname in groups.iterkeys():
		startf = float(groups[groupname][0][6]), groupname
		startlist.append(startf)
	startlist.sort()
	#print startlist
	

	#Reorder groups before printing, and add in the new group number along with genus/species designation
	newgroupnum = 0
	for startcoord, grouplist in startlist:
		newgroupnum = newgroupnum + 1
		for i in range(len(groups[grouplist])):
			linelist = groups[grouplist][i][:] #This makes it so that group is no longer a reference
			nnewgroupnum = 'Group' + str(newgroupnum)
			linelist.append(nnewgroupnum)
			print '\t'.join(str(linelist[j]) for j in range(len(linelist)))	 		
	
	groups = {}
	
