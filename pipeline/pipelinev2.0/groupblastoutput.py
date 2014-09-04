#!/usr/bin/python

'''
This script should create groups based on length of contig. Before the script is run, the BLAST output should be sorted first by contig name, second by length, and third by bit score. 

This script differs from 'aggregatebylen_revised.py' in that it does not merge groups as it runs. Instead, the script aggregatebylen_comboverlaps.py should be run to combine overlapping groups at the end. 
'''

#Import
import argparse
import hgtmodules

#Define arguments
parser = argparse.ArgumentParser()
parser.add_argument('--blastoutput', help = 'Location and file of sorted BLAST output')
parser.add_argument('--overlap', type = float, help = 'Percentage of overlap between two reads')
args = parser.parse_args()


#Create a dictionary which has the contigs as keys, and all the BLAST output as a value (key is in list format)
dictContigHits = hgtmodules.contigHits(open(args.blastoutput))

#Group by length. Anything that overlaps by more than 'overlap' will be grouped together.
ddictContigGroupHits = {}
for contig in dictContigHits.iterkeys():
	count = 0
	groups = {}
	for i in range(len(dictContigHits[contig])): #For each hit in this contig, record the start & end sites
		hit = dictContigHits[contig][i]
		#print hit
		start, end = float(hit[6]), float(hit[7])
		counter = 0
		if i == 0: #For the first hit, add it to a new group.
			groupname = 'Group1'
			groups.setdefault(groupname, [])
			groups[groupname].append(hit)
		else: #For all other hit, see if it goes in any of the groups.
			for group in groups.iterkeys():
				gstart = float(groups[group][0][6]) #Record the group start/end sites.
				gend  = float(groups[group][0][7])
				overlap = hgtmodules.calcOverlap(gstart, gend, start, end)
				#print group, gstart, gend, overlap
				if overlap == None: #If the two do not overlap at all, count number of times this happens
					counter += 1
				else: #If there is any overlap, calculate the %overlap.	
					if overlap > args.overlap: #Add to group if it overlaps with one, then continue to compare groups. 
						groups[group].append(hit)
						#print group, 'append b/c overlap'
					else: #If they do not overlap by enough, count and continue to the next group. 
						counter += 1
		if counter == len(groups): #If we've bypassed all the groups, and none overlap or overlap by enough, add it as a new group.
			groupnum = len(groups) + 1
			groupname = 'Group' + str(groupnum)
			groups.setdefault(groupname, [])
			groups[groupname].append(hit)
			#print groupname, 'append b/c end + \n'
				
	ddictContigGroupHits[contig] = groups

hgtmodules.printBlast_output(ddictContigGroupHits)
'''			
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
'''	
