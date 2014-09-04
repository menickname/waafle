#!/usr/bin/python

'''
This script should merge the groups that overlap by the same threshold as aggregatebylen.py. It will:
1. For each contig, it will take the minimum start and maximum end of each group.
2. It will then look at the overlap of each group.
3. Groups that need to be merged will then be merged.
4. Output new BLAST formatting.
'''

#Import
import argparse
import re
import hgtmodules

#Define arguments 
parser = argparse.ArgumentParser()
parser.add_argument('--grouped', help = 'Location and file of grouped BLAST output')
parser.add_argument('--overlap', type = float, help = 'Minimum overlap to merge groups')
args = parser.parse_args()

#Organize BLAST output
dictContigHits = hgtmodules.contigHits(open(args.grouped))
ddictContigGroupHits = hgtmodules.contigGroupHits(dictContigHits)

#For each group in each contig, find the smallest and largest start/end sites among the alignments
dictContigGroupSorted = hgtmodules.sortGroups(ddictContigGroupHits)
print dictContigGroupSorted
mod_ddictContigGroupHits = {}
for contig in dictContigGroupSorted.keys():
	grouplist = []
	for i in range(len(dictContigGroupSorted[contig])):	
		group = dictContigGroupSorted[contig][i]
		oldname, newname = group[0], 'Group'+str(i+1)
		start, end = float(group[1]), float(group[2])
		if i == 0: #For the first group
			startf, endf = start, end
			print contig, oldname, newname
			grouplist.append(newname)
		else: #For the next group, see if it overlaps
			prevgroup = dictContigGroupSorted[contig][i-1]
			overlap = hgtmodules.calcOverlap(startf, endf, start, end)
			if overlap > args.overlap: #If there is enough overlap, merge the groups
				if len(grouplist[len(grouplist)-1]) > 1:
					print grouplist
					grouplist[len(grouplist)-1].append(newname)
				else: 
					grouplist.append(['Group'+str(i), newname])
				print contig, prevgroup[0], 'Group'+str(i), oldname, newname, '0'	
			else: #Otherwise,
				grouplist.append(newname)
				print contig, oldname, newname, 'notoverlapping'
		print grouplist
''' 
mod_ddictContigGroupHits = {}
for contig in ddictContigGroupHits.iterkeys():
	overlapgroups = []
	for j in range(len(grouplist_sorted)):
		start = grouplist_sorted[j][1]
		end = grouplist_sorted[j][2]
		if j == 0:
			startf, endf = start, end
		else:
			overlap = hgtmodules.calcOverlap(startf, endf, start, end)
			if overlap == None: #If there is no overlap
				pass	
			else: #If there is overlap
				if overlap > args.overlap: #Keep a list of lists for the groups that should be merged. We cannot merge them now in case there a >2 groups that need to be merged in a row. 
					fgroup = grouplist_sorted[j-1][0]
					sgroup = grouplist_sorted[j][0]
					if len(overlapgroups) == 0:
						newline = [fgroup, sgroup] 
						overlapgroups.append(newline)
					elif overlapgroups[len(overlapgroups)-1][1] == fgroup: #If previous group should be merged to current group
						overlapgroups[len(overlapgroups)-1].append(sgroup)
					else:
						overlapgroups.append([fgroup, sgroup])
		startf, endf = start, end

	#After going through every group, merge the ones that should be merged. 
	for groupings in overlapgroups:
		mergedGroupName = 'merged' 
		mergedAlignments = [] #Merge those groups
		for k in range(len(groupings)):
			mergedGroupName += groupings[k] #Create the new group names
			mergedAlignments += ddictContigGroupHits[contig][groupings[k]] #Adding list of lists to create a new set of BLAST hits
			del ddictContigGroupHits[contig][groupings[k]] #Delete the original groups
		for l in range(len(mergedAlignments)):
			mergedAlignments[l][14] = mergedGroupName #Replace the group annotations
		ddictContigGroupHits[contig][mergedGroupName] = mergedAlignments #Add to the contig
	
	
	#Now that this contig's groups has been fixed, sort the groups.
	groupnumlist = []
	for groups in ddictContigGroupHits[contig].iterkeys():
		groupnum = re.search(r'[0-9]+', groups).group()
		groupnumber = [float(groupnum), groups]
		groupnumlist.append(groupnumber)
	groupnumlist.sort(key = lambda y: y[0])

	#Print the results.	
	for number, groupname in groupnumlist:
		alignments = ddictContigGroupHits[contig][groupname]
		for m in range(len(alignments)):
			newLine = alignments[m]
			print '\t'.join(str(newLine[n]) for n in range(len(newLine)))	
	'''
