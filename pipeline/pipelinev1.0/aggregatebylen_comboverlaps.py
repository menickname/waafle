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
import sys
import re
import collections

#Define arguments 
parser = argparse.ArgumentParser()
parser.add_argument('blast_output', help = 'Location and file of BLAST output from aggregate.py')
parser.add_argument('overlap', type = float, help = 'Minimum overlap to merge groups')
args = parser.parse_args()


#Define function for calculating overlap
def calcOverlap(iStart1, iEnd1, iStart2, iEnd2):
        listcoord = [float(iStart1), float(iEnd1), float(iStart2), float(iEnd2)]
        listcoord.sort()
        contig1len = iEnd1 - iStart1 + 1
        contig2len = iEnd2 - iStart2 + 1
        divisor = min(contig1len, contig2len)
        overlap = (listcoord[2]-listcoord[1] + 1)/divisor
        return overlap


#Create a dictionary of contigs and its alignments
dictContigHits = {}
for astrline in open(args.blast_output):
	aItems = astrline.strip().split('\t')
	dictContigHits.setdefault(aItems[0], []).append(aItems)


#Create a dictionary of dictionaries, with format {contig: {group1:[alignments], group2:[alignments]}, contig2...}
ddictContigGroups = {}
for bstrline in dictContigHits.iterkeys():
	dictGroupHits = {}
	for i in range(len(dictContigHits[bstrline])):
		group = dictContigHits[bstrline][i][15]
		dictGroupHits.setdefault(group,[]).append(dictContigHits[bstrline][i])
	ddictContigGroups.setdefault(bstrline, {}).update(dictGroupHits)


#For each group in each contig, find the smallest and largest start/end sites among the alignments
for contig in ddictContigGroups.iterkeys():
	group_startendlist = []
	for group in ddictContigGroups[contig].iterkeys():
		startlist = []
		endlist = []
		for i in range(len(ddictContigGroups[contig][group])):
			startlist.append(float(ddictContigGroups[contig][group][i][7]))
			endlist.append(float(ddictContigGroups[contig][group][i][8]))
		group_startendlist.append([group, min(startlist), max(endlist)])
	group_startendlist.sort(key = lambda x: float(x[1])) #Sort the group start end sites by start site


	#Look for overlap between groups
	overlapgroups = []
	for j in range(len(group_startendlist)):
		#print j, len(group_startendlist)
		start = group_startendlist[j][1]
		end = group_startendlist[j][2]
		#print contig, group_startendlist[j][0], start, end
		if j == 0:
			startf, endf = start, end
		else:
			#print j, startf, endf, start, end, group_startendlist[j][0], contig
			if start - endf > 0 or startf - end > 0: #If the two do not overlap at all
				pass
			else:
				overlap = calcOverlap(startf, endf, start, end)
				if overlap > args.overlap: #Keep a list of lists for the groups that should be merged. We cannot merge them now in case there a >2 groups that need to be merged in a row. 
					fgroup = group_startendlist[j-1][0]
					sgroup = group_startendlist[j][0]
					#print startf, endf, start, end, fgroup, sgroup, contig
					if len(overlapgroups) == 0:
						newline = [fgroup, sgroup]
						overlapgroups.append(newline)
					elif overlapgroups[len(overlapgroups)-1][1] == fgroup:
						overlapgroups[len(overlapgroups)-1].append(sgroup)
					else:
						overlapgroups.append([fgroup, sgroup])
					mergedGroup = 'merged' + fgroup + sgroup
				else:
					pass
		startf, endf = start, end
	

	#After going through every group, merge the ones that should be merged. 
	for groupings in overlapgroups:
		mergedGroupName = 'merged' #Create the new group names
		mergedAlignments = [] #Merge those groups
		for k in range(len(groupings)):
			mergedGroupName = mergedGroupName + groupings[k]
			mergedAlignments = mergedAlignments + ddictContigGroups[contig][groupings[k]]
			#print contig, mergedAlignments
			del ddictContigGroups[contig][groupings[k]] #Delete the original groups
		for l in range(len(mergedAlignments)):
			mergedAlignments[l][15] = mergedGroupName #Replace the group annotations
		ddictContigGroups[contig][mergedGroupName] = mergedAlignments #Add to the contig
		#print mergedGroupName

	
	#Print out results. Now that this contig's groups has been fixed, sort the groups.
	groupnumlist = []
	for groups in ddictContigGroups[contig].iterkeys():
		groupnum = re.search(r'[0-9]+', groups).group()
		groupnumber = [float(groupnum), groups]
		groupnumlist.append(groupnumber)
	groupnumlist.sort(key = lambda y: y[0])
	
	for number, groupname in groupnumlist:
		alignments = ddictContigGroups[contig][groupname]
		#print alignments
		for m in range(len(alignments)):
			newLine = alignments[m]
			print '\t'.join(str(newLine[n]) for n in range(len(newLine)))
		

