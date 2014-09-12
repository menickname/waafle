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

mod_ddictContigGroupHits = {}
for contig in dictContigGroupSorted.keys():
	grouplist = []
	for i in range(len(dictContigGroupSorted[contig])):	
		currgroup = dictContigGroupSorted[contig][i]
		curroldname, currnewname = currgroup[0], 'Group'+str(i+1)
		currstart, currend = float(currgroup[1]), float(currgroup[2])
		if i == 0: #For the first group, assign the current names
			info = curroldname, currnewname
			grouplist.append(info)
			#print contig, info, 'first'
		else: #For all other groups, see if they overlap with the previous group
			prevgroup = dictContigGroupSorted[contig][i-1]
			prevoldname, prevnewname = prevgroup[0], 'Group'+str(i)
			overlap = hgtmodules.calcOverlap(prevstart, prevend, currstart, currend)
			if overlap > args.overlap: #If there is enough overlap, merge the groups.
				if len(grouplist[len(grouplist)-1][1].split(',')) == 1: #If previous group is not merged, create a new merged group.	
					oldmergedname = prevoldname + ',' + curroldname
					newmergedname = prevnewname + ',' + currnewname
				elif len(grouplist[len(grouplist)-1][1].split(',')) > 1: #If the previous group is merged, check for multiple merges.
					if prevnewname in grouplist[len(grouplist)-1][1].split(','): #Merge 3+ groups if the previous merged groups are consecutive to the current one.
						oldmergedname = grouplist[len(grouplist)-1][0] + ',' + curroldname
						newmergedname = grouplist[len(grouplist)-1][1] + ',' + currnewname
					else: #Add as a new merged group if the previous merged groups are NOT consecutive.
						oldmergedname = prevoldname + ',' + curroldname
                                        	newmergedname = prevnewname + ',' + currnewname
				info = oldmergedname, newmergedname
				grouplist.pop() #Since we are adding a merged group, the previous entry should be deleted.
				grouplist.append(info)
				#print contig, info, 'overlap'
			else: #If the groups do not overlap, assign the current names
				info = curroldname, currnewname
				grouplist.append(info)
				#print contig, info, 'notoverlapping'
		prevstart, prevend = currstart, currend #Assign the current groups as previous groups, so that we can continue comparisons


	#With the newly ordered groups, create a new dictionary that has the correct group names
	GroupHits = {}
	for oldnames, newnames in grouplist:
		if len(oldnames.split(',')) == 1:
			hits = ddictContigGroupHits[contig][oldnames]
			for j in range(len(hits)):
				#print hits[j], 'old'
				hits[j][14] = newnames
				#print hits[j], 'new', '\n'
			GroupHits[newnames] = hits
		else:
			mergednewname = 'merged' #Create the merged name first
			for k in range(len(newnames.split(','))):
				mergednewname += newnames.split(',')[k]
			newHits = [] #Put all the revised hits into a list
			for l in range(len(oldnames.split(','))): #Get all the hits for the groups to be merged
				hits = ddictContigGroupHits[contig][oldnames.split(',')[l]]
				for m in range(len(hits)):
					hits[m][14] = mergednewname
					newHits.append(hits[m])
			GroupHits[mergednewname] = newHits
	mod_ddictContigGroupHits[contig] = GroupHits #Create the new final dictionary

hgtmodules.printBlast_output(mod_ddictContigGroupHits)#This prints it out of order, could potentially add order to it by combining 2 modules...do I wnat this...


#Could use this to get all group information
#mod_ddict_sorted = hgtmodules.sortGroups(mod_ddictContigGroupHits)
#for contig in mod_ddictContigGroupHits:
#	for info in mod_ddict_sorted[contig]:
#		groupname = info[0]
#		for 

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
