#!/usr/bin/python

'''
This script will take the BLAST output and group them into genes. It will then remove all groups less than some specified length.
'''

#Import
import argparse
import hgtmodules

#Define arguments
parser = argparse.ArgumentParser()
parser.add_argument('--blastoutput', help = 'Location and file of sorted BLAST output')
parser.add_argument('--hitoverlap', type = float, default = 0.5, help = 'Percentage that BLAST hits should overlap to be combined into a group. Enter 0.5 for 50%%.')
parser.add_argument('--groupoverlap', type = float, default = 0.5, help = 'Percentage that two groups should overlap to be merged into a larger group. Enter 0.5 for 50%%.')
parser.add_argument('--length', type = float, default = 100, help = 'Groups below this length should be removed.')
args = parser.parse_args()

#Create a dictionary which has the contigs as keys, and all the BLAST output as a value (key is in list format)
dictContigHits = hgtmodules.contigHits(open(args.blastoutput))

#Group overlapping hits
ddictContigGroupHits = {}
for contig in dictContigHits.iterkeys():
        count = 0
        dictGroupHits = {}
        for i in range(len(dictContigHits[contig])): #For each hit in this contig, record the start & end sites
                hit = dictContigHits[contig][i]
                #print contig, hit
                start, end = float(hit[6]), float(hit[7])
                counter = 0
                if i == 0: #For the first hit, add it to a new group.
                        groupname = 'Group1'
                        dictGroupHits.setdefault(groupname, [])
			#print contig, groupname, hit, 'first'
			dictGroupHits[groupname].append(hit)
                else: #For all other hit, see if it goes in any of the groups.
                        for group in dictGroupHits.iterkeys():
                                gstart = float(dictGroupHits[group][0][6]) #Record the group start/end sites.
                                gend  = float(dictGroupHits[group][0][7])
                                overlap = hgtmodules.calcOverlap(gstart, gend, start, end)
                                #print group, gstart, gend, overlap
                                if overlap == None: #If the two do not overlap at all, count number of times this happens
                                        counter += 1
                                else: #If there is any overlap, calculate the %overlap.
                                        if overlap > args.hitoverlap: #Add to group if it overlaps with one, then continue to compare groups.
						dictGroupHits[group].append(hit)
                                        else: #If they do not overlap by enough, count and continue to the next group.
                                                counter += 1
                if counter == len(dictGroupHits): #If we've bypassed all the groups, and none overlap or overlap by enough, add it as a new group.
                        groupnum = len(dictGroupHits) + 1
                        groupname = 'Group' + str(groupnum)
                        dictGroupHits.setdefault(groupname, [])
                        dictGroupHits[groupname].append(hit)
	ddictContigGroupHits[contig] = dictGroupHits

#Merge overlapping groups
dictContigGroupSorted = hgtmodules.sortGroups(ddictContigGroupHits) #First, sort the groups by start site.
mod_ddictContigGroupHits = {}
for contig in dictContigGroupSorted.keys(): #Continue by deciding in whether each group should be merged to the previous.
        grouplist = []
        for i in range(len(dictContigGroupSorted[contig])):
                currgroup = dictContigGroupSorted[contig][i][0]
                currstart, currend = float(dictContigGroupSorted[contig][i][1]), float(dictContigGroupSorted[contig][i][2])
                if i == 0: #For the first group, assign the current names
                        grouplist.append(currgroup)
                else: #For all other groups, see if they overlap with the previous group
                        prevgroup = dictContigGroupSorted[contig][i-1][0]
                        overlap = hgtmodules.calcOverlap(prevstart, prevend, currstart, currend)
                        if overlap > args.groupoverlap: #If there is enough overlap, merge the groups.
                                if len(grouplist[len(grouplist)-1].split(',')) == 1: #If previous group is not merged, create a new merged group.
                                        mergedname = prevgroup + ',' + currgroup
                                elif len(grouplist[len(grouplist)-1].split(',')) > 1: #If the previous group is merged, check for multiple merges.
                                        if prevgroup in grouplist[len(grouplist)-1].split(','): #Merge 3+ groups if the previous merged groups are consecutive to the current one.
                                                mergedname = grouplist[len(grouplist)-1] + ',' + currgroup
                                        else: #Add as a new merged group if the previous merged groups are NOT consecutive.
                                                mergedname = prevgroup + ',' + currgroup
                                grouplist.pop() #Since we are adding a merged group, the previous entry should be deleted.
                                grouplist.append(mergedname)
                        else: #If the groups do not overlap, assign the current names
                                grouplist.append(currgroup)
                prevstart, prevend = currstart, currend #Assign the current groups as previous groups, so that we can continue comparisons

        #With the newly ordered groups, create a new dictionary. These still use the 'oldnames' that may be out of order. 
        dictGroupHits2 = {}
        for names in grouplist:
		if len(names.split(',')) == 1:
			hits = ddictContigGroupHits[contig][names]
			dictGroupHits2[names] = hits
                else: #For names that need to be merged, create new name and then replace with new merged name
                        mergednewname = 'merged' #Create the merged name first
                        for k in range(len(names.split(','))):
                                mergednewname += names.split(',')[k]
                        newHits = [] #Put all the revised hits into a list
                        for l in range(len(names.split(','))): #Get all the hits for the groups to be merged
                                hits = ddictContigGroupHits[contig][names.split(',')[l]]
				for m in range(len(hits)):
					newHits.append(hits[m])
                        dictGroupHits2[mergednewname] = newHits
        mod_ddictContigGroupHits[contig] = dictGroupHits2 #Create the new final dictionary

contiggroupcoords = open('contiggroupcoord.txt', 'w')
#Remove groups that are below the length threshold
mod_ddict_sorted = hgtmodules.sortGroups(mod_ddictContigGroupHits)
for contig in mod_ddict_sorted.iterkeys():
	cnt = 0
	for n in range(len(mod_ddict_sorted[contig])):
		info = mod_ddict_sorted[contig][n]	
		groupname, start, end = info[0], info[1], info[2]
		contiglen = float(end) - float(start) + 1
		if contiglen > args.length:
			cnt += 1	
			#Renumber and print out the new group names
			newgroupname = 'Group' + str(cnt)
			contiggroupcoords.write(contig + '\t' + newgroupname + '\t' + str(start) + '\t' + str(end) + '\n')
			#print contig, newgroupname, start, end #Use this to see how many groups there are, and their lengths
			for o in range(len(mod_ddictContigGroupHits[contig][groupname])):
				hit = mod_ddictContigGroupHits[contig][groupname][o]
				hit.append(newgroupname)
				print '\t'.join(str(hit[p]) for p in range(len(hit)))

contiggroupcoords.close()
