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

# group overlapping hits
ddictContigGroupHits = {}
for contig in dictContigHits.iterkeys():
        count = 0
        dictGroupHits = {}
        for i in range(len(dictContigHits[contig])): #For each hit in this contig, record the start & end sites
                hit = dictContigHits[contig][i]
                
		# calculate scoverage
		# record query start and end sites
		slen, sstart, send = float(hit[3]), float(hit[7]), float(hit[8])
		qlen, qstart, qend = float(hit[2]), float(hit[5]), float(hit[6])
		if send < sstart:
			send = abs(send - slen) + 1
			sstart = abs(sstart - slen) + 1
		scov = (abs(send - sstart) + 1)/slen
		ltrim = max(sstart - qstart, 0)
		rtrim = max(slen-(sstart + qlen - qstart), 0)
		scov2 = (send - sstart + 1)/(slen - ltrim - rtrim)
		hit.append(scov)
		hit.append(scov2)
		
                counter = 0
		
		# begin grouping hits
		# for the first hit, add to new 'Group1'
                if i == 0: 
                        groupname = 'Group1'
                        dictGroupHits.setdefault(groupname, [])
			dictGroupHits[groupname].append(hit)

		# if not the first hit, check for overlap with current groups
                else:
                        for group in dictGroupHits.iterkeys():
                                gstart = float(dictGroupHits[group][0][5]) #Record the group start/end sites.
                                gend  = float(dictGroupHits[group][0][6])
                                overlap = hgtmodules.calcOverlap(gstart, gend, qstart, qend)
                                #print group, gstart, gend, overlap
                                if overlap == None: #If the hit does not overlap with the group, count number of times this happens
                                        counter += 1
                                else: #If the hit overlaps with the group, calculate the %overlap.
                                        if overlap > args.hitoverlap: #Add to group if it overlaps with one, then continue to compare groups.
						dictGroupHits[group].append(hit)
                                        else: #If they do not overlap by enough, count and continue to the next group.
                                                counter += 1
                if counter == len(dictGroupHits): #If the hit matches none of the groups, add it as a new group.
                        groupnum = len(dictGroupHits) + 1
                        groupname = 'Group' + str(groupnum)
                        dictGroupHits.setdefault(groupname, [])
                        dictGroupHits[groupname].append(hit)
	ddictContigGroupHits[contig] = dictGroupHits


# merge overlapping groups
# first, sort groups by start site rather than length
dictContigGroupSorted = hgtmodules.sortGroups(ddictContigGroupHits) 

# second, create new dictionary with merged groups (2 single groups that overlap by at least _%)
mod_ddictContigGroupHits = {}
for contig in dictContigGroupSorted.keys():
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
                                if len(grouplist[len(grouplist)-1].split(',')) == 1: #If previous group is not merged, create a merged group.
                                        mergedname = prevgroup + ',' + currgroup
                                elif len(grouplist[len(grouplist)-1].split(',')) > 1: #If the previous group is merged, check to see if new group should be merged with it.
                                        if prevgroup in grouplist[len(grouplist)-1].split(','): #New "merged" group is actually part of previous merged group.
                                                mergedname = grouplist[len(grouplist)-1] + ',' + currgroup
                                        else: #New mreged group is separate from previous merged group.
                                                mergedname = prevgroup + ',' + currgroup
                                grouplist.pop() #Since we are adding a merged group, the previous entry should be deleted.
                                grouplist.append(mergedname)
                        else: #If the groups do not overlap, assign the current names
                                grouplist.append(currgroup)
                prevstart, prevend = currstart, currend #Assign the current groups as previous groups, so that we can continue comparisons

	
        # create a new dictionary with contigs, groups, and hits with the newly merged hits. At this point, we will only rename the mergedgroups. 
        dictGroupHits2 = {}
        for names in grouplist:
		if len(names.split(',')) == 1:
			hits = ddictContigGroupHits[contig][names]
			dictGroupHits2[names] = hits
                else: # for names that need to be merged, create new name and then replace with new merged name
                        mergednewname = 'merged' #Create the merged name first
                        for k in range(len(names.split(','))):
                                mergednewname += names.split(',')[k]
                        newHits = [] # put all the revised hits into a list
                        for l in range(len(names.split(','))): # get all the hits for the groups to be merged
                                hits = ddictContigGroupHits[contig][names.split(',')[l]]
				for m in range(len(hits)):
					newHits.append(hits[m])
                        dictGroupHits2[mergednewname] = newHits
        mod_ddictContigGroupHits[contig] = dictGroupHits2 # create the new final dictionary


# create file for contigs, groups, and their coordinates.
contiggroupcoords = open('contiggroupcoord.txt', 'w')


# remove groups that are below the length threshold
mod_ddict_sorted = hgtmodules.sortGroups(mod_ddictContigGroupHits)
for contig in mod_ddict_sorted.iterkeys():
	cnt = 0
	for n in range(len(mod_ddict_sorted[contig])):
		info = mod_ddict_sorted[contig][n]	
		groupname, start, end = info[0], info[1], info[2]
		grouplen = float(end) - float(start) + 1
		if grouplen > args.length:
			cnt += 1	
			#Renumber and print out the new group names
			newgroupname = 'Group' + str(cnt)
			contiggroupcoords.write(contig + '\t' + str(grouplen) + '\t' + newgroupname + '\t' + str(start) + '\t' + str(end) + '\n')
			#print contig, newgroupname, start, end #Use this to see how many groups there are, and their lengths
			for o in range(len(mod_ddictContigGroupHits[contig][groupname])):
				hit = mod_ddictContigGroupHits[contig][groupname][o]
				hit.append(newgroupname)
				print '\t'.join(str(hit[p]) for p in range(len(hit)))
contiggroupcoords.close()

#Append to info sheet
infosheet = open('info.txt', 'w')
infosheet.write('hitoverlap' + '\t' + str(args.hitoverlap) + '\n')
infosheet.write('groupoverlap' + '\t' + str(args.groupoverlap) + '\n')
infosheet.write('length cutoff' + '\t' + str(args.length) + '\n')
infosheet.close()

