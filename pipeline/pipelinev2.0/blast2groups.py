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
                #print hit
                start, end = float(hit[6]), float(hit[7])
                counter = 0
                if i == 0: #For the first hit, add it to a new group.
                        groupname = 'Group1'
                        dictGroupHits.setdefault(groupname, [])
			hit.append(groupname)
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
						if len(hit) == 15:
							hit.pop()
						hit.append(group)
						dictGroupHits[group].append(hit)
                                                #print group, 'append b/c overlap'
                                        else: #If they do not overlap by enough, count and continue to the next group.
                                                counter += 1
                if counter == len(dictGroupHits): #If we've bypassed all the groups, and none overlap or overlap by enough, add it as a new group.
                        groupnum = len(dictGroupHits) + 1
                        groupname = 'Group' + str(groupnum)
                        dictGroupHits.setdefault(groupname, [])
			hit.append(groupname)
                        dictGroupHits[groupname].append(hit)
	ddictContigGroupHits[contig] = dictGroupHits

#Merge overlapping groups
dictContigGroupSorted = hgtmodules.sortGroups(ddictContigGroupHits) #First, sort the groups by start site.
mod_ddictContigGroupHits = {}
for contig in dictContigGroupSorted.keys(): #Continue by deciding in whether each group should be merged to the previous.
        grouplist = []
        for i in range(len(dictContigGroupSorted[contig])):
                currgroup = dictContigGroupSorted[contig][i]
                curroldname, currnewname = currgroup[0], 'Group'+str(i+1)
                currstart, currend = float(currgroup[1]), float(currgroup[2])
                if i == 0: #For the first group, assign the current names
                        info = curroldname, currnewname
                        grouplist.append(info)
                else: #For all other groups, see if they overlap with the previous group
                        prevgroup = dictContigGroupSorted[contig][i-1]
                        prevoldname, prevnewname = prevgroup[0], 'Group'+str(i)
                        overlap = hgtmodules.calcOverlap(prevstart, prevend, currstart, currend)
                        if overlap > args.groupoverlap: #If there is enough overlap, merge the groups.
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
                        else: #If the groups do not overlap, assign the current names
                                info = curroldname, currnewname
                                grouplist.append(info)
                prevstart, prevend = currstart, currend #Assign the current groups as previous groups, so that we can continue comparisons

        #With the newly ordered groups, create a new dictionary that has the correct group names
        dictGroupHits2 = {}
        for oldnames, newnames in grouplist:
                if len(oldnames.split(',')) == 1: #For names that do not have to be merged, replace with correctly numbered group name
                        hits = ddictContigGroupHits[contig][oldnames]
			#print contig, hits
                        for j in range(len(hits)):
                                hits[j][14] = newnames
                        dictGroupHits2[newnames] = hits
                else: #For names that need to be merged, create new name and then replace with new merged name
                        mergednewname = 'merged' #Create the merged name first
                        for k in range(len(newnames.split(','))):
                                mergednewname += newnames.split(',')[k]
                        newHits = [] #Put all the revised hits into a list
                        for l in range(len(oldnames.split(','))): #Get all the hits for the groups to be merged
                                hits = ddictContigGroupHits[contig][oldnames.split(',')[l]]
                                for m in range(len(hits)):
                                        hits[m][14] = mergednewname
                                        newHits.append(hits[m])
			#print contig, oldnames, newnames
                        dictGroupHits2[mergednewname] = newHits
        mod_ddictContigGroupHits[contig] = dictGroupHits2 #Create the new final dictionary

#for contig in mod_ddictContigGroupHits.iterkeys():
#	for group in mod_ddictContigGroupHits[contig].iterkeys():
#		for hit in mod_ddictContigGroupHits[contig][group]:
			#print contig, group, hit

#Remove groups that are below the length threshold
mod_ddict_sorted = hgtmodules.sortGroups(mod_ddictContigGroupHits)
contiggroup_filtered = []
for contig in mod_ddict_sorted.iterkeys():
	for n in range(len(mod_ddict_sorted[contig])):
		info = mod_ddict_sorted[contig][n]	
		groupname, start, end = info[0], info[1], info[2]
		contiglen = float(end) - float(start) + 1
		if contiglen > args.length:
			contiggroup_filtered.append([contig, groupname, start, end])
			#Renumber and print out the new group names
                        #newgroupname = 'Group' + str(n+1)
                        #for o in range(len(mod_ddictContigGroupHits[contig][groupname])):
                        #        hit = mod_ddictContigGroupHits[contig][groupname][o]
                        #        hit[14] = newgroupname
                        #        print '\t'.join(str(hit[p]) for p in range(len(hit)))




#Renumber the groups that remain and create new, final dictionary
mod_ddictContigGroupHits_len = {}
counter = 0 
GroupHits = {}

for o in range(len(contiggroup_filtered)):
	contigname, groupname, start, end = contiggroup_filtered[o][0], contiggroup_filtered[o][1], contiggroup_filtered[o][2], contiggroup_filtered[o][3]
	hits = mod_ddictContigGroupHits[contigname][groupname]
	if o == 0:
		newname = 'Group' + str(counter)
		prevcontigname = contigname
	#print contigname, prevcontigname, groupname, hits
	if prevcontigname == contigname:
		counter += 1
		newname = 'Group' + str(counter)
	        if groupname == newname:
			#print contigname, 'adding because names same'
	                GroupHits[newname] = hits
       		else:
			#print contigname, 'adding because names diff'
	                for p in range(len(hits)):
				hits[p][14] = newname
			GroupHits[newname] = hits
	else: #New contig has appeared
		counter = 1 
		newname = 'Group' + str(counter)
		#print contigname, contiggroup_filtered[o-1][0], 'newcontighasappeared'
		mod_ddictContigGroupHits_len[contiggroup_filtered[o-1][0]] = GroupHits
		GroupHits = {}
		if groupname == newname:
			GroupHits[newname] = hits
		else:
			for p in range(len(hits)):
				hits[p][14] = newname
			GroupHits[newname] = hits
	prevcontigname = contigname
	#print prevcontigname,  contiggroup_filtered[o-1][0], newname
mod_ddictContigGroupHits_len[contigname] = GroupHits


#Print the final dictionary in sorted order by start site.
finaldict = hgtmodules.sortGroups(mod_ddictContigGroupHits_len)
#for contig in mod_ddictContigGroupHits_len.iterkeys():
	#for group in mod_ddictContigGroupHits_len[contig].iterkeys():
		#print contig, group, mod_ddictContigGroupHits_len[contig][group] 

for contig in finaldict.iterkeys():
	#print contig, finaldict[contig]
	for q in range(len(finaldict[contig])):
		groupname = finaldict[contig][q][0]
		hits = mod_ddictContigGroupHits_len[contig][groupname]
		for r in range(len(hits)):
			print '\t'.join(str(hits[r][s]) for s in range(len(hits[r])-1)) + '\t' + groupname

