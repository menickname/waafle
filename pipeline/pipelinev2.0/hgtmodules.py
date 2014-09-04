#!/usr/bin/python

'''
This contains all my functions for the hgt scripts.
This will simplify the actual code.
'''

#This function calculates the overlap between 2 BLAST hit alignments using their respective start and end coordinates. It returns 'None' if there is no overlap, and the percent overlap if there is an overlap.
def calcOverlap(iStart1, iEnd1, iStart2, iEnd2):
        if iStart1 - iEnd2 > 0 or iStart2 - iEnd1 > 0:
                return None
        else:
                status = True
                listcoord = [float(iStart1), float(iEnd1), float(iStart2), float(iEnd2)]
                listcoord.sort()
                contig1len = iEnd1 - iStart1 + 1
                contig2len = iEnd2 - iStart2 + 1
                divisor = min(contig1len, contig2len)
                overlap = (listcoord[2]-listcoord[1] + 1)/divisor
                return overlap

#This function generates a dictionary with {contig: [alignments]}.
def contigHits(blastoutput):
	dictContigHits = {}
	for astrline in blastoutput:
        	aastrline = astrline.strip().split('\t')
        	dictContigHits.setdefault(aastrline[0],[]).append(aastrline)
	return dictContigHits

#This function generates a dictionary of dictionaries, with {contig: {group: [alignments]}}.
def contigGroupHits(dictContigHits):
	ddictContigGroupHits = {}
	for bstrline in dictContigHits.iterkeys():
        	dictGroupHits = {}
        	for i in range(len(dictContigHits[bstrline])):
                	group = dictContigHits[bstrline][i][14]
	                dictGroupHits.setdefault(group,[]).append(dictContigHits[bstrline][i])
        	ddictContigGroupHits.setdefault(bstrline, {}).update(dictGroupHits)
	return ddictContigGroupHits

#This function takes a dictioanry of dictionaries and sorts the groups based on start site.
def sortGroups(ddictContigGroupHits2, grouped=True):
	dictContigGroupOrder = {}
	for contig in ddictContigGroupHits2.iterkeys():
		groupcoordlist = []
		for group in ddictContigGroupHits2[contig].iterkeys():
			startlist = []
			endlist = []
			for i in range(len(ddictContigGroupHits2[contig][group])):
				startlist.append(float(ddictContigGroupHits2[contig][group][i][6]))
				endlist.append(float(ddictContigGroupHits2[contig][group][i][7]))
			minStart = min(startlist)
			maxStart = max(endlist)
			info = group, minStart, maxStart
			groupcoordlist.append(info)
		groupcoordlist.sort(key=lambda x: x[1])
		dictContigGroupOrder[contig] = groupcoordlist
	return dictContigGroupOrder
	'''	
		dictGroupHits = {}
		for j in range(len(groupcoordlist)):
			oldgroupname = groupcoordlist[j][1]
			newgroupname = 'Group' + str(j+1)
			#if grouped == True: 
			#	if oldgroupname != newgroupname:
			#		ddictContigGroupHits2[contig][newgroupname] = ddictContigGroupHits2[contig].pop(oldgroupname)
			#		for k in range(len(ddictContigGroupHits2[contig][newgroupname])):
			#			ddictContigGroupHits2[contig][newgroupname][k][14] = newgroupname
			for k in range(len(ddictContigGroupHits2[contig][oldgroupname])):
				if grouped == True:
					if oldgroupname != newgroupname:
						ddictContigGroupHits2[contig][oldgroupname][k][14] = newgroupname
					else:
						pass
				else:
					#ddictContigGroupHits2[contig][oldgroupname][k].append(newgroupname)
					linelist = ddictContigGroupHits2[contig][oldgroupname][k][:]
					linelist.append(newgroupname)
					print '\t'.join(str(linelist[l]) for l in range(len(linelist)))
					#print newline
				#dictGroupHits.setdefault(newgroupname, []).append(newline)
		#newddictContigGroupHits[contig] = dictGroupHits
	#return newddictContigGroupHits
	'''

def printBlast_output(ddictContigGroupHits3):
	for contig in ddictContigGroupHits3.iterkeys():
		for group in ddictContigGroupHits3[contig].iterkeys():
			for l in range(len(ddictContigGroupHits3[contig][group])):
				print '\t'.join(ddictContigGroupHits3[contig][group][l]) + '\t' + group	
