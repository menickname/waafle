#!/usr/bin/env python

"""
aafTable refers to a 2D numpy array that represents a contig 
with N bugs (rows) and M groups (columns); each row of the table 
contains the group scores for one bug, as in

		group1	group2	group3
bug1  [[0.7,	0.8,	0.0],
bug2   [0.1, 	0.9,	0.8],
bug3   [0.1,	0.2,	0.3]]

This script generates random tables and then applies an lgt
scoring function to find tables that look like they contain lgt.
Best/worst cases are plotted for visual inspection.
-----------------------------------------------
Author: Eric Franzosa (eric.franzosa@gmail.com)
"""

# import
import argparse
import numpy as np
import matplotlib
matplotlib.use( "Agg" )
import matplotlib.pyplot as plt
import json
import re
from operator import itemgetter, attrgetter, methodcaller

# arguments
parser = argparse.ArgumentParser()
parser.add_argument( '--dddictCOGS', help='Location and file of json file with {contig:{org:group}}.' )
parser.add_argument( '--dddictCGOS', help='Location and file of json file with {contig:{group:org}}.' )
parser.add_argument( '--realization', help='Value above which the organism likely participates in LGT.' )
parser.add_argument( '--score', type=float, help='Value above which we call LGT, value below we call nonLGT' )
args = parser.parse_args()

# generate synthetic data
"""
# constants
iBugs = 10
iGroups = 3
iTrials = 1
iTopPlots = 5

# generate a bunch of random bug x gene tables
aContigTables = [np.random.rand( iBugs, iGroups ) for iTrial in range( iTrials )]
"""

# read in real data
# read in dictionaries for scored organisms
with open(args.dddictCOGS) as infile1:
        dddictCOGS = json.load(infile1)
with open(args.dddictCGOS) as infile2:
        dddictCGOS = json.load(infile2)

# generate bug x gene tables
contigtable = {}
contigorgindex = {}
for contig in dddictCOGS.keys():
        orgscorelist = []
        orglist = []
        for org in dddictCOGS[contig].keys():

                # sort the groups by number
                groupnumlist = []
                for group in dddictCGOS[contig].keys():
                        groupnum = re.search('[0-9]+', group).group()
                        groupnumlist.append([groupnum, group])
                grouplist_sorted = sorted(groupnumlist)

                # assign scores for each group
                scorelist = []
                for groupnums, groups in grouplist_sorted:
                        if groups not in dddictCOGS[contig][org].keys():
                                score = float(0)
                        else:
                                score = dddictCOGS[contig][org][groups][0]
                        scorelist.append(score)
                orgscorelist.append(scorelist)
		orglist.append(org)
	
	# get scores for unknown organism	
	table = np.array(orgscorelist)
	unknownorg_scorelist = []
	for i in range(table.shape[1]):
		unknownscore = 1 - np.amax(table[:,i])
		unknownorg_scorelist.append(unknownscore)
	orgscorelist.append(unknownorg_scorelist)
	orglist.append('unknown')

	# generate np array for unknown and all else
        contigtable[contig] = np.array(orgscorelist)
	contigorgindex[contig] = orglist

#This function will divide the organisms into donors and recipients
def finddonorrecip ( scorelist, numgenes ):
	donorlist = []
	reciplist = []
	for score, donorscore, recipscore, perc50, realization, index in scorelist:
		if donorscore < 0.5 and recipscore > 0.5:
			reciplist.append(index)
		if donorscore >= 0.5 and recipscore <= 0.5:
			donorlist.append(index)
	return [donorlist, reciplist]

#This function should calculate the complementarity between pairs
def calcComplement ( scorelist, aafTable ):
	complementlist = []
	for i in range(len(scorelist)-1):
		for j in range(i+1, len(scorelist)):
			org1 = aafTable[i]
			org2 = aafTable[j]
			newarray = np.array([org1, org2])
			complementarity = np.average(np.amax(newarray, 0))
			donorscore_1, recipscore_1 = scorelist[i][1], scorelist[i][2]
			donorscore_2, recipscore_2 = scorelist[j][1], scorelist[j][2]
			complementlist.append([complementarity, donorscore_1, recipscore_1, donorscore_2, recipscore_2, i, j])
	return sorted(complementlist, reverse=True)

# Remove one gene only contigs, classify remainder.
def classifycontigs ( aafTable, contig ):
	bugs = aafTable.shape[0] #number of bugs
	genes = aafTable.shape[1] #number of genes
	
	# from here, determine whether the bug has lgtpotential or not.
	LGTcandidates = []
	recipcandidates = []
	donorcandidates = []
	nonLGTcandidates = []
	
	# calculate 1 bugness score (in case)r does not get rid of everything
	maxmin = max(np.amin(aafTable, 1))

	# calculate the perc50score and realization for each bug in the contig
	scorelist = []
	for iBug in range( len( aafTable ) ):
                bugMed = np.median(aafTable[iBug])
                bugAvg = np.average(aafTable[iBug])
		bugName = contigorgindex[contig][iBug]

                #Calculate the top 50%
                count = 0
                percscore = 0
                for score in aafTable[iBug]:
                        if score >= bugMed:
                                count += 1
                                percscore += score
		perc50score = percscore/float(count)
		
                #Calculate the realization
                realization = max(aafTable[iBug]) - min(aafTable[iBug])
		
		#Calculate the score
		score = perc50score - realization
		donorscore = realization - perc50score
		recipscore = min(realization, perc50score)
		scorelist.append([score, donorscore, recipscore, perc50score, realization, iBug])

	scoresonly = np.array(scorelist)[:,0]
	
	#Determine if there is LGT
	if genes == 1: #If there is only 1 gene to begin with, no chance of LGT
		pass
	elif np.amax(scoresonly) > 0.70:
		pass
	else:
		#print contig
		#print aafTable
		complementlist = calcComplement(scorelist, aafTable)
		
		candidatelist = []
		for i in range(len(complementlist)):
			complementarity = complementlist[i][0]
			donorscore1, donorscore2 = complementlist[i][1], complementlist[i][3]
			recipscore1, recipscore2 = complementlist[i][2], complementlist[i][4]
			orgindex1, orgindex2 = complementlist[i][5], complementlist[i][6]
			if complementarity > 0.5: 
				if genes > 3:
					if donorscore1 > 0.5 and recipscore2 > 0.5:
						#print complementlist[i]
						candidatelist.append(complementlist[i])
					elif donorscore2 > 0.5 and recipscore1 > 0.5:
						#print complementlist[i]
						candidatelist.append(complementlist[i])
				else:
					scorelist = complementlist[i][1:5]
					newscorelist = [j for j in complementlist[i][1:5] if j > 0.5]
					if len(newscorelist) >= 2:
						candidatelist.append(complementlist[i])
		if len(candidatelist) == 0:
			print contig, 'No complementation'
		else:
			print contig, 'LGT'
		
		"""
		# get the contigcov for this 
		grouplen = 0
		coveredgrouplen = 0
		for groups in dddictCGOS[contig].keys():
			if bugName in dddictCGOS[contig][groups].keys():
				info = dddictCOGS[contig][bugName][groups]
				grouplen += float(info[7]) #contigcov*(end-start+1)
				coveredgrouplen += float(info[5])-float(info[4])+1
			else:
				firstorg = dddictCGOS[contig][groups].keys()[0]
				grouplen = float(dddictCGOS[contig][groups][firstorg][7])
				coveredgrouplen += float(0)
			totalgrouplen += grouplen
		totalcovered += coveredgrouplen
		totalgroupcov = totalcovered/totalgrouplen
		
	
	
			if realization > float(args.realization):
				LGTcandidates.append([perc50score, realization, totalgroupcov, maxmin, bugAvg, iBug])
				if totalgroupcov >= 0.5:
					recipcandidates.append([perc50score, realization, totalgroupcov, maxmin, bugAvg, iBug])
				if totalgroupcov <= 0.5:
					donorcandidates.append([perc50score, realization, totalgroupcov, maxmin, bugAvg, iBug])
			else:
				nonLGTcandidates.append([perc50score, realization, totalgroupcov, maxmin, bugAvg, iBug])
		
		recipcandidates_sorted = sorted(recipcandidates, key=lambda numbers: numbers[0], reverse=True)
		donorcandidates_sorted = sorted(donorcandidates, key=lambda numbers: numbers[0], reverse=True)
		LGTcandidates_sorted = sorted(LGTcandidates, key=itemgetter(2,1), reverse=True) #lambda numbers: numbers[0], reverse=True)
		nonLGTcandidates_sorted = sorted(nonLGTcandidates, key=lambda numbers: numbers[4], reverse=True)
		
		if len(LGTcandidates) == 0:
			return ['1+orgonly', [0]]
		else:
			# these are potential LGT
			return ['potentialLGT', [recipcandidates_sorted, donorcandidates_sorted, LGTcandidates_sorted, nonLGTcandidates_sorted]]
			"""
"""
def funcLGTscore ( aafTable, penalty ):
	infoList = []
        for iBug1 in range( len( aafTable ) ):
                for iBug2 in range( iBug1+1, len( aafTable ) ):
                        afScores1, afScores2 = aafTable[iBug1], aafTable[iBug2]
                        afDiff = afScores1 - afScores2
                        # lgt score ("0.5" factor keeps this in [0,1])
                        fTempScore = 0.5 * ( max( afDiff ) - min( afDiff ) ) * min( abs( afDiff ) )
                        info = [iBug1, iBug2, fTempScore, afDiff]
                        infoList.append(info)
        # determine which 2 organisms are most likely to have LGT between them
        infoList_sorted = sorted(infoList, key=lambda Info: Info[2], reverse=True)
        bug1, bug2, fMaxScore = dddictCOGS[contig].keys()[infoList_sorted[0][0]], dddictCOGS[contig].keys()[infoList_sorted[0][1]], infoList_sorted[0][2]
        # determine the order of the sequence for the organisms
        sequencestring = ''
        for score in infoList_sorted[0][3]:
                if score > 0:
                        sequencestring += 'A'
                elif score < 0:
                        sequencestring += 'B'
                elif score == 0:
                        sequencestring += '-'
        # penalize if one bug covers whole contig well
        # Note: squaring penalty punishes ambiguous cases less
        # fPenalty = max( [min( afScores ) for afScores in aafTable] )
       	fPenalty = penalty 
	fScore = fMaxScore * ( 1 - fPenalty**2 )

        # determine if this is a 1+orghighonly contig
        # this should catch 1 group only, and 1 org above all
        #if sum(np.sign(infoList_sorted[0][3])) == len(infoList_sorted[0][3]):
        #        return [contig, bug1, '-', sequencestring, fMaxScore, fPenalty, fScore, '1+orghigh']
        #else: # the remainder are all potential LGT
        return [contig, bug1, bug2, sequencestring, fMaxScore, fPenalty, fScore]
"""

for contig in contigtable:
	contiginfo = classifycontigs(contigtable[contig], contig)
	"""status = contiginfo[0]
	if status == '1+orgonly':
		pass #print '1+orgonly'
	if status == 'potentialLGT':
		recipcandidates, donorcandidates, LGTcandidates, nonLGTcandidates = contiginfo[1][0], contiginfo[1][1], contiginfo[1][2], contiginfo[1][3]
		#print contig, contigtable[contig]
		LGTarray = np.zeros(shape=(len(LGTcandidates), np.shape(contigtable[contig])[1])) 
		for i in range(len(LGTcandidates)):
			bugindex = LGTcandidates[i][5]
			LGTarray[i] = contigtable[contig][bugindex]
		maxLGTarray_avg = np.average(np.amax(LGTarray, 0)) 
		
		toprecip_top50score = LGTcandidates[0][0]
		penalty = LGTcandidates[0][3]

		if len(nonLGTcandidates) > 0:
			nonLGTcandidate_topavg = nonLGTcandidates[0][4]
		else:
			nonLGTcandidate_topavg = 0
		
		# Start binning results
		if toprecip_top50score > nonLGTcandidate_topavg:
			if np.shape(LGTarray)[0] > 1: # and maxLGTarray_avg > nonLGTcandidate_topavg:
				pass #print contig, contigtable[contig], '\n'
				#print 'LGTarray', LGTarray, '\n'
				#print np.amax(LGTarray,0)
				score = funcLGTscore(LGTarray, penalty)[6]
				if score > args.score:
					pass #print 'LGT'
				else:
					print 'unknown'
					print contig, contigtable[contig], '\n'
	                                print 'LGTarray', LGTarray, '\n'
        	                        print 'nonLGT'
					for entry in nonLGTcandidates:
						print entry
                	               
					print '\n', 'maxLGT combination', np.amax(LGTarray,0)
                        	        print contig, 'potentialLGTscore', funcLGTscore(LGTarray, penalty)[6], '\n', '\n'

			else:
				#print contig, contigtable[contig], '\n'
                                #print 'LGTarray', LGTarray, '\n'
				#print 'nonLGT', nonLGTcandidates
                                #print np.amax(LGTarray,0)
                                #print contig, 'potentialLGT', funcLGTscore(LGTarray, penalty)[6]
				print 'unknown'
		else:
			pass
			#print contig, contigtable[contig]
			#print LGTarray # from here, we need to distinguish whether there is true LGT
			#print '1+orgonly'
			"""
# Use synthetic tables
"""		
for aafTable in aContigTables:
	print aafTable
	funcRecipScore( aafTable )
"""

