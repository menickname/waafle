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
	
        contigtable[contig] = np.array(orgscorelist)
	contigorgindex[contig] = orglist

# determine one bug only and 1+orghigh contigs
def classifycontigs ( aafTable, contig ):
	bugs = aafTable.shape[0]
	genes = aafTable.shape[1]
	
	# remove contigs with only 1 gene
	if genes == 1:
		return ['1+orgonly', [0]]
	else:

		# from here, determine whether the bug has lgtpotential or not.
		LGTcandidates = []
		recipcandidates = []
		donorcandidates = []
		nonLGTcandidates = []
	
		# calculate 1 bugness score (in case)r does not get rid of everything
		maxmin = max(np.amin(aafTable, 1))
		
		
		# calculate the perc50score and realization for each bug in the contig
		for iBug in range( len( aafTable ) ):
	                bugMed = np.median(aafTable[iBug])
        	        bugMin = min(aafTable[iBug])
                	bugMax = max(aafTable[iBug])
	                bugAvg = np.average(aafTable[iBug])
			bugName = contigorgindex[contig][iBug]
	
			# get the contigcov for this 
			totalgrouplen = 0
			totalcovered = 0
			for groups in dddictCGOS[contig].keys():
				if bugName in dddictCGOS[contig][groups].keys():
					info = dddictCOGS[contig][bugName][groups]
					grouplen = float(info[7]) #contigcov*(end-start+1)
					coveredgrouplen = float(info[5])-float(info[4])+1
				else:
					firstorg = dddictCGOS[contig][groups].keys()[0]
					grouplen = float(dddictCGOS[contig][groups][firstorg][7])
					coveredgrouplen = float(0)
				totalgrouplen += grouplen
				totalcovered += coveredgrouplen
			totalgroupcov = totalcovered/totalgrouplen
			
			#Calculate the top 50%
                	count = 0
	                percscore = 0
        	        for score in aafTable[iBug]:
				if score >= bugMed:
                			count += 1
                			percscore += score
	                perc50score = percscore/float(count)

			#Calculate the realization
        		realization = bugMax - bugMin
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
                #else:
                        #return [contig, bug1, bug2, sequencestring, fMaxScore, fPenalty, fScore, '1+orghigh']


for contig in contigtable:
	contiginfo = classifycontigs(contigtable[contig], contig)
	status = contiginfo[0]
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
# Use synthetic tables
"""		
for aafTable in aContigTables:
	print aafTable
	funcRecipScore( aafTable )
"""

