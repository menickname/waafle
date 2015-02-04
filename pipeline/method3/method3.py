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
import math

# arguments
parser = argparse.ArgumentParser()
parser.add_argument( '--dddictCOGS', help='Location and file of json file with {contig:{org:group}}.' )
parser.add_argument( '--dddictCGOS', help='Location and file of json file with {contig:{group:org}}.' )
parser.add_argument( '--LGTscore', type=float, help='Value below which the contig likely has LGT.' )
parser.add_argument( '--complement', type=float, help='Value above which we call LGT, value below we call no LGT' )
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
			drscore2_1, drscore2_2 = scorelist[i][7], scorelist[j][7]
			complementlist.append([complementarity, donorscore_1, recipscore_1, donorscore_2, recipscore_2, i, j, drscore2_1, drscore2_2])
	return sorted(complementlist, reverse=True)

def callDonorRecip ( LGTpair, numgenes ):
	donorindex, recipindex = 0, 0
	if numgenes >= 3:
		pair1donorscore, pair1recipscore = LGTpair[1], LGTpair[4]
		pair2donorscore, pair2recipscore = LGTpair[2], LGTpair[3]
		if pair1donorscore > 0.5 and pair1recipscore > 0.5: 
			donorindex, recipindex = LGTpair[5], LGTpair[6]
		elif pair2donorscore > 0.5 and pair2recipscore > 0.5:
			donorindex, recipindex = LGTpair[6], LGTpair[5]
		else:
			donorindex, recipindex = 'unknown', 'unknown'
	else:
		donorindex, recipindex = 'unknown', 'unknown'
	return [donorindex, recipindex]

# Remove one gene only contigs, classify remainder.
def classifycontigs ( aafTable, contig ):
	bugs = aafTable.shape[0] #number of bugs
	genes = aafTable.shape[1] #number of genes
	
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
                drscore2 = bugAvg-bugMed #Test Method 1
                drscore3 = (perc50score-bugMed)/(perc50score+0.0000001) #Test Method 2

		scorelist.append([score, donorscore, recipscore, perc50score, realization, iBug, bugAvg, drscore2, drscore3])

		#Print all the different types of scores for comparison
		#print contig, iBug, bugs, genes, score, donorscore, recipscore, perc50score, realization, bugMed, bugAvg, bugAvg-bugMed, (perc50score-bugMed)/(perc50score+0.0000001)

	scorelist_sort = sorted(scorelist, reverse=True)
	scoresonly = np.array(scorelist)[:,0]
	
	answerlist = []
	#Determine if there is LGT
	if genes == 1: #If there is only 1 gene to begin with, no chance of LGT
		answerlist.append(['1orgcontig', contig, bugs, genes, scorelist_sort[0][6], '1geneonly', contigorgindex[contig][scorelist_sort[0][5]], scorelist_sort[0][0], scorelist_sort[0][3], scorelist_sort[0][4]])
	elif np.amax(scoresonly) > args.LGTscore:
		answerlist.append(['1orgcontig', contig, bugs, genes, scorelist_sort[0][6], '1+genes', contigorgindex[contig][scorelist_sort[0][5]], scorelist_sort[0][0], scorelist_sort[0][3], scorelist_sort[0][4]])

	else:
		# calculate complement for all pairs of bugs
		complementlist = calcComplement(scorelist, aafTable)
		
		# top of the list is the top LGT candidate
		topLGTpair = complementlist[0]
		complementscore = complementlist[0][0]
		org1index, org2index = complementlist[0][5], complementlist[0][6]
		org1info = scorelist[org1index]
                org1score, org1_perc50, org1_real = org1info[0], org1info[3], org1info[4]
                org2info = scorelist[org2index]
                org2score, org2_perc50, org2_real = org2info[0], org2info[3], org2info[4]
		
		if complementscore > args.complement:
			print aafTable
				
			# identify the donor and recipient - NEW METHOD #1
			org1, org2 = topLGTpair[5], topLGTpair[6]
			org1_drscore2, org2_drscore2 = topLGTpair[7], topLGTpair[8]
			print contig, org1, org1_drscore2
			print contig, org2, org2_drscore2

			# identify the donor and recipient - NEW METHOD #2
			#if drscore3 >
			
			# identify the donor and recipient if possible - OLD METHOD
			"""
			donorindex, recipindex = callDonorRecip(topLGTpair, genes)[0], callDonorRecip(topLGTpair, genes)[1]
			if donorindex != 'unknown':
				donorinfo = scorelist[donorindex]
				donorscore, donor_perc50, donor_real = donorinfo[0], donorinfo[3], donorinfo[4]
				recipinfo = scorelist[recipindex]
				recipscore, recip_perc50, recip_real = recipinfo[0], recipinfo[3], recipinfo[4]
				answerlist.append(['LGT', contig, bugs, genes, complementscore, 'Donor', contigorgindex[contig][donorindex], donorscore, donor_perc50, donor_real])
				answerlist.append(['LGT', contig, bugs, genes, complementscore, 'Recipient', contigorgindex[contig][recipindex], recipscore, recip_perc50, recip_real])
			else:
				pass
				answerlist.append(['LGT', contig, bugs, genes, complementscore, 'Unknown', contigorgindex[contig][org1index], org1score, org1_perc50, org1_real])
				answerlist.append(['LGT', contig, bugs, genes, complementscore, 'Unknown', contigorgindex[contig][org2index], org2score, org2_perc50, org2_real])
		else:
			pass#no LGT
			answerlist.append(['No LGT', contig, bugs, genes, complementscore, 'Unknown', contigorgindex[contig][org1index], org1score, org1_perc50, org1_real])
			answerlist.append(['No LGT', contig, bugs, genes, complementscore, 'Unknown', contigorgindex[contig][org2index], org2score, org2_perc50, org2_real])
			answerlist.append(['No LGT', contig, bugs, genes, scorelist_sort[0][6], 'Best_one_bug', contigorgindex[contig][scorelist_sort[0][5]], scorelist_sort[0][0], scorelist_sort[0][3], scorelist_sort[0][4]])

	return answerlist
	"""

"""
		candidatelist = []
		for i in range(len(complementlist)):
			complementarity = complementlist[i][0]
			donorscore1, donorscore2 = complementlist[i][1], complementlist[i][3]
			recipscore1, recipscore2 = complementlist[i][2], complementlist[i][4]
			orgindex1, orgindex2 = complementlist[i][5], complementlist[i][6]
			if complementarity > 0.5: 
				if genes > 2:
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
"""	

#Script begins

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
		
		# rename org for "R"
		orgname = re.search('g__\w*', org).group()
                orglist.append(orgname)

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


for contig in contigtable:
	contiginfo = classifycontigs(contigtable[contig], contig)
	#for astrline in contiginfo:
	#	print '\t'.join(str(i) for i in astrline)

#Use synthetic tables		
"""		
for aafTable in aContigTables:
	print aafTable
	funcRecipScore( aafTable )
"""

