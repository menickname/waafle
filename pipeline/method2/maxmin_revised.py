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
Modified by Tiffany Hsu (tyhsu389@gmail.com)
"""

# import
import argparse
import numpy as np
import json
import re
#import la
import matplotlib
matplotlib.use( "Agg" )
import matplotlib.pyplot as plt

# define args
parser = argparse.ArgumentParser()
parser.add_argument('--dict1', help='Dictionary of {contig: {org: {group: info}}}')
parser.add_argument('--dict2', help='Dictionary of {contig: {group: {org: info}}}')
parser.add_argument('--LGTcut', type=float, help='Cutoff above which is considered LGT, below is considered 1+orghigh')
args = parser.parse_args() 

# constants
#iBugs = 3
#iGroups = 3
#iTrials = 25000
#iTopPlots = 5

# define functions
# define the lgt scoring function
def funcLGTScore ( aafTable, contig ):
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
        fPenalty = max( [min( afScores ) for afScores in aafTable] )
	fScore = fMaxScore * ( 1 - fPenalty**2 )
	
	# determine if this is a 1+orghighonly contig
	# this should catch 1 group only, and 1 org above all
	if sum(np.sign(infoList_sorted[0][3])) == len(infoList_sorted[0][3]):
		return [contig, bug1, '-', sequencestring, fMaxScore, fPenalty, fScore, '1+orghigh']
	else: # the remainder are all potential LGT
		if fScore > args.LGTcut:
        		return [contig, bug1, bug2, sequencestring, fMaxScore, fPenalty, fScore, 'potentialLGT']
		else:
			return [contig, bug1, bug2, sequencestring, fMaxScore, fPenalty, fScore, '1+orghigh']

# read in dictionaries for scored organisms
with open(args.dict1) as infile1:
        dddictCOGS = json.load(infile1)
with open(args.dict2) as infile2:
        dddictCGOS = json.load(infile2)

# generate bug x gene tables
contigtable = {}
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
	contigtable[contig] = np.array(orgscorelist)

# score each contig in terms of "likelihood of LGT"
resultslist = []
for contig in contigtable.keys():
	contigarray = contigtable[contig]
	if len(contigarray) == 1:
		#This indicates only 1 organism annotated the contig
		bug = dddictCOGS[contig].keys()[0]
		sequencestring = ''
		for group in dddictCGOS[contig].keys():
			sequencestring += 'A'
		fAverage = np.mean(contigarray) 
		fPenalty = min(min(contigarray))
		fScore = fAverage * (1 - fPenalty**2)
		results = [contig, bug, '-', sequencestring, fAverage, fPenalty, fScore, '1orgonly'] 
	else:
		results = funcLGTScore(contigarray, contig)
	resultslist.append(results)

# print results in a table
for line in resultslist:
	print '\t'.join(str(line[i]) for i in range(len(line)))
		

"""
# score and save each table; sort so most lgt-like at the end
aResults = [( funcLGTScore( contig ), contig ) for contig in myarray]
aResults.sort()

# plot some of the results
def funcMakePlot( fHGTScore, aafTable, strName ):
	fig = plt.figure()
	fig.set_size_inches( 3, 3 )
	ax = plt.subplot( 111 )
	for afScores in aafTable:
		ax.plot( range( iGroups ), afScores )
	ax.set_title( "LGT Potential = %.3f" % ( fHGTScore ) )
	ax.set_xlabel( "Group" )
	ax.set_ylabel( "Score" )
	ax.set_xticks( range( iGroups ) )
	ax.set_xticklabels( range( 1, iGroups + 1 ) )
	ax.set_xlim( [-0.5, iGroups - 0.5] )
	ax.set_ylim( [-0.1, 1.1] )
	plt.tight_layout()
	plt.savefig( strName+".pdf" )
	
# make plots for the best scoring tables
for iPlot in range( iTopPlots ):
	fHGTScore, aafTable = aResults[-(1+iPlot)]
	funcMakePlot( fHGTScore, aafTable, "BEST%03d" % ( iPlot+1 ) )
	
# make plots for the worst scoring tables
for iPlot in range( iTopPlots ):
	fHGTScore, aafTable = aResults[iPlot]
	funcMakePlot( fHGTScore, aafTable, "WORST%03d" % ( iPlot+1 ) )
"""
