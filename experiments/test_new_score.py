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

import numpy as np
import matplotlib.pyplot as plt

# constants
iBugs = 3
iGroups = 3
iTrials = 25000
iTopPlots = 5

# generate a bunch of random bug x gene tables
aContigTables = [np.random.rand( iBugs, iGroups ) for iTrial in range( iTrials )]

# define the lgt scoring function
def funcLGTScore ( aafTable ):
	fMaxScore = 0
	for iBug1 in range( len( aafTable ) ):
		for iBug2 in range( iBug1+1, len( aafTable ) ):
			afScores1, afScores2 = aafTable[iBug1], aafTable[iBug2]
			afDiff = afScores1 - afScores2
			# **** new lgt score formula from email ****
			fTempScore = ( max( afDiff ) - min( afDiff ) ) * min( abs( afDiff ) )
			if fTempScore > fMaxScore:
				fMaxScore = fTempScore
	# **** penalize the best score if one bug covers the whole contig ****
	fPenalty = max( [min( afScores ) for afScores in aafTable] )
	return fMaxScore - fPenalty

# score and save each table; sort so most lgt-like at the end
aResults = [( funcLGTScore( aafTable ), aafTable ) for aafTable in aContigTables]
aResults.sort()

# plot some of the results
def funcMakePlot( fHGTScore, aafTable, strName ):
	fig = plt.figure()
	fig.set_size_inches( 3, 3 )
	ax = plt.subplot( 111 )
	for afScores in aafTable:
		ax.plot( range( iGroups ), afScores )
	ax.set_title( "LGT Potential = %.2f" % ( fHGTScore ) )
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
