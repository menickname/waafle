#!/usr/bin/python


'''
01/07/2015:
This script will classify the answerkey as what they should be and match the results from method3-1.
'''

#Import
import argparse
import re

#Define arguments
parser = argparse.ArgumentParser()
parser.add_argument('--answerkey', help = 'Location and file of answerkey')
parser.add_argument('--taxa', help='Level to detect LGT at. All above this level is LGT, all below is NOT LGT.')
parser.add_argument('--hgtresults', help='Location and file for all results')
args = parser.parse_args()

#Define functions
def standardOrgs( orgset, taxalevel ):
	searchstring = taxalevel + '__\w*'
	newset = set()
	for org in orgset:
		neworg = re.search(searchstring, org).group()
		newset.add(neworg)
	return newset

def checkOrgs ( resultset, answerset, newstatus ):
	comparisonset = resultset & answerset
	if newstatus == 'TP':
		if len(comparisonset) == 2:
			pass #perfect match
			stringentstatus = 'TP'
			orgmatch = 'full'
		elif len(comparisonset) == 1:
			pass #1/2 match
			stringentstatus = 'FP'
			orgmatch = 'partial'
		else:
			pass #no match
			stringentstatus = 'FP'
			orgmatch = 'none'
	elif newstatus == 'FP':
		stringentstatus = 'FP'
		if len(comparisonset) == 1:
			pass 
			orgmatch = 'partial'
		else:
			orgmatch = 'none'
			pass #print resultset, answerset, 'wrong'
	elif newstatus == 'TN':
		if len(comparisonset) == 1:
			pass 
			stringentstatus = 'TN'
			if len(resultset) == len(answerset):
				# perfect match
				orgmatch = 'full'
			else:
				# match if recipient in resultset is the same as the 1 org in answerset
				orgmatch = 'partial'
		else: 
			pass
			stringentstatus = 'FN'
			orgmatch = 'none'
	elif newstatus == 'FN':
		stringentstatus = 'FN'
		if len(comparisonset) == 2:
			orgmatch = 'full'
		elif len(comparisonset) == 1:
			orgmatch = 'partial'
		elif len(comparisonset) == 0:
			orgmatch = 'none'
	return stringentstatus, orgmatch


#Set taxalevels
leveldiff = ['t', 's', 'g', 'f', 'o', 'c', 'p', 'k']
levelnum = leveldiff.index(args.taxa)

#Read in answerkey
#The current goal is to see if the organisms match, and whether they are properly call LGT or not.
dictanswerkey = {}
for astrline in open(args.answerkey):
	aastrline = astrline.strip().split('\t')
	contigname, donorrecipinfo, donortaxa, reciptaxa = aastrline[0], aastrline[1], aastrline[2], aastrline[3]
	taxadiff, contiglen, numgenes  = aastrline[4], aastrline[5], aastrline[6]
	if leveldiff.index(taxadiff) >= leveldiff.index(args.taxa):
		dictanswerkey[contigname] = ['LGT', donortaxa, reciptaxa]
	elif leveldiff.index(taxadiff) < leveldiff.index(args.taxa):
		dictanswerkey[contigname] = ['noLGT', donortaxa, reciptaxa]

#This tells me how many P and N from the answerkey
truepositives = 0
truenegatives = 0
for contig in dictanswerkey.keys():
	if dictanswerkey[contig][0] == 'LGT':
		truepositives += 1
	else:
		truenegatives += 1

# Of the total number of TP, FP, TN, FN as called by the pipeline
pipelineP = []
TPcount = 0
FPcount = 0
TNcount = 0
FNcount = 0
for astrline in open(args.hgtresults):
	aastrline = astrline.strip().split(' ')
	contig, label, score, numbugs, numgenes = aastrline[0], aastrline[1], aastrline[2], aastrline[3], aastrline[4]

	# Call TP and FP
	if label == 'LGT':
		
		if dictanswerkey[contig][0] == 'LGT':
			newstatus = 'TP'
			TPcount += 1
		else:
			newstatus = 'FP'
			FPcount += 1
	else:
		if dictanswerkey[contig][0] == 'LGT':
			newstatus = 'FN'
			FNcount += 1
		else:
			newstatus = 'TN'
			TNcount += 1

TPR = TPcount/float(truepositives)
FPR = FPcount/float(truenegatives)

print truepositives, truenegatives, TPcount, FPcount, TNcount, FNcount, TPR, FPR	

