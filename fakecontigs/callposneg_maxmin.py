#!/usr/bin/python


'''
This script will classify the answerkey as what they should be.
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
def standardOrgs(orgset, taxalevel):
	searchstring = taxalevel + '__\w*'
	newset = set()
	for org in orgset:
		neworg = re.search(searchstring, org).group()
		newset.add(neworg)
	return newset

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
	aastrline = astrline.strip().split('\t')
	contig, status, contigstring = aastrline[0], aastrline[7], aastrline[3]
	
	# Call TP and FP
	if status == 'potentialLGT':
		pipelineorgs = set([aastrline[1], aastrline[2]])
	        pipelineorgset = standardOrgs(pipelineorgs, args.taxa)
		if dictanswerkey[contig][0] == 'LGT':
			newstatus = 'TP'
			TPcount += 1
			answerorgs = set(dictanswerkey[contig][1:])
        		answerorgset = standardOrgs(answerorgs, args.taxa)
		else:
			newstatus = 'FP'
			FPcount += 1
			answerorgs = set([dictanswerkey[contig][2]]) #only the recipient should be counted
			answerorgset = standardOrgs(answerorgs, args.taxa)
		orgdiff = pipelineorgset - answerorgset
		orgdiffnum = len(orgdiff)
	else:
		if status == '1orgonly':
			pipelineorgs = set([aastrline[1]])
			pipelineorgset = standardOrgs(pipelineorgs, args.taxa)
		elif status == '1+orghigh':
			if contigstring.count('A') > contigstring.count('B'):
				pipelineorgs = set([aastrline[1]])
			elif contigstring.count('A') < contigstring.count('B'):
				pipelineorgs = set([aastrline[2]])
			elif contigstring.count('A') == contigstring.count('B'):
				pipelineorgs = set([aastrline[1], aastrline[2]])
		pipelineorgset = standardOrgs(pipelineorgs, args.taxa)
		
		if dictanswerkey[contig][0] == 'noLGT':
			newstatus = 'TN'
			TNcount += 1
			answerorgs = set([dictanswerkey[contig][2]])
			answerorgset = standardOrgs(answerorgs, args.taxa)
				
		else:
			newstatus = 'FN'
			FNcount += 1
			answerorgs = set(dictanswerkey[contig][1:])
			answerorgset = standardOrgs(answerorgs, args.taxa)
		orgdiff = pipelineorgset - answerorgset
		orgdiffnum = len(orgdiff)
	
	print contig, newstatus, orgdiff, orgdiffnum, answerorgset, pipelineorgset
	#Use this to determine whether organisms are correct

TPR = TPcount/float(truepositives)
FPR = FPcount/float(truenegatives)

#print truepositives, truenegatives, TPcount, FPcount, TNcount, FNcount, TPR, FPR	


#Output results
#print "There are " + str(truepositives) + " with LGT and " + str(truenegatives) + " without LGT in this fake contig set." 
#print "There are " + str(pipelinepositives) + " called LGT and " + str(pipelinenegatives) + " called noLGT from the pipeline."
#print "There are " + str(TPcount) + " true positives, " + str(TNcount) + " true negatives, " + str(FPcount) + " false positives, and " + str(FNcount) + " false negatives."
#print "The TP rate is " + str(TPcount/float(truepositives)) + ". The FP rate is " + str(FPcount/float(truenegatives)) + "."
