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
parser.add_argument('--genetable', help='Location and file for genetable')
parser.add_argument('--hgtresults', help='Location and file for LGT results')
args = parser.parse_args()

#Set taxalevels
leveldiff = ['t', 's', 'g', 'f', 'o', 'c', 'p']
levelnum = leveldiff.index(args.taxa)

#Read in answerkey
#The current goal is to see if the organisms match, and whether they are properly call LGT or not.
dictanswerkey = {}
for astrline in open(args.answerkey):
	aastrline = astrline.strip().split('\t')
	contigname, donorrecipinfo, donortaxa, reciptaxa = aastrline[0], aastrline[1], aastrline[2], aastrline[3]
	splitdonortaxa, splitreciptaxa = donortaxa.split('|'), reciptaxa.split('|')

	for i in range(len(splitreciptaxa)):
		if splitreciptaxa[i] == splitdonortaxa[i] and i == len(splitreciptaxa)-1:
			dictanswerkey[contigname] = ['noLGT', donortaxa, reciptaxa]
			print contigname, 'noLGT'
		elif splitreciptaxa[i] != splitdonortaxa[i]:
			#print contigname, splitreciptaxa[i], splitdonortaxa[i]
			taxadiff = re.search('.+__', splitreciptaxa[i]).group()[0]
			currlevelnum = leveldiff.index(taxadiff)
			if currlevelnum >= levelnum:
				dictanswerkey[contigname] = ['LGT', donortaxa, reciptaxa]
				print contigname, 'LGT'
			else:
				dictanswerkey[contigname] = ['noLGT', donortaxa, reciptaxa]
				print contigname, 'noLGT'
			break

"""
highconfLGT = open(args.hgtresults).readlines()
highconfLGTlist = []
for i in range(len(highconfLGT)):
	contigname = highconfLGT[i].split(' ')[0]
	highconfLGTlist.append(contigname)

#Read in answers from pipeline
for bstrline in open(args.genetable):
	bbstrline = bstrline.strip().split(' ')
	contigname, groupnum, taxa = bbstrline[0], bbstrline[1], bbstrline[2]
	finalscore, finalpercID, finalgroupcov, contigcov = bbstrline[3], bbstrline[4], bbstrline[5], bbstrline[6]
	start, end, combhitlen, grouplen, contiglen, status = bbstrline[7], bbstrline[8], bbstrline[9], bbstrline[10], bbstrline[11], bbstrline[12]
	if status == '1orgonly':
		answer_status = dictanswerkey[contigname][0]
		if answer_status == 'noLGT':
			print contigname, 'TN' #Supposed to be noLGT, and algorithm gave 1orgonly
		else:
			print contigname, 'FN' #Supposed to be LGT, and algorithm gave 1orgonly
	elif status == '1+orghigh':
		answer_status = dictanswerkey[contigname][0]
		if answer_status == 'noLGT':
			print contigname, 'TN' #Supposed to be noLGT, and algorithm gave 1+orghigh
		else:
			print contigname, 'FN' #Supposed to be LGT, and algorithm gave 1+orghigh
	elif status == 'potentialLGT':
		if contigname in highconfLGTlist:
			print contigname, 'TP'
		else:
			if answer_status == 'noLGT':
				print contigname, 'FP'  #Supposed to be noLGT, but called LGT
			else:
				print contigname, 'FP-notreally'
"""
"""
#Read in answers from hgt results
for cstrline in open(args.hgtresults):
	ccstrline = cstrline.strip().split(' ')
	contigname, order, status, completecoverage = ccstrline[0], ccstrline[1], ccstrline[2], ccstrline[3]
	answer_status = dictanswerkey[contigname][0]
	if answer_status == 'LGT': 
		print contigname, 'TP', status #Supposed to be LGT, and algorithm gave LGT
	else:
		print contigname, 'FP', status
	#	for bstrline in open(args.genetable):
	#		bbstrline = bstrline.strip().split(' ')
	#		contigname, status = bbstrline[0], bbstrline[12]
	#		if status == 'potentialLGT':
	#			print contigname, 'FN', 'notdetected_highconf' #Supposed to be no LGT, and algorithm gave LGT
"""	

