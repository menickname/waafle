#!/usr/bin/python

'''
This script detects high confidence LGT events.
'''

# Import
import argparse
import json
import re
import collections

# Define arguments
parser = argparse.ArgumentParser()
parser.add_argument('--dict1', help='Dictionary of {contig: {org: {group: info}}}')
parser.add_argument('--dict2', help='Dictionary of {contig: {group: {org: info}}}')
parser.add_argument('--delta', type=float, help='Upper bound for high confidence scoring BLAST hits')
parser.add_argument('--epsilon', type=float, help='Lower bound for low confidence hits')
args = parser.parse_args()

#Read in dictionaries from previous script
with open(args.dict1) as infile1:
        dddictCOGS = json.load(infile1)
with open(args.dict2) as infile2:
        dddictCGOS = json.load(infile2)

#Categorize organisms within each contig as mid, low, highmid, highlow, highmidlow, midlow.
dictContigOrgLabel = {}
for contig in dddictCOGS.keys():
	dictOrgLabel = {}
	for taxa in dddictCOGS[contig].keys():
		dictLabel = {}
		dictQuality = {}
		for group in dddictCGOS[contig].keys():
			if group not in dddictCOGS[contig][taxa].keys():
				dictQuality.setdefault('low', []).append(group)
			else:
				info = dddictCOGS[contig][taxa][group]
				finalscore, finalpercID, finalgroupcov, contigcov = info[0], info[1], info[2], info[3]
				newstart, newend, combhitlen = info[4], info[5], info[6]
				grouplen, contiglen = info[7], info[8]
				if finalscore > args.delta:
					dictQuality.setdefault('high', []).append(group)
				elif finalscore <= args.delta and finalscore >= args.epsilon:
					dictQuality.setdefault('mid', []).append(group)
				elif finalscore < args.epsilon:
					dictQuality.setdefault('low', []).append(group)
		label=set(dictQuality.keys())
		labelname=''.join(sorted(label))
		dictLabel.setdefault(labelname, {}).update(dictQuality)
		dictOrgLabel[taxa] = dictLabel
	dictContigOrgLabel[contig] = dictOrgLabel

#Count the number of times each category comes up. Contigs with at least 2 highlow annotations are considered "more potential LGT." 
contiglist = []
for contig in dictContigOrgLabel.keys():
	labelcnt = collections.Counter()
	labellist = []
	for org in dictContigOrgLabel[contig].keys():
		for label in dictContigOrgLabel[contig][org].keys():
			labellist.append(label)
	for labels in labellist:
		labelcnt[labels] += 1
	if labelcnt['highlow'] >= 2:
		contiglist.append(contig)

#We will rank each bug in the contig using Curtis's scoring idea
#This means that we will sum the score*contigcov across the contig, and the highest ranking one covers most of the contig. This weights by how large the group is.
dictContigOrgRank = {}
for contig in contiglist:
	ranklist = []
	for org in dddictCOGS[contig].keys():
		rankbug = 0
		for group in dddictCOGS[contig][org].keys():
			info = dddictCOGS[contig][org][group]
			finalscore, finalpercID, finalgroupcov, contigcov = info[0], info[1], info[2], info[3]
                        newstart, newend, bp, grouplen, contiglen = info[4], info[5], info[6], int(info[7]), info[8]
			#finalscore, finalpercID, finalgroupcov, contigcov = info[0], info[1], info[2], info[3]
			#newstart, newend, combhitlen = info[4], info[5], info[6]
			#grouplen, contiglen = info[7], info[8]
			rankbug += float(contigcov)*float(finalscore)
		ranklist.append([org, rankbug])
	dictContigOrgRank[contig] = ranklist

#We will then see how many organisms with highlow annotations it takes to fill a contig.
lgtlist = []
for contig in dictContigOrgRank.keys():
	list_sorted = sorted(dictContigOrgRank[contig], key=lambda x:x[1], reverse=True)
	totalgroupsfilled = set()
	orglist = []
	for i in range(len(list_sorted)):
		org, score = list_sorted[i][0], list_sorted[i][1]
		label = dictContigOrgLabel[contig][org].keys()[0]
		if 'high' in label:
			groups = set(dictContigOrgLabel[contig][org][label]['high'])
			oldlen = len(totalgroupsfilled)
			totalgroupsfilled = totalgroupsfilled | groups
			newlen = len(totalgroupsfilled)
			if newlen > oldlen:
				orglist.append(org)
		if len(totalgroupsfilled) == len(dddictCGOS[contig].keys()):
			lgtlist.append([contig, orglist, 'complete'])
			break
		elif i == len(list_sorted)-1 and len(totalgroupsfilled) < len(dddictCGOS[contig].keys()):
			groupsleft = set(dddictCGOS[contig].keys()) - totalgroupsfilled
			lgtlist.append([contig, orglist, 'incomplete'])

#Print out the contig information/answer sheet
#Calculate order of organism
for contig, orglist, status in lgtlist:
	import string
	allLetters = string.uppercase
	order = ''		
	groupnumlist = []
	for group in dddictCGOS[contig].keys():
		groupnum = int(group.split('p')[1])
		groupnumlist.append(groupnum)
	groupssorted = sorted(groupnumlist)

	for number in groupssorted:
		groupname = 'Group' + str(number)
		for i in range(len(orglist)):
			org = orglist[i]
			label = dictContigOrgLabel[contig][org].keys()[0]
			if groupname in dictContigOrgLabel[contig][org][label]['high']:
				order += allLetters[i]
				break
			elif i == len(orglist)-1:
				order += '-'

	#Determine contig coverage and average scores for these contigs
	orginfolist = []
	completecoverage = 0
	for taxa in orglist:
		label = dictContigOrgLabel[contig][taxa].keys()[0]
		totalcoverage, totalgroupcov, totalpercID, totalscore, totalgrouplen, totalcomblen = 0, 0, 0, 0, 0, 0
		for groups in dictContigOrgLabel[contig][taxa][label]['high']:
			info = dddictCOGS[contig][taxa][groups]
			finalscore, finalpercID, finalgroupcov, contigcov = info[0], info[1], info[2], info[3]
                        newstart, newend, bp, grouplen, contiglen = info[4], info[5], info[6], int(info[7]), info[8]
                        #finalscore, finalpercID, finalgroupcov, contigcov = info[0], info[1], info[2], info[3]
                        #newstart, newend, combhitlen = info[4], info[5], info[6]
                        #grouplen, contiglen = info[7], info[8]
			
			totalcoverage += contigcov
			totalgroupcov += finalgroupcov #combhitlen
			totalgrouplen += grouplen
			totalcomblen += bp #combhitlen
			totalpercID += finalpercID*bp #combhitlen
			totalscore += finalscore*bp #combhitlen
		orginfo = [taxa, totalcoverage, str(totalgroupcov/len(dictContigOrgLabel[contig][taxa][label]['high'])), str(totalpercID/totalgrouplen), str(totalscore/totalgrouplen)]
		orginfolist.append(orginfo)
		completecoverage += totalcoverage
	print contig, order, status, completecoverage, ' '.join(orglist), ' '.join(str(i) for i in orginfolist)  
