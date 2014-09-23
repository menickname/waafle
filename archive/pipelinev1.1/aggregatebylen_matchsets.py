#!/usr/bin/python

'''
This script will 
'''

#Import
import argparse
import re
import json
import collections
import copy

#Define arguments
parser = argparse.ArgumentParser()
#parser.add_argument('dict_output1', help = 'Location and file of dictionary of contigs, groups, organisms, and scores')
parser.add_argument('dict1', help = 'JSON file for dictionary containing {Contigs: {Org: {Group: Score}}}')
parser.add_argument('dict2', help = 'JSON file for dictionary containing {Contigs: {Group: {Org: Score}}}')
parser.add_argument('delta', type = float, help = 'Determine the higher threshold score. Those above this score are high confidence BLAST hits')
parser.add_argument('epsilon',type = float, help = 'Determine the lower threshold score. Those below this score are low confidence BLAST hits')
parser.add_argument('length', type = float, help = 'Length cutoff for group to be included')
args = parser.parse_args()


#Read in dictionaries, define delta and epsilon
with open(args.dict1) as infile1:
	dddictOrgGroupScore = json.load(infile1)
with open(args.dict2) as infile2:
	dddictGroupOrgScore = json.load(infile2)
with open('dictContigLength.json') as infile3:
	dictContigLength = json.load(infile3)
with open('ddictContigGroupLen.json') as infile4:
	ddictContigGroupLen = json.load(infile4)
delta, epsilon = args.delta, args.epsilon


#Classify each organism within each contig into categories of annotation
mod_dddictOGS = copy.deepcopy(dddictOrgGroupScore) #Create copies of each dictionary
mod_dddictGOS = copy.deepcopy(dddictGroupOrgScore)
#mod_ddictContigGroupLen = copy.deepcopy(ddictContigGroupLen)


#Remove all groups that are below a certain length
for contig in dddictOrgGroupScore.keys():
	#print contig
	for org in dddictOrgGroupScore[contig].keys():
		for group in dddictOrgGroupScore[contig][org].keys():
			#info = dddictOrgGroupScore[contig][org][group]
			#score, grouplen, start, end = info[0], float(info[1]), info[2], info[3]
			if float(ddictContigGroupLen[contig][group]) < args.length:
				del mod_dddictOGS[contig][org][group] #This deletes the groups within each organism with len < args.length
	for org2 in dddictOrgGroupScore[contig].keys(): #Removes the organisms that no longer are in any groups
		if org2 in mod_dddictOGS[contig].keys():
			if len(mod_dddictOGS[contig][org2].keys()) == 0:
				del mod_dddictOGS[contig][org2]
	for group2 in dddictGroupOrgScore[contig].keys():
		if float(ddictContigGroupLen[contig][group2]) < args.length:
			del mod_dddictGOS[contig][group2]

	
dictCutLen = {} #This dictionary keeps track of the number of base pairs covered by all groups >args.length across the contig. This should be the maximum number of bp that can be covered. Note that overlapping groups are not accounted for, so it could potentially be just >100% coverage.
dictCoverage = {} #This dictionary will keep track of the total coverage all groups that are >args.length make across the contig
for contig in ddictContigGroupLen.keys():
	totallen = 0
	cutlen = 0
	for group in ddictContigGroupLen[contig].keys():
		#totallen += float(ddictContigGroupLen[contig][group])
		if float(ddictContigGroupLen[contig][group]) >= args.length:
			cutlen += float(ddictContigGroupLen[contig][group])
	totalcoverage = cutlen/dictContigLength[contig]
	#print contig, cutlen, dictContigLength[contig], totalcoverage #Use this to print overall contig coverage
	dictCoverage[contig] = totalcoverage
	dictCutLen[contig] = cutlen


#Use new dictionary, but create another copy
mod_dddictOGS2 = copy.deepcopy(mod_dddictOGS)
mod_dddictGOS2 = copy.deepcopy(mod_dddictGOS)

myfile1 = open('hgtnegative_highonly.txt', 'w')
#Assign scores to all alignments, continue modifying the dictionary copies
dictContigOrgLabel = {}
for contig in mod_dddictOGS.keys():
	dictOrgLabel = {}
	highorglist = []
	#print contig
	for org in mod_dddictOGS[contig].keys():
		#print contig, org
		dictLabel = {}
		dictQuality = {}
		for group in mod_dddictGOS[contig].keys(): #We still want to assign every single group an annotation
			#print contig, org, mod_dddictOGS[contig][org].keys()
			#print mod_dddictGOS[contig].keys()
			if group not in mod_dddictOGS[contig][org]:
				dictQuality.setdefault('low', []).append(group)
			else: 
				info = mod_dddictOGS[contig][org][group]
				score, length, start, end = float(info[0]), info[1], info[2], info[3]
				if score > delta: #These are high quality hits
					dictQuality.setdefault('high', []).append(group) 
				elif epsilon <= score <= delta: #These are medium quality hits
					dictQuality.setdefault('med', []).append(group)
				elif score < args.epsilon: #These are low quality hits
					dictQuality.setdefault('low', []).append(group)
		#print contig, org, dictQuality
		label = ""
		highonly, medonly, lowonly = set('high'), set('med'), set('low')
		highlow, highmed, medlow, highmedlow = set(['high', 'low']), set(['high', 'med']), set(['med', 'low']), set(['high', 'med', 'low'])
		if len(set(dictQuality.keys())) == 1:
			if 'high' in dictQuality:
				#print contig, org
				label = 'highonly'
				orglen = 0
				orgscore = 0
				groupnum = ''
				#highorglist = []
				for group3 in mod_dddictOGS[contig][org].keys():
					highonlyinfo = mod_dddictOGS[contig][org][group3]
					orglen += float(highonlyinfo[1])
					orgscore += float(highonlyinfo[0])
					groupnum += 'A'
				highonlycov = orglen/dictContigLength[contig]
				highonlyscore = orgscore/len(mod_dddictOGS[contig][org].keys())
				overallscore = highonlycov*highonlyscore
				highorginfo = overallscore, org, highonlyscore, highonlycov, groupnum
				highorglist.append(highorginfo)
				#myfile1.write(contig + '\t' + str(dictContigLength[contig]) + '\t' + org + '\t' + str(highonlyscore) + '\t' + str(highonlycov) + '\t' + 'nohgt' + '\t' + 'highonly' + '\n')
				#print contig#, dictQuality #Tells us all the contigs that are multiple genes, 1 org high_only
				if contig in mod_dddictOGS2:
					del mod_dddictOGS2[contig] #Remove contigs associated with 'highonly'
				if contig in mod_dddictGOS2:
					del mod_dddictGOS2[contig]
				#continue
			elif 'med' in dictQuality:
				label = 'medonly'
			elif 'low' in dictQuality:
				label = 'lowonly'
				if contig in mod_dddictOGS2:
					del mod_dddictOGS2[contig][org]  #Remove organisms associated with 'lowonly'
					#print contig, org, mod_dddictOGS[contig][org].keys(), '\n'
				if contig in mod_dddictGOS2:
					for group2 in mod_dddictGOS[contig].keys():
						if group2 in mod_dddictGOS2[contig].keys():
							if org in mod_dddictGOS2[contig][group2].keys():
								del mod_dddictGOS2[contig][group2][org]
		else:
			if set(dictQuality.keys()).issubset(highmed):
				label = 'highmed'
			elif set(dictQuality.keys()).issubset(highlow):
				label = 'highlow'
			elif set(dictQuality.keys()).issubset(medlow):
				label = 'medlow'
			elif set(dictQuality.keys()).issubset(highmedlow):
				label = 'highmedlow'
	
		#Assign each contig into a set label, and include which groups are in which labels
		dictLabel.setdefault(label, {}).update(dictQuality) #This dictionary, as well as dictOrgLabel and dictContigOrgLabel will contain everything. Use these for reference, but loop through dddictOGS or dddictGOS instead.
		dictOrgLabel[org] = dictLabel
		#print contig, org, label, dictLabel
	dictContigOrgLabel[contig] = dictOrgLabel
	#print contig, highorglist
	if len(highorglist) != 0:
		highorglist.sort()
		myfile1.write(contig + '\t' +  str(dictContigLength[contig]) + '\t' + highorglist[0][1] + '\t' + str(highorglist[0][2]) + '\t' + str(highorglist[0][3]) + '\t' + highorglist[0][4] + '\t' + 'highonly' + '\t' + 'nohgt' + '\n')
	
		
	
	#Remove all groups that do not have any more organisms
	if contig in mod_dddictGOS2:
		for group3 in mod_dddictGOS2[contig].keys():
			if group3 in mod_dddictGOS2[contig]:
				if len(mod_dddictGOS2[contig][group3].keys()) == 0:
					del mod_dddictGOS2[contig][group3]


#Find all contigs with highlow annotations
contigSet = set()
for contig in mod_dddictOGS2.keys():
	count = 0
	for org in mod_dddictOGS2[contig].keys():
		for label in dictContigOrgLabel[contig][org].keys():
			if label == 'highlow': #or label == 'highmedlow':
				count = count + 1
	if count >= 2:
		contigSet.add(contig)


#Decide which bugs will represent the contig before...
#...calculating the coverage of all contigs within the contigSet
dictRank = {}
for contig in contigSet:
	for org in mod_dddictOGS2[contig].keys():
		rankbug, scorebug1, scorebug2, coveragebug = 0, 0, 0, 0
		dictInfo = {}
		count = 0
		for group in mod_dddictOGS2[contig][org].keys():
			info = mod_dddictOGS2[contig][org][group]
			score, org_len, start, end = float(info[0]), float(info[1]), float(info[2]), float(info[3])
			label = dictContigOrgLabel[contig][org].keys()[0]
			delta_change = score
			#if score > delta:
			#	delta_change = 1
			#	count += 1
			#else:
			#	delta_change = 0
			rankbug += org_len*delta_change #This is weighted for %ID and length. Otherwise we end up ignoring the medonly.
			#scorebug1 += score*delta_change
			scorebug2 += score
			coveragebug += org_len
			#print contig, org, group, org_len, scorebug2, scorebug2
		#if count > 0:
		#	avgscorebug1 = scorebug1/count
		#else:
		#	avgscorebug1 = 0
		avgscorebug2 = scorebug2/len(mod_dddictOGS2[contig][org].keys())
		#highcoveragebug = rankbug/dictContigLength[contig]
		avgcoveragebug = coveragebug/dictContigLength[contig]
		#print contig, org, rankbug, avgscorebug1, avgscorebug2, avgcoveragebug, dictCutLen[contig], dictContigLength[contig], dictCoverage[contig]
		info2  = org, rankbug, avgscorebug2, avgcoveragebug, dictCoverage[contig]
		if rankbug != 0:
			#print contig, org, rankbug, scorebug1, avgscorebug, avgcoveragebug, dictCoverage[contig]
			dictRank.setdefault(contig, []).append(info2)

from operator import itemgetter, attrgetter
finalL = []
for contig in dictRank:
	nextContig = False
	setUnion = set()
	orgHGT= []
	orgList_sorted = sorted(dictRank[contig], key=itemgetter(1,2), reverse=True) #Essentially, this tells us what is covered across in base pairs
	
	#This tells us which organisms are needed to "fill the contig". It starts from the most prevalent and moves to the least
	groupone = set()
	grouptwo = set() 
	for i in range(len(orgList_sorted)):
		if nextContig == True:
			continue
		org1 = orgList_sorted[i][0]
		allGroups = set(mod_dddictGOS2[contig].keys())
		group1 = set(mod_dddictOGS2[contig][org1].keys())
		label = dictContigOrgLabel[contig][org1].keys()[0]
		'''if 'high' in label:
			for j in range(len(dictContigOrgLabel[contig][org1][label]['high'])):
				groupone.add(dictContigOrgLabel[contig][org1][label]['high'][j])
		if 'med' in label:
			for j in range(len(dictContigOrgLabel[contig][org1][label]['med'])):
				groupone.add(dictContigOrgLabel[contig][org1][label]['med'][j])'''
		if i == 0:
			setUnion = group1
			orgHGT.append(org1)
			#setUnion = groupone
		if i < len(orgList_sorted)-1:
			org2 = orgList_sorted[i+1][0]
			group2 = set(mod_dddictOGS2[contig][org2].keys())
			label2 = dictContigOrgLabel[contig][org2].keys()[0]
			'''if 'high' in label2:
				for j in range(len(dictContigOrgLabel[contig][org2][label2]['high'])):
					grouptwo.add(dictContigOrgLabel[contig][org2][label2]['high'][j])
			if 'med' in label2:
				for j in range(len(dictContigOrgLabel[contig][org2][label2]['med'])):
					grouptwo.add(dictContigOrgLabel[contig][org2][label2]['med'][j])'''
			oldLength = len(setUnion)
			setUnion = setUnion | group2
			#setUnion = setUnion | grouptwo
			newLength = len(setUnion)
			#print newLength - oldLength
			if newLength - oldLength >= 1:
				orgHGT.append(org2)
			newSet = allGroups - setUnion
			if len(newSet) == 0:
				nextContig = True
	
	if len(orgHGT) > 1:
		#Using the list of most prevalent, we will determine the organism order of ABA, or AAB, etc.
		allGroups =  mod_dddictGOS2[contig].keys()
		groupList = []
		import string
		allLetters = string.uppercase
		for group in allGroups:
			num = float(re.search('[0-9]+', group).group()) 
			info3 = num, group
			groupList.append(info3)
		groupList_sorted = sorted(groupList)
		orgOrder = ''
		for group2 in groupList_sorted:
			orgList = mod_dddictGOS2[contig][group2[1]].keys()
			for i in range(len(orgHGT)):
				orgLetter = allLetters[i]
				if orgHGT[i] in orgList:
					orgOrder = orgOrder + orgLetter
					#print orgHGT, contig, group2, orgHGT[i]
					break
		orgStatus = '' #This tells us where the LGT is within the contig
		if orgOrder[0] == orgOrder[len(orgOrder)-1]:
			orgStatus = 'mid'
		else:
			orgStatus = 'end'
		
		#Format for results	
		finalList = []
		finalList.append(contig)
		finalList.append(dictContigLength[contig])
		finalList.append(len(orgHGT))
		for finalorg in orgHGT:
			for j in range(len(orgList_sorted)):
				#print orgList
				if finalorg == orgList_sorted[j][0]:
					finalList.append(finalorg) 
					finalList.append(float(orgList_sorted[j][2])) #Change these according to info2
					finalList.append(float(orgList_sorted[j][3]))
		finalList.append(dictCoverage[contig])
		finalList.append(orgOrder)
		finalList.append(orgStatus)
		finalList.append('hgt')
		finalL.append(finalList)
	finalL.sort(key=lambda x: x[len(x)-4], reverse=True)

for astrline in finalL:
	print '\t'.join(str(i) for i in astrline)

myfile1.close()
'''
#The below is to print out graph formats
#To print for up/down graphs
for contig in contigSet:
	maxSum = 0
	orgSet = set()
	for org in mod_dddictOGS[contig].keys():
		for group in dddictGroupOrgScore[contig]:
			if group not in mod_dddictOGS[contig][org]:
				print contig, group, org, 0, 0
			else:
				print contig, group, org, mod_dddictGOS[contig][group][org][0], mod_dddictGOS[contig][group][org][1]
	print contig, len(mod_dddictOGS[contig].keys()), mod_dddictOGS[contig]

'''
'''
#To print for line graphs
for contig in contigSet:
	maxSum = 0
	orgSet = set()
	for group in mod_dddictGOS[contig].keys():
		for org in mod_dddictGOS[contig][group]:
			print contig, org, group, mod_dddictGOS[contig][group][org][0], mod_dddictGOS[contig][group][org][2], mod_dddictGOS[contig][group][org][3]
'''
