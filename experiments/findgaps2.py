#!/usr/bin/env python

"""
This script contains a function that will determine the 'dippiness' of each organism in each contig.
It will also identify the regions that are not part of the dips, and quantify the average of those scores.
It should also determine what proportion of the contig is part of the rises and the dips: 
1. If more of the contig is part of "rises", it is likely a recipient.
2. If more of the contig is part of "dips", it is likely a donor.
"""

# import
import argparse
import numpy as np
import json
import re
import collections

# define args
parser = argparse.ArgumentParser()
parser.add_argument('--dict1', help='Dictionary of {contig: {org: {group: info}}}')
parser.add_argument('--dict2', help='Dictionary of {contig: {group: {org: info}}}')
parser.add_argument('--taxa', type=str, help='Taxa level at which to do this.')
parser.add_argument('--dropthresh', type=float, help='Threshhold above which to call a rise or drop transition.')
args = parser.parse_args()

# read in dictionaries for scored organisms
with open(args.dict1) as infile1:
        dddictCOGS = json.load(infile1)
with open(args.dict2) as infile2:
        dddictCGOS = json.load(infile2)

# generate bug x gene tables per contig
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

# function for calling gaps
def scorerecip(org_array):
	transitions = org_array[1:] - org_array[0:len(org_array)-1]
	dropstring = ''
	for i in range(len(transitions)):
		value = transitions[i]
		if value > args.dropthresh:
			dropstring += 'R'
		elif value < -args.dropthresh:
			dropstring += 'D'
		else:
			dropstring += 'C'
	if re.search('[RD]', dropstring):
		firsttrans = re.search('[RD]', dropstring).group()
		hllist = ['-']*(len(org_array))
		
		# Determine recipient score based on consecutive low and high regions
		if firsttrans == 'R':
			firsttransindex = dropstring.find('R')
			for i in range(firsttransindex + 1):
				hllist[i] = 'L'
		else:
			firsttransindex = dropstring.find('D')
			for i in range(firsttransindex + 1):
				hllist[i] = 'H'	
		for i in range(len(hllist)): 
			if re.search('[HL]', hllist[i]):
				continue
			else:
				transition = transitions[i-1]
				if transition > args.dropthresh:
					hllist[i] = 'H'
				elif transition < -args.dropthresh:
					hllist[i] = 'L'
				else:
					hllist[i] = hllist[i-1]
		recipscore = 0
		maxscore = 1
		for i in range(len(hllist)):
			status = hllist[i]
			if i == 0:
				recipscore += org_array[i]
			else:
				prevstatus = hllist[i-1]
				if status == prevstatus:
					recipscore += (org_array[i] + 0.3)
				else:
					recipscore += org_array[i]
				maxscore += 1.3
		print recipscore/maxscore
	else:
		recipscore = 0
		maxscore = 1
		for i in range(len(org_array)):
			score = org_array[i]
			recipscore += (score + 0.3)
			if i > 0:
				maxscore += 1.3
		print recipscore/maxscore, 'no transitions'
		"""
		for j in range(len(transitions)):
			prevtransition = transitionstring[len(transitionstring)-1]
			value = transitions[j]
			if value > args.dropthresh:
				transitionstring += 'r'
				recipscore += 1
			elif value < -args.dropthresh:
				transitionstring += 'd'
				recipscore += -1
			else:
				transitionstring += prevtransition
				if prevtransition == 'd': 
					recipscore += -0.5
				else:
					recipscore += 1.5
			print value, transitionstring, recipscore, 'adding'
		maxrecipscore = (len(transitionstring)-1)*1.5 + 1
		print transitionstring, recipscore, maxrecipscore, recipscore/maxrecipscore
	else:
		print dropstring, 'no transitions'
		"""

for contig in contigtable.keys():
	LGTcandidates = []
	nonLGTcandidates = []
	for i in range(len(contigtable[contig])):
		orgarray = contigtable[contig][i]
		orgname = dddictCOGS[contig].keys()[i]
		shortorgname = re.search(args.taxa + '__\w*', orgname).group()
		print contig, shortorgname, scorerecip(orgarray)
