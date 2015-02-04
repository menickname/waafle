#!/usr/bin/env python

"""
This script calls LGT from contigs based on several sets of scores:
1) The "onebugscore" determines which contigs are likely explained by 1 bug
2) The "complementscore" determins which contigs are likely explained by 2 bugs.
3) For all 2 bug contigs, we try to call the donor and recipient.

We then print out:
contigname, onebugscore/complementscore, bugName, donorrecipient_status
-----------------------------------------------
Author: Tiffany Hsu (tyhsu389@gmail.com)
"""

# import
import argparse
import numpy as np
import json
import re


# arguments
parser = argparse.ArgumentParser()
parser.add_argument( '--dddictCOGS', help='Location and file of json file with {contig:{org:group}}.' )
parser.add_argument( '--dddictCGOS', help='Location and file of json file with {contig:{group:org}}.' )
parser.add_argument( '--onebug', type=float, help='Value above which the contig is likely explained by one bug.' )
parser.add_argument( '--complement', type=float, help='Value above which we call LGT, value below we call no LGT, or unexplainable contig.' )
parser.add_argument( '--taxa', help='Level taxa we are doing this at.' )
args = parser.parse_args()


# functions
def calcOneBug ( aafTable, index ):
	onebugscore = np.max(np.amin(aafTable, axis=1))
        bugindex = [i for i,j in enumerate(np.amin(aafTable, axis=1)) if j==onebugscore]
	bugnames_list = []
	for i in bugindex:
		bugName = index[i]
		bugnames_list.append(bugName)
	return onebugscore, bugnames_list, bugindex

def calcComplement ( aafTable, index ):
        complementlist = []
	bugnames_list = []
	bugnames_index = []
	bugs = aafTable.shape[0]
        for i in range(bugs-1):
                for j in range(i+1, bugs):
                        org1 = aafTable[i]
                        org2 = aafTable[j]
                        twobugarray = np.array([org1, org2])
		        complementarity = np.min(np.amax(twobugarray, axis=0))
                        complementlist.append([complementarity, i, j])
	complementlist_sorted = sorted(complementlist, reverse=True)
        complementscore = complementlist_sorted[0][0]
	org1index, org2index = complementlist_sorted[0][1], complementlist_sorted[0][2]
	org1name, org2name = index[org1index], index[org2index]
	return complementscore, [org1name, org2name], [org1index, org2index]

def findRankPerc ( aafTable, bugindex ):
	# this is using the rank percentile
	bugarray = aafTable[bugindex]
	bugarray_sorted = sorted(aafTable[bugindex])
	buggaparray = np.array(bugarray_sorted[1:]) - np.array(bugarray_sorted[:len(bugarray_sorted)-1])
	percentile = (list(buggaparray).index(max(buggaparray))+1)/float(len(bugarray))
	if percentile < 0.5:
		# call recipient
		return 'recipient'
	elif percentile > 0.5: 
		# call donor
		return 'donor'
	else:
		# call ambiguous
	 	return 'unclear'

def callDonorRecip ( aafTable, complementbugindex):
	complementbugstatus = []
	for i in complementbugindex:
		status = findRankPerc(contigarray, i)
                complementbugstatus.append(status)
	return complementbugstatus


def outputInfo ( contig, score1, score2, nameslist, indexlist, label, bugs, genes):
	# this function should give us the info we need to print out to screen
	print contig, label, score1, score2, bugs, genes, ' '.join(nameslist), ' '.join(str(i) for i in indexlist)


# Script for detecting LGT
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
                        groupnum = int(re.search('[0-9]+', group).group())
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
		searchstring = args.taxa + '__\w*'
                orgname = re.search(searchstring, org).group()
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


# decide which contigs:
# 1) Have only 1 gene 
# 2) Are likely to be covered by 1 bug
# 3) Are likely to be covered by 2 bugs (calculate complementarity scores)
# 4) And then, for those with LGT, determine the donor and recipient
for contig in dddictCOGS.keys():
	contigarray = contigtable[contig]
	contigorglist = contigorgindex[contig]
	
	bugs = contigarray.shape[0]
	genes = contigarray.shape[1]

	# remove those with only 1 gene
	if genes == 1:
		onegenescore = np.max(contigarray)
		onegeneindex = [i for i,j in enumerate(contigarray) if j==onegenescore]
	        onegenenames_list = []
        	for i in onegeneindex:
                	onegeneName = contigorglist[i]
	                onegenenames_list.append(onegeneName)
		outputInfo(contig, onegenescore, 0, onegenenames_list, onegeneindex, 'onegene', bugs, genes)
	else:

		# remove those likely to be one bug organisms
		onebugscore, oneorgname, onebugindex = calcOneBug(contigarray, contigorglist)

                # determine whether remaining contigs are better explained as pairs
                if bugs > 1:
	                complementscore, complementbugnames, complementbugindex = calcComplement(contigarray, contigorglist)
                else:
                        complementscore = 0

		if onebugscore > args.onebug:
			pass # call onebug contig
			outputInfo(contig, onebugscore, complementscore, oneorgname, onebugindex, 'oneorg', bugs, genes)
		else:
			# See if  LGT
			if complementscore > args.complement:
				pass # call LGT
				complementbugstatus = callDonorRecip(contigarray, complementbugindex)
				outputInfo(contig, onebugscore, complementscore, complementbugnames, complementbugstatus, 'LGT', bugs, genes)
			else:
				pass # unexplained contig
				if onebugscore > complementscore:
					outputInfo(contig, onebugscore, complementscore, oneorgname, onebugindex, 'ambiguous-onebug', bugs, genes)
				else:
					complementbugstatus = callDonorRecip(contigarray, complementbugindex)
					outputInfo(contig, onebugscore, complementscore, complementbugnames, complementbugstatus, 'ambiguous-LGT', bugs, genes)
