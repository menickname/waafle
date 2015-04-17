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
	bugnames_abb = []
	searchstring = args.taxa + '__\w*'
	for i in bugindex:
		bugName = index[i]
		bugnames_list.append(bugName)
		if bugName != 'unknown':
			bugnames_abb.append(re.search(searchstring, bugName).group())
		else:
			bugnames_abb.append('unknown')
	return onebugscore, bugnames_list, bugnames_abb, bugindex

def calcComplement ( aafTable, index ):
        complementlist = []
	bugnames_list = []
	bugnames_index = []
	searchstring = args.taxa + '__\w*'
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
	if org1name != 'unknown' and org2name != 'unknown':
		org1nameabb, org2nameabb = re.search(searchstring, org1name).group(), re.search(searchstring, org2name).group()
	elif org1name == 'unknown' and org2name != 'unknown':
		org1nameabb, org2nameabb = 'unknown', re.search(searchstring, org2name).group()
	else:
		org1nameabb, org2nameabb = re.search(searchstring, org1name).group(), 'unknown'
	return complementscore, [org1name, org2name], [org1nameabb, org2nameabb], [org1index, org2index]

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

def callDonorRecip ( aafTable, complementbugindex, orgindex):
	complementbugstatus = {}
	searchstring = args.taxa + '__\w*'
	recipient, donor, unknown = [], [], []
	for i in complementbugindex:
		status = findRankPerc(contigarray, i)
		if orgindex[i] != 'unknown':
			org = re.search(searchstring, orgindex[i]).group()
		else:
			org = 'unknown'
		complementbugstatus.setdefault(status, []).append(org)
	if 'recipient' in complementbugstatus:
		for r in complementbugstatus['recipient']:
			recipient.append(r)
	else: 
		recipient.append('NA')
	if 'donor' in complementbugstatus:
		for d in complementbugstatus['donor']:
			donor.append(d)
	else: 
		donor.append('NA')
	if 'unclear' in complementbugstatus:
		for u in complementbugstatus['unclear']:
			unknown.append(u)
	else:
		unknown.append('NA')
	return ';'.join(recipient), ';'.join(donor), ';'.join(unknown)
	

def getTopUniref ( unireflist ):
	reformat_list = []
	for element in unireflist.split(';'):
		uniref, value = element.split('-')[0], element.split('-')[1]
		reformat_list.append([uniref, value])
	topuniref = sorted(reformat_list, key=lambda score: score[1], reverse=True)[0]
	topuniref_final = str('-'.join(topuniref))
	return topuniref_final

def getGroupInfo ( dddictCOGS, contig, orglist, numgenes ):
	dictOrg_Info = {}
	for org in orglist:
		scovline, scov_modline, uniref90_line, uniref90_allline, uniref50_line, uniref50_allline = [], [], [], [], [], []
		for i in range(numgenes):
			group = "Group" + str(i+1)
			if group in dddictCOGS[contig][org]:
				info = dddictCOGS[contig][org][group]
				scov, scov_mod = str(info[3]), str(info[4])
				uniref90, uniref90_all, uniref50, uniref50_all = getTopUniref(info[11]), getTopUniref(info[12]), getTopUniref(info[13]), getTopUniref(info[14])
				scovline += [scov]
				scov_modline += [scov_mod]
				uniref90_line += [uniref90]
				uniref90_allline += [uniref90_all]
				uniref50_line += [uniref50]
				uniref50_allline += [uniref50_all]
			else:
				scovline += ['NA']
                                scov_modline += ['NA']
                                uniref90_line += ['NA']
                                uniref90_allline += ['NA']
                                uniref50_line += ['NA']
                                uniref50_allline += ['NA']
		dictOrg_Info[org] = [';'.join(scovline), ';'.join(scov_modline), ';'.join(uniref90_line), ';'.join(uniref90_allline), ';'.join(uniref50_line), ';'.join(uniref50_allline)]
	return dictOrg_Info
					

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
                orglist.append(org)

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

# print header
print '\t'.join(['contig', 'status', 'numbugs', 'numgenes', 'onebugscore', 'complementarity', 'organism', 'recipient', 'donor', 'unclear'])


for contig in dddictCOGS.keys():
	contigarray = contigtable[contig]
	contigorglist = contigorgindex[contig]
	
	bugs = contigarray.shape[0] #number of unique organisms
	genes = contigarray.shape[1] #number of genes
	
	searchstring = args.taxa + '__\w*'

	# calculate onebug scores
	onebugscore, oneorgname, oneorgname_abb, onebugindex = calcOneBug(contigarray, contigorglist)

	# remove those with only 1 gene
	if genes == 1:
		onegenescore = np.max(contigarray)
		onegeneindex = [i for i,j in enumerate(contigarray) if j==onegenescore]
	        onegenenames_list = []
		onegenenames_abb = []
        	for i in onegeneindex:
                	onegeneName = contigorglist[i]
	                onegenenames_list.append(onegeneName)
			if onegeneName != 'unknown':
				onegenenames_abb.append(re.search(searchstring, onegeneName).group())
			else:
				onegenenames_abb.append('unknown')
		print '\t'.join([contig, 'onegene', str(bugs), str(genes), str(onebugscore), 'NA', ';'.join(onegenenames_abb), 'NA', 'NA', 'NA'])
	else:
                # determine whether remaining contigs are better explained as pairs
                if bugs > 1:
	                complementscore, complementbugnames, complementbugnames_abb, complementbugindex = calcComplement(contigarray, contigorglist)
                else:
                        complementscore = 'NA'

		if onebugscore > args.onebug:
			pass # call onebug contig
			print '\t'.join([contig, 'oneorg', str(bugs), str(genes), str(onebugscore), str(complementscore), ';'.join(oneorgname_abb), 'NA', 'NA', 'NA'])
		else:
	
			# determine whether remaining contigs are better explained as pairs
			if onebugscore >= complementscore:
				pass # call onebug contig
				print '\t'.join([contig, 'bad-oneorg', str(bugs), str(genes), str(onebugscore), str(complementscore), ';'.join(oneorgname_abb), 'NA', 'NA', 'NA'])
				
			else:
				if complementscore > args.complement:
					pass # call LGT
					recipient, donor, unknown = callDonorRecip(contigarray, complementbugindex, contigorgindex[contig])
					print '\t'.join([contig, 'LGT', str(bugs), str(genes), str(onebugscore), str(complementscore), ';'.join(complementbugnames_abb), recipient, donor, unknown])
				else:
					pass # unexplained contig
					if onebugscore > complementscore:
						pass
						recipient, donor, unknown = callDonorRecip(contigarray, complementbugindex, contigorgindex[contig])
						print '\t'.join([contig, 'ambiguous-onebug', str(bugs), str(genes), str(onebugscore), str(complementscore), ';'.join(oneorgname_abb), recipient, donor, unknown])
					else:
						recipient, donor, unknown = callDonorRecip(contigarray, complementbugindex, contigorgindex[contig])
						print '\t'.join([contig, 'ambiguous-LGT', str(bugs), str(genes), str(onebugscore), str(complementscore), ';'.join(complementbugnames_abb), recipient, donor, unknown])
