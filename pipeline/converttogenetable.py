#!/usr/bin/python

'''
This script converts the json files to the gene table which can be plotted in R.
'''

# import
import argparse
import json
import re

# arguments
parser = argparse.ArgumentParser()
parser.add_argument( '--dddictCGOS', help='Json file that consists of a three tiered dictionary.' )
args = parser.parse_args()

# get taxa level
taxa = args.dddictCGOS.split('_')[1].replace('.json', '')
searchstring = taxa + '__\w*'

# read in dictionary
with open(args.dddictCGOS) as infile:
        dddictCGOS = json.load(infile)

# print it out
for contig in dddictCGOS.keys():
	#for group in dddictCGOS[contig].keys():
	for i in range(len(dddictCGOS[contig].keys())):
		groupname = 'Group' + str(i+1)
		ranklist = []
		for org in dddictCGOS[contig][groupname]:
			finalscore = dddictCGOS[contig][groupname][org][0]
			ranklist.append([finalscore, org])
			ranklist_sort = sorted(ranklist, key=lambda organism: organism[0], reverse=True)
		for score, organism in ranklist_sort:
			neworg = re.search(searchstring, organism).group()
			info = dddictCGOS[contig][groupname][organism]
			finalscore, finalpercID, finalgroupcov, contigcov = info[0], info[1], info[2], info[3]
			newstart, newend, bp, grouplen, contiglen = info[4], info[5], info[6], int(info[7]), info[8]
			print contig, contiglen, grouplen, bp, newstart, newend, groupname, neworg, finalscore, finalpercID, contigcov, finalgroupcov
			#print contig, group, neworg, ' '.join(str(info[i]) for i in range(len(info)))
