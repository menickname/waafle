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
parser.add_argument( '--taxa', help='Level of taxa pipeline has been run at.')
args = parser.parse_args()

# read in dictionary
with open(args.dddictCGOS) as infile:
        dddictCGOS = json.load(infile)

searchstring = args.taxa + '__\w*'

# print it out
for contig in dddictCGOS.keys():
	for group in dddictCGOS[contig].keys():
		for org in dddictCGOS[contig][group]:
			neworg = re.search(searchstring, org).group()
			info = dddictCGOS[contig][group][org]
			print contig, group, neworg, ' '.join(str(info[i]) for i in range(len(info)))
