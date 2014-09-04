#!/bin/python

#Import 
import json
import argparse

#Filename
parser = argparse.ArgumentParser()
parser.add_argument('taxalevel', help = 'Taxon level you would like to print the score information')
args = parser.parse_args()
filename = 'dddictOrgGroupScore_' + args.taxalevel + '.json'

#Read in json for scores
with open(filename) as infile:
        dddictOGS = json.load(infile)

myfilename = 'scores_' + args.taxalevel + '.txt'
myfile = open(myfilename, 'w')
#Print contig, org, group, information
for contig in dddictOGS.keys():
	for org in dddictOGS[contig].keys():
		for group in dddictOGS[contig][org].keys():
			info = dddictOGS[contig][org][group]
			infoline = ' '.join(str(i) for i in info)
			newline = contig + ' ' + org + ' ' + group + ' ' + infoline + '\n'
			myfile.write(newline)

myfile.close()
