#!/usr/bin/python

"""
"""

# import 
import argparse
import re

# arguments
parser = argparse.ArgumentParser()
parser.add_argument( '--answerkey', help='Answer key for synthetic contigs.' )
parser.add_argument( '--taxa', help='Taxa which we want to visualize contigs at.' )
args = parser.parse_args()

# taxa levels
leveldiff = ['t', 's', 'g', 'f', 'o', 'c', 'p', 'k']
levelnum = leveldiff.index(args.taxa)

# print header
print 'contig', 'contiglen', 'grouplen', 'groupstart', 'groupend', 'genelen', 'genestart', 'geneend', 'groupname', 'taxa', 'score', 'percID', 'contigcov', 'groupcov', 'status'

# read in answerkey
for astrline in open(args.answerkey):
	aastrline = astrline.strip().split('\t')
	contig, recipient, donor, taxadiff, contiglen, numgenes = aastrline[0], aastrline[2], aastrline[3], aastrline[4], aastrline[5], aastrline[6]
	geneinfo = aastrline[7]
	count = 0
	searchstring = args.taxa + '__\w*'

	status = ''
	answerlevel = leveldiff.index(taxadiff)
	if answerlevel >= levelnum:
		status = 'P'
	else:
		status = 'N' 

	for info in geneinfo.split('],'):
		myinfo = info.split(',')
		if re.search('R', myinfo[0]):
			genename = 'R' + taxadiff + status
			org = re.search(searchstring, recipient).group()
		else:
			genename = 'D' + taxadiff + status
			org = re.search(searchstring, donor).group()
		genelen, start, end = myinfo[1].strip(), myinfo[5].strip(), myinfo[6].strip(']').strip()
		print contig, contiglen, genelen, start, end, genelen, start, end, genename, org, 1, 1, 1, 1, 'answerkey'  
