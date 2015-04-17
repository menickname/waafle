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
parser.add_argument( '--contiggroupcoord', help='contiggroupcoord file created in blast2groups2.py' )
args = parser.parse_args()

# get taxa level
taxa = args.dddictCGOS.split('_')[1].replace('.json', '')
searchstring = taxa + '__\w*'

# read in dictionary
with open(args.dddictCGOS) as infile:
        dddictCGOS = json.load(infile)

# read in contig group coordinates
dict_contiggroupcoord = {}
dict_groupcoord = {}
fcontig = ''
contiggroupcoord = open(args.contiggroupcoord, 'r').readlines()

for i in range(len(contiggroupcoord)): 
	aastrline = contiggroupcoord[i].strip().split('\t')
	scontig, grouplen, groupname, groupstart, groupend = aastrline[0], aastrline[1], aastrline[2], float(aastrline[3]), float(aastrline[4])
	if scontig != fcontig:
		if len(dict_groupcoord) != 0:
			dict_contiggroupcoord[fcontig] = dict_groupcoord
		dict_groupcoord = {}
	dict_groupcoord[groupname] = [groupstart, groupend]
	fcontig = scontig
	if i == len(contiggroupcoord)-1:
		dict_contiggroupcoord[fcontig] = dict_groupcoord

# print header
print '\t'.join(['contig', 'contiglen', 'grouplen', 'groupstart', 'groupend', 'genelen', 'genestart', 'geneend', 'groupname', 'taxa', 'uniref90', 'uniref90_all', 'uniref50', 'uniref50_all', 'score', 'percID', 'contigcov', 'scov', 'scov_mod', 'group_cov', 'status'])

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
			#print contig, info
			finalscore, finalpercID, finalgroupcov, scov1, scov2, contigcov = info[0], info[1], info[2], info[3], info[4], info[5]
			newstart, newend, bp, grouplen, contiglen = info[6], info[7], info[8], int(info[9]), info[10]
			groupstart, groupend = int(dict_contiggroupcoord[contig][groupname][0]), int(dict_contiggroupcoord[contig][groupname][1])
			uniref90, uniref90_all, uniref50, uniref50_all = info[11], info[12], info[13], info[14]
			print '\t'.join([contig, str(contiglen), str(grouplen), str(groupstart), str(groupend), str(bp), str(newstart), str(newend), str(groupname), neworg, uniref90, uniref90_all, uniref50, uniref50_all, str(finalscore), str(finalpercID), str(contigcov), str(scov1), str(scov2), str(finalgroupcov), 'waafle'])
			#print contig, group, neworg, ' '.join(str(info[i]) for i in range(len(info)))
