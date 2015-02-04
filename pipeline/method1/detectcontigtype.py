#!/usr/bin/python

"""
This script will split all the scored contigs into:
1) 1 organism only contigs
2) 
"""

# import
import argparse
import json
import re

# arguments
parser = argparse.ArgumentParser()
parser.add_argument( '--dddictCGOS', help='Location and file of dddictCGOS json file.' )
parser.add_argument( '--dddictCOGS', help='Location and file of dddictCOGS json file.' )
parser.add_argument( '--delta', type=float, help='Upper level of confidence for BLAST hit annotation.' )
args = parser.parse_args()

# functions
def getgenetable(dictselectcontigs, dictallcontigs, taxalevel, phrase):
	genelist = []
	for contig in dictselectcontigs:
		for i in range(len(dictallcontigs[contig].keys())):
			groupname = 'Group' + str(i+1)
                	ranklist = []
                	for org in dddictCGOS[contig][groupname]:
                        	finalscore = dddictCGOS[contig][groupname][org][0]
                        	ranklist.append([finalscore, org])
                        	ranklist_sort = sorted(ranklist, key=lambda organism: organism[0], reverse=True)
                	for score, organism in ranklist_sort:
                        	info = dddictCGOS[contig][groupname][organism]
				searchstring = taxalevel + '__\w*'
				neworg = re.search(searchstring, organism).group() 
                        	finalscore, finalpercID, finalgroupcov, contigcov = info[0], info[1], info[2], info[3]
                        	newstart, newend, bp, grouplen, contiglen = info[4], info[5], info[6], int(info[7]), info[8]
				genelist.append([contig, contiglen, grouplen, bp, newstart, newend, groupname, neworg, finalscore, finalpercID, contigcov, finalgroupcov, phrase])
	return genelist


# import dddictCGOS and dddictCOGS
with open(args.dddictCGOS) as infile1:
        dddictCGOS = json.load(infile1)
with open(args.dddictCOGS) as infile2:
        dddictCOGS = json.load(infile2)

# detect taxa level
taxalevel = args.dddictCGOS.split('_')[1].replace('.json', '')

# Detect 1) One organism only 2) Multiple organisms across the whole contig 3) Potential LGT
all_contigset = set(dddictCGOS.keys())
oneorg_contigset = set()
highconforg_contigset = set()
dict_highconforg = {}
lgt_contigset = set()

genetable_annot = open('genetable' + taxalevel + '_annotated.txt', 'w')
infosheet = open('info.txt', 'a')

for contig in dddictCGOS.iterkeys():
        totalgroups = len(dddictCGOS[contig].keys())
        totalorgs = len(dddictCOGS[contig].keys())
        if totalorgs == 1: #These are contigs with only 1 organism explanation, which could have multiple taxon annotations
		oneorg_contigset.add(contig)
        else:
		highconf_orglist = []
		for org in dddictCOGS[contig].keys():
			scorelist = []
			for group in dddictCOGS[contig][org].keys():
				score = dddictCOGS[contig][org][group][0]
		                if score >= args.delta:
					scorelist.append(score)
			if len(scorelist) == totalgroups:
				highconforg_contigset.add(contig)
				highconf_orglist.append(org)
		if len(highconf_orglist) >= 1:
			dict_highconforg[contig] = highconf_orglist

#Remainder should be potential LGT
lgt_contigset = all_contigset - (oneorg_contigset | highconforg_contigset)

#Include new information in infosheet
infosheet.write('taxalevel' + '\t' + taxalevel + '\n')
infosheet.write('delta' + '\t' + str(args.delta) + '\n')
infosheet.write('number of contigs with BLAST annotations' + '\t' + str(len(all_contigset)) + '\t' + str(len(dddictCGOS)) + '\t' + str(len(dddictCOGS)) + '\n')
infosheet.write('oneorgonly' + '\t' + str(len(oneorg_contigset)) + '\n')
infosheet.write('highconforg' + '\t' + str(len(highconforg_contigset)) + '\n')
infosheet.write('potentiallgt' + '\t' + str(len(lgt_contigset)) + '\n')

#Print new genetable
oneorg_list = getgenetable(oneorg_contigset, dddictCGOS, taxalevel, 'oneorgonly')
highconforg_list = getgenetable(highconforg_contigset, dddictCGOS, taxalevel, 'highconforg')
lgt_list = getgenetable(lgt_contigset, dddictCGOS, taxalevel, 'lgt')

#Decide how to call organism annotations here for oneorgonly and highconforg
# oneorgonly - nothing really need be done
# highconforg, ideas include: 1) Summing across all groups to see which gives the highest set of scores for those that cover all
"""
for contig in dict_highconforg:
	orgscore_list = []
	for org in dict_highconforg[contig]:
		sumscore = 0
		for group in dddictCOGS[contig][org]:
			score = dddictCOGS[contig][org][group][0]
			sumscore += score
		avgscore = sumscore/len(dddictCOGS[contig][org].keys())
		orgscore_list.append([org, avgscore])
	orgscore_sorted = sorted(orgscore_list, key=lambda score: score[1], reverse=True)
	print contig, orgscore_sorted[0]
"""

#Currently, this genetable outputs everything
for entry in oneorg_list:
	genetable_annot.write('\t'.join(str(i) for i in entry ) + '\n')
for entry2 in highconforg_list:
	genetable_annot.write('\t'.join(str(j) for j in entry2) + '\n')
for entry3 in lgt_list:
	genetable_annot.write('\t'.join(str(k) for k in entry3) + '\n')

genetable_annot.close()
infosheet.close()

#Generate new dictionaries for detectLGT step
dddictCOGS_lgt = {}
dddictCGOS_lgt = {}
for contig in lgt_contigset:
        dddictCOGS_lgt[contig] = dddictCOGS[contig]
        dddictCGOS_lgt[contig] = dddictCGOS[contig]

filename = 'dddictCOGS_lgt_' + taxalevel + '.json'
with open(filename, "w") as outfile:
        json.dump(dddictCOGS_lgt, outfile)

filename2 = 'dddictCGOS_lgt_' + taxalevel + '.json'
with open(filename2, "w") as outfile2:
        json.dump(dddictCGOS_lgt, outfile2)

