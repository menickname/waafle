#!/usr/bin/python

'''
This code will:
1a) Create a 3-tiered dictionary of {contig:{group:{org:score}}} and {contig{org:{group:score}}}. 
1b) Scores will be calculated based on percID*groupcoverage over that particular group. BLAST hits from the same taxon level will be aggregated. The combined percent ID will be calculated by averaging the percent IDs from the BLAST hit with the highest overall score at each base pair. The combined length will be calculated by the number of base pairs covered by all BLAST hits divided by the group length. 
2) We will then output a gene table that categorizes contigs by those with  A) 1 gene only B) Multiple genes covered by a single taxon C) Multiple genes completely covered by a high-scoring taxon D) Potential LGT events.
3) Those in category D above will continue on to the LGT algorithm.
'''

#Import
import argparse
import hgtmodules
import re
from Bio import SeqIO
import json

#Define arguments
parser = argparse.ArgumentParser()
parser.add_argument('--fasta', help = 'Location and file of contig fasta file.')
parser.add_argument('--blastoutput', help = 'Location and file of grouped and length-filtered BLAST output.')
parser.add_argument('--taxa', help = "Determine what taxonomy level to do this at; it should be one of the letters: 'k', 'p', 'c', 'o', 'f', 'g', 's'")
args = parser.parse_args()

#Create a dictionary consisting of contig total lengths.
dictContigLength = {}
for seq_record in SeqIO.parse(args.fasta, "fasta"):
         dictContigLength[seq_record.id] = len(seq_record.seq)

#Organize BLAST output into dictionaries
dictContigHits = hgtmodules.contigHits(open(args.blastoutput)) #Dictionary of {contig: [hits]}
ddictContigGroupHits = hgtmodules.contigGroupHits(dictContigHits) #Dictionary of {contig: {group: [hits]}}

#Create a dictionary consisting of contig group lengths.
dictContigGroupOrder = hgtmodules.sortGroups(ddictContigGroupHits)
ddictContigGroupLen = {}
for contig in dictContigGroupOrder.iterkeys():
        dictGroupLen = {}
        for i in range(len(dictContigGroupOrder[contig])):
                info = dictContigGroupOrder[contig][i]
                group, start, end = info[0], info[1], info[2]
                grouplen = end - start + 1
                dictGroupLen[group] = grouplen
        ddictContigGroupLen[contig] = dictGroupLen

#Generate 3-tiered dictionary which scores organisms within each group
dddictContigOrgGroupHits = {}
dddictContigGroupOrgHits = {}
for contig in ddictContigGroupHits.iterkeys():
	ddictGroupOrgHits = {}
	ddictOrgGroupHits = {}
	dictGroupHits = {}
	for group in ddictContigGroupHits[contig].iterkeys():
		dictOrgHits = {}
		for hit in ddictContigGroupHits[contig][group]: 
			organism = re.search(str('.*' + args.taxa + '__\w*'), hit[1].strip().split('|')[5]).group()
			dictOrgHits.setdefault(organism, []).append(hit)
			dictGroupHits.setdefault(group, []).append(hit)
			ddictOrgGroupHits[organism] = dictGroupHits 
		ddictGroupOrgHits[group] = dictOrgHits
	dddictContigGroupOrgHits[contig] = ddictGroupOrgHits
	dddictContigOrgGroupHits[contig] = ddictOrgGroupHits

#Score each unique taxon within each group
dddictContigOrgGroupScores = {}
dddictContigGroupOrgScores = {}			
for contig in dddictContigGroupOrgHits.iterkeys():
	ddictGroupOrgScores = {}
	ddictOrgGroupScores = {}
	for group in dddictContigGroupOrgHits[contig].iterkeys():
		dictOrgScores = {}
		dictGroupScores = {}
		for org in dddictContigGroupOrgHits[contig][group].iterkeys():
			dictIndexScore = {}
			for hit in dddictContigGroupOrgHits[contig][group][org]:
				start, end = int(hit[6]), int(hit[7])
				percID = float(hit[2])/100
				groupcov = (end - start + 1)/float(ddictContigGroupLen[contig][group])
				hitscore = percID*groupcov
				info = [hitscore, percID]

				#Use vote method to aggregate BLAST hits to the same taxon
				#At each base, add the percent ID of the hit with the highest overall score. 
				for iRange in range(start, end+1): 
					if dictIndexScore.get(iRange, 0) != 0:	
						oldscore = dictIndexScore[iRange][0]
						if oldscore < hitscore:
							dictIndexScore[iRange] = info
					else:
						dictIndexScore[iRange] = info
				
			#Calculate overall score for the taxon
			sumpercID = 0
			for base in dictIndexScore:
				sumpercID += dictIndexScore[base][1]
			finalpercID = sumpercID/len(dictIndexScore)
			finalgroupcov = len(dictIndexScore)/float(ddictContigGroupLen[contig][group])
			finalscore = finalpercID*finalgroupcov
			newstart, newend = sorted(dictIndexScore.keys())[0], sorted(dictIndexScore.keys())[len(dictIndexScore)-1]
			contigcov = len(dictIndexScore)/float(dictContigLength[contig])

			#Add to dictionaries
			dictOrgScores[org] = [finalscore, finalpercID, finalgroupcov, contigcov, newstart, newend, len(dictIndexScore), ddictContigGroupLen[contig][group], dictContigLength[contig]]
			dictGroupScores[group] = [finalscore, finalpercID, finalgroupcov, contigcov, newstart, newend, len(dictIndexScore), ddictContigGroupLen[contig][group], dictContigLength[contig]]
			ddictOrgGroupScores.setdefault(org, {}).update(dictGroupScores)
		ddictGroupOrgScores[group] = dictOrgScores
	dddictContigOrgGroupScores[contig] = ddictOrgGroupScores
	dddictContigGroupOrgScores[contig] = ddictGroupOrgScores

#temporarily disable for Eric's scoring method
"""
#Determine which contigs are:
#1) One organism only 2) Multiple organisms across the whole contig 3) Potential LGT
all_contigset = set(dddictContigOrgGroupScores.keys())
oneorg_contigset = set()
highconforg_contigset = set()
lgt_contigset = set()

genetable = open('genetable' + args.taxa + '.txt', 'w')

for contig in dictContigGroupOrder.iterkeys(): 
	totalgroups = len(dddictContigGroupOrgScores[contig].keys())
	totalorgs = len(dddictContigOrgGroupScores[contig].keys())
	if totalorgs == 1: #These are contigs with only 1 organism explanation, which could have multiple taxon annotations
		for i in range(len(dictContigGroupOrder[contig])):
			group = dictContigGroupOrder[contig][i][0]
			for org in dddictContigGroupOrgScores[contig][group].keys():
				oneorg_contigset.add(contig)
				newinfo = [contig, group, org, ' '.join(str(j) for j in dddictContigGroupOrgScores[contig][group][org]), '1orgonly', '\n']
				genetable.write(' '.join(newinfo))
	else:
		#Determine which contigs are explained fully by 1 or more organisms at high confidence
		for org in dddictContigOrgGroupScores[contig].keys():
			scorelist = []
			for group in dddictContigOrgGroupScores[contig][org].keys():
				score = dddictContigOrgGroupScores[contig][org][group][0]
				if score >= args.delta:
					scorelist.append(score) 
			if len(scorelist) == totalgroups:
				#print contig, org
				#print dddictContigOrgGroupScores[contig][org].keys()
				#print dddictContigGroupOrgScores[contig].keys()
				highconforg_contigset.add(contig)

#Print out the contigs which are explained at high confidence fully by 1 or more organisms
for contig in dictContigGroupOrder.iterkeys():
	if contig in highconforg_contigset:
		for i in range(len(dictContigGroupOrder[contig])):
			group = dictContigGroupOrder[contig][i][0]
			for org in dddictContigGroupOrgScores[contig][group].keys():
				newinfo = [contig, group, org, ' '.join(str(j) for j in dddictContigGroupOrgScores[contig][group][org]), '1+orghigh', '\n']
				genetable.write(' '.join(newinfo))

#Annotate remainder as potential LGT.
lgt_contigset = all_contigset - (oneorg_contigset | highconforg_contigset)
dddictCOGS = {}
dddictCGOS = {}
for contig in lgt_contigset:
	dddictCOGS[contig] = dddictContigOrgGroupScores[contig]
	dddictCGOS[contig] = dddictContigGroupOrgScores[contig]
	for i in range(len(dictContigGroupOrder[contig])):
		group = dictContigGroupOrder[contig][i][0]
		for org in dddictContigGroupOrgScores[contig][group].keys():
			newinfo = [contig, group, org, ' '.join(str(j) for j in dddictContigGroupOrgScores[contig][group][org]), 'potentialLGT', '\n']
			genetable.write(' '.join(newinfo))

genetable.close()
"""

#This is a temporary fix so I can use the old code to detect LGT; currently changed to original dictionaries for now
filename = 'dddictCOGS_' + args.taxa + '.json'
with open(filename, "w") as outfile:
        json.dump(dddictContigOrgGroupScores, outfile)

filename2 = 'dddictCGOS_' + args.taxa + '.json'
with open(filename2, "w") as outfile2:
	json.dump(dddictContigGroupOrgScores, outfile2)
