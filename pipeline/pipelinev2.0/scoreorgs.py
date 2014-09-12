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
import json

#Define arguments
parser = argparse.ArgumentParser()
parser.add_argument('--blastoutput', help = 'Location and file of grouped and length-filtered BLAST output.')
parser.add_argument('--taxa', help = "Determine what taxonomy level to do this at; it should be one of the letters: 'k', 'p', 'c', 'o', 'f', 'g', 's'")
#parser.add_argument('delta', type = float, help = 'Determine the higher threshold score. Those above this score are high confidence BLAST hits')
#parser.add_argument('epsilon', type = float, help = 'Determine the lower threshold score. Those below this score are low confidence BLAST hits')
args = parser.parse_args()

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
	dictGroupScores = {}
	for group in dddictContigGroupOrgHits[contig].iterkeys():
		dictOrgScores = {}
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
			combhitlen = len(dictIndexScore)
			for base in dictIndexScore:
				sumpercID += dictIndexScore[base][0]
			finalpercID = sumpercID/combhitlen
			finalgroupcov = combhitlen/float(ddictContigGroupLen[contig][group])
			finalscore = finalpercID*finalgroupcov
			newstart, newend = sorted(dictIndexScore.keys())[0], sorted(dictIndexScore.keys())[combhitlen - 1]
		
			#Add to dictionaries
			dictOrgScores[org] = [finalscore, finalgroupcov, finalpercID, newstart, newend, combhitlen]
			dictGroupScores[group] = [finalscore, finalgroupcov, finalpercID, newstart, newend, combhitlen]
			ddictOrgGroupScores[org] = dictGroupScores
		ddictGroupOrgScores[group] = dictOrgScores
	dddictContigOrgGroupScores[contig] = ddictOrgGroupScores
	dddictContigGroupOrgScores[contig] = ddictGroupOrgScores

#Determine which contigs are:
#1) One gene only 2) Explained by one organism at high confidence
for contig in dddictContigGroupOrgScores.iterkeys():
	totalgroups = len(dddictContigGroupOrgScores[contig].keys())
	totalorgs = len(dddictContigOrgGroupScores[contig].keys())
	if totalgroups == 1: #These are contigs with only 1 group, which could have multiple taxon annotations
		print contig, hgtmodules.scoreOrgs(dddictContigOrgGroupScores, ddictContigGroupLen, contig)
		for org in dddictContigGroupOrgScores[contig]['Group1'].iterkeys():
			print contig, org, ' '.join(str(dddictContigGroupOrgScores[contig]['Group1'][org][i]) for i in range(len(dddictContigGroupOrgScores[contig]['Group1'][org]))), '1GroupOnly'
			#pass
	#else:
		#if totalorgs == 1: #These are contigs explained by 1 taxon, or strongly explained by at least 1 taxa across all genes.
			#print contig, hgtmodules.scoreOrgs(dddictContigOrgGroupScores, ddictContigGroupLen, contig)
			
		'''
			sumscore_len = 0
			totallen = 0
			for org in dddictContigOrgGroupScores[contig].iterkeys():
				numgroups = len(dddictContigOrgGroupScores[contig][org].keys())
				for group in dddictContigOrgGroupScores[contig][org].iterkeys():
					info = dddictContigOrgGroupScores[contig][org][group]
					sumscore_len += info[0]*info[5]
					totallen += info[5]
				score_groups = sumscore_len/float(totallen)
				print contig, org, score_groups
					

			
			pass
		else:
			for org in dddictContigOrgGroupScores[contig].iterkeys(): 
				count = 0
				if totalgroups == len(dddictContigOrgGroupScores[contig][org].keys()):
					for group in dddictContigOrgGroupScores[contig][org].iterkeys():
						score = dddictContigOrgGroupScores[contig][org][group]
						if score > 0.75:
							count += 1
				if totalgroups == count:
					#print contig, org, totalgroups, count
					pass
'''				
#This is for printing out information for R
#for contig in dddictOrgGroupScore.keys():
#	for org in dddictOrgGroupScore[contig].keys():
#		for group in dddictOrgGroupScore[contig][org].keys():
#			score = dddictOrgGroupScore[contig][org][group]
#			print contig, org, group, score[0], score[1]
