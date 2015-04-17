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
import collections

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
			organism = re.search(str('.*' + args.taxa + '__\w*'), hit[1].strip().split('|')[6]).group()
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

		#Add Uniref annotations for only BLAST hits that are aggregated
		uniref50list = []
		uniref90list = []
		#Add Uniref annotations for all BLAST hits in group
		uniref50_total = [] 
		uniref90_total = []
		#Track Uniref annotation per organism
		dict_OrgUniref50 = {}		
		dict_OrgUniref90 = {}

		#Score each organism by multiplying %ID*coverage_over_group
		for org in dddictContigGroupOrgHits[contig][group].iterkeys():
			dictIndexScore = {}
			for hit in dddictContigGroupOrgHits[contig][group][org]:
				start, end = int(hit[5]), int(hit[6])
				percID = float(hit[9])/100
				groupcov = (end - start + 1)/float(ddictContigGroupLen[contig][group])
				scov1, scov2 = float(hit[14]), float(hit[15])
				hitscore = percID*groupcov
				info = [hitscore, percID, scov1, scov2]
				uniref90, uniref50 = hit[1].split('|')[7], hit[1].split('|')[8]
				uniref90_total.append(uniref90)
				uniref50_total.append(uniref50)				
				dict_OrgUniref90.setdefault(org, set()).add(uniref90)
				dict_OrgUniref50.setdefault(org, set()).add(uniref50)	

				#Use vote method to aggregate BLAST hits to the same taxon
				#At each base, add the percent ID of the hit with the highest overall score. 
				add_button = False
				for iRange in range(start, end+1): 
					if dictIndexScore.get(iRange, 0) != 0:	
						oldscore = dictIndexScore[iRange][0]
						if oldscore < hitscore:
							dictIndexScore[iRange] = info
							if add_button == False:
								uniref90list.append(uniref90)
								uniref50list.append(uniref50)
								add_button = True
					else:
						dictIndexScore[iRange] = info
						if add_button == False:
							uniref90list.append(uniref90)
							uniref50list.append(uniref50)
							add_button = True
			
			#Calculate overall score for the taxon
			sumpercID, sumscov1, sumscov2 = 0, 0, 0
			for base in dictIndexScore:
				sumpercID += dictIndexScore[base][1]
				sumscov1 += dictIndexScore[base][2]
				sumscov2 += dictIndexScore[base][3]
			finalpercID = sumpercID/len(dictIndexScore)
			finalscov1 = sumscov1/len(dictIndexScore)
			finalscov2 = sumscov2/len(dictIndexScore)
			finalgroupcov = len(dictIndexScore)/float(ddictContigGroupLen[contig][group])
			finalscore = finalpercID*finalgroupcov
			newstart, newend = sorted(dictIndexScore.keys())[0], sorted(dictIndexScore.keys())[len(dictIndexScore)-1]
			contigcov = len(dictIndexScore)/float(dictContigLength[contig])

			#Add scores to dictionaries
			info = [finalscore, finalpercID, finalgroupcov, finalscov1, finalscov2, contigcov, newstart, newend, len(dictIndexScore), ddictContigGroupLen[contig][group], dictContigLength[contig], '', '', '', '']
			dictOrgScores[org] = info
			dictGroupScores[group] = info
			ddictOrgGroupScores.setdefault(org, {}).update(dictGroupScores)
		
		#Calculate Uniref scores per group
                dict_Uniref90scores = {}
                dict_Uniref90total_scores = {}
                dict_Uniref50scores = {}
                dict_Uniref50total_scores = {}
                uniref90c = collections.Counter(uniref90list)
                uniref90total_c = collections.Counter(uniref90_total)
                for uniprot90id in uniref90total_c:
                        if uniprot90id in uniref90c:
                                uniprot90score = uniref90c[uniprot90id]/float(len(uniref90list))
                                dict_Uniref90scores[uniprot90id] = uniprot90score
                        uniprot90_totalscore = uniref90total_c[uniprot90id]/float(len(uniref90_total))
                        dict_Uniref90total_scores[uniprot90id] = uniprot90_totalscore
                uniref50c = collections.Counter(uniref50list)
                uniref50total_c = collections.Counter(uniref50_total)
                for uniprot50id in uniref50total_c:
                        if uniprot50id in uniref50c:
                                uniprot50score = uniref50c[uniprot50id]/float(len(uniref50list))
                                dict_Uniref50scores[uniprot50id] = uniprot50score
                        uniprot50_totalscore = uniref50total_c[uniprot50id]/float(len(uniref50_total))
                        dict_Uniref50total_scores[uniprot50id] = uniprot50_totalscore

                #Assign Uniref scores to the right organism
                for orgs in dddictContigGroupOrgHits[contig][group].keys():
			u90list = []
	                u90_alllist = []
        	        u50list = []
                	u50_alllist = []
                        for id90 in dict_OrgUniref90[orgs]:
                                if id90 in dict_Uniref90scores:
					newinfo = str(id90) + '-' + str(dict_Uniref90scores[id90])
					u90list.append(newinfo)
				newinfo_all = str(id90) + '-' + str(dict_Uniref90total_scores[id90])
				u90_alllist.append(newinfo_all)
                        for id50 in dict_OrgUniref50[orgs]:
                                if id50 in dict_Uniref50scores:
					newinfo = str(id50) + '-' + str(dict_Uniref50scores[id50])
					u50list.append(newinfo)
				newinfo_all = str(id50) + '-' + str(dict_Uniref50total_scores[id50])
				u50_alllist.append(newinfo_all)
			
			# Add Uniref scores	
			ddictOrgGroupScores[orgs][group][len(info)-4] = ';'.join(u90list)
			ddictOrgGroupScores[orgs][group][len(info)-3] = ';'.join(u90_alllist)
			ddictOrgGroupScores[orgs][group][len(info)-2] = ';'.join(u50list)
			ddictOrgGroupScores[orgs][group][len(info)-1] = ';'.join(u50_alllist)
			
			dictOrgScores[orgs][len(info)-4] = ';'.join(u90list)
			dictOrgScores[orgs][len(info)-3] = ';'.join(u90_alllist)
			dictOrgScores[orgs][len(info)-2] = ';'.join(u50list)
			dictOrgScores[orgs][len(info)-1] = ';'.join(u50_alllist)

		#Add to remaining dictionaries
                ddictGroupOrgScores[group] = dictOrgScores
        dddictContigOrgGroupScores[contig] = ddictOrgGroupScores
        dddictContigGroupOrgScores[contig] = ddictGroupOrgScores
	

#This is a temporary fix so I can use the old code to detect LGT; currently changed to original dictionaries for now
filename = 'dddictCOGS_' + args.taxa + '.json'
with open(filename, "w") as outfile:
        json.dump(dddictContigOrgGroupScores, outfile)

filename2 = 'dddictCGOS_' + args.taxa + '.json'
with open(filename2, "w") as outfile2:
	json.dump(dddictContigGroupOrgScores, outfile2)

