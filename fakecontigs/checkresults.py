#!/usr/bin/python

#Import
import argparse
import re
import collections

#Define arguments
parser = argparse.ArgumentParser()
parser.add_argument('--allresults', help = 'Location and file of results')
parser.add_argument('--answerkey', help = 'Location and file of contig answer key')
parser.add_argument('--taxalevel', help = 'Taxonomy level you wanted to search for in LGT.')
args = parser.parse_args()

#Read the answerkey into a dictionary
dictResults = {}
for astrline in open(args.answerkey):
	aastrline = astrline.strip().split('\t')
	contig = aastrline[0]
	donortaxa, recipienttaxa = aastrline[2], aastrline[3]
	matchlevel = args.taxalevel + '__\w+'
	donor, recipient = re.search(matchlevel, donortaxa).group(), re.search(matchlevel, recipienttaxa).group()
	orgSet = set([donor,recipient])
	contiglen = aastrline[7]
	taxadiff = aastrline[4]
	hgtstatus = aastrline[5]
	info = hgtstatus, taxadiff, contiglen, donor, recipient
	dictResults[contig] = info
	#print hgtstatus

#Indicate the taxa levels: this is so we can compare and re-classify what is and is not LGT
taxalevels = {'p': 0, 'c': 1, 'o': 2, 'f': 3, 'g': 4, 's':5}


cnt = collections.Counter()
cnt2 = collections.Counter()
cnt3 = collections.Counter()
listresults = []
list_answerkeyresults = []
list_algorithmresults = []
myendresults = open('answerkeycomparison.txt', 'w')
myendresults.write('\t'.join(['contig', 'contiglen', 'phyldist', 'result', 'answerkey_result', 'algorithm_result', 'organism', 'answerkey_organism', 'algorithm_organism', 'algorithm_call']) + '\n')

#Iterate through the results and determine whether each contig is a TN, TP, FN, and FP
for bstrline in open(args.allresults):
	result = ''
	organism = ''
	bbstrline = bstrline.strip().split('\t')
	contig, contiglen, hgt_status = bbstrline[0], bbstrline[1], bbstrline[len(bbstrline)-1]

	#Get answerkey for that contig
	result_status, phyldist, result_len, result_donor, result_recipient = dictResults[contig][0], dictResults[contig][1], dictResults[contig][2], dictResults[contig][3], dictResults[contig][4]
	if hgt_status == 'nohgt':
		org = bbstrline[2]
		algorithmcall = bbstrline[5]
		if hgt_status == result_status:
			result, answerkey_result, algorithm_result = 'TN', result_status, hgt_status
			if result_recipient == org:
				organism, answerkey_organism, algorithm_organism = '1/1', result_recipient, org
			else:
				organism, answerkey_organism, algorithm_organism = '0/1', result_recipient, org
			
		else:
			result, answerkey_result, algorithm_result = 'FN', result_status, hgt_status
			resultorgSet = set([result_recipient, result_donor])
			if len(resultorgSet & set([org])) == 1:
				organism, answerkey_organism, algorithm_organism = '1/2', result_recipient + ';' + result_donor, org
			else:
				organism, answerkey_organism, algorithm_organism = '0/2', result_recipient + ';' + result_donor, org
	else:
		#First, reclassify those that are 'HGT' as not based on taxa levels.
		if taxalevels[args.taxalevel] < phyldist:
			#This contig becomes nohgt
		
		#Second, those that are no longer HGT due to the length requirement should be reclassified.
		elif 
		orgList = []
		algorithmcall = bbstrline[len(bbstrline)-3]
		startcoord = 0
		for i in range(int(bbstrline[2])):
			startcoord += 3
			orgList.append(bbstrline[startcoord])
		orgSet = set(orgList)
		orgLine = ';'.join(orgList)
		if hgt_status == result_status:
			result, answerkey_result, algorithm_result = 'TP', result_status, hgt_status
			if len(orgSet & set([result_recipient, result_donor])) >= 2:
				organism, answerkey_organism, algorithm_organism = '2/2', result_recipient + ';' + result_donor, orgLine
			if len(orgSet & set([result_recipient, result_donor])) == 1:
				organism, answerkey_organism, algorithm_organism = '1/2', result_recipient + ';' + result_donor, orgLine
			if len(orgSet & set([result_recipient, result_donor])) == 0:
				organism, answerkey_organism, algorithm_organism = '0/2', result_recipient + ';' + result_donor, orgLine
		else:	
			result, answerkey_result, algorithm_result = 'FP', result_status, hgt_status
                        if len(orgSet & set([result_recipient, result_donor])) == 1:
                                organism, answerkey_organism, algorithm_organism = '1/1', result_recipient + ';' + result_donor, orgLine
                        if len(orgSet & set([result_recipient, result_donor])) == 0:
                                organism, answerkey_organism, algorithm_organism = '0/1', result_recipient + ';' + result_donor, orgLine
	
	endresults = [contig, contiglen, phyldist, result, answerkey_result, algorithm_result, organism, answerkey_organism, algorithm_organism, algorithmcall]
	myendresults.write('\t'.join(endresults) + '\n')
	listresults.append(result)
	list_answerkeyresults.append(answerkey_result)
	list_algorithmresults.append(algorithm_result)

myendresults.close()

for word in listresults:
	cnt[word] += 1
for word2 in list_answerkeyresults:
	cnt2[word2] += 1 
for word3 in list_algorithmresults:
	cnt3[word3] += 1

print 'The TP number is: ' + str(cnt['TP'])
print 'The TN number is: ' + str(cnt['TN'])
print 'The FP number is: ' + str(cnt['FP'])
print 'The FN number is: ' + str(cnt['FN']) 
print 'There were ' + str(len(listresults)) + ' contigs. There were ' + str(cnt2['hgt']) + ' with HGT, and ' + str(cnt2['nohgt']) + ' with no HGT.'
print 'The scripts called ' + str(cnt3['hgt']) + ' with HGT, and ' + str(cnt3['nohgt']) + ' with no HGT.' 

