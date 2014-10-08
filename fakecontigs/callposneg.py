#!/usr/bin/python


'''
This script will classify the answerkey as what they should be.
'''

#Import
import argparse
import re

#Define arguments
parser = argparse.ArgumentParser()
parser.add_argument('--answerkey', help = 'Location and file of answerkey')
parser.add_argument('--taxa', help='Level to detect LGT at. All above this level is LGT, all below is NOT LGT.')
parser.add_argument('--genetable', help='Location and file for genetable')
parser.add_argument('--hgtresults', help='Location and file for LGT results')
args = parser.parse_args()

#Set taxalevels
leveldiff = ['t', 's', 'g', 'f', 'o', 'c', 'p', 'k']
levelnum = leveldiff.index(args.taxa)

#Read in answerkey
#The current goal is to see if the organisms match, and whether they are properly call LGT or not.
dictanswerkey = {}
for astrline in open(args.answerkey):
	aastrline = astrline.strip().split('\t')
	contigname, donorrecipinfo, donortaxa, reciptaxa = aastrline[0], aastrline[1], aastrline[2], aastrline[3]
	taxadiff, contiglen, numgenes  = aastrline[4], aastrline[5], aastrline[6]
	if leveldiff.index(taxadiff) >= leveldiff.index(args.taxa):
		dictanswerkey[contigname] = ['LGT', donortaxa, reciptaxa]
	elif leveldiff.index(taxadiff) < leveldiff.index(args.taxa):
		dictanswerkey[contigname] = ['noLGT', donortaxa, reciptaxa]

#This tells me how many P and N from the answerkey
truepositives = 0
truenegatives = 0
for contig in dictanswerkey.keys():
	if dictanswerkey[contig][0] == 'LGT':
		truepositives += 1
	else:
		truenegatives += 1

falsecontigs = open('FPFN.txt', 'w')
# Of the total number of positives, these are the positives from the pipeline
pipelineP = []
TPcount = 0
FPcount = 0
highconfLGT = open(args.hgtresults).readlines()
for i in range(len(highconfLGT)):
	contigname = highconfLGT[i].split(' ')[0]
	orgset = set(highconfLGT[i].split('[')[0].strip().split(' ')[4:])
	answer_status = dictanswerkey[contigname][0]
	answer_orgset = set([dictanswerkey[contigname][1], dictanswerkey[contigname][2]])
	answerabbr_orgset = set()
	for org in answer_orgset:
		searchstring = '.*' + args.taxa + '__\w*'
		neworg = re.search(searchstring, org).group()
		neworg2 = neworg.replace('|', '.')
		answerabbr_orgset.add(neworg2)
	numdiff = len(orgset - answerabbr_orgset)
	orgset_diff = orgset-answerabbr_orgset
	if answer_status == 'LGT' and numdiff == 0:
		TPcount += 1
		pipelineP.append(contigname)
	else:
		FPcount += 1
		falsecontigs.write(contigname + '\t' + 'FP' + '\t' + answer_status + '\t' + str(len(orgset)) + '\t' + str(numdiff) + '\t' + str(orgset_diff) + '\n')
		pipelineP.append(contigname)
pipelinepositives = len(pipelineP)

#Of all the negatives, these are the negatives from the pipeline
pipelineN = set()
for bstrline in open(args.genetable):
	bbstrline = bstrline.strip().split(' ')
	contigname, groupnum, taxa = bbstrline[0], bbstrline[1], bbstrline[2]
	finalscore, finalpercID, finalgroupcov, contigcov = bbstrline[3], bbstrline[4], bbstrline[5], bbstrline[6]
	start, end, combhitlen, grouplen, contiglen, status = bbstrline[7], bbstrline[8], bbstrline[9], bbstrline[10], bbstrline[11], bbstrline[12]
	pipelineN.add(contigname + ':' + status)


TNcount = 0
FNcount = 0
pipelinenegatives = 0
for contigstatus in pipelineN:
	print contigstatus
	contigname, status = contigstatus.split(':')[0], contigstatus.split(':')[1]
	if status == '1orgonly':
		answer_status = dictanswerkey[contigname][0]
		answer_orgset = set([dictanswerkey[contigname][1], dictanswerkey[contigname][2]])
		if answer_status == 'noLGT':
			TNcount += 1
			pipelinenegatives += 1
			#print contigname, 'TN-1orgonly' #Supposed to be noLGT, and algorithm gave 1orgonly
		else:
			FNcount += 1
			pipelinenegatives += 1
			falsecontigs.write(contigname + '\t' + 'FN' + '\t' + '1orgonly' + '\n')
			#print contigname, 'FN-1orgonly' #Supposed to be LGT, and algorithm gave 1orgonly
	elif status == '1+orghigh':
		answer_status = dictanswerkey[contigname][0]
		if answer_status == 'noLGT':
			TNcount += 1
			pipelinenegatives += 1
			#print contigname, 'TN-1+orghigh' #Supposed to be noLGT, and algorithm gave 1+orghigh
		else:
			FNcount += 1
			pipelinenegatives += 1
			falsecontigs.write(contigname + '\t' + 'FN' + '\t' + '1org+high' + '\n')
			#print contigname, 'FN-1+orghigh' #Supposed to be LGT, and algorithm gave 1+orghigh
	elif status == 'potentialLGT':
		answer_status = dictanswerkey[contigname][0]
		if contigname in pipelineP:
			pass
		else:
			if answer_status == 'noLGT':
				TNcount += 1
				pipelinenegatives += 1
				#print contigname, 'TN-potentiallgt'  #Supposed to be noLGT, and didn't make it through pipeline
			else:
				FNcount += 1
				pipelinenegatives += 1
				falsecontigs.write(contigname + '\t' + 'FN' + '\t' + 'potentialLGT' + '\n')
				#print contigname, 'FN-potentiallgt' #Supposed to be LGT, and didn't make it through the pipeline

falsecontigs.close()

#Output results
#print "There are " + str(truepositives) + " with LGT and " + str(truenegatives) + " without LGT in this fake contig set." 
#print "There are " + str(pipelinepositives) + " called LGT and " + str(pipelinenegatives) + " called noLGT from the pipeline."
#print "There are " + str(TPcount) + " true positives, " + str(TNcount) + " true negatives, " + str(FPcount) + " false positives, and " + str(FNcount) + " false negatives."
#print "The TP rate is " + str(TPcount/float(truepositives)) + ". The FP rate is " + str(FPcount/float(truenegatives)) + "."
print ' '.join([str(truepositives), str(truenegatives), str(pipelinepositives), str(pipelinenegatives), str(TPcount), str(TNcount), str(FPcount), str(FNcount), str(TPcount/float(truepositives)), str(FPcount/float(truenegatives))])
