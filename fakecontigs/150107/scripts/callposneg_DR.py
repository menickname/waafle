#!/usr/bin/python


'''
01/07/2015:
This script will classify the answerkey as what they should be and match the results from method3-1.
'''

#Import
import argparse
import re

#Define arguments
parser = argparse.ArgumentParser()
parser.add_argument('--answerkey', help = 'Location and file of answerkey')
parser.add_argument('--taxa', help='Level to detect LGT at. All above this level is LGT, all below is NOT LGT.')
parser.add_argument('--hgtresults', help='Location and file for all results')
args = parser.parse_args()

#Define functions
def standardOrgs( orgset, taxalevel ):
	searchstring = taxalevel + '__\w*'
	newset = set()
	for org in orgset:
		neworg = re.search(searchstring, org).group()
		newset.add(neworg)
	return newset

#Check whether TPs and TNs got it correct or not
def checkOrgs ( resultset, answerset, newstatus ):
	comparisonset = resultset & answerset
	if newstatus == 'TP':
		if len(comparisonset) == 2:
			pass #perfect match
			orgmatch = 'full'
		elif len(comparisonset) == 1:
			pass #1/2 match
			orgmatch = 'partial'
		else:
			pass #no match
			orgmatch = 'none'
	elif newstatus == 'TN':
		if len(comparisonset) >= 1:
			pass 
			orgmatch = 'full'
		else: 
			pass
			orgmatch = 'none'
	return orgmatch

def callDR ( orglist ):
	statusset = set(orglist[2:4])
	if len(statusset) == 1:
		pass# Could be donor-donor, recip-recip, or unclear-unclear
		if 'unclear' in statusset:
			pass
			status = 'UU'
			org1, org2 = orglist[0], orglist[1]
		elif 'donor' in statusset:
			pass #print 'donoronly'
			status = 'DD'
			org1, org2 = orglist[0], orglist[1]
		else:
			pass #print 'reciponly'
			status = 'RR' 
			org1, org2 = orglist[0], orglist[1]
	else:
		pass # Could be donor-recip, unclear-donor, unclear-recip
		if 'unclear' in statusset and 'donor' in statusset:
			status = 'DU'
			if orglist[2] == 'donor':
				org1 = orglist[0]
				org2 = orglist[1]
			elif orglist[3] == 'donor':
				org1 = orglist[1]
				org2 = orglist[0]
		elif 'unclear' in statusset and 'recipient' in statusset:
			status = 'RU'
			if orglist[2] == 'recipient':
				org1 = orglist[0]
				org2 = orglist[1]
			elif orglist[3] == 'recipient':
				org1 = orglist[1]
				org2 = orglist[0]
		else:
			status = 'DR'
			if orglist[2] == 'donor':
				org1  = orglist[0]
				org2 = orglist[1]
			else:
				org1 = orglist[1]
				org2 = orglist[0]
	return status, org1, org2 
		

def checkDR ( DRstatus, orgstatus, w_donor, w_recip, ans_donor, ans_recip ):
	if DRstatus == 'DR':
		if w_donor == ans_donor and w_recip == ans_recip:
			return 'DR_correct'
		else:
			if orgstatus == 'full':
				return 'DR_incorrect'
			elif orgstatus == 'partial':
				return 'DR_incorrectorgp'
			else:
				return 'DR_incorrectorgn'
	else:
		if orgstatus == 'full':
			return 'DR_uncalledf'
		elif orgstatus == 'partial':
			return 'DR_uncalledp'
		else:
			return 'DR_uncalledn'

#Set taxalevels
leveldiff = ['t', 's', 'g', 'f', 'o', 'c', 'p', 'k']
levelnum = leveldiff.index(args.taxa)

#Read in answerkey
#The current goal is to see if the organisms match, and whether they are properly call LGT or not.
dictanswerkey = {}
for astrline in open(args.answerkey):
	aastrline = astrline.strip().split('\t')
	contigname, donorrecipinfo, reciptaxa, donortaxa = aastrline[0], aastrline[1], aastrline[2], aastrline[3]
	taxadiff, contiglen, numgenes  = aastrline[4], aastrline[5], aastrline[6]
	if leveldiff.index(taxadiff) >= leveldiff.index(args.taxa):
		dictanswerkey[contigname] = ['LGT', reciptaxa, donortaxa]
	elif leveldiff.index(taxadiff) < leveldiff.index(args.taxa):
		dictanswerkey[contigname] = ['noLGT', reciptaxa, donortaxa]

#This tells me how many P and N from the answerkey
truepositives = 0
truenegatives = 0
for contig in dictanswerkey.keys():
	if dictanswerkey[contig][0] == 'LGT':
		truepositives += 1
	else:
		truenegatives += 1

# Of the total number of TP, FP, TN, FN as called by the pipeline
pipelineP = []

#For TP/TN/FP/FN
TPcount = 0
FPcount = 0
TNcount = 0
FNcount = 0

#Keep count for whether orgs were correctly called
LGT_full = 0
LGT_partial = 0 
LGT_none = 0
oneorg_full = 0
oneorg_none = 0

#Keep count for whether DR were correctly called
LGT_drcorrect = 0
LGT_drincorrect = 0
LGT_drincorrectop = 0
LGT_drincorrecton = 0
LGT_druncalledf = 0
LGT_druncalledp = 0
LGT_druncalledn = 0

for astrline in open(args.hgtresults):
	aastrline = astrline.strip().split(' ')
	contig, label, score1, score2, numbugs, numgenes = aastrline[0], aastrline[1], aastrline[2], aastrline[3], aastrline[4], aastrline[5]
	
	# Define organism sets and DR status for WAAFLE Results
	# If there are 2 orgs involved:
	if label == 'LGT' or label == 'ambiguous-LGT':
		bug1name, bug2name = aastrline[6], aastrline[7]
		bug1status, bug2status = aastrline[8], aastrline[9]

		#Define sets (to see if orgs are correct)
		resultorgset = set([bug1name, bug2name])
		resultstatus = aastrline[6:10]

		#Define DR (to see if it called DR correctly)
		w_status, w_donor, w_recipient = callDR(resultstatus)

	# If there is only 1 org involved:
	else:
		bug1name = aastrline[6]
		resultorgset = set([bug1name])
		w_status, w_org = 'oneorg', bug1name
	
	# Define organism sets and DR status for answerkey
	# Define sets (for answerkey) to see if orgs are correct (this can be done at once for LGT and non-LGT candidates)	
	answerorgs = set(dictanswerkey[contig][1:])	
	answerorgset = standardOrgs(answerorgs, args.taxa)
	ans_label = dictanswerkey[contig][0]
	searchstring = args.taxa + '__\w*'	
	# Define DR (for answerkey)
	if ans_label == 'LGT':
		answerstatus = [re.search(searchstring, dictanswerkey[contig][1]).group(), re.search(searchstring, dictanswerkey[contig][2]).group(), 'recipient', 'donor']
		ans_status, ans_donor, ans_recipient = callDR(answerstatus)
	else:
		ans_status, ans_org = 'oneorg', re.search(searchstring, dictanswerkey[contig][1]).group()

	# Call TP and FP
	if label == 'LGT':
		
		if dictanswerkey[contig][0] == 'LGT':
			newstatus = 'TP'
			TPcount += 1
			orgmatch = checkOrgs(resultorgset, answerorgset, newstatus)
			if orgmatch == 'full':
				LGT_full += 1
			elif orgmatch == 'partial':
				LGT_partial += 1
			else:
				LGT_none += 1
			drmatch = checkDR(w_status, orgmatch, w_donor, w_recipient, ans_donor, ans_recipient)
			if drmatch == 'DR_correct':
				LGT_drcorrect += 1
			elif drmatch == 'DR_incorrect':
				LGT_drincorrect += 1
			elif drmatch == 'DR_incorrectorgp':
				LGT_drincorrectop += 1
			elif drmatch == 'DR_incorrectorgn':
				LGT_drincorrecton += 1
			elif drmatch == 'DR_uncalledf':
				LGT_druncalledf += 1
			elif drmatch == 'DR_uncalledp':
				LGT_druncalledp += 1
			else:
				LGT_druncalledn += 1			
			print contig, 'TP', numbugs, numgenes, orgmatch, drmatch, w_status, label, 'waafle:', w_donor, w_recipient, ans_status, 'answer:', ans_donor, ans_recipient
		else:
			newstatus = 'FP'
			FPcount += 1
			print contig, 'FP', numbugs, numgenes, 'NA', 'NA', 'NA', label, 'waafle:', w_donor, w_recipient, ans_status, 'answer:', ans_org
	else:
		if dictanswerkey[contig][0] == 'LGT':
			newstatus = 'FN'
			FNcount += 1
			if label == 'ambiguous-LGT':
				print contig, 'FN', numbugs, numgenes, 'NA', 'NA', 'NA', label, 'waafle:', w_donor, w_recipient, ans_status, 'answer:', ans_donor, ans_recipient  
			else:
				print contig, 'FN', numbugs, numgenes, 'NA', 'NA', 'NA', label, 'waafle:', w_org, ans_status, 'answer:', ans_donor, ans_recipient
		else:
			newstatus = 'TN'
			TNcount += 1
			orgmatch = checkOrgs(resultorgset, answerorgset, newstatus)		
			if orgmatch == 'full':
				oneorg_full += 1
			else:
				oneorg_none += 1
			print contig, 'TN', numbugs, numgenes, orgmatch, 'NA', w_status, label, 'waafle:', w_org, ans_status, 'answer:', ans_org
				
	#print contig, newstatus, dictanswerkey[contig][0], label, answerorgset, resultstatus, orgmatch, check
	#Use this to determine whether organisms are correct

TPR = TPcount/float(truepositives)
FPR = FPcount/float(truenegatives)

#Orgs correct
LGTfull = LGT_full/float(TPcount)
LGTpartial = LGT_partial/float(TPcount)
LGTnone = LGT_none/float(TPcount)
oneorgfull = oneorg_full/float(TNcount)
oneorgnone = oneorg_none/float(TNcount)

#DR correct
LGTdr = LGT_drcorrect/float(TPcount)
LGTdrincorrect = LGT_drincorrect/float(TPcount)
LGT_drincorrecto = LGT_drincorrectop + LGT_drincorrecton
LGTdrio = LGT_drincorrecto/float(TPcount)
LGTdriop = LGT_drincorrectop/float(TPcount)
LGTdrion = LGT_drincorrecton/float(TPcount)
LGT_druncalled = LGT_druncalledf + LGT_druncalledp + LGT_druncalledn
LGTu = LGT_druncalled/float(TPcount)
LGTuf = LGT_druncalledf/float(TPcount)
LGTup = LGT_druncalledp/float(TPcount)
LGTun = LGT_druncalledn/float(TPcount)


print truepositives, truenegatives, TPcount, FPcount, TNcount, FNcount, LGT_full, LGT_partial, LGT_none, oneorg_full, oneorg_none, LGT_drcorrect, LGT_drincorrect, LGT_drincorrecto, LGT_drincorrectop, LGT_drincorrecton, LGT_druncalled, LGT_druncalledf, LGT_druncalledp, LGT_druncalledn, TPR, FPR, LGTfull, LGTpartial, LGTnone, oneorgfull, oneorgnone, LGTdr, LGTdrincorrect, LGTdrio, LGTdriop, LGTdrion, LGTu, LGTuf, LGTup, LGTun 
