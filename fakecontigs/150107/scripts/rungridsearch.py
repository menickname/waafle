#!/bin/python

'''
This will run the pipeline to generate a grid search.
'''

# import
import subprocess
import argparse
import numpy

# arguments
parser = argparse.ArgumentParser()
parser.add_argument( '--fasta', help='Fasta file of contigs.' )
args = parser.parse_args()

# fixed values
pipelineloc = '/n/home05/thsu/bitbucket/hgt_project/pipeline'
currentloc = subprocess.check_output(["pwd"])
fakecontigloc = '/n/home05/thsu/bitbucket/hgt_project/fakecontigs/150107'
taxa = ['p', 'c', 'o', 'f', 'g', 's']
#taxa = ['g','s']
complement_range = numpy.arange(0.8, 1.0, 0.05)
onebug_range = numpy.arange(0.8, 1.0, 0.05)

myresults = open('gridsearch.txt', 'w')

# read in the parameters
for level in taxa:
	# generate scores at that level
	scriptloc = pipelineloc + '/scoreorgs.py'
	groupedcombhit_loc = fakecontigloc + '/groupedcombhits.out'
	#subprocess.call(["python", scriptloc, "--fasta", args.fasta, "--blastoutput", groupedcombhit_loc, "--taxa", level])
	
	for i in range(len(onebug_range)):
		onebug = onebug_range[i]
		for j in range(len(complement_range)):
			complement = complement_range[j]
			print level, onebug, complement
			skip = False
			if onebug > min(complement_range):
				for obscore in onebug_range[0:i]:
					if numpy.isclose(obscore, complement):
						skip = True
						break
			if skip == False:
				hgtfilename = 'hgt_results_' + level + '.txt'
				hgtresults = open(hgtfilename, 'w')
				scriptloc = pipelineloc + '/method3/method3-1.py'
				dictnamelgt_1 = 'dddictCGOS_' + level + '.json'				
				dictnamelgt_2 = 'dddictCOGS_' + level + '.json'
				subprocess.call(["python", scriptloc, "--dddictCGOS", dictnamelgt_1, "--dddictCOGS", dictnamelgt_2, "--onebug", str(onebug), "--complement", str(complement), "--taxa", level], stdout=hgtresults)
				hgtresults.close()
				#print 'python', dictnamelgt_1, dictnamelgt_2, onebug, complement, level
				#determine TP/TN/FP/FN rates
				scriptloc = fakecontigloc + '/scripts/callposneg_DR.py'
				answerkeyloc = fakecontigloc + '/answerkey.txt'
				partialresults = subprocess.check_output(["python", scriptloc, "--taxa", level, "--hgtresults", hgtfilename, "--answerkey", answerkeyloc])
				myresults.write(' '.join([level, str(onebug), str(complement), partialresults.strip()]) + '\n')
				#print 'python', 'checkresults', answerkeyloc
myresults.close()

