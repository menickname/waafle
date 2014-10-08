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
#parser.add_argument( '--params', help='Upper threshold for calling high confidence BLAST hits' )
args = parser.parse_args()

# fixed values
scriptloc = '/n/home05/thsu/bitbucket/hgt_project/fakecontigs/runpipeline.py'
taxa = ['p', 'c', 'o', 'f', 'g', 's']
zerotoone = numpy.arange(0, 1, 0.1)
#epsilon = 0.25
#delta = 0.75

myresults = open('gridsearch.txt', 'w')
resultlist = []
# read in the parameters
for level in taxa:
	for delta in zerotoone:
		for epsilon in zerotoone:
			if epsilon <= delta:
				result = subprocess.check_output(['python', scriptloc, '--taxa', str(level), '--delta', str(delta), '--epsilon', str(epsilon), '--hitoverlap', '0.5', '--groupoverlap', '0.5'])
				results = result.strip()
				resultlist.append(results)

for line in resultlist:
	myresults.write(line + '\n')
myresults.close()
