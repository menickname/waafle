#!/bin/python

# import 
import subprocess
import argparse

# arguments to change
parser = argparse.ArgumentParser()
parser.add_argument( '--delta', type=str, help='Upper threshold for calling high confidence BLAST hits' )
parser.add_argument( '--epsilon', type=str, help='Lower threshold for calling low confidence BLAST hits' )
parser.add_argument( '--length', default=100, help='Length of groups to exclude.' )
parser.add_argument( '--taxa', help='Taxa level to detect LGT at.' )
parser.add_argument( '--hitoverlap', default=0.5, help='Amount hits should overlap to join a group.')
parser.add_argument( '--groupoverlap', default=0.5, help='Amount groups should overlap to merge.')
args = parser.parse_args()

# file locations
pipelineloc = '/n/home05/thsu/bitbucket/hgt_project/pipeline'
fakecontigloc = '/n/home05/thsu/bitbucket/hgt_project/fakecontigs' 
currentloc = subprocess.check_output(["pwd"])

# run BLAST
# ignore this for now

# sort BLAST hits by name, length, and bitscore
# ignore this for now

# identify gene groups from the BLAST hits based on group length and overlap
groupedresults = open('groupedcombhits2.out', 'w')
scriptloc = pipelineloc + '/blast2groups2.py'
blastloc = fakecontigloc + '/141006/gridsearch/blastn.sorted'
subprocess.call(["python", scriptloc, "--blastoutput", blastloc, "--hitoverlap", args.hitoverlap, "--groupoverlap", args.groupoverlap, "--length", "100"], stdout=groupedresults)
groupedresults.close()

# score all taxa within the BLAST results
scriptloc = pipelineloc + '/scoreorgs.py'
fastaloc = fakecontigloc + '/141006/gridsearch/fakecontigs4.fasta'
subprocess.call(["python", scriptloc, "--fasta", fastaloc, "--blastoutput", "groupedcombhits2.out", "--taxa", args.taxa, "--delta", args.delta])
subprocess.call(["rm", "groupedcombhits2.out"])

# detect LGT
hgtfilename = 'hgt_' + args.taxa + '2.txt'
hgtresults = open(hgtfilename, 'w')
scriptloc = pipelineloc + '/detectlgt.py'
dict1 = "dddictCOGS_" + args.taxa + ".json"
dict2 = "dddictCGOS_" + args.taxa + ".json"
subprocess.call(["python", scriptloc, "--dict1", dict1, "--dict2", dict2, "--delta", args.delta, "--epsilon", args.epsilon], stdout=hgtresults)
hgtresults.close()

# determine FP, FN, TP, TN results
scriptloc = fakecontigloc + '/callposneg.py'
answerkey = fakecontigloc + '/141006/gridsearch/answerkey4.txt'
genetable = fakecontigloc + '/141006/gridsearch/genetable' + args.taxa + '.txt'
info = subprocess.check_output(["python", scriptloc, "--answerkey", answerkey, "--taxa", args.taxa, "--genetable", genetable, "--hgtresults", hgtfilename])
print ' '.join([args.taxa, args.delta, args.epsilon, args.hitoverlap, args.groupoverlap]) + info

