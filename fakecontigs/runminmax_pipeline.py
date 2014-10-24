#!/bin/python

# import 
import subprocess
import argparse

# arguments to change
parser = argparse.ArgumentParser()
parser.add_argument( '--length', default=100, help='Length of groups to keep.')
parser.add_argument( '--taxa', help='Taxa level to detect LGT at.' )
parser.add_argument( '--hitoverlap', default=0.5, help='Amount hits should overlap to join a group.')
parser.add_argument( '--groupoverlap', default=0.5, help='Amount groups should overlap to merge.')
parser.add_argument( '--LGTthresh', help='Where to determine cutoff for LGT-ness score.')
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
groupedresults = open('groupedcombhits.out', 'w')
scriptloc = pipelineloc + '/blast2groups2.py'
blastloc = fakecontigloc + '/141020/blastn.sorted'
subprocess.call(["python", scriptloc, "--blastoutput", blastloc, "--hitoverlap", args.hitoverlap, "--groupoverlap", args.groupoverlap, "--length", args.length], stdout=groupedresults)
groupedresults.close()

# score all taxa within the BLAST results
scriptloc = pipelineloc + '/scoreorgs.py'
fastaloc = fakecontigloc + '/141020/fakecontigs4.fasta'
subprocess.call(["python", scriptloc, "--fasta", fastaloc, "--blastoutput", "groupedcombhits.out", "--taxa", args.taxa])
subprocess.call(["rm", "groupedcombhits.out"])

# detect LGT
hgtfilename = 'hgt_' + args.taxa + args.LGTthresh + '.txt'
hgtresults = open(hgtfilename, 'w')
scriptloc = pipelineloc + '/maxmin_method.py'
dict1 = "dddictCOGS_" + args.taxa + ".json"
dict2 = "dddictCGOS_" + args.taxa + ".json"
subprocess.call(["python", scriptloc, "--dict1", dict1, "--dict2", dict2, "--LGTcut", args.LGTthresh], stdout=hgtresults)
hgtresults.close()

# determine FP, FN, TP, TN results
scriptloc = fakecontigloc + '/callposneg_maxmin.py'
answerkey = fakecontigloc + '/141020/answerkey4.txt'
info = subprocess.check_output(["python", scriptloc, "--answerkey", answerkey, "--taxa", args.taxa, "--hgtresults", hgtfilename])
subprocess.call(["rm", dict1])
print ' '.join([args.taxa, args.LGTthresh, args.hitoverlap, args.groupoverlap]) + ' ' + info

