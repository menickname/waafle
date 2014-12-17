#!/bin/python

"""
This script will take us through our first method for detecting LGT. To do this, it will:
1) BLAST contigs against Repophlan 2013. 
"""

# import 
import subprocess
import argparse

# arguments to change
parser = argparse.ArgumentParser()
parser.add_argument( '--fasta', help='Location and file of assembled contigs. This should be in fasta format.' )
parser.add_argument( '--hitoverlap', default=0.5, help='Amount hits should overlap to join a group.')
parser.add_argument( '--groupoverlap', default=0.5, help='Amount groups should overlap to merge.')
parser.add_argument( '--delta', type=str, help='Upper threshold for calling high confidence BLAST hits' )
parser.add_argument( '--epsilon', type=str, help='Lower threshold for calling low confidence BLAST hits' )
parser.add_argument( '--length', default=100, help='Length of groups to exclude.' )
parser.add_argument( '--taxa', type=str, help='Taxa level to detect LGT at.' )
args = parser.parse_args()

# file locations
pipelineloc = '/n/home05/thsu/bitbucket/hgt_project/pipeline'
currentloc = subprocess.check_output(["pwd"])

# run BLAST
#database = '/n/huttenhower_lab_nobackup/data/hgt/blast/blast_db_updated/repophlan_31122013_speciescentroids.db'
#input = args.fasta
output = 'blastn_1000.out'
#subprocess.call(["blast", "-db", database, "-query", input, "-out", output, "-outfmt", "'6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp'"])

# sort BLAST hits by name, length, and bitscore
blastoutput_sorted = open('blastn_1000.sorted', 'w')
subprocess.call(["sort", "-k1,1", "-k4,4nr", "-k12,12nr", output], stdout=blastoutput_sorted)
blastoutput_sorted.close()

# identify gene groups from the BLAST hits based on group length and overlap
groupedresults = open('groupedcombhits.out', 'w')
scriptloc = pipelineloc + '/blast2groups2.py'
subprocess.call(["python", scriptloc, "--blastoutput", "blastn_1000.sorted", "--hitoverlap", args.hitoverlap, "--groupoverlap", args.groupoverlap, "--length", "100"], stdout=groupedresults)
groupedresults.close()

# score all taxa within the BLAST results
scriptloc = pipelineloc + '/scoreorgs.py'
subprocess.call(["python", scriptloc, "--fasta", args.fasta, "--blastoutput", "groupedcombhits.out", "--taxa", args.taxa])
#subprocess.call(["rm", "groupedcombhits.out"])

# generate the gene table
scriptloc = pipelineloc + '/converttogenetable.py'
dictname1 = 'dddictCGOS_' + args.taxa + '.json'
dictname2 = 'dddictCOGS_' + args.taxa + '.json'
genetable = open('genetable.txt', 'w')
subprocess.call(["python", scriptloc, "--dddictCGOS", dictname1], stdout=genetable)
genetable.close()

### Up until this section the pipeline is the same for all methods
# split contigs into oneorgonly, highconforg, and potentialLGT
scriptloc = pipelineloc + '/method1/detectcontigtype.py'
subprocess.call(["python", scriptloc, "--dddictCGOS", dictname1, "--dddictCOGS", dictname2, "--delta", args.delta])

# detect high confidence lgt
hgtfilename = 'hgt_results_' + args.taxa + '.txt'
hgtresults = open(hgtfilename, 'w')
scriptloc = pipelineloc + '/method1/detectlgt.py'
dictnamelgt_1 = 'dddictCGOS_lgt_' + args.taxa + '.json'
dictnamelgt_2 = 'dddictCOGS_lgt_' + args.taxa + '.json'
subprocess.call(["python", scriptloc, "--dict2", dictnamelgt_1, "--dict1", dictnamelgt_2, "--delta", args.delta, "--epsilon", args.epsilon], stdout=hgtresults)
hgtresults.close()

"""
# determine FP, FN, TP, TN results
scriptloc = fakecontigloc + '/callposneg.py'
answerkey = fakecontigloc + '/141006/gridsearch/answerkey4.txt'
genetable = fakecontigloc + '/141006/gridsearch/genetable' + args.taxa + '.txt'
info = subprocess.check_output(["python", scriptloc, "--answerkey", answerkey, "--taxa", args.taxa, "--genetable", genetable, "--hgtresults", hgtfilename])
print ' '.join([args.taxa, args.delta, args.epsilon, args.hitoverlap, args.groupoverlap]) + info
"""
