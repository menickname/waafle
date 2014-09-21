#!/usr/bin/python

'''
'''

# Import
import argparse
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pickle

# Global constants


# Define arguments
parser = argparse.ArgumentParser()
parser.add_argument( '--drlist' )
parser.add_argument( '--ncontigs', type=int )
args = parser.parse_args()


# Utilities


# Determine minimum number of contigs that will be made
# This section only works occassionallyL I think this is because 'generatelist.py' sometimes fails.
"""
contignum = 0
while contignum < args.ncontigs:
	drtext = subprocess.check_output(["python", "generatelist.py"])
	drlist = drtext.strip().split('\n')
	for drpair in drlist:
		drinfo = drpair.strip().split('\t')
		taxadiff, donorGCF, recipGCF, donortaxa, reciptaxa = drinfo[0], drinfo[1], drinfo[2], drinfo[3], drinfo[4]
		contignum += 1
"""

# Temporary solution is to read in pre-generated files:
# Until then, the number of contigs will be determined by the number of lines in this file.
contignum = 0
fastalist = []
infolist = []

for astrline in open(args.drlist):
	aastrline = astrline.strip().split('\t')
	taxadiff, donorGCF, recipGCF, donortaxa, reciptaxa = aastrline[0], aastrline[1], aastrline[2], aastrline[3], aastrline[4]	
	subprocess.call(["python", "fakemake.py", "--recipient", ''+recipGCF+'', "--donor", ''+donorGCF+'', "--reciptaxa", ''+reciptaxa+'', "--donortaxa", ''+donortaxa+''])
		
	#Generate the fasta file
	seqrecord = pickle.load(open('fastaresult'))
	while seqrecord == None: #If makefake.py fails, we need to try again. This needs to be minimized to avoid infinite looping though...
		#print taxadiff, donorGCF, recipGCF
		subprocess.call(["python", "fakemake.py", "--recipient", ''+recipGCF+'', "--donor", ''+donorGCF+'', "--reciptaxa", ''+reciptaxa+'', "--donortaxa", ''+donortaxa+''])
		seqrecord = pickle.load(open('fastaresult'))
	else: #If we get a contig, we should number it and add to the info sheet.
		contignum += 1
		seqrecord.id = 'contig' + str(contignum)
		fastalist.append(seqrecord)
	
		#Generate the answer sheet
		contiglen = len(seqrecord.seq) 
		info = [seqrecord.id, contiglen, taxadiff, seqrecord.description, seqrecord.name, donortaxa, reciptaxa] #Only other thing would be to consider outputting # of genes in recip
		infolist.append(info)
	
	#Remove pickle
	subprocess.call(["rm", "fastaresult"])

SeqIO.write(fastalist, 'fakecontigs.fasta', 'fasta')

infosheet = open('infosheet.txt', 'w')
for list in infolist:
	infosheet.write('\t'.join(str(i) for i in list) + '\n')
	#infosheet.write('\t'.join(list))
infosheet.close()
	
