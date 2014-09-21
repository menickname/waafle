#!/usr/bin/python

'''
This script will take a recipient GCF, donor GCF, and create ncontigs for them.
'''

# Import
import argparse
import random
import re
import math
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pickle

# any global constants here
mingenelength=100
mincontiglen=300
mediancontiglen=905
seed=1

# Define arguments
parser = argparse.ArgumentParser()
parser.add_argument( '--recipient' )
parser.add_argument( '--donor' )
parser.add_argument( '--reciptaxa' )
parser.add_argument( '--donortaxa' )
parser.add_argument( '--ngenes' ) #I have not implemented this, but it could be a solution for generating the same donor/recip with multiple genes to see if the donor/recip pair is the problem, or not
args = parser.parse_args()

# utilities
def load_scaffolds ( GCF ):
	""" returns dict with key=scaffold GI, item = SeqRecord """
	fileloc = '/n/huttenhower_lab/data/2014_repophlan/repophlan_microbes/fna/'
	genomefileloc = fileloc  + GCF + '.fna'
        dictScaffolds = SeqIO.index(genomefileloc, 'fasta')
        return dictScaffolds 

def load_bog ( GCF ):
	""" returns dict with key=gene GI, item=SeqRecord """
	fileloc = '/n/huttenhower_lab/data/2014_repophlan/repophlan_microbes/ffn/'
	bogfileloc = fileloc + GCF + '.ffn'
	dictBog = SeqIO.index(bogfileloc, 'fasta')
	return dictBog

def getcoords ( GI ):
	""" return coordinates in sorted order """
	coordinates = GI.split(':')[1]
	if re.search('c', coordinates):
		end, start = int(coordinates.split('-')[0][1:]), int(coordinates.split('-')[1])
	else:
		start, end = int(coordinates.split('-')[0]), int(coordinates.split('-')[1])
	return start, end

def bog2coords ( dictBog ):
	""" return dict with key=scaffold name, value=[ start, end, geneGI ] """
	coords = {}
	for geneGI in dictBog:
		start, end = getcoords(geneGI)
		scaffoldname = geneGI.split(':')[0]
		coords.setdefault(scaffoldname, []).append([start, end, geneGI])
	return coords

# actually load the scaffolds and the genes
scaffolds = load_scaffolds( args.recipient )
rbog = load_bog( args.recipient )
dbog = load_bog( args.donor )

# get coords
rcoords = bog2coords( rbog )

# make a random contig
def make1contig( dbog, rbog, rcoords):
	
	#random.seed(seed)
	dgenename = random.choice(dbog.keys()) #choose a donor gene
	rgenename = random.choice(rbog.keys()) #choose a recipient gene
	dGI, rGI = dgenename.split(':')[0], rgenename.split(':')[0]
	rstart, rend = getcoords(rgenename) #get coords of recipient gene

	#random.seed()
	contiglen = max(int(random.expovariate(1/float(mediancontiglen))), mincontiglen) #determine contig length
	contigstart, contigend = rstart-contiglen/2, rstart+contiglen/2 #get start and end sites for the hybrid sequence
	newseq = scaffolds[rGI].seq[0: rstart - 1] + dbog[dgenename].seq + scaffolds[rGI].seq[rend:]
	if contigstart < 0:
		contigstart = 0
	if contigend > len(newseq):
		contigend = len(newseq)-1

	#determine if recipient genes will be in contig
	surroundgenes = []
	for start, end, genename in sorted(rcoords[rGI]):
		if start < contigstart and end > contigstart:
			if end <= rstart:
				genelen = end - contigstart + 1
                                if genelen > mingenelength:
                                        surroundgenes.append([genelen, genename, contigstart, end])
			elif end > rstart:
				genelen = rstart - contigstart + 1
                                if genelen > mingenelength:
                                        surroundgenes.append([genelen, genename, contigstart, rstart])
		elif start >= contigstart and start <= rstart:
			if end <= rstart:
				genelen = end - start + 1
				if genelen > mingenelength:
					surroundgenes.append([genelen, genename, start, end])
			elif end > rstart:
				genelen = rstart - start + 1
				if genelen > mingenelength:
					surroundgenes.append([genelen, genename, start, rstart])
		elif start > contigstart:
			break
	
	#determine if there are any problems with detecting hgt in this contig
	if len(surroundgenes) == 0 or len(dbog[dgenename].seq) < mingenelength:
		return None
	else:
		newcontiglen = len(newseq[contigstart:contigend])
		# Generate the new seqRecord
		new_SeqRec = SeqRecord(newseq[contigstart:contigend])
        	new_SeqRec.description = 'donor:' + args.donor + '|' + dgenename + '||recipient:' + args.recipient + '|' + rgenename
        	new_SeqRec.annotations['source'] = args.reciptaxa
        	new_SeqRec.annotations['taxonomy'] = args.donortaxa
		new_SeqRec.name = surroundgenes
		return new_SeqRec	
	 
	"""
	* determine its scaff and coords from the rcoords dict
	* pick random dgene (could be from rec. if making a TN)
	* do the swap
	* tests...
	** is the hgt'ed piece too short?
	** do we not have SOME recipient gene to the left?
	"""

# do that N times
"""
This does not work because each time, it creates only one contig. There are several issues with this:
1. The contig counter needs to be outside of this script, and it doesn't seem to number all the contigs correctly
2. Writing the the fasta file overwrites each previous entry. Biopython seems to only allow for writing a whole list, or writing a one sequence.

9/21/14
Current solution is to pickle the SeqRecord and pass it back to the wrapper script.
"""
pickleresults = open('fastaresult', 'w')
result = make1contig(dbog, rbog, rcoords)	
if result is not None:
	pickle.dump(result, pickleresults, 2)
else:
	#print 'noresults'
	pickle.dump(None, pickleresults, 2)
