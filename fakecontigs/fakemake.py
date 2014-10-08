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

# testing
random.seed( 1 )

# any global constants here
mingenelength=100
mincontiglen=300
mediancontiglen=905

# Define arguments
parser = argparse.ArgumentParser()
parser.add_argument( '--recipient', help='This is the recipient GCF number.' )
parser.add_argument( '--donor', help='This is the donor GCF number.' )
parser.add_argument( '--reciptaxa', help='This is the complete recipient taxonomy.' )
parser.add_argument( '--donortaxa', help='This is the complete donor taxonomy.' )
parser.add_argument( '--taxadiff', help='This is the phylogenetic difference between the 2.')
parser.add_argument( '--ngenes', type=int, default=1, help='This is the number of contigs we would like to generate for this donor-recipient pair.')
#I have not implemented this, but it could be a solution for generating the same donor/recip with multiple genes to see if the donor/recip pair is the problem, or not
args = parser.parse_args()

outname= "-".join( [args.recipient, args.donor] )

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
	# gi like ">bug|:1-1213,1213-154242 annotations" 
	# first string annotation
	# then isolate coordinates
	# then, in rare case of ,-separate coords, take the first pair
	coordinates = GI.split( " " )[0].split(':')[1].split( "," )[0]
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
		# These should be calculated, since we will always have rgenes on the left side
		if start < contigstart and end > contigstart:
			if end <= rstart:
				genelen = end - contigstart + 1
                                if genelen > mingenelength:
                                        surroundgenes.append(["R", genelen, genename, contigstart, end, 0, genelen-1])
			elif end > rstart:
				genelen = rstart - contigstart + 1
                                if genelen > mingenelength:
                                        surroundgenes.append(["R", genelen, genename, contigstart, rstart, 0, genelen-1])
		elif start >= contigstart and start <= rstart:
			if end <= rstart:
				genelen = end - start + 1
				if genelen > mingenelength:
					surroundgenes.append(["R", genelen, genename, start, end, start-contigstart, start-contigstart+genelen-1])
			elif end > rstart:
				genelen = rstart - start + 1
				if genelen > mingenelength:
					surroundgenes.append(["R", genelen, genename, start, rstart, start-contigstart, start-contigstart+genelen-1])
		#These should be on the other side of the donor gene
		leftover = contiglen/2 - len(dbog[dgenename].seq)
		rcontigend = leftover + rend
		if leftover > 0:
			if start >= rend and start <= rcontigend:
				if end <= rcontigend:
					genelen = end - start + 1
					if genelen > mingenelength:
						surroundgenes.append(["Rafter", genelen, genename, start, end, rstart-contigstart+len(dbog[dgenename].seq)+start-rend, rstart-contigstart+len(dbog[dgenename].seq)+start-rend+genelen])
				elif end > rcontigend:
					genelen = rcontigend - start + 1
					if genelen > mingenelength:
						surroundgenes.append(["Rafter", genelen, genename, start, rcontigend, rstart-contigstart+len(dbog[dgenename].seq)+start-rend, rstart-contigstart+len(dbog[dgenename].seq)+start-rend+genelen])
		if leftover < 0 and start > rstart:
			break
		elif leftover > 0 and start > rend + leftover:
			break
	
	#determine if there are any problems with detecting hgt in this contig
	if len(surroundgenes) == 0 or len(dbog[dgenename].seq[0:contiglen/2]) < mingenelength:
		return None
	else:
		leftover = contiglen/2 - len(dbog[dgenename].seq)
		if leftover < 0:
			surroundgenes.append(["D", len(dbog[dgenename].seq[0:contiglen/2]), dgenename, 0, contiglen/2, rstart-contigstart, rstart-contigstart+contiglen/2])
		else:
			surroundgenes.append(["D", len(dbog[dgenename].seq), dgenename, 0, len(dbog[dgenename].seq), rstart-contigstart, rstart-contigstart+len(dbog[dgenename].seq)])
		newcontiglen = len(newseq[contigstart:contigend])
		# Generate the new seqRecord
		new_SeqRec = SeqRecord(newseq[contigstart:contigend])
		new_SeqRec.features = surroundgenes
        	new_SeqRec.description = 'donor:' + args.donor + '|' + dgenename + '||recipient:' + args.recipient + '|' + rgenename
        	new_SeqRec.annotations['source'] = args.reciptaxa
        	new_SeqRec.annotations['taxonomy'] = args.donortaxa
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

fh_fnt = open( outname+".fnt", "w" )
fh_ans = open( outname+"-ans.txt", "w" )

escape_counter = 0
successes=0
while successes < args.ngenes:
	escape_counter += 1
	if escape_counter >= 25 * args.ngenes:
		sys.exit( "got stuck in the loop!" )
	result = make1contig(dbog, rbog, rcoords)	
	if result is not None:
		successes += 1
		# write the sequence
		contig_name = "%s|contig%05d|" % ( outname, successes )
		print >>fh_fnt, ">"+contig_name 
		print >>fh_fnt, str( result.seq )
		# write the answer 
		
		print >>fh_ans, "\t".join( 
			[contig_name, 
			 result.description, 
			 args.reciptaxa, 
			 args.donortaxa,
			 args.taxadiff,
			 str(len(result.seq)),
			 str(len(result.features)),
			 str(result.features),
			 ]
			  )

# cleanup
fh_fnt.close()
fh_ans.close()
