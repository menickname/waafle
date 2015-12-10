#!/usr/bin/python

"""
WAAFLE VALIDATION: fakemake.py

Authors:
Tiffany Hsu
Eric Franzosa

This script generates the synthetic contigs from a list of GCFs.
Taxa pairs are drawn from /n/huttenhower_lab/data/2014_repophlan/repophlan_microbes.txt
Run with the "-h" flag for usage help.
"""

from __future__ import print_function # python 2.7+ required
import argparse, random, re, math, pickle
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# ---------------------------------------------------------------
# constants
# ---------------------------------------------------------------
# testing
#random.seed( 1 )
mingenelength=100
mincontiglen=300
mediancontiglen=905


# ---------------------------------------------------------------
# functions
# ---------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument( '--recipient', help='This is the recipient GCF number.' )
parser.add_argument( '--donor', help='This is the donor GCF number.' )
parser.add_argument( '--reciptaxa', help='This is the complete recipient taxonomy.' )
parser.add_argument( '--donortaxa', help='This is the complete donor taxonomy.' )
parser.add_argument( '--taxadiff', help='This is the phylogenetic difference between the 2.')
parser.add_argument( '--ngenes', type=int, default=1, help='This is the number of contigs we would like to generate for this donor-recipient pair.')
#I have not implemented this, but it could be a solution for generating the same donor/recip with multiple genes to see if the donor/recip pair is the problem, or not
args = parser.parse_args()

# ---------------------------------------------------------------
# functions
# ---------------------------------------------------------------
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

def calc_overlap( onestart, oneend, twostart, twoend ):
    """ Calculate overlap between two hits or genes. """
    if onestart - twoend > 0 or twostart - oneend > 0:
        return 0
    else:
        if onestart < twostart and twoend < oneend or onestart > twostart and oneend < twoend:
            divisor = float( min( oneend - onestart, twoend - twostart ) )
	else:
            divisor = float( min( oneend - onestart + 1, twoend - twostart + 1 ) )
        coord_sorted = sorted( [onestart, oneend, twostart, twoend] )
        overlap = float( ( coord_sorted[2] - coord_sorted[1] )/divisor )
        return overlap

def make_contig( donor, recipient, dbog, rbog, rcoords, scaffolds ):
    adjcontiglen, dgenelen = 0, 1
    dgenename, rgenename, dGI, rGI, newseq = "", "", "", "", ""
    start, rend, contigstart, contigend = 0, 0, 0, 0
    while adjcontiglen < dgenelen:
        contiglen = max( int( random.expovariate( 1/float( mediancontiglen ) ) ), mincontiglen ) #determine contig length
        dgenename = random.choice(dbog.keys()) #choose a donor gene
        rgenename = random.choice(rbog.keys()) #choose a recipient gene
        dGI, rGI = dgenename.split(':')[0], rgenename.split(':')[0]
        rstart, rend = getcoords(rgenename) #get coords of recipient gene
        newseq = scaffolds[rGI].seq[0: rstart - 1] + dbog[dgenename].seq + scaffolds[rGI].seq[rend:]
        contigstart, contigend = max( rstart-contiglen/2, 0), min( rstart+contiglen/2, len(newseq) - 1 ) #get start and end sites for the hybrid sequenc
        adjcontiglen = contigend - contigstart
        dgenelen = len( dbog[dgenename].seq )

    genelist = []
    for start, end, genename in sorted( rcoords[rGI] ):
        overlap = 0
	if contigstart < start and start < rstart:
            overlap = calc_overlap( contigstart, rstart, start, end )
	elif rstart+dgenelen < start and start < contigend:
	    overlap = calc_overlap( rstart+dgenelen, contigend, start, end )
	if overlap > 0.8:
	    genelist.append( [start, end, genename] )
    
    new_SeqRec = SeqRecord( newseq[contigstart:contigend] )
    new_SeqRec.description = 'donor:' + donor + '==' + dgenename + '||recipient:' + recipient + '==' + rgenename
    if len(genelist) > 0:
        status = True
    else:
        status = False
    
    return status, new_SeqRec, genelist

# ---------------------------------------------------------------
# main
# ---------------------------------------------------------------

def main():
	outname= "-".join( [args.recipient, args.donor] )
	escape_counter = 0
        successes = 0

	fh_fnt = open( outname+".fnt", "w" )
        fh_ans = open( outname+"-ans.txt", "w" )

	while successes < args.ngenes:
		escape_counter += 1
		if escape_counter >= 25*args.ngenes:
			sys.exit( "got stuck in the loop!" )

		# actually load the scaffolds and the genes
		scaffolds = load_scaffolds( args.recipient )
		rbog = load_bog( args.recipient )
		dbog = load_bog( args.donor )
		rcoords = bog2coords( rbog )
	
		# pick a donor and recipient genes and make a contig that contains surrounding genes
		status, seqrec, genelist = False, SeqRecord( '' ), []
		while status == False:
			status, seqrec, genelist = make_contig( args.donor, args.recipient, dbog, rbog, rcoords, scaffolds )
		successes += 1
		donorinfo, recipinfo = seqrec.description.split('||')[0], seqrec.description.split('||')[1]
		dGI, rGI = donorinfo.split('==')[1], recipinfo.split('==')[1]

		# generatefasta
		seqrec.id = "%s|contig%05d|" % ( outname, successes )
		fh_fnt.write( ">" + seqrec.id + '\n')
		fh_fnt.write( str( seqrec.seq ) + '\n' )

	        #generate answerkey
		fh_ans.write( '\t'.join( [str(x) for x in [args.taxadiff, seqrec.id, args.reciptaxa, args.donortaxa, len(seqrec.seq), dGI, len(genelist), genelist]]) + '\n' )

	fh_fnt.close()
	fh_ans.close()


if __name__ == "__main__":
    main()
