#!/usr/bin/env python

"""
WAAFLE MODULE: waafle_orgscorer.py

Authors:
Tiffany Hsu
Eric Franzosa

This script combines gff output with BLAST hits and annotates genes with microbial taxa.
Run with the "-h" flag for usage help.

"""

from __future__ import print_function # python 2.7+ required
import os, sys, csv, argparse
import waafle_utils as wu
import numpy as np

# ---------------------------------------------------------------
# constants
# ---------------------------------------------------------------

c__dict_taxa = {
    "k": 0,
    "p": 1,
    "c": 2,
    "o": 3,
    "f": 4,
    "g": 5,
    "s": 6,
}

# ---------------------------------------------------------------
# functions
# ---------------------------------------------------------------
def get_args():
    """
    Get arguments passed to script
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument(
        "-g", "--gff",
        required=True,
        help="output from waafle_genecaller or user supplied gff"
        )
    parser.add_argument(
	"-b", "--blast",
	required=True,
	help="output from waafle_search"
	)
    parser.add_argument(
	"-t", "--taxa",
	help="level of taxa to score"
	)
    args = parser.parse_args()
    return args

def score_taxa( hitlist, genelen ):
    """
    Loop through hits that correspond to a single taxon and gene.
    Calculate a score for each hit. Score = Percent identity * Coverage.
    Assign the highest score (from any hit) to each base position covered by all hits.
    Calculate a final score by averaging these scores. 
    """
    dict_indexscore = {}
    for hit in hitlist:
	groupcov = hit.length/float( genelen )
        percid = hit.pident/float( 100 )
	score = float( groupcov ) * percid
	info = [score, percid, groupcov]
	for coordinate in range( hit.qstart, hit.qend + 1 ):
	    if dict_indexscore.get( coordinate, [0, 0, 0] )[0] != 0:
		oldscore = dict_indexscore[coordinate][0]
		if oldscore < score:
			dict_indexscore[coordinate] = info
	    else:
		dict_indexscore[coordinate] = info

    indexlist, scorelist, percidlist, groupcovlist = [], [], [], []
    for index in dict_indexscore:
        indexlist.append( index )
	scorelist.append( dict_indexscore[index][0] )
	percidlist.append( dict_indexscore[index][1] )
        groupcovlist.append( dict_indexscore[index][2] )

    finalscore = np.mean( scorelist )
    finalpercid = np.mean( percidlist )
    finalgroupcov = np.mean( groupcovlist )
    index_sort = sorted( indexlist )
    finalstart = index_sort[0]
    finalend = index_sort[len(index_sort)-1]
    return finalscore, finalpercid, finalgroupcov, finalstart, finalend
	
def hits2genes( hitlist, gene ):
    """
    Assign hits to genes based on coordinates and strandedness.
    """
    genehits = []
    for hit in hitlist:
        myhit = wu.Hit( hit )
        overlap = wu.calc_overlap( myhit.qstart, myhit.qend, gene.start, gene.end )
        if gene.strand == "-":
            genestrand = "minus"
        else:
            genestrand = "plus"
        if myhit.sstrand == genestrand and overlap >= 0.5:
            genehits.append( myhit )
    return genehits
        
def hits2orgs( genelist, hitlist, taxalevel, contig ):
    """
    Loop through genes in order and assign hits to each gene.
    Assign hits (within a gene) to each organism present.
    For each organism, generate a class Taxa object.
    """
    taxalist = []
    for gene in genelist:
        mygene = wu.GFF( gene )
        genelen = mygene.end - mygene.start + 1
	genehits = hits2genes( hitlist, mygene )
	
        dict_orghits = {}
	for myhit in genehits:
            org = myhit.taxonomy[taxalevel]
            dict_orghits.setdefault( org, [] ).append( myhit )
    
    	for org in dict_orghits:
            taxa = wu.Taxa( [] )
	    taxa.taxa = org
            taxa.strand = mygene.strand
            taxa.score, taxa.percid, taxa.genecov, taxa.start, taxa.end = score_taxa( dict_orghits[org], genelen )
            taxa.gene = mygene.genenum
            taxa.contig = contig
            taxalist.append( taxa )
    return taxalist
   
 	    
def print_taxa( taxalist ):
    for taxa in taxalist:
        orderedlist = [ str( taxa.contig ),
                        str( taxa.gene ),
                        str( taxa.strand ),
                        str( taxa.start ),
                        str( taxa.end ),
                        str( taxa.taxa ),
                        str( taxa.score ),
                        str( taxa.percid ),
                        str( taxa.genecov ),
                        ]
        print( '\t'.join( orderedlist ) )

# ---------------------------------------------------------------
# main
# ---------------------------------------------------------------

def main():
    args = get_args()
    taxalevel = c__dict_taxa[args.taxa]

    dict_contigitems = {}
    for contig, genelist in wu.iter_contigs( args.gff ):	
	mycontig = wu.Contig( [] )
	mycontig.genes = genelist
	mycontig.numgenes = len( genelist )
	dict_contigitems[contig] = mycontig
    
    for contig, hitlist in wu.iter_contigs( args.blast ):
        if contig in dict_contigitems:
            mycontig = dict_contigitems[contig]
            mycontig.hits = hitlist
            mycontig.numhits = len( hitlist )
            mycontig.length = wu.Hit( hitlist[0] ).qlen 

    for contig in dict_contigitems:
        mycontig = dict_contigitems[contig]
        taxalist = hits2orgs( mycontig.genes, mycontig.hits, taxalevel, contig )
        mycontig.taxa = taxalist
        print_taxa( taxalist )	 


if __name__ == "__main__":
    main()


