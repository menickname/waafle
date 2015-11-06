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
        help="output from waafle_genecaller or user supplied gff",
        )
    parser.add_argument(
	"-b", "--blast",
	required=True,
	help="output from waafle_search",
	)
    parser.add_argument(
	"-t", "--taxa",
	help="level of taxa to score",
	)
    parser.add_argument(
        "-o", "--out",
        default="waafle-scoredorgs.tsv",
        help="output for scored taxa",
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

    finalscore = np.mean( [x[1][0] for x in dict_indexscore.items()] )
    finalpercid = np.mean( [x[1][1] for x in dict_indexscore.items()] )
    finalgroupcov = np.mean( [x[1][2] for x in dict_indexscore.items()] )
    index_sort = sorted( dict_indexscore.keys() )
    finalstart = index_sort[0]
    finalend = index_sort[len(index_sort)-1]
    return finalscore, finalpercid, finalgroupcov, finalstart, finalend

	
def hits2orgs( contig, gene, hitlist, taxalevel ):
    """
    Assign hits to a gene based on coordinates and strandedness.
    Assign hits within that gene to taxa.
    Generate taxa classes and return the list.
    """
    genelen = gene.end - gene.start + 1
    genehits = []
    for hit in hitlist:
        overlap = wu.calc_overlap( hit.qstart, hit.qend, gene.start, gene.end )
        if gene.strand == "-":
            genestrand = "minus"
        else:
            genestrand = "plus"
        if hit.sstrand == genestrand and overlap >= 0.5:
            genehits.append( hit ) 
    
    dict_orghits = {}
    for hit in genehits:
        org = hit.taxonomy[taxalevel]
        dict_orghits.setdefault( org, [] ).append( hit )
    
    taxalist = []
    for org in dict_orghits:
        taxa = wu.Taxa( [] )
        taxa.taxa = org
        taxa.strand = gene.strand
        taxa.score, taxa.percid, taxa.genecov, taxa.start, taxa.end = score_taxa( dict_orghits[org], genelen )
        taxa.gene = gene.genenum
        taxa.contig = contig
        taxalist.append( taxa )
    if len( taxalist ) == 0:
        taxa = wu.Taxa( [] )
        taxa.taxa = "unknown"
        taxa.strand = gene.strand
        taxa.score, taxa.percid, taxa.genecov, taxa.start, taxa.end = 0, 0, 0, 0, 0
        taxa.gene = gene.genenum
        taxa.contig = contig
        taxalist.append( taxa )
    return taxalist
   
 	    
def print_taxa( taxa ):
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
    return orderedlist

# ---------------------------------------------------------------
# main
# ---------------------------------------------------------------

def main():
    args = get_args()
    taxalevel = c__dict_taxa[args.taxa]

    dict_contiggenes = {}
    for contig, genelist in wu.iter_contig_genes( args.gff ):
        dict_contiggenes[contig] = genelist

    with wu.try_open( args.out, "w" ) as fh:
        writer = csv.writer( fh, dialect="excel-tab" )
        for contig, hitlist in wu.iter_contig_hits( args.blast ):
            if contig in dict_contiggenes:
                genes = dict_contiggenes[contig]
                for gene in genes:
                    taxalist = hits2orgs( contig, gene, hitlist, taxalevel )
                    for taxa in taxalist:
                        taxaline = print_taxa( taxa )
                        writer.writerow( taxaline )
    fh.close()

if __name__ == "__main__":
    main()


