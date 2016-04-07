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
from collections import Counter

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

c__list_taxa = ["k", "p", "c", "o", "f", "g", "s"]

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
        "-o", "--out",
        default="waafle-scoredorgs.tsv",
        help="output for scored taxa",
        )
    parser.add_argument(
        "-lap", "--overlap_hits",
        help="amount of overlap to include hit in a gene",
        default=0.5,
        type=float,
        )
    parser.add_argument(
        "-scov", "--scov_hits",
        help="cutoff for gene coverage or subject coverage when grouping hits",
        default=0,
        type=float,
        )
    parser.add_argument(
        "-t", "--taxalevel",
        help="level of taxa to score",
        )
    parser.add_argument(
        "-s", "--strand",
        help="strand-specific splitting, default=False",
        default=False,
        type=bool,
        )
    args = parser.parse_args()
    return args

def calc_overlap( a1, b1, a2, b2 ):
    """ compute overlap between two intervals """
    if b1 < a2 or b2 < a1:
        return 0
    else:
        outleft, inleft, inright, outright = sorted( [a1, b1, a2, b2] )
        denom = min( ( b1 - a1 + 1 ), ( b2 - a2 + 1 ) )
        return ( inright - inleft + 1 ) / float( denom )

def hits2genes( gene, hits, strand_specific, lap, scov, taxalevel ):
    genehits = []
    taxaset = set()
    uniref50, uniref90 = [], []
    for hit in hits:
        hitstrand = wu.convert_strand( hit.sstrand )
        if ( hitstrand == gene.strand or not strand_specific ) and hit.scov_modified > scov:
            overlap = calc_overlap( hit.qstart, hit.qend, gene.start, gene.end )    
            if overlap > lap:
                genehits.append( hit )
                taxaset.add( hit.taxonomy[taxalevel] )
                uniref50.append( hit.uniref50 )
                uniref90.append( hit.uniref90 )
    uniref50_c = ','.join( [ str(x) + ':' + str(y) for x, y in Counter( uniref50 ).most_common( 3 )] )
    uniref90_c = ','.join( [ str(x) + ':' + str(y) for x, y in Counter( uniref90 ).most_common( 3 )] )
    info = [genehits, taxaset, uniref50_c, uniref90_c, taxalevel ]
    return info

def score_taxa( gene, info, contiglen ):
    """
    Group hits that correspond to specific taxa.
    If the gene did not have hits to begin with, output an "unknown" taxa with score 0.
    If the gene has hits that corresponds to orgs, score the orgs and output the taxa annotations.
    """
    genehits, taxaset, uniref50_c, uniref90_c, taxalevel = info
    gene.sepattr( gene.attribute )
    genelen = gene.end - gene.start + 1
    taxalist = []
    if len( genehits ) == 0:
        #set unknown taxa or output
        taxa = c__list_taxa[taxalevel] + '__Unknown'
        start_list, end_list, score = gene.start, gene.end, 1 
        taxa = wu.Taxa( [gene.seqname, contiglen, gene.genenum, gene.strand, gene.start, gene.end, taxa, start_list, end_list, score, 'No_Uniref50', 'No_Uniref90', 0] )
        taxalist.append( taxa )
    else:
        dict_orgscores = {}
        for taxa in taxaset:
            dict_orgscores.setdefault( taxa, np.zeros( genelen ) )
            numhits = 0
            for hit in genehits:
                if hit.taxonomy[taxalevel] == taxa:
                    numhits += 1
                    orgarray = dict_orgscores.get( taxa )
                    hitarray = np.zeros( genelen )
                    arraystart = max( hit.qstart - gene.start, 0 )
                    arrayend = min( hit.qend - gene.start + 1, genelen )
                    hitarray[arraystart: arrayend] = (hit.pident/float(100))*hit.scov_modified
                    dict_orgscores[taxa] = np.maximum( hitarray, orgarray )
            score = np.mean( dict_orgscores[taxa] )
            startstop = np.split(np.argwhere( dict_orgscores[taxa] ), np.where(np.diff(np.argwhere( dict_orgscores[taxa] ), axis=0)!= 1)[0]+1)
            start_list = ','.join( [str(element[0][0] + gene.start) for element in startstop] )
            end_list = ','.join( [str(element[-1][0] + gene.start) for element in startstop] ) #check this
            taxa = wu.Taxa( [gene.seqname, contiglen, gene.genenum, gene.strand, gene.start, gene.end, taxa, start_list, end_list, score, uniref50_c, uniref90_c, numhits] )
            taxalist.append( taxa )
    return taxalist

# ---------------------------------------------------------------
# main
# ---------------------------------------------------------------

def main():
    """
    This script does:
    1) Groups hits into orgs per gene.
    2) Score each org.
    4) Prints out the taxa into a new file.
    """
    args = get_args()
    taxalevel = c__dict_taxa[args.taxalevel]

    fh = wu.try_open( args.out, "w" )
    writertaxa = csv.writer( fh, dialect="excel-tab" )
    writertaxa.writerow( ["contig", "contiglen", "genenum", "strand", "genestart", "geneend", "taxa", "taxastart", "taxaend", "score", "uniref50", "uniref90", "orghitnum"] )
       
    #Build dictionary of genes
    dict_genes = {}
    for contig, genelist in wu.iter_contig_genes( args.gff ):
        dict_genes[contig] = genelist

    #Group taxa based on new genes
    for contig, hitlist in wu.iter_contig_hits( args.blast ):
        if contig in dict_genes:
            genelist = dict_genes[contig]
            contiglen = hitlist[0].qlen
        
            for gene in genelist:
                info = hits2genes( gene, hitlist, args.strand, args.overlap_hits, args.scov_hits, taxalevel )
                taxalist = score_taxa( gene, info, contiglen )
           
                for taxa in taxalist:
                    writertaxa.writerow( wu.print_taxa( taxa ) )
    fh.close()

if __name__ == "__main__":
    main()


