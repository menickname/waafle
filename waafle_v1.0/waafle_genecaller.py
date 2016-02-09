#!/usr/bin/env python

"""
WAAFLE MODULE: waafle_genecaller.py

Authors:
Tiffany Hsu
Eric Franzosa

This script interprets waafle blast output to make gene calls.
Run with the "-h" flag for usage help.

Details of the GTF/GFF file format:

Fields must be tab-separated. 
Empty columns should be denoted with a '.'.

0: seqname - name of the chromosome or scaffold
1: source - name of the program that generated this feature
2: feature - feature type name, e.g. Gene, Variation, Similarity
3: start - Start position of the feature, with sequence numbering starting at 1.
4: end - End position of the feature, with sequence numbering starting at 1.
5: score - A floating point value.
6: strand - defined as + (forward) or - (reverse).
7: frame - One of '0', '1' or '2'.
8: attribute - A semicolon-separated list of tag-value pairs.
"""

from __future__ import print_function # python 2.7+ required
import os, sys, csv, argparse
import waafle_utils as wu
from operator import itemgetter, attrgetter, methodcaller
from collections import Counter


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
        "-i", "--input",
        required=True,
        help="output from waafle_search"
        )
    parser.add_argument( 
        "-o", "--out",
        default="waafle-genes.gff",
        help="waafle gene calls",
        )
    parser.add_argument(
        "-lap_h", "--overlap_hits",
        default=0.5,
        type=float,
        help="overlap at which to merge hits into groups",
        )
    parser.add_argument(
        "-lap_g", "--overlap_genes",
        default=0.5,
        type=float,
        help="overlap at which to merge groups into genes",
        )
    parser.add_argument(
        "-l", "--length",
        default=0,
        type=float,
        help="length constraint for genes",
        )
    parser.add_argument(
        "-scov_h", "--scov_hits",
        default=0,
        type=float,
        help="scoverage filter for hits",
        )
    parser.add_argument(
        "-scov_g", "--scov_genes",
        default=0,
        type=float,
        help="scoverage filter for genes",
        )
    args = parser.parse_args()
    return args

def hits2coords( hitlist, scov ):
    """
    For a list of hits, filter by some metric (current method: scoverage )
    Get start/end sites for positive and negative strands (two lists).
    """
    poscoordslist, negcoordslist = [], []
    for hit in hitlist:
        if hit.scov_modified >= scov:
            if hit.sstrand == 'minus':
                negcoordslist.append( [hit.qstart, hit.qend, '-'] )
            else:
                poscoordslist.append( [hit.qstart, hit.qend, '+'] )
    return poscoordslist, negcoordslist
"""
def coords2groups( coordslist, overlap_thresh ):
    #Bin a list of hits into groups by start/end coordinates.
    groups = []
    for start, end, strand in coordslist:
        group_add = False 
        if len( groups )== 0:
            groups.append( [start, end, strand] )
	    continue
        else:
            for i in range( len( groups ) ):
                overlap = wu.calc_overlap( groups[i][0], groups[i][1], start, end )
                if overlap >= overlap_thresh:
                   coord_sorted = sorted( [groups[i][0], groups[i][1], start, end] )
                   groups[i] = [coord_sorted[0], coord_sorted[3], strand]
                   group_add = True    
                if i == len( groups ) - 1 and group_add == False:
                   groups.append( [start, end, strand] )
    return groups

def groups2genes( groups, overlap_thresh ):
    #Sort groups by start/end coordinates.
    #Bin a list of groups into genes by start/end coordinates (of groups).
    genes = []
    #Sort groups by end and start coordinates
    for gene in groups:
        length = int(gene[1]) - int(gene[0])
        gene.append( length )
    groups.sort( key=itemgetter( 3 ), reverse=True )
    group_sorted = [x[0:3] for x in groups]
    if len(groups) == 1:
        genes = group_sorted
    else:
        genes = coords2groups( group_sorted, overlap_thresh )
    return genes
"""     

# ---------------------------------------------------------------
# main
# ---------------------------------------------------------------

def main():
    """
    This script does the following:
    1) Sorts hits per contig by length and bitscore.
    2) Separate hits by strand and potentially filter by scoverage_modified.
    3) Form genes by grouping overlapping hits, and then grouping overlapping groups.
    4) Group hits by genes to calculate a scores and annotate unirefs.
    5) Filter genes by length, score, or number of hits.
    6) Print out gff.
    """
    args = get_args()
    with wu.try_open( args.out, "w" ) as fh:
            writer = csv.writer( fh, dialect="excel-tab" )

            for contig, hitlist in wu.iter_contig_hits( args.input ):
                hitlist_sorted = sorted( hitlist, key=attrgetter( 'length', 'bitscore' ), reverse=True )
                poscoordlist, negcoordlist = hits2coords( hitlist_sorted, args.scov_hits )

                posgrouplist = wu.coords2groups( poscoordlist, args.overlap_hits )
                neggrouplist = wu.coords2groups( negcoordlist, args.overlap_hits )

                posgenelist = wu.groups2genes( posgrouplist, args.overlap_genes )
                neggenelist = wu.groups2genes( neggrouplist, args.overlap_genes )
                pos_filtered = wu.filter_genes( posgenelist, hitlist, args.overlap_hits, args.length, args.scov_genes, args.scov_hits )
                neg_filtered = wu.filter_genes( neggenelist, hitlist, args.overlap_hits, args.length, args.scov_genes, args.scov_hits )
                
                allgenelist = pos_filtered + neg_filtered
                allgenelist.sort( key=itemgetter( 0 ) )
                gfflist = wu.write_gff( contig, allgenelist )
                for gff in gfflist:
                    writer.writerow( gff )
    fh.close()

if __name__ == "__main__":
    main()
