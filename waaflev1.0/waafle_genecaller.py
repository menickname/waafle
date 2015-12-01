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
        "-lap", "--overlap",
        default=0.5,
        type=float,
        help="overlap at which to merge hits into genes",
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
                negcoordslist.append( [hit.qstart, hit.qend] )
            else:
                poscoordslist.append( [hit.qstart, hit.qend] )
    return poscoordslist, negcoordslist

def coords2groups( coordslist, overlap_thresh ):
    """
    Bin a list of hits into groups by start/end coordinates.
    """
    groups = []
    for start, end in coordslist:
        group_add = False 
        if len( groups )== 0:
            groups.append( [start, end] )
	    continue
        else:
            for i in range( len( groups ) ):
                overlap = wu.calc_overlap( groups[i][0], groups[i][1], start, end )
                if overlap >= overlap_thresh:
                   coord_sorted = sorted( [groups[i][0], groups[i][1], start, end] )
                   groups[i] = [coord_sorted[0], coord_sorted[3]]
                   group_add = True    
                if i == len( groups ) - 1 and group_add == False:
                   groups.append( [start, end] )
    return groups

def groups2genes( groups, overlap_thresh ):
    """
    Sort groups by start/end coordinates.
    Bin a list of groups into genes by start/end coordinates (of groups).
    """
    genes = []
    #Sort groups by end and start coordinates
    groups.sort( key=itemgetter( 1 ), reverse=True )
    groups.sort( key=itemgetter( 0 ) )
    if len(groups) == 1:
        genes = groups
    else:
        genes = coords2groups( groups, overlap_thresh )
    return genes

def filter_genes( genelist, sign, hitlist, overlap, length, scov_genes, scov_hits ):
    """
    First, group hits into genes.
    Second, filter genes based on length, scoverage, or # of BLAST hits.
    """
    genelist_filtered = []
    for start, end in genelist:
        gene_hits, gene_info = wu.hits2genes( start, end, sign, hitlist, overlap, scov_hits )
	genelen = end - start + 1
	gene_score, gene_percid, gene_cov, gene_starts, gene_ends = wu.score_hits( gene_hits, start, end )
            
        #if genelen >= length:
        #    genelist_filtered.append( [start, end, sign, gene_score, gene_info[1], gene_info[2]] )
        if gene_cov >= scov_genes:
            genelist_filtered.append( [start, end, sign, gene_score, gene_info[1], gene_info[2]] )
        #if gene_info[1] > some_count:
        #    genelist_filtered.append( [start, end, sign, gene_score, gene_info[1], gene_info[2]] )
    return genelist_filtered

def printgff( gffrow ):
    """
    Format the gff class into a ordered list for printing.
    """
    gfflist = [ gffrow.seqname,
                gffrow.source,
                gffrow.feature,
                gffrow.start,
                gffrow.end,
                gffrow.score,
                gffrow.strand,
                gffrow.frame,
                gffrow.attribute
                ]
    return gfflist

def writegff( contig, posgenelist, neggenelist ):
    """
    For each contig, sort all genes by start coordinate regardless of strand.
    Output each gene in gff format to a new file.
    """
    allgenelist = posgenelist + neggenelist
    allgenelist.sort( key=itemgetter( 0 ) )
    counter = 1
    gfflist = []
    for gene in allgenelist: 
        gffrow = wu.GFF( [] )
        gffrow.seqname = contig
        gffrow.source = 'WAAFLE'
        gffrow.feature = 'CDS'
        gffrow.start = gene[0]
        gffrow.end = gene[1]
        gffrow.score = gene[3]
        gffrow.strand = gene[2]
        gffrow.frame = '0'
        ID = 'ID=' + contig + '_' + str(counter)
        Uniref50 = 'Uniref50=' + gene[4]
        Uniref90 = 'Uniref90=' + gene[5]
        gffrow.attribute = ';'.join( [ID, Uniref50, Uniref90] )
        counter += 1
        gffline = printgff( gffrow )
        gfflist.append( gffline )
    return gfflist
        

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

                posgrouplist = coords2groups( poscoordlist, args.overlap )
                neggrouplist = coords2groups( negcoordlist, args.overlap )

                posgenelist = groups2genes( posgrouplist, args.overlap )
                neggenelist = groups2genes( neggrouplist, args.overlap )
		
                pos_filtered = filter_genes( posgenelist, '+', hitlist, args.overlap, args.length, args.scov_genes, args.scov_hits )
                neg_filtered = filter_genes( neggenelist, '-', hitlist, args.overlap, args.length, args.scov_genes, args.scov_hits )
                               
                gfflist = writegff( contig, pos_filtered, neg_filtered )
                for gff in gfflist:
                    writer.writerow( gff )
    
    fh.close()

if __name__ == "__main__":
    main()
