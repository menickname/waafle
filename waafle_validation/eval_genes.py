#!/usr/bin/env python

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
        "-w_gff", "--waaflegff",
        required=True,
        help="output from waafle_genecaller"
        )
    parser.add_argument(
        "-o_gff", "--othergff",
        required=True,
        help="gff to compare"
        )
    args = parser.parse_args()
    return args

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

# ---------------------------------------------------------------
# main
# ---------------------------------------------------------------

def main():
    args = get_args()
    
    fh_missed = open( "missed_contigs.tsv", "w" )
    dict_contigogenes, dict_contigwgenes = {}, {} 
    
    for contig, ogenes in wu.iter_contig_genes( args.othergff ):
        dict_contigogenes[contig] = ogenes

    for contig, wgenes in wu.iter_contig_genes( args.waaflegff ):
        dict_contigwgenes[contig] = wgenes

    missed_genes = []
    for contig in dict_contigwgenes.keys(): 
        wgenes = dict_contigwgenes[contig]
        
        if contig in dict_contigogenes:
            ogenes = dict_contigogenes[contig]
            for wgene in wgenes: #list of genes in waafle
                wgene_found = False
                for i in range( len( ogenes ) ): #compare to list of genes in reference
                    ogene = ogenes[i]
                    if wgene_found == True:
                        break
                    
                    if wgene.strand == ogene.strand:
                        overlap = wu.calc_overlap( wgene.start, wgene.end, ogene.start, ogene.end )
                        if overlap > 0.5:
                            wgene_found = True 
                    if i == len( ogenes ) - 1 and wgene_found == False:
                        missed_genes.append( printgff( wgene ) )

    for gene in missed_genes:
        fh_missed.write( '\t'.join( [str(x) for x in gene] ) )

    fh_missed.close()


if __name__ == "__main__":
    main()
