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

def calc_overlap( w_start, w_end, ans_start, ans_end ):
    """
    Calculate overlap between two hits or genes.
    """
    if w_start - ans_end > 0 or ans_start - w_end > 0:
        status = 'no_call'
    else:
        w_divisor = float( w_end - w_start + 1 )
        ans_divisor = float( ans_end - ans_start + 1 )
        coord_sorted = sorted( [w_start, w_end, ans_start, ans_end] )
        ans_overlap = float( ( coord_sorted[2] - coord_sorted[1] )/ans_divisor )
        w_overlap = float( ( coord_sorted[2] - coord_sorted[1] )/w_divisor )
        if w_overlap >= 0.9 and ans_overlap >= 0.5:
            status = 'call'
        elif w_overlap >= 0.9 and ans_overlap < 0.5:
            status = 'partial_call'
        else:
            status = 'no_call'
    return status

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
    dict_contigogenes, dict_contigwgenes = {}, {} 
    
    for contig, ogenes in wu.iter_contig_genes( args.othergff ):
        dict_contigogenes[contig] = ogenes

    for contig, wgenes in wu.iter_contig_genes( args.waaflegff ):
        dict_contigwgenes[contig] = wgenes

    genes_ans = []
    for contig in dict_contigwgenes.keys(): 
        wgenes = dict_contigwgenes[contig]
        
        if contig in dict_contigogenes:
            ogenes = dict_contigogenes[contig]
            for wgene in wgenes: #list of genes to annotate
                wgene_found = False
                genes_missed = []
                for i in range( len( ogenes ) ): #list of genes to compare to 
                    ogene = ogenes[i]
                    strand_counter = 0
    
                    if wgene_found == True:
                        break
                    
                    if wgene.strand == ogene.strand:
                        strand_counter += 1
                        status = calc_overlap( wgene.start, wgene.end, ogene.start, ogene.end )
                        
                        if status == 'call':
                            wgene_found = True
                            line = '\t'.join( str(x) for x in [contig, wgene.genenum, wgene.start, wgene.end, ogene.start, ogene.end, 'called'] )
                        if status == 'partial_call':
                            genes_missed.append( [wgene.genenum, wgene.start, wgene.end, ogene.start, ogene.end] )

                    if i == len( ogenes ) - 1 and wgene_found == False:
                        if len( genes_missed ) > 0:
                            newname = ''
                            for x in genes_missed:
                                newname = newname + '\t'.join( [str(y) for y in x] ) + ';'
                            line = '\t'.join( [contig, newname, 'partial'] ) 
                            wgene_found = True
                        else:
                            if strand_counter == 0:
                                line = '\t'.join( str(x) for x in [contig, wgene.genenum, wgene.start, wgene.end, ogene.start, ogene.end, 'missed_strand'] )
                            else:
                                line = '\t'.join( str(x) for x in [contig, wgene.genenum, wgene.start, wgene.end, ogene.start, ogene.end, 'missed_overlap'] )
                print( line )
        else:
            for wgene in wgenes: #genes not in comparison because no genes in contig
                line = '\t'.join( str(x) for x in [contig, wgene.genenum, wgene.start, wgene.end, 'NA', 'NA', 'missed_nogene'] )
                print( line )


if __name__ == "__main__":
    main()
