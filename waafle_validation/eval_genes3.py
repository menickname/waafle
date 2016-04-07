#!/usr/bin/env python

from __future__ import print_function # python 2.7+ required
import os, sys, csv, argparse
import waafle_utils as wu
from operator import itemgetter, attrgetter, methodcaller
from collections import Counter
import numpy as np

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
        "-gff", "--gff",
        required=True,
        help="gff"
        )
    parser.add_argument(
        "-rgff", "--refgff",
        required=True,
        help="reference gff"
        )
    parser.add_argument(
        "-lap", "--lap",
        help="lower overlap threshold",
        type=float,
        )
    parser.add_argument(
        "-rlap", "--rlap",
        help="higher overlap threshold",
        type=float,
        )
    parser.add_argument(
        "-s", "--strand",
        help="overlap only genes of same strand, default:False",
        default=False,
        type=bool
        )        
    args = parser.parse_args()
    return args

def calc_overlap( start, end, refstart, refend ):
    """
    Calculate overlap between two hits or genes.
    """
    a1, b1 = start, end
    a2, b2 = refstart, refend
    if b1 < a2 or b2 < a1:
        return 0, 0
    else:
        outleft, inleft, inright, outright = sorted( [a1, b1, a2, b2] )
        denom = end-start+1
        refdenom = refend-refstart+1
        overlap = ( inright - inleft + 1 ) / float( denom )
        refoverlap = ( inright - inleft + 1 ) / float( refdenom )
        return overlap, refoverlap

def det_status( overlap, refoverlap, arglap ):
    if overlap >= arglap and refoverlap >= arglap:
        status = 'call'
    else:
        status = 'no_call'
    return status
    
def partial_to_call( one_array, notone_array, startarr, endarr, strandarr, genestart, geneend, lap ):
    where = np.where( one_array >= lap )
    if len( one_array[where] ) > 1:
        newstart = min( startarr[where] )
        newend =  max( endarr[where] )
        newoverlap, newrefoverlap = calc_overlap( genestart, geneend, newstart, newend )
        newstatus = det_status( newoverlap, newrefoverlap, lap )
        if newstatus == 'call':
            return ['partial-call', newoverlap, newrefoverlap, startarr[where], endarr[where], strandarr[where]]
        else:
            return ['partial-multiple', newoverlap, newrefoverlap, startarr[where], endarr[where], strandarr[where]]
    else:
        return ['partial-single', startarr[where], endarr[where], strandarr[where]]

def print_gff( gffrow ):
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

    #store genes and reference genes
    dict_contigrefgenes, dict_contiggenes = {}, {}     
    for contig, refgenes in wu.iter_contig_genes( args.refgff ):
        dict_contigrefgenes[contig] = refgenes
    for contig, genes in wu.iter_contig_genes( args.gff ):
        dict_contiggenes[contig] = genes

    #compare genes to reference genes
    genes_ans = []
    for contig in dict_contiggenes.keys(): 
        genes = dict_contiggenes[contig]
        
        for gene in genes: #looking at the genes to annotate
            gene.sepattr( gene.attribute )
            gene_found = False
            start_all, end_all, strand_all, lap_all, reflap_all = np.array([]), np.array([]), np.array([]), np.array([]), np.array([])
            line = [contig, gene.genenum, gene.start, gene.end, gene.strand]


            if contig in dict_contigrefgenes:
                refgenes = dict_contigrefgenes[contig]        
                for i in range( len( refgenes ) ): #looking at refgenes for comparison 
                    refgene = refgenes[i]
                    refgene.sepattr( refgene.attribute )

                    if gene_found == True: #end loop if you find a matching gene
                        break
                    
                    #calculate overlap region w/regard to gene and reference gene, respectively
                    overlap, refoverlap = calc_overlap( gene.start, gene.end, refgene.start, refgene.end )
                    status = det_status( overlap, refoverlap, args.lap )
                    if args.strand and gene.strand == refgene.strand: #do it w/regard to strand if flagged
                        overlap, refoverlap = calc_overlap( gene.start, gene.end, refgene.start, refgene.end )
                        status = det_status( overlap, refoverlap, args.lap )

                    start_all = np.append( start_all, refgene.start )
                    end_all = np.append( end_all, refgene.end )
                    strand_all = np.append( strand_all, refgene.strand )
                    lap_all = np.append( lap_all, overlap )
                    reflap_all = np.append( reflap_all, refoverlap )

                    if status == 'call': #if you get a matching gene
                        gene_found = True
                        nline = [refgene.start, refgene.end, refgene.strand, 'called', overlap, refoverlap]
                        result = '\t'.join( str(x) for x in line + nline )
                        print( result ) 

                    if i == len( refgenes )-1 and gene_found == False: #if you never get a total match...
                        zero = np.count_nonzero( lap_all )
                        refzero = np.count_nonzero( reflap_all ) 
                        if zero == 0 and refzero == 0: #no overlap
                            nline = ['NA', 'NA', 'NA', 'no-call', 0, 0]
                            result = '\t'.join( str(x) for x in line + nline )
                            print( result )

                        else: #some overlap
                            if max( lap_all ) >= args.lap and max( reflap_all ) < args.lap: #partial-product_in_reference
                                where = np.where( lap_all >= args.lap )
                                info = partial_to_call( lap_all, reflap_all, start_all, end_all, strand_all, gene.start, gene.end, args.lap )
                                if info[0] == 'partial-single':
                                    nline = [list(info[1])[0], list(info[2])[0], list(info[3])[0], 'partial-single-inref', list(lap_all[where])[0], list(reflap_all[where])[0]] 
                                    result = '\t'.join( str(x) for x in line + nline )
                                    print( result )
                                else:
                                    nline = [','.join( str(y) for y in list(info[3]) ), ','.join( str(y) for y in list(info[4]) ), ','.join( str(y) for y in list(info[5]) ), str(info[0]), info[1], info[2]]
                                    result = '\t'.join( str(x) for x in line + nline )
                                    print( result )

                            elif max( lap_all ) < args.lap and max( reflap_all ) >= args.lap: #partial-reference_in_product
                                where = np.where( reflap_all >= args.lap )
                                info = partial_to_call( reflap_all, lap_all, start_all, end_all, strand_all, gene.start, gene.end, args.lap ) 
                                if info[0] == 'partial-single':
                                    nline = [list(info[1])[0], list(info[2])[0], list(info[3])[0], 'partial-single-ingene', list(lap_all[where])[0], list(reflap_all[where])[0]]
                                    result = '\t'.join( str(x) for x in line + nline )
                                    print( result )
                                else:
                                    nline = [','.join( str(y) for y in list(info[3]) ), ','.join( str(y) for y in list(info[4]) ), ','.join( str(y) for y in list(info[5]) ), info[0], info[1], info[2] ]
                                    result = '\t'.join( str(x) for x in line + nline )
                                    print( result )

                            else:
                                where = np.where( (lap_all - reflap_all) != 0 )
                                nline = [ ','.join( str(y) for y in list(start_all[where])), ','.join( str(y) for y in list(end_all[where])), ','.join( str(y) for y in list(strand_all[where])), 'no-call', ','.join( str(y) for y in list(lap_all[where])), ','.join(list(str(y) for y in reflap_all[where]))]
                                result = '\t'.join( str(x) for x in line + nline )
                                print( result )

            else:
                nline = [ 'NA', 'NA', 'NA', 'no-contig', 0, 0 ]
                result = '\t'.join( str(x) for x in line + nline )
                print( result )            

if __name__ == "__main__":
    main()
