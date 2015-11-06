#!/usr/bin/env python

"""
WAAFLE MODULE: waafle_lgtscorer.py

Authors:
Tiffany Hsu
Eric Franzosa

This script takes the output from waafle_orgscorer.py and determines which contigs likely have LGT.
Run with the "-h" flag for usage help.

"""

from __future__ import print_function # python 2.7+ required
import os, sys, csv, argparse
import waafle_utils as wu
import numpy as np
import numpy.ma as ma

# ---------------------------------------------------------------
# constants
# ---------------------------------------------------------------

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
        help="output from waafle_orgscorer",
        )
    parser.add_argument(
        "-o", "--output",
        help="results",
        default="waafle-scoredcontigs.tsv",
        )
    parser.add_argument(
        "-1", "--onebug",
        default=0.8
        )
    args = parser.parse_args()
    return args

def count_strands( taxalist ):
    """
    Determine whether both strands in the contig have genes.
    """
    posgeneset, neggeneset = set([]), set([])
    for taxa in taxalist:
        if taxa.strand == "+":
            posgeneset.add( taxa.gene )
        else:
            neggeneset.add( taxa.gene )
    if len( posgeneset ) != 0 and len( neggeneset ) != 0:
        return 2
    else:
        return 1

def find_overlap_genes( taxalist ):
    """
    For contigs with genes on both strands, determine if genes overlap.
    For contigs without overlapping genes, return "False."
    For contigs with overlapping genes, return "True", and the set of genes overlapping.
    """
    for i in range( len( taxalist ) - 1 ):
        overlapset = set([])
        first_taxa = taxalist[i]
        next_taxa = taxalist[i + 1]
        if first_taxa.gene != next_taxa.gene and first_taxa.strand != next_taxa.strand:
            overlap = wu.calc_overlap( first_taxa.start, first_taxa.end, next_taxa.start, next_taxa.end )
            if overlap > 0.5:
                overlapset.add( (first_taxa.gene, next_taxa.gene) )
    
    if len( overlapset ) == 0:
        return "no_overlap", None
    else:
        return "overlap", overlapset


def spike_unknown( contigarray, taxaorder ):
    """
    Spike in the unknown organism.
    Assign the unknown values if it has already been called,
    or generate a new unknown.
    """
    final_contigarray = np.array( contigarray )
    final_taxaorder = [item for item in taxaorder ]
    numorgs, numgenes = np.shape( final_contigarray )
    if numgenes == 0:
        return final_contigarray, final_taxaorder
    else:
        if "unknown" in taxaorder:
            unknown_index = taxaorder.index( "unknown" )
            for i in range( numgenes ):
                unknown_score = 1 - np.max( contigarray[:, i] )
                final_contigarray[ unknown_index, i ] = unknown_score
            return final_contigarray, final_taxaorder
        else:
            unknown = []
            for i in range( numgenes ):
                unknown_score = 1 - np.max( contigarray[:, i] )
                unknown.append( unknown_score )
            final_contigarray = np.append( final_contigarray, [unknown], axis=0 )
            final_taxaorder.append( 'unknown' )
            return final_contigarray, final_taxaorder

            
def generate_array( taxalist, genelist, target_taxon ):
    genelist_sort = sorted( genelist )
    initarray = np.zeros( len( genelist_sort ) )
    for taxon in taxalist:
        if target_taxon == taxon.taxa:
            for i in range( len( genelist_sort ) ):
                if taxon.gene == genelist_sort[i]:
                    initarray[i] = taxon.score
    return initarray


def generate_tables( taxalist ):
    """
    Return all the genes in the contig in an array.
    """
    taxaset, geneset = set([]), set([])
    for taxa in taxalist:
        taxaset.add( taxa.taxa )
        geneset.add( taxa.gene )
    taxaorder, contigarray = [], []
    for taxon in taxaset:
        taxonarray = generate_array( taxalist, geneset, taxon )
        taxaorder.append( taxon )
        contigarray.append( taxonarray )
    return np.array( contigarray ), taxaorder


def account_overlap( contigarray, overlapset  ):
    newmaskarray = np.copy( contigarray )
    numbugs, numgenes = contigarray.shape
    replacement_dict = {}
    counter = 1
    for firstgene, secondgene in overlapset:
        firstindex = firstgene - 1
        newindex = firstgene - counter
        replacement_dict[ newindex ] = np.average( contigarray[:, firstindex:secondgene], axis=1 )
        counter += 1
        for i in range( numgenes ):
            if i == firstgene - 1 or i == secondgene - 1:
                newmaskarray[:, i] = 1
            else:
                newmaskarray[:, i] = 0

    if numgenes > 2:
        newcontigarray = ma.compress_cols( ma.array( np.copy( contigarray ), mask=newmaskarray ) )
        for index in replacement_dict.keys():
            numrows, numcol = newcontigarray.shape
            if index - 1 == numcol:
                newcontigarray = np.append( newcontigarray, replacement_dict[index], axis=1 )
            else:
                newcontigarray = np.insert( newcontigarray, index, replacement_dict[index], axis=1 )
    else:
        newcontigarray = replacement_dict[0]
    return newcontigarray
        

def calc_onebug( array, bugindex ):
    """
    Calculates the one bug score. High scores indicate that contig is likely explained by one taxa.
    """
    onebugscore = np.max( np.amin( array, axis=1 ) )
    indices = [i for i,j in enumerate( list(np.amin(array, axis=1) ) ) if j==onebugscore]
    onebuglist = [bugindex[index] for index in indices]
    return onebugscore, onebuglist

def calc_complement( array, bugindex ):
    """
    """
    pass

    
    

# ---------------------------------------------------------------
# main
# ---------------------------------------------------------------

def main():
    """
    1. Generate array tables for each contig.
    2. Determine if the contig has genes on opposite strands.
    3. If so, determine if these genes overlap.
    4. From 3, collapse those 2 genes into 1 column within the array. Average scores across those 2 columns.
    5. For all contigs, perform LGT algorithm.
    """
    args = get_args()
    for contig, taxalist in wu.iter_contig_taxa( args.input ):
        contigarray, taxaorder = generate_tables( taxalist )
        print( contig )
        if count_strands( taxalist ) == 2:
            gene_status, overlapgenes = find_overlap_genes( taxalist )
            if gene_status == "overlap":
                contigarray = account_overlap( contigarray, overlapgenes )
        
        spikearray, spikeorder = spike_unknown( contigarray, taxaorder )
        numgenes, numbugs = spikearray.shape
        onebugscore, onebuglist = calc_onebug( spikearray, spikeorder )
        if onebugscore >= args.onebug or numgenes == 1:
            #call as one taxa/no lgt
            pass
        else:
            calc_complement( spikearray, spikeorder )

        
if __name__ == "__main__":
    main()

