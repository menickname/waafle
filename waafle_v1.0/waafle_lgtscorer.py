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
        "-o", "--out",
        help="results",
        default="waafle-scoredcontigs.tsv",
        )
    parser.add_argument(
        "-s1", "--onebug",
        type=float,
        default=0.8
        )
    parser.add_argument(
        "-s2", "--twobug",
        type=float,
        default=0.8
        )
    parser.add_argument(
        "-u", "--unknown",
        help="By default this is 'False'. If you want to output unknowns, type 'True'.",
        type=bool,
        default=False
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
            overlap = wu.calc_overlap( first_taxa.genestart, first_taxa.geneend, next_taxa.genestart, next_taxa.geneend )
            if overlap > 0.5:
                overlapset.add( (first_taxa.gene, next_taxa.gene) )
    
    if len( overlapset ) == 0:
        return "no_overlap", None
    else:
        return "overlap", overlapset


def spike_unknown( contigarray, taxaorder, unknown ):
    """
    Spike in the unknown organism.
    Assign the unknown values if it has already been called,
    or generate a new unknown.
    """
    final_contigarray = np.copy( contigarray )
    final_taxaorder = [item for item in taxaorder ]
    numorgs, numgenes = final_contigarray.shape
    if numgenes == 0 or unknown:
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
    """
    Create the row for an array, which represents all scores for 1 taxon across all genes.
    """
    genelist_sort = sorted( genelist )
    initrow = [0]*len( genelist_sort )
    for taxon in taxalist:
        if target_taxon == taxon.taxa:
            for i in range( len( genelist_sort ) ):
                if taxon.gene == genelist_sort[i]:
                    initrow[i] = taxon.score
    return initrow


def generate_tables( taxalist ):
    """
    Returns an array for all scores across all bugs and genes in the contig.
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
    """
    If genes on different strands overlap, isolate the corresponding columns and take the average.
    Mask the array as well in order to delete the columns that correspond to those columns.
    """
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
        newcontigarray = np.transpose([replacement_dict[0]])
    return newcontigarray
        

def calc_onebug( array, bugindex ):
    """
    Calculates the one bug score. Equals the maximum of the minimum of scores across all bugs per gene.
    High scores indicate that contig is likely explained by one taxa.
    """
    onebugscore = np.max( np.amin( array, axis=1 ) )
    indices = [i for i,j in enumerate( list(np.amin(array, axis=1) ) ) if j==onebugscore]
    onebuglist = [bugindex[index] for index in indices]
    return onebugscore, onebuglist


def calc_twobug( array, bugindex ):
    """
    Calculates the complement score. Equals the minimum of the maximum of scores across all genes per pair of bugs.
    High scores indicate that contig is likely explained by two taxa.
    """
    complementlist = []
    for i in range( array.shape[0]-1 ):
        bug1 = array[i, :]
        for j in range( i+1, array.shape[0] ):
            bug2 = array[j, :]
            twobugarray = np.array( [bug1, bug2] )
            complement = np.min(np.amax(twobugarray, axis=0))
            bugnames = str(bugindex[i]) + '-' + str(bugindex[j])
            complementlist.append( [ bugnames, complement ] )
    complement_sorted = sorted( complementlist, key=lambda x: x[1], reverse=True )
    twobugscore = complement_sorted[0][1]
    twobuglist = []
    for pair in complement_sorted:
        if pair[1] == twobugscore:
            twobuglist.append( pair[0] )
        elif pair[1] < twobugscore:
            break
    return twobugscore, twobuglist


def calc_donorrecip( contigarray, taxaorder, orgpair ):
    """
    Determines donor and recipient if there are 2 transitions, and the outer genes are from the same bug. (For eg, "ABA").
    The bug in between (eg. "B") is the donor. The bug that explains the outside genes (eg. "A") is the recipient.
    If there is only 1 transition, we call it "unknown". If there is more than 1 transition, we call it a "hybrid."
    """
    top_one, top_two = orgpair[0].split('-')[0], orgpair[0].split('-')[1]
    index_one, index_two = taxaorder.index( top_one ), taxaorder.index( top_two )
    array_one, array_two = contigarray[index_one, :], contigarray[index_two, :]
    twobugmax = np.amax( np.array( [array_one, array_two] ), axis=0 )
    element_list, shift_counter = [], 0
    for i in range( len( twobugmax ) ):
        element = twobugmax[i]
        one_element, two_element = array_one[i], array_two[i]
        if element == one_element:
            element_list.append( top_one )
            state = top_one
        else:
            element_list.append( top_two )
            state = top_two
        if i != 0:
            if element_list[i-1] != state:
                shift_counter += 1
    if shift_counter == 2 and element_list[0] == element_list[len( element_list ) - 1]:
        recipient = element_list[0]
        if recipient == top_one:
            donor = top_two
        else:
            donor = top_one
    elif shift_counter > 2:
        recipient, donor = 'hybrid', 'hybrid' 
    elif shift_counter == 1:
        recipient, donor = 'unknown', 'unknown'
    return recipient, donor
    

def rank_ambiguous( onebugscore, onebuglist, twobugscore, orgpair ):
    """
    If the contig cannot be called "noLGT" or "LGT", decide which one it is closer to and call ambiguous.
    """
    if onebugscore >= twobugscore:
        score = onebugscore
        org = onebuglist
        status = "ambiguous-NoLGT"
    else:
        score = twobugscore
        org = orgpair
        status = "ambiguous-LGT"
    return status, score, org


def print_result( contig, status, score, orgs, recipient, donor ):
    """
    Print out the results in a neat list.
    """
    orglist = []
    if status == "ambiguous-NoLGT" or status == "NoLGT":
        if len(orgs) > 1:
            if len(orgs) > 3:
                orglist = orgs[0:3]
                orglist.append( str(len( orgs) - 1 -3 ) + ' other taxa with same score' )
                orgs = orglist
    finalorgs = ';'.join( orgs )
    return [str(x) for x in [contig, status, score, finalorgs, recipient, donor] ]
    

# ---------------------------------------------------------------
# main
# ---------------------------------------------------------------

def main():
    """
    1. Generate array tables for each contig. If the contig has genes on opposite strands, check if they overlap.
    2. If so, generate a new array table where those genes are merged.
    3. Separate out contigs with only 1 taxa or gene.
    4. For remaining contigs, perform LGT algorithm.
    """
    args = get_args()
    with wu.try_open( args.out, "w" ) as fh:
        writer = csv.writer( fh, dialect="excel-tab" )
        writer.writerow( ['contig', 'status', 'score', 'taxa', 'receipient', 'donor'] )
        for contig, taxalist in wu.iter_contig_taxa( args.input ):
            contigarray, taxaorder = generate_tables( taxalist )
            recipient, donor = "NA", "NA"

            if count_strands( taxalist ) == 2:
                gene_status, overlapgenes = find_overlap_genes( taxalist )
                if gene_status == "overlap":
                    contigarray = account_overlap( contigarray, overlapgenes )
            
            spikearray, spikeorder = spike_unknown( contigarray, taxaorder, args.unknown )
            onebugscore, onebuglist = calc_onebug( spikearray, spikeorder )
            numbugs, numgenes = spikearray.shape
                    
            if onebugscore >= args.onebug or numgenes == 1 or numbugs == 1:
                score, taxa = onebugscore, onebuglist
                status = "NoLGT"
            else:
                twobugscore, orgpair = calc_twobug( spikearray, spikeorder )
                if twobugscore >= args.twobug:
                    score, taxa = twobugscore, orgpair
                    recipient, donor = calc_donorrecip( spikearray, spikeorder, orgpair )
                    status = "LGT"
                else:
                    status, score, taxa = rank_ambiguous( onebugscore, onebuglist, twobugscore, orgpair )
        
            result_line = print_result( contig, status, score, taxa, recipient, donor )
            writer.writerow( result_line )
            

if __name__ == "__main__":
    main()
