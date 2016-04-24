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
from operator import itemgetter, attrgetter, methodcaller
import waafle_utils as wu
import numpy as np
import numpy.ma as ma

# ---------------------------------------------------------------
# constants
# ---------------------------------------------------------------

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
        default=0.5
        )
    parser.add_argument(
        "-u", "--unknown",
        help="By default this is 'False'. If you want to output unknowns, type 'True'.",
        type=bool,
        default=False
        )
    args = parser.parse_args()
    return args


def spike_unknown( contigarray, taxaorder, unknown ):
    """
    Spike in the unknown organism.
    Assign the unknown values if it has already been called,
    or generate a new unknown.
    """
    final_contigarray = np.copy( contigarray )
    final_taxaorder = [item for item in taxaorder ]
    numorgs, numgenes = final_contigarray.shape
    taxalevel = taxaorder[0].split('__')[0]
    unknown_name = taxalevel + "__Unknown"
    
    if numgenes == 0:
        return final_contigarray, final_taxaorder
    if unknown_name in taxaorder:
        unknown_index = taxaorder.index( unknown_name )
        for i in range( numgenes ):
            unknown_score = 1 - np.max( contigarray[:, i] )
            if contigarray[unknown_index, i] == 1:
                unknown_score = 1
            final_contigarray[ unknown_index, i ] = unknown_score
        return final_contigarray, final_taxaorder
    elif not unknown:
        return final_contigarray, final_taxaorder
    else:
        unknown = []
        for i in range( numgenes ):
            unknown_score = 1 - np.max( contigarray[:, i] )
            unknown.append( unknown_score )
        final_contigarray = np.append( final_contigarray, [unknown], axis=0 )
        final_taxaorder.append( unknown_name )
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
    uniref50, uniref90 = '', ''
    for taxa in taxalist:
        taxaset.add( taxa.taxa )
        geneset.add( taxa.gene )
    taxaorder, contigarray = [], []
    for taxon in taxaset:
        taxonarray = generate_array( taxalist, geneset, taxon )
        taxaorder.append( taxon )
        contigarray.append( taxonarray )
    return np.array( contigarray ), taxaorder


def calc_onebug( array, bugindex ):
    """
    Calculates the one bug score. Equals the maximum of the minimum of scores across all bugs per gene.
    High scores indicate that contig is likely explained by one taxa.
    """
    onebugscore = np.max( np.amin( array, axis=1 ) )
    indices = [i for i,j in enumerate( list(np.amin(array, axis=1) ) ) if j==onebugscore]
    bugmeanlist = []
    for index in indices:
        bugname = bugindex[index]
        bugmean = np.mean(array[index])
        bugmeanlist.append( [bugname, bugmean] )
    bugmeanlist_sorted = sorted( bugmeanlist, key=itemgetter(1), reverse=True )
    topmean = bugmeanlist_sorted[0][1]
    onebuglist = []
    for bug, avg in bugmeanlist_sorted:
        if avg >= topmean:
            onebuglist.append( bug )
    return onebugscore, onebuglist


def find_taxa_order( contigarray, taxaorder, orgpair ):
    """
    Determines the order of taxa by gene number, aka ['TaxaA', 'TaxaB', 'TaxaA'] for [Gene1, Gene2, Gene3].
    """
    top_one, top_two = orgpair.split('-')[0], orgpair.split('-')[1]
    index_one, index_two = taxaorder.index( top_one ), taxaorder.index( top_two )
    array_one, array_two = contigarray[index_one, :], contigarray[index_two, :]
    twobugmax = np.amax( np.array( [array_one, array_two] ), axis=0 )
    element_list, taxastring, shift_counter = [], '', 0
    for i in range( len( twobugmax ) ):
        element = twobugmax[i]
        one_element, two_element = array_one[i], array_two[i]
        if element == one_element:
            element_list.append( top_one )
            taxastring += 'A'
            state = top_one
        else:
            element_list.append( top_two )
            taxastring += 'B'
            state = top_two
        if i != 0:
            if element_list[i-1] != state:
                shift_counter += 1
    return element_list, taxastring, shift_counter


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
            average = np.mean(np.amax(twobugarray, axis=0))
            #indices = np.argmax(twobugarray, axis=0)
            bugsort = sorted( [str(bugindex[i]), str(bugindex[j])] )
            bugnames = bugsort[0] + '-' + bugsort[1]
            complementlist.append( [ bugnames, complement, average ] )
    complement_sorted = sorted( complementlist, key=itemgetter(1,2), reverse=True )

    twobugscore, twobugavg = complement_sorted[0][1], complement_sorted[0][2]
    twobuglist = []
    for info in complement_sorted:
        element_list, taxastring, shift = find_taxa_order( array, bugindex, info[0] )
        if info[1] == twobugscore and info[2] >= twobugavg and len(set(element_list)) > 1:
            twobuglist.append( info[0] )
        elif info[1] == twobugscore and info[2] >= twobugavg and len(set(element_list)) <= 1:
            twobugscore = 0
            twobuglist = ['None-None']
            break
        elif info[1] < twobugscore or (info[1] == twobugscore and info[2] < twobugavg):
            break
    return twobugscore, twobuglist


def calc_donorrecip( contigarray, taxaorder, taxa ):
    """
    Determines donor and recipient if there are 2 transitions, and the outer genes are from the same bug. (For eg, "ABA").
    The bug in between (eg. "B") is the donor. The bug that explains the outside genes (eg. "A") is the recipient.
    If there is only 1 transition, we call it "unknown". If there is more than 1 transition, we call it a "hybrid."
    """
    ans = []
    for orgs in taxa:
        recipient, donor, status = '', '', ''
        top_one, top_two = orgs.split('-')[0], orgs.split('-')[1]
        element_list, taxastring, shift_counter = find_taxa_order( contigarray, taxaorder, orgs )
        if shift_counter == 2 and element_list[0] == element_list[len( element_list ) - 1]:
            recipient = element_list[0]
            if recipient == top_one:
                donor = top_two
            else:
                donor = top_one
            status = recipient + '-' + donor
        elif shift_counter > 2:
            status = 'NA' 
        elif shift_counter == 1:
            status = 'NA'
        ans.append( status )
    return ans
    

def rank_ambiguous( onebugscore, onebuglist, twobugscore, twobuglist ):
    """
    If the contig cannot be called "noLGT" or "LGT", decide which one it is closer to and call ambiguous.
    """
    if onebugscore >= twobugscore:
        score = onebugscore
        org = onebuglist
        status = "ambiguous-NoLGT"
    else:
        score = twobugscore
        org = twobuglist
        status = "ambiguous-LGT"
    return status, score, org


def find_uniref( status, contigarray, taxaorder, taxa, taxalist ):
    taxastringlist, uniref50list, uniref90list = [], [], [] #single list of taxastrings per org, single list of Uniref by gene (for all orgs is same)
    if status == 'LGT' or status == 'ambiguous-LGT':
        for orgs in taxa:
            element_list, taxastring, shift_counter = find_taxa_order( contigarray, taxaorder, orgs ) 
            uniref50list = [ 'NoUniref' ]*len( element_list )
            uniref90list = [ 'NoUniref' ]*len( element_list )
            for i in range(len(element_list)):
                for info in taxalist:
                    if info.gene == i+1 and info.taxa == element_list[i] and uniref50list[i] == 'NoUniref':
                        uniref50list[i] = info.uniref50
                        uniref90list[i] = info.uniref90
                    else:
                        continue
            taxastringlist.append(taxastring)
    else:
        for orgs in taxa:
            numbugs, numgenes = contigarray.shape
            uniref50list, uniref90list = [ 'NoUniref' ]*numgenes, [ 'NoUniref' ]*numgenes
            taxastring = ''
            for i in range(numgenes):
                taxastring += 'A'
                for info in taxalist:
                    if info.gene == i+1 and info.taxa == orgs and uniref50list[i] == 'NoUniref':
                        uniref50list[i] = info.uniref50
                        uniref90list[i] = info.uniref90
                    else:
                        continue
            taxastringlist.append(taxastring)
    return taxastringlist, uniref50list, uniref90list


def print_result( contig, length, status, score, orgs, rdstatus_list, taxaannot, uniref50, uniref90 ):
    """
    Print out the results in a neat list.
    """
    orglist = []
    if len(orgs) > 1:
        if len(orgs) > 3:
            orglist = orgs[0:3]
            orglist.append( str(len( orgs) - 1 - 3 ) + ' other taxa with same score' )
            orgs = orglist
    finalorgs = ';'.join( orgs )
    rdstatus = ';'.join(rdstatus_list)
    finaltaxaannot = ';'.join( taxaannot )
    finaluniref50 = '|'.join(uniref50)
    finaluniref90 = '|'.join(uniref90)
    return [str(x) for x in [contig, length, status, score, finalorgs, rdstatus, finaltaxaannot, finaluniref50, finaluniref90] ]
    

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
        writer.writerow( ['contig', 'contiglen', 'status', 'score', 'taxa', 'recipient-donor', 'taxaorder', 'uniref50', 'uniref90'] )
        for contig, taxalist in wu.iter_contig_taxa( args.input ):
            contiglen = taxalist[0].length
            contigarray, taxaorder = generate_tables( taxalist )
            recipient, donor = "NA", "NA"
           
            spikearray, spikeorder = spike_unknown( contigarray, taxaorder, args.unknown )
            onebugscore, onebuglist = calc_onebug( spikearray, spikeorder )
            numbugs, numgenes = spikearray.shape
            rdstatus_list = ['NA']
                    
            if onebugscore >= args.onebug or numgenes == 1 or numbugs == 1:
                score, taxa = onebugscore, onebuglist
                status = "NoLGT"
            else:
                twobugscore, twobuglist = calc_twobug( spikearray, spikeorder )
                if twobugscore >= args.twobug:
                    score, taxa = twobugscore, twobuglist
                    rdstatus_list = calc_donorrecip( spikearray, spikeorder, twobuglist )
                    status = "LGT"
                else:
                    status, score, taxa = rank_ambiguous( onebugscore, onebuglist, twobugscore, twobuglist )
                    if status == 'ambiguous-LGT':
                        rdstatus_list = calc_donorrecip( spikearray, spikeorder, twobuglist )
            taxaannot, uniref50, uniref90 = find_uniref( status, spikearray, spikeorder, taxa, taxalist ) 
            result_line = print_result( contig, contiglen, status, score, taxa, rdstatus_list, taxaannot, uniref50, uniref90 )
            writer.writerow( result_line )
            
            
if __name__ == "__main__":
    main()
