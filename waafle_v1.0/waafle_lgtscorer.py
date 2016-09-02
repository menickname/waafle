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
import os, sys, csv, argparse, re
from operator import itemgetter, attrgetter, methodcaller
import waafle_utils as wu
import numpy as np
import numpy.ma as ma

# ---------------------------------------------------------------
# constants
# ---------------------------------------------------------------

"""
@codereview 9/2/2016
Outsource this to the utils file
"""

c__list_taxa = ["k", "p", "c", "o", "f", "g", "s"]

# ---------------------------------------------------------------
# functions
# ---------------------------------------------------------------

"""
@codereview 9/2/2016
"-s1", "--onebug",
"-s2", "--twobug",
Usually the flag with a single dash is just one character
"""

"""
@codereview 9/2/2016
Use the action="store_true" option for --unknown
"""

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


def find_taxa_order( status, contigarray, taxaorder, org ):
    """
    Returns order of taxa by gene number:
    1) In a list, aka ['TaxaA', 'TaxaB', 'TaxaA'] for [Gene1, Gene2, Gene3].
    2) As a string, aka 'ABA'
    3) As well as number of shifts from A-->B or B-->A; C<-->A|B does not count.
    """
    element_list, taxastring, shift_counter = [], '', 0
    if status == 'LGT' or status == 'ambiguous-LGT':
        top_one, top_two = re.split('[><\?-]', org )
        index_one, index_two = taxaorder.index( top_one ), taxaorder.index( top_two )
        array_one, array_two = contigarray[index_one, :], contigarray[index_two, :]
        twobugmax = np.amax( np.array( [array_one, array_two] ), axis=0 )
        for i in range( len( twobugmax ) ):
            element = twobugmax[i]
            one_element, two_element = array_one[i], array_two[i]
            if element == one_element and element == two_element:
                element_list.append( [top_one, top_two] )
                taxastring += 'C'
            elif element == one_element and element != two_element:
                element_list.append( top_one )
                taxastring += 'A'
            else:
                element_list.append( top_two )
                taxastring += 'B'
        reduced_taxastring = taxastring.replace( 'C', '')
        for j in range( len( reduced_taxastring ) -1 ):
            if reduced_taxastring[j] != reduced_taxastring[j+1]:
                shift_counter += 1
    else:
        numbugs, numgenes = contigarray.shape
        element_list, taxastring = [org]*numgenes, 'A'*numgenes
    return element_list, taxastring, shift_counter


def calc_twobug( array, bugindex, onebugscore ):
    """
    Calculates the complement score. Equals the minimum of the maximum of scores across all genes per pair of bugs.
    High scores indicate that contig is likely explained by two taxa.
    """
    complementlist = []
    for i in range( array.shape[0]-1 ):
        bug1 = array[i, :]
        if max( bug1 ) < onebugscore:
            continue
        for j in range( i+1, array.shape[0] ):
            bug2 = array[j, :]
            if max( bug2 ) < onebugscore:
                continue
            twobugarray = np.array( [bug1, bug2] )
            complement = np.min(np.amax(twobugarray, axis=0))
            average = np.mean(np.amax(twobugarray, axis=0))
            bugsort = sorted( [str(bugindex[i]), str(bugindex[j])] )
            bugnames = bugsort[0] + '-' + bugsort[1]
            complementlist.append( [ bugnames, complement, average ] )
    if len( complementlist ) == 0:
        twobugscore, twobuglist = 0, ['None-None']
    else:
        complement_sorted = sorted( complementlist, key=itemgetter(1,2), reverse=True )
        twobugscore, twobugavg = complement_sorted[0][1], complement_sorted[0][2]
        twobuglist = []
        for info in complement_sorted:
            element_list, taxastring, shift = find_taxa_order( 'LGT', array, bugindex, info[0] )
            lgt_status = False
            reduced_taxastring = taxastring.replace( 'C', '')
            if len( set( reduced_taxastring ) ) == 2:
                lgt_status = True
            if info[1] == twobugscore and info[2] >= twobugavg and lgt_status == True:
                twobuglist.append( info[0] )
            elif info[1] == twobugscore and info[2] >= twobugavg and lgt_status == False:
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
    annottaxa = []
    for orgs in taxa:
        newtaxa = ''
        top_one, top_two = orgs.split('-')[0], orgs.split('-')[1]
        element_list, taxastring, shift_counter = find_taxa_order( 'LGT', contigarray, taxaorder, orgs )
        reduced_taxastring = taxastring.replace( 'C', '' )
        if shift_counter == 2 and reduced_taxastring[0] == reduced_taxastring[-1]:
            recip_index = taxastring.index( reduced_taxastring[0] )
            recipient = element_list[ recip_index ]
            if recipient == top_one:
                newtaxa = top_one + '<' + top_two
            else:
                newtaxa = top_one + '>' + top_two
        else:
            newtaxa = top_one + '?' + top_two 
        annottaxa.append( newtaxa )
    return annottaxa
    

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
    numbugs, numgenes = contigarray.shape
    uniref50list, uniref90list = ['unknown']*numgenes, ['unknown']*numgenes
    taxaset = set([])
    for org in taxa:
        if status == 'LGT' or status == 'ambiguous-LGT':
            taxaone, taxatwo = re.split('[?><]', org )
            taxaset.add( taxaone )
            taxaset.add( taxatwo )
        else:
            taxaset.add( org )
    for i in range( numgenes ):
        for info in taxalist:
            if info.gene == i+1 and info.taxa in taxaset:
                uniref50list[i] = info.uniref50.replace(',', '|')
                uniref90list[i] = info.uniref90.replace(',', '|')
                break
    return uniref50list, uniref90list


def print_result( contig, length, status, onescore, twoscore, orgs, taxaannot, uniref50, uniref90 ):
    """
    Print out the results in a neat list.
    """
    finalorgs = '|'.join( [org.replace('|', '.') for org in orgs] ) #can't join with | because full names
    finaltaxaannot = '|'.join( taxaannot ) #should follow orgs
    finaluniref50 = ';'.join(uniref50)
    finaluniref90 = ';'.join(uniref90)
    return [str(x) for x in [contig, length, status, onescore, twoscore, finalorgs, finaltaxaannot, finaluniref50, finaluniref90] ]
    

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
        writer.writerow( ['contig', 'contiglen', 'status', 'onebug_score', 'twobug_score', 'taxa', 'recipient-donor', 'taxaorder', 'uniref50', 'uniref90'] )
        for contig, taxalist in wu.iter_contig_taxa( args.input ):
            contiglen = taxalist[0].length

            #Generate arrays
            contigarray, taxaorder = generate_tables( taxalist )
            spikearray, spikeorder = spike_unknown( contigarray, taxaorder, args.unknown )

            #Calculate onebug and twobug scores
            onebugscore, onebuglist = calc_onebug( spikearray, spikeorder )
            numbugs, numgenes = spikearray.shape
            if numgenes == 1 or numbugs == 1:
                twobugscore, twobuglist = 0, ["NA"]
            else:
                twobugscore, twobuglist = calc_twobug( spikearray, spikeorder, onebugscore )
                    
            #Apply thresholds for one/twobug scores
            if onebugscore >= args.onebug: 
                score, taxa = onebugscore, onebuglist
                status = "NoLGT"
            else:
                if twobugscore >= args.twobug:
                    score, taxa = twobugscore, twobuglist
                    taxa = calc_donorrecip( spikearray, spikeorder, twobuglist ) #assign d/r if possible
                    status = "LGT"
                else:
                    status, score, taxa = rank_ambiguous( onebugscore, onebuglist, twobugscore, twobuglist )
                    if status == 'ambiguous-LGT':
                        taxa = calc_donorrecip( spikearray, spikeorder, twobuglist )

            #Get synteny for all orgs
            taxaannot = []
            for org in taxa:
                element_list, taxastring, shift_counter = find_taxa_order( status, spikearray, spikeorder, org )
                taxaannot.append( taxastring )

            #Get unirefs for contig
            uniref50, uniref90 = find_uniref( status, spikearray, spikeorder, taxa, taxalist ) 
            result_line = print_result( contig, contiglen, status, onebugscore, twobugscore, taxa, taxaannot, uniref50, uniref90 )
            writer.writerow( result_line )
            
            
if __name__ == "__main__":
    main()
