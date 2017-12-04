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
import os, sys, csv, argparse, re
import waafle_utils as wu
import numpy as np
from collections import Counter
from operator import itemgetter, attrgetter, methodcaller

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
        "-p", "--plotfile",
        default="waafle-scoredorgs.tsv",
        help="output for scored taxa",
        )
    parser.add_argument(
        "-o", "--results",
        default="waafle-scoredcontigs.tsv",
        help="output for contig annotations",
        )
    parser.add_argument(
        "-lap", "--overlap-hits",
        help="amount of overlap to include hit in a gene",
        default=0.1,
        type=float,
        )
    parser.add_argument(
        "-scov", "--scov-hits",
        help="cutoff for gene coverage or subject coverage when grouping hits",
        default=0,
        type=float,
        )
    parser.add_argument(
        "-s", "--strand",
        action="store_true",
        help="turn on strand specific gene calling",
        )
    parser.add_argument(
        "-u", "--unknown",
        action="store_true",
        help="include unknown taxa",
        )
    parser.add_argument(
        "-s1", "--onebug",
        default=0.5,
        type=float,
        help="onebug score cutoff"
        )
    parser.add_argument(
        "-s2", "--twobug",
        default=0.8,
        type=float,
        help="twobug score cutoff"
        )
    parser.add_argument(
        "-r", "--custrange",
        default=0.1,
        type=float,
        help="Range of top bugs to consider",
        )
    parser.add_argument(
        "-a", "--adjust",
        action='store_true',
        help="Adjust k2 or not (default not)",
        )
    args = parser.parse_args()
    return args

def hits2genes( gene, hits, strand_specific, lap, scov ):
    genehits = []
    taxaset = set()
    uniref50, uniref90 = [], []
    for hit in hits:
        hitstrand = wu.convert_strand( hit.sstrand )
        if ( hitstrand == gene.strand or not strand_specific ) and hit.scov_modified > scov:
            overlap = wu.calc_overlap( hit.qstart, hit.qend, gene.start, gene.end )    
            if overlap > lap:
                genehits.append( hit )
                taxaset.add( '|'.join( hit.taxonomy ) )
                uniref50.append( hit.uniref50 )
                uniref90.append( hit.uniref90 )
    uniref50_c = ','.join( [ str(x) + ':' + str(y) for x, y in Counter( uniref50 ).most_common( 3 )] )
    uniref90_c = ','.join( [ str(x) + ':' + str(y) for x, y in Counter( uniref90 ).most_common( 3 )] )
    info = [genehits, taxaset, uniref50_c, uniref90_c]
    return info


def score_taxa( gene, info, contiglen ):
    """
    Group hits that correspond to specific taxa.
    If the gene did not have hits to begin with, output an "unknown" taxa with score 0.
    If the gene has hits that corresponds to orgs, score the orgs and output the taxa annotations.
    """
    genehits, taxaset, uniref50_c, uniref90_c = info
    genelen = gene.end - gene.start + 1
    dict_taxaarray, taxalist = {}, []
    if len( genehits ) == 0:
        #set unknown taxa or output
        taxaname = []
        for level in wu.c__list_taxa:
            taxa = level + '__Unknown'
            taxaname.append( taxa )
        taxa = '|'.join( taxaname )
        start_list, end_list, score = gene.start, gene.end, 1
        dict_taxaarray[ taxa ] = np.ones( genelen ) 
        taxa = wu.Taxa( [
                    gene.seqname, 
                    contiglen, 
                    int( gene.attributes['ID'].split('_')[-1] ), 
                    gene.strand, 
                    gene.start, 
                    gene.end, 
                    taxa, 
                    start_list, 
                    end_list, 
                    score, 
                    'UniRef50_unknown:0', 
                    'UniRef90_unknown:0', 
                    0,
                    ] )
        taxalist.append( taxa )
        #dict_taxaarray[ taxa ] = np.ones( genelen )
    else:
        for taxa in taxaset:
            dict_taxaarray.setdefault( taxa, np.zeros( genelen ) )
            numhits = 0
            for hit in genehits:
                if '|'.join( hit.taxonomy ) == taxa:
                    numhits += 1
                    orgarray = dict_taxaarray.get( taxa )
                    hitarray = np.zeros( genelen )
                    arraystart = max( hit.qstart - gene.start, 0 )
                    arrayend = min( hit.qend - gene.start + 1, genelen )
                    hitarray[arraystart: arrayend] = (hit.pident/float(100))*hit.scov_modified
                    dict_taxaarray[taxa] = np.maximum( hitarray, orgarray )
            startstop = np.split(np.argwhere( dict_taxaarray[taxa] ), np.where(np.diff(np.argwhere( dict_taxaarray[taxa] ), axis=0)!= 1)[0]+1)
            start_list = ','.join( [str(element[0][0] + gene.start) for element in startstop] )
            end_list = ','.join( [str(element[-1][0] + gene.start) for element in startstop] ) #check this
            score = np.mean( dict_taxaarray[taxa] )
            taxa = wu.Taxa( [
                    gene.seqname, 
                    contiglen,
                    int( gene.attributes['ID'].split('_')[-1] ),
                    gene.strand, 
                    gene.start, 
                    gene.end, 
                    taxa, 
                    start_list,
                    end_list, 
                    score, 
                    uniref50_c, 
                    uniref90_c, 
                    numhits,
                    ] )
            taxalist.append( taxa )
    return dict_taxaarray, taxalist


def generate_matrix( dict_genearray, taxalevel ):
    """ Generate a table in which rows are taxa and columns are genes.
        Start by aggregating score arrays by taxa (at some taxalevel).
        Calculate the new score and return table and its taxa. """
    # get taxa across all genes
    taxaindex = wu.c__list_taxa.index( taxalevel )
    complete_taxaset = set([])
    dict_genearray_revised = {}
    for gene in sorted( dict_genearray.keys() ):

        # aggregate arrays for taxa level
        dict_taxaarray = dict_genearray[gene]
        dict_taxaarray_revised = {}
        for taxa in dict_taxaarray.keys():
            taxaname = '|'.join( taxa.split('|')[:taxaindex + 1] )
            dict_taxaarray_revised.setdefault( taxaname, [] ).append( dict_taxaarray[ taxa ] )
            complete_taxaset.add( taxaname )
        # calculate scores per taxa
        dict_taxascores = {}
        for taxa in dict_taxaarray_revised.keys():
            score = np.mean( np.amax( np.array( dict_taxaarray_revised[ taxa ] ), axis=0 ) )
            dict_taxascores[taxa] = score
        dict_genearray_revised[ gene ] = dict_taxascores

    # generate matrix
    table, taxaorder = [], []
    for taxa in complete_taxaset:
        line = []
        taxaorder.append( taxa )
        for gene in sorted( dict_genearray.keys() ):
            dict_taxascores = dict_genearray_revised[gene]
            if taxa in dict_taxascores.keys():
                line.append( dict_taxascores[taxa] )
            else:
                line.append( 0 )
        table.append( line )
    return np.array( table ), taxaorder


def spike_unknown( table, taxaorder, unknown ):
    """ Spike in the unknown organism.
    Assign the unknown values if it has already been called,
    or generate a new unknown. """
    final_table = np.copy( table )
    final_taxaorder = [item for item in taxaorder ]
    numorgs, numgenes = final_table.shape
    taxaindex = wu.c__list_taxa.index( taxaorder[0].split('|')[-1].split('__')[0] )
    taxaname = []
    for i in range( len( wu.c__list_taxa ) ):
        if i <= taxaindex:
            level = wu.c__list_taxa[i]
            taxa = level + '__Unknown'
            taxaname.append( taxa )
    unknown_name = '|'.join( taxaname )
    
    # check if unknown needs to be spiked
    if numgenes == 0:
        return final_table, final_taxaorder
    if unknown_name in taxaorder:
        unknown_index = taxaorder.index( unknown_name )
        for i in range( numgenes ):
            unknown_score = 1 - np.max( table[:, i] )
            if table[unknown_index, i] == 1:
                unknown_score = 1
            final_table[ unknown_index, i ] = unknown_score
        return final_table, final_taxaorder
    elif not unknown:
        return final_table, final_taxaorder
    else:
        unknown = []
        for i in range( numgenes ):
            unknown_score = 1 - np.max( table[:, i] )
            unknown.append( unknown_score )
        final_table = np.append( final_table, [unknown], axis=0 )
        final_taxaorder.append( unknown_name )
        return final_table, final_taxaorder


def find_taxa_order( status, contigarray, taxaorder, org, threshold ):
    """
    For a NoLGT taxon, returns taxon and synteny. (no Cs)
    For a LGT pair, returns order of taxa by gene number:
        1) In a list, aka ['TaxaA', 'TaxaB', 'TaxaA'] for [Gene1, Gene2, Gene3].
        2) As a string, aka 'ABA'
        3) As well as number of shifts from A-->B or B-->A; C<-->A|B does not count.
    """
    element_list, taxastring, shift_counter = [], '', 0
    if status == 'LGT':
        top_one, top_two = re.split('[><\?-]', org )
        index_one, index_two = taxaorder.index( top_one ), taxaorder.index( top_two )
        array_one, array_two = contigarray[index_one, :], contigarray[index_two, :]
        twobugmax = np.amax( np.array( [array_one, array_two] ), axis=0 )
        for i in range( len( twobugmax ) ):
            # element = twobugmax[i]
            element = threshold
            one_element, two_element = array_one[i], array_two[i]
            if element <= one_element and element <= two_element:
                element_list.append( [top_one, top_two] )
                taxastring += 'C'
            elif element <= one_element and element > two_element:
                element_list.append( top_one )
                taxastring += 'A'
            elif element > one_element and element <= two_element:
                element_list.append( top_two )
                taxastring += 'B'
            else:
                pass #this should not be possible, taxastring += '@'
        reduced_taxastring = taxastring.replace( 'C', '')
        for j in range( len( reduced_taxastring ) -1 ):
            if reduced_taxastring[j] != reduced_taxastring[j+1]:
                shift_counter += 1
    else:
        numbugs, numgenes = contigarray.shape
        element_list, taxastring = [org]*numgenes, 'A'*numgenes
    return element_list, taxastring, shift_counter


def calc_donorrecip( orgs, element_list, taxastring, shiftcounter ):
    """
    Determines donor and recipient if there are 2 transitions, and the outer genes are from the same bug. (For eg, "ABA").
    The bug in between (eg. "B") is the donor. The bug that explains the outside genes (eg. "A") is the recipient.
    If there is only 1 transition, we call it "unknown". If there is more than 1 transition, we call it a "hybrid."
    """
    newtaxa = ''
    top_one, top_two = orgs.split('-')[0], orgs.split('-')[1]
    reduced_taxastring = taxastring.replace( 'C', '' )
    if shiftcounter == 2 and reduced_taxastring[0] == reduced_taxastring[-1]:
        recip_index = taxastring.index( reduced_taxastring[0] )
        recipient = element_list[ recip_index ]
        if recipient == top_one:
            newtaxa = top_one + '<' + top_two
        else:
            newtaxa = top_one + '>' + top_two
    else:
        newtaxa = top_one + '?' + top_two
    return newtaxa


def calc_onebug( table, buglist, onebugthresh, myrange ):
    """ Calculates the 1 bug scores by:
    1. Get top onebug score.
    2. Get all taxa within x% of top onebug score.
    3. Get all taxa within x% of top average values across genes. <-- remove
    4. Get all taxa above the onebug threshold.  """
    onebugscore = np.max( np.amin( table, axis=1 ) )
    nolgt_status = False
    # Get all taxa within x% of top onebugscore
    onebugscoremin = onebugscore - (float(onebugscore)*myrange)
    indices = [i for i,j in enumerate( list(np.amin(table, axis=1) ) ) if j>=onebugscoremin]
    bugmeanlist = []
    for index in indices:
        bugname = buglist[index]
        bugmean = np.mean(table[index])
        bugscore = np.amin(table[index])
        bugmeanlist.append( [bugname, bugmean, bugscore ] )
    # Find top average
    bugmeanlist_sorted = sorted( bugmeanlist, key=itemgetter(1), reverse=True )
    topmean = bugmeanlist_sorted[0][1]
    topmeanmin = topmean - (float(topmean)*myrange)
    # Filter for all taxa within x% of top onebugavg
    onebuginfo = []
    for bug, avg, bugscore in bugmeanlist:
        if avg >= topmeanmin:
            element_list, taxastring, shift_counter = find_taxa_order( 'NoLGT', table, buglist, bug, onebugscoremin )
            onebuginfo.append( [bug, bugscore, avg, taxastring] )
    # Filter by onebug score
    onebuginfo_filt = []
    for bug, bugscore, avg, taxastring in onebuginfo:
        if bugscore >= onebugthresh:
            onebuginfo_filt.append( [bug, bugscore, avg, taxastring] )
    if len( onebuginfo_filt ) > 0:
        nolgt_status = True
    return nolgt_status, onebuginfo, onebuginfo_filt


def filter_twobug( complement_sorted, table, buglist, myrange ):
    """ Filters complement scores by:
    1. Get top twobug score and top average.
    2. Get all taxon pairs within x% of top twobug score and x% of top average values.
    3. Get all taxon pairs with k3 >= y.  """
    twobugscore, twobugavg = complement_sorted[0][1], complement_sorted[0][2] #top scores
    twobugscore_min = twobugscore - (float(twobugscore)*myrange)
    twobugavg_min = twobugavg - (float(twobugavg)*myrange)
    # For twobug scores within range, classify as LGT or not based on synteny
    lgt_status, possibilities, non_possibilities = False, [], []
    for info in complement_sorted:
        #element_list, taxastring, shift = find_taxa_order( 'LGT', table, buglist, info[0], twobugscore_min )
        #reduced_taxastring = taxastring.replace( 'C', '')
        if info[1] >= twobugscore_min and info[2] >= twobugavg_min: #see if entry makes threshold
            element_list, taxastring, shift = find_taxa_order( 'LGT', table, buglist, info[0], twobugscore_min )
            reduced_taxastring = taxastring.replace( 'C', '')    
            if len( set( reduced_taxastring ) ) == 2: #see if two bugs together made threshold
                newname = calc_donorrecip( info[0], element_list, taxastring, shift )
                possibilities.append( info + [newname, taxastring] )
            else: #1 bug explains top scores (eg. Bug1 [0.8, 0.9, 0.8], Bug2 [0.8, 0.7, 0.8])
                non_possibilities.append( info + ['-', '-'] )
        if info[1] < twobugscore_min or (info[1] == twobugscore_min and info[2] < twobugavg_min):
            break
    # Determine whether there is possibility of LGT
    if len( non_possibilities ) >= 1: #Assume NoLGT, some "pair" above min was really onebug
        lgt_status = False
        return lgt_status, non_possibilities
    else: # Assume potential LGT, all pairs above min were two bugs
        lgt_status = True
        return lgt_status, possibilities
            

def calc_twobug( table, buglist, twobugthresh, myrange, adjust ):
    """ Calculates multiple k scores and store.
    1. For every pair, calculate k2, which equals the minimum of 
    the maximum of scores across all genes per pair of bugs.
    High scores indicate that contig is likely explained by two taxa.
    2. Store only pairs with k2 above threshold. """
    complementlist = []
    for i in range( table.shape[0]-1 ):
        bug1 = table[i, :]
        if max( bug1 ) < twobugthresh: #will never be able to form a pair
            #print( max( bug1 ), twobugthresh, buglist[i] )
            continue
        else:
            for j in range( i+1, table.shape[0] ):
                bug2 = table[j, :]
                if max( bug2 ) < twobugthresh: #will never be able to form a pair
                    continue
                else:
                    twobugarray = np.array( [bug1, bug2] )
                    # k2 score adjustments for testing
                    k2_modmean = np.mean( np.absolute( twobugarray[0] - twobugarray[1] ), axis=0 )
                    k2_modmed = np.median( np.absolute( twobugarray[0] - twobugarray[1] ), axis=0 )
                    # original k2 score is complement
                    complement = np.min(np.amax(twobugarray, axis=0))
                    if adjust: #adjust k2 score by penalty
                        complement = max( complement - max( 1 - k2_modmean, 0 ), float(0) )
                    average = np.mean(np.amax(twobugarray, axis=0))
                    # get bug names
                    bugsort = sorted( [str(buglist[i]), str(buglist[j])] )
                    bugnames = bugsort[0] + '-' + bugsort[1]
                    # place all into list after filtering for k2
                    if complement >= twobugthresh:
                        complementlist.append( [ bugnames, complement, average, k2_modmean, k2_modmed ] )
    # for all pairs, filter results
    lgt_status, list_of_possibilities = False, [] #assume all taxa pairs filtered out, no LGT possible
    if len( complementlist ) != 0: #not all pairs filtered out, find possible pairs
        complement_sorted = sorted( complementlist, key=itemgetter(1,2,3), reverse=True )
        lgt_status, list_of_possibilities = filter_twobug( complement_sorted, table, buglist, myrange )
    return lgt_status, list_of_possibilities


def flip_annot( sign, synteny ):
    newsign = sign
    if sign == '>':
        newsign = '<'
    if sign == '<':
        newsign = '>'
    syn1 = synteny.replace( 'A', 'D' )
    syn2 = syn1.replace( 'B', 'A' )
    newsynteny = syn2.replace( 'D', 'B' )
    return newsign, newsynteny


def reduce_set( myset, highestlevel ):
    newset = set([])
    i = 0
    for taxon in myset:
        taxonlist = taxon.split('|')
        if i == 0:
            newset |= set( taxonlist )
        else:
            newset &= set( taxonlist )
        i += 1
    taxaname, multiple_count = [], 0
    for level in wu.c__list_taxa:
        r = re.compile( level + '__*' )
        nextname = filter( r.match, list( newset ) )
        if len( nextname ) == 1:
            taxaname.append( nextname[0] )
        else:
            #taxaname.append( level + '__multiple' )
            multiple_count += 1
        if level == highestlevel:
            break
    finalname = '|'.join( taxaname )
    return finalname, multiple_count


def get_new_taxastring( first_taxon_sorted, second_taxon_sorted, oldtaxapair, oldtaxastring ):
    oldtaxaone, oldtaxatwo = re.split( '[><\?]', oldtaxapair )
    oldtaxasign = re.search( '[><\?]', oldtaxapair ).group()
    newtaxastring = ''
    if oldtaxaone == first_taxon_sorted:
        newtaxastring = oldtaxastring
        newtaxasign = oldtaxasign
    else:
        newtaxasign, newtaxastring = flip_annot( oldtaxasign, oldtaxastring )
    newname = first_taxon_sorted + newtaxasign + second_taxon_sorted
    return newname, newtaxastring
    

def match_synteny( synteny_check ):
    consensus_synteny = ''
    consensus_vals = []
    for i in range( len( synteny_check[0] ) ): #go through each letter
        position_vals = set( [] )
        for j in range( len( synteny_check ) ): #go through each pair
            if synteny_check[j][i] != 'C':
                position_vals.add( synteny_check[j][i] )
        if len( list( position_vals ) ) == 1:
            consensus_synteny += list( position_vals )[0]
        elif len( list( position_vals ) ) == 0:
            consensus_synteny += 'C'
        else:
            consensus_synteny += 'C' # just make them C, /'.join( list( position_vals ) )
        consensus_vals.append( len( list( position_vals ) ) )
    consensus_value = max( consensus_vals )
    return consensus_synteny, consensus_value


def det_multiple_two( list_of_possibilities ):
    """ Determine if we can annotate LGT at a lower level (do we have to walk up) """
    # Among potential pairs, see if all pairs share a common partner
    firsttaxon_set, secondtaxon_set = set([]), set([])
    array_values, synteny_check = [], []
    for info in list_of_possibilities:
        sortedname, k2, k2avg, k3avg, k3med, newname, taxastring = info
        array_values.append( [k2, k2avg, k3avg, k3med] )
        firsttaxon = ''
        for letter in taxastring:
            if letter == 'B':
                firsttaxon = 'B'
                break
            elif letter == 'A':
                firsttaxon = 'A'
                break
            else:
                continue
        taxaone, taxatwo = '', ''
        if firsttaxon == 'A':
            taxaone, taxatwo = re.split( '[><\?]', newname )
            synteny_check.append( taxastring )
        else:
            taxatwo, taxaone = re.split( '[><\?]', newname )
            new_sign, new_taxastring = flip_annot( '<', taxastring ) #ignore new_sign
            synteny_check.append( new_taxastring )
        firsttaxon_set.add( taxaone )
        secondtaxon_set.add( taxatwo )
    # Check synteny matches
    consensus_synteny, consensus_value = match_synteny( synteny_check )
    # Get new k2 and k3 values
    k2, k2avg, k3avg, k3med = np.average( np.array( array_values ), axis=0 )
    # If all pairs share a common partner, determine what level the other partner can be annotated to
    newcombined, newcombinedstring = 'Multiple', 'C'*len( taxastring )
    high_lvl = False
    first_taxon, second_taxon, multiple_count = '', '', 0
    if len( firsttaxon_set ) == 1 and len( secondtaxon_set ) > 1 and consensus_value == 1:
        first_taxon = list( firsttaxon_set )[0]
        highest_level = first_taxon.split('|')[-1].split('__')[0]
        second_taxon, multiple_count = reduce_set( secondtaxon_set, highest_level )
        last_known_lvl = second_taxon.split('|')[-1]
        if not re.search( last_known_lvl, first_taxon ):
            high_lvl = True
    elif len( firsttaxon_set ) > 1 and len( secondtaxon_set ) == 1 and consensus_value == 1:
        second_taxon = list( secondtaxon_set )[0]
        highest_level = second_taxon.split('|')[-1].split('__')[0]
        first_taxon, multiple_count = reduce_set( firsttaxon_set, highest_level )
        last_known_lvl = first_taxon.split('|')[-1]
        if not re.search( last_known_lvl, second_taxon ):
            high_lvl = True
    # If only 1 level of the other partner is annotated as "multiple", then get the annotation
    # Otherwise, return "Unknown" and continue walking up the taxonomic levels
    # Reasoning: interesting to know g__TaxaA|s__TaxaB exchanges with g__TaxaA|s__Multiple, but not interesting to know g__TaxaA|s__TaxaB exchanges with Multiple.
    if multiple_count == 1 and high_lvl:
        first_taxon_sorted = sorted( [first_taxon, second_taxon] )[0]
        second_taxon_sorted = sorted( [first_taxon, second_taxon] )[1]
        newcombined, newcombinedstring = get_new_taxastring( first_taxon_sorted, second_taxon_sorted, list_of_possibilities[0][-2], list_of_possibilities[0][-1] )
    return newcombined, newcombinedstring, k2, k2avg, k3avg, k3med
    

def annotate_onebug( onebuglist, onebuglist_filt ):
    """ Annotate the onebug with its taxon """
    k1, onebug, onesyn = 0, '-', '-'
    # multiple taxa passed the filter, call "Multiple" and average one bug score
    if len( onebuglist_filt ) > 1:
        k1 = np.average( [float( x[1] ) for x in onebuglist_filt] )
        k1_avg = np.average( [float( x[2] ) for x in onebuglist_filt] )
        onebug, onesyn = 'Multiple', 'C'*len( onebuglist_filt[0][3] )
    # one bug passed the filter, annotate the taxon
    elif len( onebuglist_filt ) == 1:
        onebug, k1, k1_avg, onesyn = onebuglist_filt[0]
    # nothing passed the filter
    else:
        # multiple taxa did not pass filter, call "Multiple" and average socres
        if len( onebuglist ) > 1:
            k1 = np.average( [float( x[1] ) for x in onebuglist] )
            k1_avg = np.average( [float( x[2] ) for x in onebuglist] )
            onebug, onesyn = 'Multiple', 'C'*len( onebuglist[0][3] )
        # one bug did not pass filter
        elif len( onebuglist ) == 1:
            onebug, k1, k1_avg, onesyn = onebuglist[0]
        # no hits to begin with
        else:
            pass
    return k1, onebug, onesyn


def annotate_twobug( lgt_status, twobuglist ):
    """ Annotate the twobugs with their taxa """
    k2, k3_avg, k3_med, twobug, twosyn = 0, 0, 0, '-', '-'
    # one or multiple pairs passed the threshold
    if lgt_status == True:
        if len( twobuglist ) == 1: #one pair passed
            twobug_sorted, k2, k2_avg, k3_avg, k3_med, twobug, twosyn = twobuglist[0]
        elif len( twobuglist ) > 1: #multiple pairs passed
            twobug, twosyn, k2, k2_avg, k3_avg, k3_med = det_multiple_two( twobuglist )
        else: # (should not happen)
            pass #use defaults from above
    # no pairs passed the threshold, or pairs that passed were actually singletons
    else:
        pass #use defaults from above
    return k2, k3_avg, k3_med, twobug, twosyn


def det_status( matrix_u, taxaorder_u, onebug_thresh, twobug_thresh, myrange, adjust ):
    """ Find whether the status is LGT or not, and annotate with taxon calls.
    We must call LGT status and taxa at the same time. This is to avoid giving an LGT call but with "multiple" taxa, which is not helpful.
    1. Define the output.
    2. Calculate 1 bug score.
    3. Determine if there is only 1 gene or bug, or not.
    4. For those with only 1 gene or bug, we annotated with "NoLGT" and annotate if 1 clear taxon, otherwise we walk up.
    5. For those with >1 gene and/or >1 bug:
        a) If pass "NoLGT" threshold and 1 clear taxon, annotate. 
        b) If pass "NoLGT" threshold and 1+ taxon, walk up. 
        c) If fail "NoLGT" threshold and make "LGT" threshold and 1 clear pair or 1 clear taxon with multiple, annotate.
        d) If fail "NoLGT" threshold and make "LGT" threshold and no clear pairs, walk up.
        e) If fail "NoLGT" threshold and fail "LGT" threshold, walk up. """
    # Define output
    status = 'Ambiguous'
    k1, onebug, onesyn, k2, k3_avg, k3_med, twobug, twosyn = 0, '-', '-', 0, 0, 0, '-', '-'
    
    # Calculate onebug score
    nolgt_status, onebuglist, onebuglist_filt = calc_onebug( matrix_u, taxaorder_u, onebug_thresh, myrange )
    # Annotate with onebug values regardless onebug or not
    k1, onebug, onesyn = annotate_onebug( onebuglist, onebuglist_filt )

    # Check number of genes and orgs 
    numbugs, numgenes = np.shape( matrix_u )
    if numgenes == 1 or numbugs == 1:
        # Call NoLGT and annotate if single organism, otherwise walk up (NoLGT status will not change)
        if nolgt_status == True: #NoLGT
            if len( onebuglist_filt ) == 1: #one bug top
                status = 'NoLGT'
    else:
        # Calculate twobug scores
        lgt_status, twobuglist = calc_twobug( matrix_u, taxaorder_u, twobug_thresh, myrange, adjust )
        # Annotate with twobug values regardless of twobug or not
        k2, k3_avg, k3_med, twobug, twosyn = annotate_twobug( lgt_status, twobuglist )

        # Determine NoLGT, LGT, or ambiguous status
        if nolgt_status == True: # Pass one bug threshold
            if len( onebuglist_filt ) == 1: #one bug dominates, call NoLGT and annotate
                status = 'NoLGT'
            else: # pass one bug threshold but unsure about taxon
                pass
        else: #Did NOT pass one bug threshold
            if lgt_status == True: #Pass two bug threshold
                if len( twobuglist ) == 1: #one pair dominates, call LGT and annotate
                    status = 'LGT'
                else: # pass two bug threshold but do not know which pairs
                    if twobug != 'Multiple': # check if we know 1 of the bugs out of the pairs
                        status = 'LGT'
                    # else: # do not know anything within pairs, walk up
                    pass
            else: # Did NOT pass two bug threshold either or have any pairs, walk up
                pass
    taxoninfo = [k1, onebug, onesyn, k2, k3_avg, k3_med, twobug, twosyn]
    return status, taxoninfo

def check_percent( twobugsyn, genelist ):
    dict_genetypes = {}
    all_genes = []
    for i in range( len( genelist ) ):
        genelen = genelist[i].end - genelist[i].start
        dict_genetypes.setdefault( twobugsyn[i], [] ).append( genelen )
        all_genes.append( genelen )
    # get smallest gene
    smallest_gene = min( all_genes )

    # calculate LGT percentages
    total = sum( dict_genetypes['A'] ) + sum( dict_genetypes['B'] )
    c_total = 0
    if 'C' in dict_genetypes:
        c_total = sum( dict_genetypes['C'] )
        total = total + sum( dict_genetypes['C'] )
    numerator = min( [sum( dict_genetypes['A'] ), sum( dict_genetypes['B'] )] )
    lgt_perc = numerator/float( total )
    c_perc = c_total/float( total )
    return numerator, total, lgt_perc, smallest_gene, c_perc, c_total

def find_top_uniref( unireflist ):
    finalunireflist = []
    for geneuniref in unireflist:
        unirefs = geneuniref.split(',')
        topscore = 0
        # establish top score for uniref that is not 0
        for i in range( len( unirefs ) ):
            unirefname, unirefnum = unirefs[i].split(':')
            if unirefname != 'unknown':
                topscore = int( unirefnum )
                break
        finaluniref = []
        for i in range( len( unirefs ) ):
            unirefname, unirefnum = unirefs[i].split(':')
            if unirefname != 'unknown' and int( unirefnum ) >= topscore:
                finaluniref.append( unirefname )
            if i == len( unirefs ) - 1 and len( finaluniref ) == 0:
                finaluniref.append( 'unknown' )
        finalunireflist.append( '|'.join( finaluniref ) )
    return ';'.join( finalunireflist )


def find_lgt_level( taxapair ):
    """ Find level at which we have lgt """
    taxa1, taxa2 = re.split('[><\?-]', taxapair )
    taxa1list, taxa2list = taxa1.split('|'), taxa2.split('|')
    lowest_known = min( len( taxa1list ), len( taxa2list ) )
    found = False
    for i in range( lowest_known ):
        level = wu.c__list_taxa[i]
        if taxa1list[i] != taxa2list[i]:
            found = True
            break
        if i == ( lowest_known - 1 ) and found == False:
            print( 'looking farther' )
            level = wu.c__list_taxa[i+1]
    return level

# ---------------------------------------------------------------
# main
# ---------------------------------------------------------------

def main():
    """
    This script does:
    1) Groups hits into species per gene.
    2) Generates score arrays per species.
    3) Prints out the taxa into a new file for plotting.
    4) Iterate through taxonomic levels (species to kingdom).
    5) Aggregate score arrays per taxa level.
    6) Generate table with final scores per taxa per gene.
    7) Calculate one and twobug scores.
    8) If status + taxa annotation is confidently assigned, output.
    """
    args = get_args()
    #Open and write result headers
    fh_taxa = wu.try_open( args.plotfile, "w" )
    writertaxa = csv.writer( fh_taxa, dialect="excel-tab" )
    writertaxa.writerow( [
                        "contig", 
                        "contiglen", 
                        "genenum", 
                        "strand", 
                        "genestart", 
                        "geneend", 
                        "taxa", 
                        "taxastart", 
                        "taxaend", 
                        "score", 
                        "uniref50", 
                        "uniref90", 
                        "orghitnum",
                        ] )
    fh_results = wu.try_open( args.results, "w" )
    writerresults = csv.writer( fh_results, dialect="excel-tab" )
    writerresults.writerow( [
                        "contig",
                        "contiglen",
                        "status",
                        "onebugscore",
                        "twobugscore",
                        "twobugmeanadj","twobugmedadj",
                        "lgtbp", "totalbp", "lgtperc", "smallest_gene", "c_perc", "c_bp",
                        "onebug",
                        "onebugsyn",
                        "twobug",
                        "twobugsyn",
                        "uniref50",
                        "uniref90"
                        ] )
   
    #Build dictionary of genes
    dict_genes = {}
    for contig, genelist in wu.iter_contig_genes( args.gff ):
        dict_genes[contig] = genelist
    
    #Group taxa based on new genes
    for contig, hitlist in wu.iter_contig_hits( args.blast ):
        if contig in dict_genes: #some contigs have hits but not genes
            genelist = dict_genes[contig]
            contiglen = hitlist[0].qlen
            #Generate score arrays per species
            dict_genearray, dict_taxalist = {}, {}
            for gene in genelist:
                info = hits2genes( gene, hitlist, args.strand, args.overlap_hits, args.scov_hits )
                dict_taxaarray, taxalist = score_taxa( gene, info, contiglen )
                dict_genearray[ int( gene.attributes.get('ID').split('_')[-1] ) ] = dict_taxaarray
                dict_taxalist[ int( gene.attributes.get('ID').split('_')[-1] ) ] = taxalist

            #Output plotting file and format unirefs
            uniref50list, uniref90list = ['']*len(dict_genearray.keys()), ['']*len(dict_genearray.keys())
            for gene in sorted( dict_taxalist.keys() ):
                taxalist_inorder = dict_taxalist[gene]
                for taxa in taxalist_inorder:
                    writertaxa.writerow( wu.order_taxa( taxa ) )
                    uniref50list[taxa.gene-1] = taxa.uniref50
                    uniref90list[taxa.gene-1] = taxa.uniref90
            #Format uniref
            uniref90, uniref50 = find_top_uniref( uniref90list), find_top_uniref( uniref50list )
            
            #Scoring onebug/twobug scores at the species level.
            #Do not stop iteration until
            #1: Clear annotation (LGT/NoLGT)
            #2: One taxon/taxa-pair
            #3: Clear filter for percentage
            for taxalevel in reversed( wu.c__list_taxa ): 
                #Generate table for contig
                matrix, taxaorder = generate_matrix( dict_genearray, taxalevel )
                matrix_u, taxaorder_u = spike_unknown( matrix, taxaorder, args.unknown )
                #Calculate one/two bug scores and find taxon annotations
                status, taxon_list = det_status( matrix_u, taxaorder_u, args.onebug, args.twobug, args.custrange, args.adjust )
                onebugscore, onebug, onebugsyn, twobugscore, k3avg, k3med, twobug, twobugsyn = taxon_list
                #Filter for percentage of taxon with LGT
                perc_status = True
                if status == 'LGT':
                    lgtbp, totalbp, lgtperc, smallest_gene, c_perc, c_bp = check_percent( twobugsyn, genelist )
                #Stop loop if status found
                if status == 'NoLGT':
                    lgtbp, totalbp, lgtperc, smallest_gene, c_perc, c_bp = 0, 0, 0, 0, 0, 0
                    break
                elif status == 'LGT' and perc_status == True:
                    break

            #Print results
            if status == 'LGT':
                lgtlevel = find_lgt_level( twobug )
                status = status + ':' + lgtlevel
            writerresults.writerow( [
                            contig, 
                            contiglen, 
                            status, 
                            onebugscore, 
                            twobugscore, 
                            k3avg,k3med,
                            lgtbp, totalbp, lgtperc, smallest_gene, c_perc, c_bp,
                            onebug, 
                            onebugsyn, 
                            twobug, 
                            twobugsyn, 
                            uniref50, 
                            uniref90
                            ] )
    
    fh_taxa.close()
    fh_results.close()

if __name__ == "__main__":
    main()


