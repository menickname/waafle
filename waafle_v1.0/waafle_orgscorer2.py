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
        help="turn on strand specific hit assignment",
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
        "-cperc", "--cpercent",
        default=0.5,
        type=float,
        help="Contig called as LGT if percentage of C genes less than this value",
        )
    args = parser.parse_args()
    return args

def get_genenum( gene ):
    """ Added this on 12/7/17 because getting the gene number is different for WAAFLE gffs vs. answer gffs"""
    genenum = ''
    if 'ID' in gene.attributes:
        genenum = gene.attributes['ID'].split('_')[-1]
    else:
        genenum = gene.attributes['gene_number']
    return genenum

def hits2genes( gene, hits, strand_specific, lap, scov ):
    """ Group hits to genes """
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
                    int( get_genenum( gene ) ), 
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
                    int( get_genenum( gene ) ),
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
    if unknown_name in taxaorder: # if there is an unknown in the table, get values
        unknown_index = taxaorder.index( unknown_name )
        for i in range( numgenes ):
            unknown_score = 1 - np.max( table[:, i] )
            if table[unknown_index, i] == 1:
                unknown_score = 1
            final_table[ unknown_index, i ] = unknown_score
        return final_table, final_taxaorder
    elif not unknown: # if there is no unknown, and arg against spike
        return final_table, final_taxaorder
    else: # if there is no unknown, and arg for spike
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
        top_one, top_two = re.split('[><\?-]', org ) # get both taxa 
        index_one, index_two = taxaorder.index( top_one ), taxaorder.index( top_two ) # get indices for taxa
        array_one, array_two = contigarray[index_one, :], contigarray[index_two, :] # get values for taxa
        twobugmax = np.amax( np.array( [array_one, array_two] ), axis=0 ) # determine max for the pair
        for i in range( len( twobugmax ) ): # loop through each gene
            one_element, two_element = array_one[i], array_two[i] # scores for each taxon
            if threshold <= one_element and  threshold <= two_element: # if both taxa >= threshold, call 'C'
                element_list.append( [top_one, top_two] )
                taxastring += 'C'
            elif threshold <= one_element and threshold > two_element: # if first taxa >= threshold but not second taxa, call 'A'
                element_list.append( top_one )
                taxastring += 'A'
            elif threshold > one_element and threshold <= two_element: # if second taxa >= threshold but not first taxa, call 'B'
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
    1. Filter for taxa with onebug scores >= k1 and onebug scores >= within x range of top onebugscore
    2. Get top average from step 1 and retain those within x range of top average.
    4. Call NoLGT status if 1+ taxa make all filters.  """
    nolgt_status = False
    # Get all taxa within x% of top onebugscore and above onebugthresh
    onebugscore = np.max( np.amin( table, axis=1 ) )
    onebugscoremin = onebugscore - myrange
    threshold = max( [onebugthresh, onebugscoremin] )
    print( 'Searching for onebug:\n', 'onebugscore: ', onebugscore, ', ', 'onebugscoremin: ', onebugscoremin )
    indices = [i for i,j in enumerate( list(np.amin(table, axis=1) ) ) if j>=onebugscoremin and j>=onebugthresh ]
    onebuginfo = []
    if len( indices ) > 0: #Some taxa survived filter
        bugmeanlist = []
        for index in indices:
            bugname = buglist[index]
            bugmean = np.mean(table[index])
            bugscore = np.amin(table[index])
            bugmeanlist.append( [bugname, bugmean, bugscore ] )
        # Find top average and filter for all taxa within x% of top onebugavg
        bugmeanlist_sorted = sorted( bugmeanlist, key=itemgetter(1), reverse=True )
        topmean = bugmeanlist_sorted[0][1]
        topmeanmin = topmean - myrange
        print( 'Getting onebug top mean:' )
        for x,y,z in bugmeanlist:
            print( x.split('|')[-1], y, z )
        print( 'topmean: ', topmean, ' topmeanmin: ', topmeanmin )
        for bug, avg, bugscore in bugmeanlist:
            if avg >= topmeanmin:
                # Find synteny if make threshold
                element_list, taxastring, shift_counter = find_taxa_order( 'NoLGT', table, buglist, bug, threshold )
                onebuginfo.append( [bug, bugscore, avg, taxastring] )
    print( 'After filtering by k1:\n', onebuginfo )
    if len( onebuginfo ) > 0:
        nolgt_status = True
    return nolgt_status, onebuginfo


def filter_c( twobugsyn, genelist ):
    dict_genetypes = {}
    all_genes = []
    for i in range( len( genelist ) ):
        genelen = genelist[i].end - genelist[i].start
        dict_genetypes.setdefault( twobugsyn[i], [] ).append( genelen )
        all_genes.append( genelen )
    # get smallest gene
    smallest_gene = min( all_genes )
    # calculate LGT percentages
    A_type, B_type, C_type = 0, 0, 0
    if 'A' in dict_genetypes:
        A_type = sum( dict_genetypes['A'] )
    if 'B' in dict_genetypes:
        B_type = sum( dict_genetypes['B'] )
    if 'C' in dict_genetypes:
        C_type = sum( dict_genetypes['C'] )
    total = A_type + B_type + C_type
    numerator = min( A_type, B_type )
    lgt_perc = numerator/float( total )
    c_perc = C_type/float( total )
    return c_perc, C_type, numerator, total


def filter_twobug( genelist, complement, table, buglist, twobugthresh, myrange, cthresh ):
    """ Filters complement scores by:
    1. Get top twobug score and retain those within x% (we have already kept only those >=k2).
    2. Get top average from step 1 and retain those within x%.
    3. Remove pairs with C% > c%. 
    5. Call LGT status if 1+ taxa pairs make all filters. """
    lgt_status = False
    # Get all taxa within x% of top twobugscore
    twobugscore = sorted( complement, key=itemgetter(1), reverse=True )[0][1]
    twobugscoremin = twobugscore - myrange
    threshold = max( [twobugthresh, twobugscoremin] )
    complement_scorefilt = [j for i, j in enumerate( complement ) if j[1] >= twobugscoremin ]
    print( 'Filter by top twobugscore range:\n', 'twobugscore: ', twobugscore, 'twobugscoremin: ', twobugscoremin )
    for x,y,z,aa in complement_scorefilt:
        print( x.split('-')[0].split('|')[-1], x.split('-')[1].split('|')[-1],y,z,aa )
    # Find top average and filter for all taxa within x% of top twobugavg
    topmean = sorted( complement_scorefilt, key=itemgetter(2), reverse=True )[0][2]
    topmeanmin = topmean - myrange
    possibilities = []
    print( 'Filter by top twobugavg range:\n', 'topmean: ', topmean, 'topmeanin: ', topmeanmin )
    for bugnames, twobugscore, twobugavg, penalty in complement_scorefilt:
        if twobugavg >= topmeanmin:
            # Find synteny and percentage of Cs
            element_list, taxastring, shift = find_taxa_order( 'LGT', table, buglist, bugnames, threshold )
            c_perc, c_bp, lgt_bp, total_bp = filter_c( taxastring, genelist )
            numtaxaleft = len( set( taxastring.replace( 'C', '' ) ) )
            print( 'Filter by C%:')
            print( bugnames.split('-')[0].split('|')[-1], bugnames.split('-')[1].split('|')[-1], 'C%: ', c_perc, ', C_bp: ', c_bp, ', LGT_bp: ', lgt_bp, ', Total_bp: ', total_bp )
            if c_perc <= cthresh and numtaxaleft == 2: # Include pair if lower than C% threshold
                newname = calc_donorrecip( bugnames, element_list, taxastring, shift )
                possibilities.append( [bugnames, twobugscore, twobugavg, penalty, newname, taxastring] )
    # Filter by twobug score
    print( 'Remaining possibilities:' )
    for info in possibilities:
        print( re.split('[><\?]', info[4] )[0].split('|')[-1], re.split('[><\?]', info[4] )[1].split('|')[-1], info[5], info[1], info[2], info[3] )
    if len( possibilities ) >= 1:
        lgt_status = True
    return lgt_status, possibilities


def calc_twobug( genelist, table, buglist, twobugthresh, myrange, cperc ):
    """ Calculates multiple k scores and store.
    1. For every pair, calculate twobug score, twobug avg, and penalty.
    2. Filter results through filter_twobug.
    """
    complementlist = []
    for i in range( table.shape[0]-1 ):
        bug1 = table[i, :]
        if max( bug1 ) < twobugthresh: #will never be able to form a pair
            continue
        else:
            for j in range( i+1, table.shape[0] ):
                bug2 = table[j, :]
                bugdiff = bug1 - bug2
                if ( bugdiff < 0 ).all() or ( bugdiff > 0 ).all(): # one taxon consistently higher than the other
                    continue
                elif max( bug2 ) < twobugthresh: #will never be able to form a pair
                    continue
                else:
                    twobugarray = np.array( [bug1, bug2] )
                    # calculate values
                    twobugscore = np.min(np.amax(twobugarray, axis=0)) #twobug score
                    twobugaverage = np.mean(np.amax(twobugarray, axis=0)) #twobug avg
                    penalty = np.mean( np.absolute( twobugarray[0] - twobugarray[1] ), axis=0 ) #penalty value
                    # get bug names
                    bugsort = sorted( [str(buglist[i]), str(buglist[j])] )
                    bugnames = bugsort[0] + '-' + bugsort[1]
                    # place those greater than k2 into the list
                    if twobugscore >= twobugthresh:
                        complementlist.append( [ bugnames, twobugscore, twobugaverage, penalty ] )
    # for all pairs, filter results
    print( 'Searching for LGT:' )
    for x,y,z,aa in complementlist:
        print( x.split('-')[0].split('|')[-1], x.split('-')[1].split('|')[-1], y, z, aa )
    lgt_status, possibilities = False, [] #assume all taxa pairs filtered out, no LGT possible
    if len( complementlist ) != 0: # if we have pairs, filter
        lgt_status, possibilities = filter_twobug( genelist, complementlist, table, buglist, twobugthresh, myrange, cperc )
    return lgt_status, possibilities


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
    """ From a set of taxa, determine the LCA. """
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


def get_new_taxastring( first_taxon_sorted, second_taxon_sorted, oldtaxapair, oldtaxastring, consensus_synteny ):
    oldtaxaone, oldtaxatwo = re.split( '[><\?]', oldtaxapair )
    oldtaxasign = re.search( '[><\?]', oldtaxapair ).group()
    newtaxastring = ''
    if oldtaxaone == first_taxon_sorted or oldtaxatwo == second_taxon_sorted:
        newtaxastring = oldtaxastring
        newtaxasign = oldtaxasign
    elif oldtaxaone == second_taxon_sorted or oldtaxatwo == first_taxon_sorted:
        newtaxasign, newtaxastring = flip_annot( oldtaxasign, oldtaxastring )
    newname = first_taxon_sorted + newtaxasign + second_taxon_sorted
    return newname, newtaxastring
    

def organize_synteny( list_of_possibilities ):
    """ Loop through each pair and:
    1. Reorganize synteny so that 'A' is first.
    2) Store synteny as a list.
    3) Bin first and second taxa into sets."""
    firsttaxon_set, secondtaxon_set = set([]), set([])
    array_values, synteny_check, stop = [], [], False
    for info in list_of_possibilities:
        sortedname, twobugscore, twobugavg, penalty, newname, taxastring = info
        array_values.append( [twobugscore, twobugavg, penalty] )
        firsttaxon = '' # Find the first taxon letter
        for letter in taxastring:
            if letter == 'B':
                firsttaxon = 'B'
                break
            elif letter == 'A':
                firsttaxon = 'A'
                break
            else:
                continue
        taxaone, taxatwo = re.split( '[><\?]', newname )
        if firsttaxon == 'A':
            synteny_check.append( taxastring )
            firsttaxon_set.add( taxaone )
            secondtaxon_set.add( taxatwo )
        elif firsttaxon == 'B':
            taxatwo, taxaone = re.split( '[><\?]', newname )
            new_sign, new_taxastring = flip_annot( '<', taxastring ) #ignore new_sign, just inverting taxastring
            synteny_check.append( new_taxastring )
            firsttaxon_set.add( taxaone )
            secondtaxon_set.add( taxatwo )
        else: #firsttaxon=='', since all are Cs, should not happen since filter 50% Cs
            synteny_check.append( taxastring )
            firsttaxon_set.add( taxaone )
            firsttaxon_set.add( taxatwo )
            secondtaxon_set.add( taxaone )
            secondtaxon_set.add( taxatwo )
    return firsttaxon_set, secondtaxon_set, array_values, synteny_check 


def match_synteny( synteny_check ):
    """ Loop through each gene and determine consensus synteny. 
    Example 1: 
    Pair 1: AABB, Pair 2: AABB
    consensus_synteny = 'AABB', consensus_vals = [1,1,1,1], consensus_value = 1

    Example 2:
    Pair 1: AABB, Pair 2: AAAB
    consensus_synteny = 'AACB', consensus_vals = [1,1,2,1], consensus_value = 2

    Example 3:
    Pair 1: AABB, Pair 2: AACB
    consensus_synteny = 'AACB', consensus_vals = [1,1,1,1], consensus_value = 1
     """
    consensus_synteny = '' #store consensus synteny
    consensus_vals = [] #store number of possibilities per gene
    for i in range( len( synteny_check[0] ) ): #go through each letter
        position_vals = set( [] )
        for j in range( len( synteny_check ) ): #go through each pair
            position_vals.add( synteny_check[j][i] )
        if len( list( position_vals ) ) == 1:
            consensus_synteny += list( position_vals )[0]
        else:
            consensus_synteny += 'C'
        consensus_vals.append( len( list( position_vals ) ) )
    consensus_value = max( consensus_vals )
    return consensus_synteny, consensus_value


def det_multiple_two( list_of_possibilities, genelist, cthresh ):
    """ Determine if we can annotate LGT at a lower level (do we have to walk up)
    1. For each LGT pair: 
        a) Get the synteny to start with A. 
        b) Get the taxa assigned to 'A' into one set, and the taxa assigned to 'B' into a second set.
    2. Get a consensus synteny and value from a)
    3. Average the twobugscore, twobugavg, and penalty across all pairs
     """
    # Among potential pairs, see if they share synteny and partners
    firsttaxon_set, secondtaxon_set, array_values, synteny_check = organize_synteny( list_of_possibilities )
    consensus_synteny, consensus_value = match_synteny( synteny_check )
    consensus_cperc = filter_c( consensus_synteny, genelist )[0]
    twobugscore, twobugavg, penalty = np.average( np.array( array_values ), axis=0 )
    print( "Consider multiple taxa:\n", 'consensus_synteny: ', consensus_synteny, 'consensus_value: ', consensus_value, 'consensus_cperc: ', consensus_cperc )
    print( firsttaxon_set, secondtaxon_set )
    # Determine if we have an LGT pair in which one taxon is known, and the second taxon is semi-annotated
    # To do this, we check 1) consensus_value and the 2) second taxon (list of individual taxa)
    # consensus_value checks that all the synteny of the pairs match
    # multiple_count and high_lvl check the second taxon annotation
        # multiple_count checks the number of taxonomic levels at which we cannot annotate the second taxon 
            # Eg 1: If the second taxon could be "k__A|p__B|c__C|o__D|f__E" or "k__A|p__B|c__C|o__D|f__F", it is unknown at 1 level, "f"
            # Eg 2: If the second taxon could be "k__A|p__B|c__C|o__D|f__E" or "k__A|p__B|c__C|o__F|f__G", it is unknown at 2 levels, "o" and "f"
            # You can only have 1 level where you do NOT know the taxon for it to be semi-annotated
        # high_lvl checks for where the "multiple" LGT is:
            # We do want to know if g__TaxaA|s__Multiple exchanges with g__TaxaB|s__TaxaC.
            # We do not care if g__TaxaA|s__Multiple exchanges with g__TaxaA|s__TaxaD.
    newcombined, newcombinedstring = 'Multiple', 'C'*len( synteny_check[0] )
    high_lvl, multiple_count = False, 0
    first_taxon, second_taxon = '', ''
    numtaxaleft = len( set( consensus_synteny.replace( 'C', '' ) ) )
    print( 'numtaxaleft: ', numtaxaleft )
    if consensus_cperc <= cthresh and numtaxaleft == 2: # check final contig is still < cthresh and contains more than 1 taxon
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
        print( 'multiple_count:', multiple_count, 'high_lvl:', high_lvl )
        # If we have semi-annotation and know the LGT level, output the information
        if multiple_count == 1 and high_lvl:
            first_taxon_sorted, second_taxon_sorted = sorted( [first_taxon, second_taxon] )
            newcombined, newcombinedstring = get_new_taxastring( first_taxon_sorted, 
                                                second_taxon_sorted, 
                                                list_of_possibilities[0][-2], 
                                                list_of_possibilities[0][-1],
                                                consensus_synteny )
    return newcombined, newcombinedstring, twobugscore, twobugavg, penalty
    

def annotate_onebug( onebuglist ):
    """ Annotate the onebug with its taxon, if nothing passed annotate with default """
    onebugscore, onebug, onesyn = 0, '-', '-'
    # multiple taxa passed the filter, call "Multiple" and average one bug score
    if len( onebuglist ) > 1:
        onebugscore = np.average( [float( x[1] ) for x in onebuglist] )
        onebugavg = np.average( [float( x[2] ) for x in onebuglist] )
        onebug, onesyn = 'Multiple', 'C'*len( onebuglist[0][3] )
    # one bug passed the filter, annotate the taxon
    elif len( onebuglist ) == 1:
        onebug, onebugscore, onebugavg, onesyn = onebuglist[0]
    # nothing passed the filter
    else:
        pass
    return onebugscore, onebug, onesyn


def annotate_twobug( lgt_status, twobuglist, genelist, cperc ):
    """ Annotate the twobugs with their taxa, if nothing passed annotate with default """
    twobugscore, twobugavg, penalty, twobug, twosyn = 0, 0, 0, '-', '-'
    # one or multiple pairs passed the threshold
    if lgt_status == True:
        if len( twobuglist ) == 1: #one pair passed
            twobug_sorted, twobugscore, twobugavg, penalty, twobug, twosyn = twobuglist[0]
        elif len( twobuglist ) > 1: #multiple pairs passed
            twobug, twosyn, twobugscore, twobugavg, penalty = det_multiple_two( twobuglist, genelist, cperc )
        else: # (should not happen)
            pass #use defaults from above
    else: # no pairs passed the threshold
        pass
    return twobugscore, penalty, twobug, twosyn


def assign_status( nolgt_status, onebuglist_filt, lgt_status, twobuglist_filt, twobug ):
    """ Assign LGT status:
    1) If pass k1 and find 1 taxon, annotate NoLGT.
    2) If not pass k1 and pass k2 and find 1 pair OR multiple pair (1 taxon known), annotate LGT. 
    """
    status = "Ambiguous"
    if nolgt_status == True: # 1+ taxa passed one bug threshold
        if len( onebuglist_filt ) == 1: #one bug dominates, call NoLGT and annotate
           status = 'NoLGT'
    else: # 0 taxa passed one bug threshold
        if lgt_status == True: # 1+ taxa pairs passed two bug threshold
            if len( twobuglist_filt ) == 1: #one pair dominates, call LGT and annotate
                status = 'LGT'
            else: # pass two bug threshold but do not know which pairs
                if twobug != 'Multiple': # check if we know 1 of the bugs out of the pairs
                    status = 'LGT'
    return status


def det_status( genelist, matrix_u, taxaorder_u, onebug_thresh, twobug_thresh, myrange, cperc ):
    """ Find whether the status is LGT or not, and annotate with taxon calls.
    We must call LGT status and taxa at the same time. This is to avoid giving an LGT call but with "multiple" taxa, which is not helpful.
    1. Define the output.
    2. Calculate 1 bug score.
    3. Determine if there is only 1 gene or bug, or not.
    4. For those with only 1 gene or bug, we annotated with "NoLGT" and annotate if 1 clear taxon, otherwise we walk up.
    5. For those with >1 gene and/or >1 bug:
        a) If pass "NoLGT" threshold and 1 clear taxon, annotate. 
        c) If fail "NoLGT" threshold and make "LGT" threshold and 1 clear pair or 1 clear taxon with multiple, annotate.
        e) All other cases, status: "Ambiguous", multiple taxa annotation, walk up. """
    # Define output
    status = 'Ambiguous'
    onebugscore, onebug, onesyn, twobugscore, penalty, twobug, twosyn = 0, '-', '-', 0, 0, '-', '-'
    # Calculate onebug score and annotate with onebug values regardless of status
    nolgt_status, onebuglist = calc_onebug( matrix_u, taxaorder_u, onebug_thresh, myrange )
    onebugscore, onebug, onesyn = annotate_onebug( onebuglist )
    # Check number of genes and orgs 
    numbugs, numgenes = np.shape( matrix_u )
    # If only 1 gene or 1 bug, call NoLGT and annotate if 1 taxon passes threshold, otherwise walk up
    if numgenes == 1 or numbugs == 1:
        if nolgt_status == True: #NoLGT
            if len( onebuglist ) == 1:
                status = 'NoLGT'
    # If >1 gene or bug, calculate twobug score and annotate with twobug values regardless of status
    else:
        lgt_status, twobuglist = calc_twobug( genelist, matrix_u, taxaorder_u, twobug_thresh, myrange, cperc )
        twobugscore, penalty, twobug, twosyn = annotate_twobug( lgt_status, twobuglist, genelist, cperc )
        # Determine NoLGT, LGT, or ambiguous status
        status = assign_status( nolgt_status, onebuglist, lgt_status, twobuglist, twobug )
    taxoninfo = [onebugscore, onebug, onesyn, twobugscore, penalty, twobug, twosyn]
    return status, taxoninfo


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
                dict_genearray[ int( get_genenum( gene ) ) ] = dict_taxaarray
                dict_taxalist[ int( get_genenum( gene ) ) ] = taxalist

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
            for taxalevel in reversed( wu.c__list_taxa ): 
                #Generate table for contig
                print( '=====================================' )
                print( 'Start analysis for contig: ', contig )
                print( 'Looping through ', taxalevel )
                matrix, taxaorder = generate_matrix( dict_genearray, taxalevel )
                print( 'Gene table:' )
                for i in range( len( taxaorder ) ):
                    print( matrix[i,:], taxaorder[i].split('|')[-1] )
                matrix_u, taxaorder_u = spike_unknown( matrix, taxaorder, args.unknown )
                print( 'Spike in unknown:' )
                for i in range( len( taxaorder_u ) ):
                    print( matrix_u[i,:], taxaorder_u[i].split('|')[-1] )
                #Calculate one/two bug scores and find taxon annotations
                status, taxon_list = det_status( genelist, matrix_u, taxaorder_u, args.onebug, args.twobug, args.custrange, args.cpercent )
                onebugscore, onebug, onebugsyn, twobugscore, penalty, twobug, twobugsyn = taxon_list
                print( '\nResults: ' )
                print( taxon_list ) 
                # Break loop if found status
                if status == 'LGT' or status == 'NoLGT':
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


