#!/usr/bin/env python

"""
WAAFLE MODULE: waafle_aggregator.py

Authors:
Tiffany Hsu
Eric Franzosa

This script takes the output from waafle_lgtscorer.py and determines at what level contigs have LGT.
Run with the "-h" flag for usage help.

"""

from __future__ import print_function # python 2.7+ required
import os, sys, csv, argparse
from operator import itemgetter, attrgetter, methodcaller
import waafle_utils as wu
import numpy as np
import numpy.ma as ma
import re

# ---------------------------------------------------------------
# constants
# ---------------------------------------------------------------

c__list_taxa = ["k", "p", "c", "o", "f", "g", "s"]
cr__list_taxa = ["s", "g", "f", "o", "c", "p", "k"]


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
        help="concatenated output from waafle_lgtscorer",
        )
    parser.add_argument(
        "-o", "--out",
        help="results",
        default="final-scoredcontigs.tsv",
        )
    args = parser.parse_args()
    return args

def get_taxalevel( status, taxonomy ):
    """
    Get taxa level (last taxa in the taxonomy will be the level).
    """
    bugone = taxonomy.split('|')[0]
    if status == 'LGT' or status == 'ambiguous-LGT':
        pair_one, pair_two = re.split('[><?]', bugone )
        taxalevel = c__list_taxa[ len( pair_one.split('.') ) - 1 ]
    else:
        taxalevel = c__list_taxa[ len( bugone.split('.') ) - 1 ]
    return taxalevel

def sort_by_taxalevels( infolist ):
    """
    Get taxa for a single contig in sorted order.
    """
    sortedlist = [None]*len(infolist)
    for info in infolist:
        taxalevel = get_taxalevel( info.status, info.taxa )
        index = c__list_taxa.index( taxalevel )
        sortedlist[index] = info
    return sortedlist

def process_lgtstatus( statuslist ):
    """
    This function calculates how certain we are about the LGT status.
    In theory, all contigs, when progressing from kingdom to species, should go from NoLGT to LGT.
    We assume only the NoLGT and LGT states are certain, while ambiguous states actually go either way. Thus,
        1. A contig will be assigned LGT if "LGT" appears.
        2. A contig will be assigned NoLGT if "LGT" never appears. It does not matter if "ambiguous" appears.
        3. If "ambiguous-LGT" appears, it is assigned that status.
        4. Otherwise, it is assigned "ambiguous-NoLGT".
    """
    lgtlevel = []
    for i in range( len( statuslist ) ):
        level = c__list_taxa[i]
        status = statuslist[i]
        if "LGT" in statuslist:
            indices = [i for i, x in enumerate(statuslist) if x == "LGT"]
            start = c__list_taxa[ indices[0] ]
            end = c__list_taxa[ indices[ len(indices)-1 ] ]
            endstatus = "LGT"
            if "NoLGT" in statuslist[ indices[ len(indices)-1 ]: ]:
                endstatus = "Uncertain_LGT" #maybe just call it NoLGT, also switch order
        elif "NoLGT" in statuslist:
            indices = [i for i, x in enumerate(statuslist) if x == "NoLGT"]
            start = c__list_taxa[ indices[0] ]
            end = c__list_taxa[ indices[ len(indices)-1 ] ]
            endstatus = "NoLGT"
        elif "ambiguous-LGT" in statuslist:
            indices = [i for i, x in enumerate(statuslist) if x == "ambiguous-LGT"]
            start = c__list_taxa[ indices[0] ]
            end = c__list_taxa[ indices[ len(indices)-1 ] ]
            endstatus = "ambiguous-LGT"
        else:
            indices = [i for i, x in enumerate(statuslist) if x == "ambiguous-NoLGT"]
            start = c__list_taxa[ indices[0] ]
            end = c__list_taxa[ indices[ len(indices) - 1 ] ]
            endstatus = "ambiguous-NoLGT"
    return start, end, endstatus

def split_taxa( taxa, status ):
    taxalist = []
    if len( taxa.split('|') ) > 1: #if there are multiple taxa
        multtaxa = taxa.split('|')
        if status == 'ambiguous-LGT' or status == 'LGT':
            for orgpair in multtaxa:
                taxaone, taxatwo = re.split('[><?]', orgpair )
                taxalist.append( (taxaone, taxatwo) )
        else:
            for org in multtaxa:
                taxalist.append( org )
    else:
        if status == 'ambiguous-LGT' or status == 'LGT':
            taxaone, taxatwo = re.split('[><?]', taxa )
            taxalist = [ (taxaone, taxatwo) ] 
        else:
            taxalist = [ taxa ]
    return taxalist

def get_taxaset( taxalist, level ):
    """
    Take taxa and create a set based on its status.
    Return both the set and whether the taxa has multiple pairs/taxa with LGT/NoLGT, respectively.
    """
    taxaset = set( [] )
    for orgs in taxalist:
        if type( orgs ) is tuple:
            orgone, orgtwo = orgs[0], orgs[1]
            orgone_abbr, orgtwo_abbr = orgone.split('.')[c__list_taxa.index( level )], orgtwo.split('.')[c__list_taxa.index( level )]
            taxaset.add( orgone_abbr )
            taxaset.add( orgtwo_abbr )
        else:
            orgs_abbr = orgs.split('.')[c__list_taxa.index( level )]
            taxaset.add( orgs_abbr )
    return taxaset

def consistent_taxa( statuslist, taxalist ):
    """
    This function will check whether the taxa status is correct by:
        1. Checking if each subsequent taxa above fits in the previous one.
            a. If not, the previous taxa will be flagged for its final level.
    """
    # start at the kingdom level and go down the list
    inconsistent_list = []
    higher_taxa = split_taxa( taxalist[0], statuslist[0] )
    higher_level = get_taxalevel( statuslist[0], taxalist[0] )
    higher_set = get_taxaset( higher_taxa, higher_level )
    consistent_taxaset = set([])
    for i in range( len( taxalist ) ):
        if i > 0:
            lower_taxa = split_taxa( taxalist[i], statuslist[i] )
            lower_set = get_taxaset( lower_taxa, higher_level )
            if len( higher_set & lower_set ) == 0:
                inconsistent_list.append( higher_level )
            else:
                consistent_taxaset |= ( higher_set & lower_set )
            higher_taxa = lower_taxa
            higher_level = get_taxalevel( statuslist[i], taxalist[i] )
            higher_set = get_taxaset( higher_taxa, higher_level )
    return inconsistent_list, consistent_taxaset

def known_taxa( statuslist, taxalist ):
    unknown_list, abbrtaxa_list = [], []
    for i in range( len( taxalist ) ):
        taxa = taxalist[i]
        if len( taxa.split('|') ) > 1:
            taxalevel = get_taxalevel( statuslist[i], taxa )
            unknown_list.append( taxalevel )
    return unknown_list
    
def get_direction( knownone, dirlist, taxalist, syntenylist ):
    if '>' in dirlist:
        first, second = re.split( '>', taxalist[ dirlist.index('>') ] )
        synteny = syntenylist[ dirlist.index('>') ]
        if first in knownone:
            direction = '>'
        else:
            direction = '<'
            synteny = synteny.replace( 'A', 'D' ).replace( 'B', 'A' ).replace( 'D', 'B' )
    elif '<' in dirlist:
        first, second = re.split( '<', taxalist[ dirlist.index('<') ] )
        synteny = syntenylist[ dirlist.index('<') ]
        if first in knownone:
            direction = '<'
        else:
            direction = '>'
            synteny = synteny.replace( 'A', 'D' ).replace( 'B', 'A' ).replace( 'D', 'B' )
    else:
        direction = '?'
        first, second = re.split( '\?', taxalist[ dirlist.index('?') ] )
        synteny = syntenylist[ dirlist.index('?') ]
        if first not in knownone:
            synteny = synteny.replace( 'A', 'D' ).replace( 'B', 'A' ).replace( 'D', 'B' )
    return direction, synteny

def format_multiple_lgt( taxa, synteny, unknown ):
    pairset, taxaset, dirlist = set([]), set([]), []
    taxalist, syntenylist = taxa.split('|'), synteny.split('|')
    finaltaxa, newsynteny = '', ''
    for i in range( len( taxalist ) ):
        taxaone, taxatwo = re.split( '[><\?]', taxalist[i] )
        dirlist.append( re.search( '[><\?]', taxalist[i] ).group() )
        taxalevel = taxaone.split('.')[-1].split('__')[0]
        taxaset |= set( [taxaone, taxatwo] )
        if i == 0:
            pairset = set( [taxaone, taxatwo] )
        else:
            pairset &= set( [taxaone, taxatwo] )
    if len( pairset ) == 1:
        taxaknown = list( pairset )[0].split('.')[c__list_taxa.index( taxalevel )]
        direction, newsynteny = get_direction( pairset, dirlist, taxalist, syntenylist )
        if unknown == True:
            finaltaxa = taxaknown + direction + taxalevel + '__multiple'
        else:
            othertaxa = taxaset - pairset
            otherformat = [x.split('.')[-1] for x in list( othertaxa )]
            finaltaxa = taxaknown + direction + '|'.join( otherformat )
    else:
        taxalevel = re.split( '[><\?]', taxalist[0] )[0].split('.')[-1].split('__')[0]
        if unknown == True:
            finaltaxa = taxalevel + '__multiple?' + taxalevel + '__multiple'
            newsynteny = 'unknown'
        else:
            reformat_list = []
            for pairs in taxalist:
                sign = re.search( '[><\?]', pairs ).group()
                taxaone, taxatwo = re.split( '[><\?]', pairs )
                abbrone, abbrtwo = taxaone.split('.')[-1], taxatwo.split('.')[-1]
                reformat_list.append( abbrone + sign + abbrtwo ) 
            finaltaxa = '|'.join( reformat_list )
            newsynteny = synteny
    return finaltaxa, newsynteny

def make_call( sortedlist, end, unknown ):
    info = sortedlist[ c__list_taxa.index( end ) ]
    taxalist = split_taxa( info.taxa, info.status )
    finaltaxa, finalsynteny = info.taxa, info.synteny
    if len( taxalist ) > 1:
        if info.status == 'LGT' or info.status == 'ambiguous-LGT':
            finaltaxa, finalsynteny = format_multiple_lgt( info.taxa, info.synteny, unknown )
        else:
            finalsynteny = list( set( info.synteny.split('|') ) )[0]
            if unknown == True:
                finaltaxa = end + '__multiple'
            else:
                finallist = []
                for singles in info.taxa.split('|'):
                    finallist.append( singles.split('.')[c__list_taxa.index(end)] )
                finaltaxa = '|'.join( finallist )
    else:
        if info.status == 'LGT' or info.status == 'ambiguous-LGT':
            direction = re.search( '[><?]', info.taxa ).group()
            taxaone, taxatwo = [x.split('.')[c__list_taxa.index(end)] for x in re.split( '[><?]', info.taxa )]
            finaltaxa = taxaone + direction + taxatwo
        else:
            finaltaxa = info.taxa.split('.')[c__list_taxa.index(end)]
    return finaltaxa, finalsynteny

def format_call( status, onescore, twoscore, taxa, synteny, inconsistent_list ):
    f_onescore, f_twoscore = format( onescore, '.4f'), format( twoscore, '.4f' )
    call = ';'.join( [status, f_onescore, f_twoscore, taxa, synteny] )
    if len( inconsistent_list ) > 0:
        call = call + '*'
    return call   

def get_topuniref( gene_uniref ):
    """
    Rank Unirefs if more than 1 with same count.
    """
    potential_unirefs, topuniref = gene_uniref.split('|'), ''
    if len( potential_unirefs ) == 1:
        topuniref = potential_unirefs[0].strip().split(':')[0]
    else:
        unireflist = [ [unirefvalue.strip().split(':')[0], float( unirefvalue.strip().split(':')[1] )] for unirefvalue in potential_unirefs ]
        sorted_uniref = sorted( unireflist, key=lambda unirefs: unirefs[1], reverse=True )
        within_gene_list = []
        topvalue = float( sorted_uniref[0][1] )
        for i in range( len( sorted_uniref ) ):
            topuniref, freq = sorted_uniref[i][0], float( sorted_uniref[i][1] )
            if topvalue <= freq:
                within_gene_list.append( topuniref )
            else:
                if 'unknown' in within_gene_list and len( within_gene_list ) == 1:
                    within_gene_list = []
                    within_gene_list.append( topuniref )
                    break
                else:
                    break
        topuniref = '|'.join( within_gene_list )
    return topuniref

def format_uniref( myuniref ):
    """
    Split Unirefs into print-able format.
    """
    final_uniref = ''
    if len( myuniref.split(';') ) > 1: #split into genes
        topunireflist = []
        for gene_uniref in myuniref.split(';'):
            topuniref = get_topuniref( gene_uniref )
            topunireflist.append( topuniref )
        final_uniref = ';'.join( topunireflist )
    else:
        final_uniref = get_topuniref( myuniref.split(';')[0] )
    return final_uniref
 
#--------------------------------------------------------------
# main
#--------------------------------------------------------------
def main( ):

    args = get_args()

    dict_contigs = {}
    # create dictionary of contigs where key=contig and items=list of lgtscorer output
    for astrline in open( args.input ):
        aastrline = astrline.split('\t')
        contigname = aastrline[0]
        if contigname != 'contig':
            dict_contigs.setdefault( contigname, [] ).append( wu.Scores( aastrline ) )

    # print header
    print( '\t'.join( ['contig', 'length', 'uniref50', 'uniref90', 'preferred_call', 'unambiguous_call', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']) )

    # for each contig in the results, sort lgtscorer output by phylogenetic lvl
    for contig in dict_contigs.keys():
        infolist = dict_contigs[contig]
        sortedlist = sort_by_taxalevels( infolist ) 
        
        # get the status and taxa in order
        statuslist, scorelist, taxalist = [], [], []
        for info in sortedlist:
            statuslist.append( info.status )
            taxalist.append( info.taxa )

        # find lgt status, unknowns, and consistency
        start, end, status = process_lgtstatus( statuslist )
        if status == 'LGT' or status == 'ambiguous-LGT':
            call = status + ':' + start
        else:
            call = status
        unknown_list = known_taxa( statuslist, taxalist )
        inconsistent_list, consistent_taxaset = consistent_taxa( statuslist, taxalist )
        
        # make final and unambiguous calls
        finaltaxa, finalsynteny = make_call( sortedlist, end, True )
        finalinfo = sortedlist[c__list_taxa.index( end )]
        finalcall = format_call( call, finalinfo.onescore, finalinfo.twoscore, finaltaxa, finalsynteny, inconsistent_list ) 
        unambigcall = finalcall
        if len( unknown_list ) >= 1:
            lastknown_index = max( min( c__list_taxa.index( unknown_list[0] ) - 1, c__list_taxa.index( end ) ), 0 )
            lastknown = c__list_taxa[ lastknown_index ]
            unambigtaxa, unambigsynteny = make_call( sortedlist, lastknown, True )
            unambiginfo = sortedlist[lastknown_index]
            if status == 'LGT' or status == 'ambiguous-LGT':
                call = status + ':' + lastknown
            else:
                call = status
            unambigcall = format_call( call, unambiginfo.onescore, unambiginfo.twoscore, unambigtaxa, unambigsynteny, inconsistent_list )
        
        # get remaining info
        infolist = []
        for i in range( len( sortedlist ) ):
            info = sortedlist[i]
            remaintaxa, remainsynteny = make_call( sortedlist, c__list_taxa[i], False )
            call = format_call( info.status, info.onescore, info.twoscore, remaintaxa, remainsynteny, [] )
            if c__list_taxa[i] in inconsistent_list:
                call = format_call( info.status, info.onescore, info.twoscore, remaintaxa, remainsynteny, inconsistent_list ) 
            infolist.append( call )
        
        # print everything
        final_uniref50, final_uniref90 = format_uniref( info.uniref50 ), format_uniref( info.uniref90 )
        final_line = [info.contig, str( info.length ), final_uniref50, final_uniref90, finalcall, unambigcall] + infolist
        print( '\t'.join( final_line ) )
        
if __name__ == "__main__":
    main()


