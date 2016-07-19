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

def get_taxaset( taxa, status ):
    """
    Take taxa and create a set based on its status.
    Return both the set and whether the taxa has multiple pairs/taxa with LGT/NoLGT, respectively.
    """
    taxaset = set([])
    elimset = set([])
    multipletaxa = False
    if len( taxa.split(';') ) > 1:
        multipletaxa = True
        if status == 'ambiguous-LGT' or status == 'LGT':
            i = 0
            for splittaxa in taxa.split(';'):
                taxaone = splittaxa.split('-')[0]
                taxatwo = splittaxa.split('-')[1]
                taxaset |= set([ taxaone, taxatwo ] )
                if i == 0:
                    elimset = set([ taxaone, taxatwo ])
                    i += 1
                else:
                    elimset &= set([ taxaone, taxatwo ])
        else:
            for splittaxa in taxa.split(';'):
                taxaset.add( splittaxa )
    else:
        if status == 'ambiguous-LGT' or status == 'LGT':
            taxaset.add( taxa.split('-')[0] )
            taxaset.add( taxa.split('-')[1] )
        else:
            taxaset.add( taxa )
    return multipletaxa, taxaset, elimset

def compare_taxasets( original_taxaset, next_taxaset, status, next_status, index, next_multipletaxa ):
    """
    Compare two taxasets and their statuses at the appropriate taxalevel.
    We return whether there was a match between the two (counter) and whether it was status we've called
    """
    unknown_taxa = [] #this tracks the indices where there is unsure-ty about the taxa
    counter = 0 #this tracks the number of times there is a match, but we only need 1 match
    for origtaxa in original_taxaset:
        for nexttaxa in next_taxaset:
            if re.search( nexttaxa.replace('|', ';'), origtaxa.replace('|', ';') ): #for some reason '|' breaks re.search
                counter += 1
    return min( counter, 1 )

def process_taxastatus( taxastart, taxaend, taxalist, status, statuslist ):
    """
    This function will check whether the taxa status is correct by:
        1. Checking if each subsequent taxa above fits in the previous one.
        2. Checking if multiple taxa were called for that level.
    """
    # start at species level and determine if a) higher level matches lower level b) we know the single representative taxon.
    lowest_taxa, lowest_status = taxalist[ c__list_taxa.index( 's' ) ], statuslist[ c__list_taxa.index( 's' ) ]
    lowest_multiple, lowest_taxaset, lowest_elimset = get_taxaset( lowest_taxa, statuslist[c__list_taxa.index( 's' )] )
    match_list, multiple_list = [], [ lowest_multiple ]
    for i in reversed( range( 1, len( c__list_taxa ) ) ):
        next_status, next_taxa = statuslist[i], taxalist[i]
        next_multiple, next_taxaset, next_elimset = get_taxaset( next_taxa, next_status )
        match_count = compare_taxasets( lowest_taxaset, next_taxaset, lowest_status, next_status, i, next_multiple )
        match_list.append( match_count )
        multiple_list.append( next_multiple )
    # determine consistency
    unmatched_list = []
    unknown_list = []
    if sum( match_list ) == c__list_taxa.index( 's' ) + 1:
        consistency = True
    else:
        consistency = False
        for i in reversed( range( len( match_list ) ) ):
            callvalue, calllevel = match_list[i], cr__list_taxa[i]
            if callvalue != 1:
                unmatched_list.append( calllevel )
    # determine where taxa is unknown
    if True not in multiple_list:
        unknown = False
    else:
        unknown = True
        for j in reversed( range( len( multiple_list ) ) ):
            if multiple_list[j] == True:
                unknown_list.append( cr__list_taxa[ j ] )
    return consistency, unmatched_list, unknown, unknown_list
    
def process_lgtstatus( statuslist, taxalist ):
    """
    This function calculates how certain we are about the LGT status.
    In theory, all contigs, when progressing from kingdom to species, should go from NoLGT to LGT.
    We assume only the NoLGT and LGT states are certain, while ambiguous states actually go either way. Thus,
        1. A contig will be assigned LGT if "LGT" appears.
        2. A contig will be assigned NoLGT if "LGT" never appears, but "NoLGT" does.
        3. If only ambiguous statuses appear, it will be assigned whatever one it was.
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

def get_taxalevel( status, taxonomy ):
    """
    Get taxa level (last taxa in the taxonomy will be the level).
    """
    if status == 'LGT' or status == 'ambiguous-LGT':
        taxasplit_bugone, taxasplit_bugtwo = taxonomy.split('-')[0], taxonomy.split('-')[1]
        taxasplit_bug1level, taxasplit_bug2level = taxasplit_bugone.split('|'), taxasplit_bugtwo.split('|')
        taxalevel = taxasplit_bug1level[ len(taxasplit_bug1level)-1 ].split('_')[0]
    else:
        taxasplit = taxonomy.split('|')
        taxalevel = taxasplit[len(taxasplit)-1].split('_')[0]
    return taxalevel

def sort_by_taxalevels( infolist ):
    """
    Get taxa for a single contig in sorted order.
    """
    sortedlist = [None]*len(infolist)
    for info in infolist:
        taxalevel = get_taxalevel( info[2], info[5] )
        index = c__list_taxa.index( taxalevel )
        sortedlist[index] = info
    return sortedlist

def get_topuniref( gene_uniref ):
    """
    Rank Unirefs if more than 1 with same count.
    """
    potential_unirefs = gene_uniref.split(',')
    if len( potential_unirefs ) == 1:
        topuniref = potential_unirefs[0].strip().split(':')[0]
    else:
        within_gene_list = []
        for i in range( len( potential_unirefs )-1 ):
            curr_unirefname, curr_unirefvalue = potential_unirefs[i].strip().split(':')
            next_unirefname, next_unirefvalue = potential_unirefs[i+1].strip().split(':')
            if float(next_unirefvalue) < float( curr_unirefvalue ):
                if curr_unirefname == 'unknown' and i == len( potential_unirefs )-2:
                    within_gene_list.append( next_unirefname )
                else:
                    if curr_unirefname not in within_gene_list:
                        within_gene_list.append( curr_unirefname )
                break
            else:
                if curr_unirefname not in within_gene_list:
                    within_gene_list.append( curr_unirefname )
                within_gene_list.append( next_unirefname )
        topuniref = '|'.join( within_gene_list )
    return topuniref

def format_uniref( myuniref ):
    """
    Split Unirefs into print-able format.
    """
    final_uniref = ''
    if len( myuniref.split('|') ) > 1:
        # if there is more than 1 gene
        topunireflist = []
        for gene_uniref in myuniref.split('|'):
            topuniref = get_topuniref( gene_uniref )
            topunireflist.append( topuniref )
        final_uniref = ';'.join( topunireflist )
    else:
        final_uniref = get_topuniref( myuniref.split('|')[0] )
    return final_uniref

def get_direction( taxaone, taxatwo, rd, taxalevel ):
    """
    Get directionality between two taxa if recipient/donor was called.
    """
    direction = ''
    if rd == 'NA':
        direction = '?'
    else:
        recipient, donor = rd.split('-')
        recip_short, donor_short = rd.split('|')[ c__list_taxa.index( taxalevel ) ], rd.split('|')[ c__list_taxa.index( taxalevel ) ]
        rd_list = [recip_short, donor_short]
        if ( len( set([ taxaone, taxatwo ] ) & set( rd_list ) ) == 2 ) or ( taxaone in rd_list and re.search( '\|', taxatwo ) ):
            recipient, donor = rd.split('-')
            if taxaone == recipient:
                direction = '<'
            else:
                direction = '>'
        else:
            direction = '?'
    return direction

def format_output( status, taxa, rd, synteny, multtaxa ):
    newsynteny = ''
    format_taxa = ''
    if status == 'LGT' or status == 'ambiguous-LGT':
        multiple_taxa, taxaset, elimset = get_taxaset( taxa, status )
        taxalevel = get_taxalevel( status, taxa.split(';')[0] )
        if multiple_taxa == True: #multiple taxa
            taxapair_list, rd_list, synteny_list = taxa.split(';'), rd.split(';'), synteny.split(';')
            direction = ''
            if len( list( elimset ) ) == 1: #we know one taxa
                taxaone = list( elimset )[0].split('|')[ c__list_taxa.index( taxalevel ) ]
                taxatwo = taxalevel + '__multiple'
                if multtaxa == True: #determine taxatwo
                    remain_taxa = list( taxaset - elimset )
                    remain_list = []
                    for r_taxa in remain_taxa:
                        r2_taxa = r_taxa.split('|')[ c__list_taxa.index( taxalevel ) ]
                        remain_list.append( r2_taxa )
                    taxatwo = '|'.join( remain_list )                
                for i in range( len( taxapair_list ) ): #get direction
                    taxapair = taxapair_list[i]
                    taxapairs = taxapair.split('-')
                    if list( elimset )[0] in taxapairs: 
                        direction = get_direction( taxaone, taxatwo, rd_list[i], taxalevel )
                        if taxapairs[0].split('|')[c__list_taxa.index( taxalevel )] == taxaone:
                            newsynteny = synteny_list[i]
                        else:
                            newsynteny_one = synteny_list[i].replace('A', 'C')
                            newsynteny_two = synteny_list[i].replace('B', 'A')
                            newsynteny = synteny_list[i].replace('C', 'B')
                        if direction != '?':
                            break
                format_taxa = taxaone + direction + taxatwo
                
            else: #we don't know either taxa
                taxaone = taxalevel + '__multiple'
                taxatwo = taxalevel + '__multiple'
                format_taxa = taxaone + direction + taxatwo
                newsynteny = 'unknown'
                direction = '?'
                if multtaxa == True:
                    pairlist = []
                    for i in range( len( taxapair_list ) ):
                        tpair = taxapair_list[i]
                        tpair_one, tpair_two = tpair.split('-')
                        tpair_oneshort, tpair_twoshort = tpair_one.split('|')[c__list_taxa.index( taxalevel )], tpair_two.split('|')[c__list_taxa.index( taxalevel )]
                        direction = get_direction( tpair_oneshort, tpair_twoshort, rd_list[i], taxalevel )
                        pairs = direction.join( [tpair_oneshort, tpair_twoshort] )
                        pairlist.append( pairs )
                    format_taxa = '|'.join( pairlist )
                    newsynteny = '|'.join( synteny_list )
        else: #we know both taxa in lgt
            taxaone, taxatwo = taxa.split('-')[0], taxa.split('-')[1]
            taxaonelist, taxatwolist = taxaone.split('|'), taxatwo.split('|')
            known_taxaone, known_taxatwo = taxaonelist[ len( taxaonelist )- 1 ], taxatwolist[ len(taxatwolist)-1 ]
            if rd != 'NA':
                recipient, donor = rd.split('-')[0], rd.split('-')[1]
                if taxaone == recipient:
                    format_taxa = known_taxaone + '<' + known_taxatwo
                else:
                    format_taxa = known_taxaone + '>' + known_taxatwo
            else:
                format_taxa = known_taxaone + '?' + known_taxatwo
            newsynteny = synteny
    else: #no lgt
        multiple_taxa, taxaset, elimset = get_taxaset( taxa, status )
        taxalevel = get_taxalevel( status, taxa.split(';')[0] )
        if multiple_taxa == True:
            taxalist = list( taxaset )[0].split('|')
            format_taxa = taxalevel + '__multiple'
            newsynteny = synteny.split(';')[0]
            if multtaxa == True:
                singlelist = []
                for singletaxa in list( taxaset ):
                    single_formattaxa = singletaxa.split('|')[c__list_taxa.index( taxalevel )]
                    singlelist.append( single_formattaxa )
                format_taxa = '|'.join( singlelist )
        else:
            taxalist = taxa.split('|')
            format_taxa = taxalist[ c__list_taxa.index( taxalevel ) ]
            newsynteny = synteny
    return format_taxa, newsynteny  
    
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
            dict_contigs.setdefault( contigname, [] ).append( aastrline )

    # print header
    print( '\t'.join( ['contig', 'length', 'uniref50', 'uniref90', 'preferred_call', 'unambiguous_call', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']) )

    # for each contig in the results, sort lgtscorer output by phylogenetic lvl
    for contig in dict_contigs.keys():
        infolist = dict_contigs[contig]
        sortedlist = sort_by_taxalevels( infolist ) 

        # get the status and taxa in order
        statuslist, scorelist, taxalist = [], [], []
        for info in sortedlist:
            scoreinfo = wu.Scores( info )
            statuslist.append( scoreinfo.status )
            taxalist.append( scoreinfo.taxa )

        # find lgt status
        start, end, status = process_lgtstatus( statuslist, taxalist )
        if status == 'LGT' or status == 'ambiguous-LGT':
            finalstatus = status + ':' + start
        else:
            finalstatus = status

        # find taxa consistency
        consistency, c_taxalist, unknown, u_taxalist = process_taxastatus( start, end, taxalist, status, statuslist )
        
        # make final and unambiguous calls
        finalinfo = wu.Scores( sortedlist[ c__list_taxa.index( end ) ] )
        final_taxa,final_synteny = format_output( finalinfo.status, finalinfo.taxa, finalinfo.rd, finalinfo.synteny, False )
        
        final_annot = ';'.join( str(x) for x in [finalstatus, format( finalinfo.onescore, '.4f'), format( finalinfo.twoscore, '.4f' ), final_taxa, final_synteny] )
        unambig_annot = ''
        if unknown == False: # final and unambiguous calls are the same
            unambig_annot = final_annot
        else: # unambiguous call should contain last known taxa
            lastknown_taxa = min( c__list_taxa.index( u_taxalist[0] ), c__list_taxa.index( end ) ) 
            knowninfo = wu.Scores( sortedlist[ lastknown_taxa ] )
            if knowninfo.status == 'LGT':
                modstatus = knowninfo.status + ':' + c__list_taxa[lastknown_taxa]
            else:
                modstatus = knowninfo.status
            unambig_taxa, unambig_synteny = format_output( knowninfo.status, knowninfo.taxa, knowninfo.rd, knowninfo.synteny, False )
            unambig_annot = ';'.join( str(x) for x in [modstatus, format( knowninfo.onescore, '.4f'), format( knowninfo.twoscore, '.4f'), unambig_taxa, unambig_synteny] )

        # get remaining info
        infolist = []
        for i in range( len( sortedlist ) ):
            info = wu.Scores( sortedlist[i] )
            info_taxa, info_synteny = format_output( info.status, info.taxa, info.rd, info.synteny, True )
            info_annot = ';'.join( str(x) for x in [info.status, format( info.onescore, '.4f'), format( info.twoscore, '.4f'), info_taxa, info_synteny] )
            if consistency == False and (c__list_taxa[i] in c_taxalist):
                info_annot = info_annot + '*'
            infolist.append( info_annot )

        # print everything
        final_uniref50, final_uniref90 = format_uniref( scoreinfo.uniref50 ), format_uniref( scoreinfo.uniref90 )
        final_line = [scoreinfo.contig, str( scoreinfo.length ), final_uniref50, final_uniref90, final_annot, unambig_annot] + infolist
        print( '\t'.join( final_line ) )
        
if __name__ == "__main__":
    main()


