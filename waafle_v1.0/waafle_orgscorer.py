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
        default=0.8,
        type=float,
        help="onebug score cutoff"
        )
    parser.add_argument(
        "-s2", "--twobug",
        default=0.8,
        type=float,
        help="twobug score cutoff"
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


def calc_onebug( table, buglist ):
    """ Calculates the one bug score, which equals the maximum of 
    the minimum of scores across all bugs per gene.
    High scores indicate that contig is likely explained by one taxa. """
    onebugscore = np.max( np.amin( table, axis=1 ) )
    indices = [i for i,j in enumerate( list(np.amin(table, axis=1) ) ) if j==onebugscore]
    bugmeanlist = []
    for index in indices:
        bugname = buglist[index]
        bugmean = np.mean(table[index])
        bugmeanlist.append( [bugname, bugmean] )
    bugmeanlist_sorted = sorted( bugmeanlist, key=itemgetter(1), reverse=True )
    topmean = bugmeanlist_sorted[0][1]
    onebuglist, onebugsyn = [], []
    for bug, avg in bugmeanlist_sorted:
        if avg >= topmean:
            onebuglist.append( bug )
            element_list, taxastring, shift_counter = find_taxa_order( 'NoLGT', table, buglist, bug )
            onebugsyn.append( taxastring )
    return onebugscore, onebuglist, onebugsyn


def calc_twobug( table, buglist, twobugthresh ):
    """ Calculates the complement score, which equals the minimum of 
    the maximum of scores across all genes per pair of bugs.
    High scores indicate that contig is likely explained by two taxa. """
    complementlist = [] 
    twobugscore, twobuglist, twobugsyn = 0, [], []
    for i in range( table.shape[0]-1 ):
        bug1 = table[i, :]
        if max( bug1 ) < twobugthresh:
            continue
        for j in range( i+1, table.shape[0] ):
            bug2 = table[j, :]
            if max( bug2 ) < twobugthresh:
                continue
            twobugarray = np.array( [bug1, bug2] )
            complement = np.min(np.amax(twobugarray, axis=0))
            average = np.mean(np.amax(twobugarray, axis=0))
            bugsort = sorted( [str(buglist[i]), str(buglist[j])] )
            bugnames = bugsort[0] + '-' + bugsort[1]
            complementlist.append( [ bugnames, complement, average ] )
    if len( complementlist ) == 0: #all taxa pairs filtered out, no lgt possible
        twobugscore, twobuglist, twobugsyn = 0, ['-'], ['-']
    else:
        complement_sorted = sorted( complementlist, key=itemgetter(1,2), reverse=True )
        twobugscore, twobugavg = complement_sorted[0][1], complement_sorted[0][2]
        twobuglist, twobugsyn = [], []
        for info in complement_sorted:
            element_list, taxastring, shift = find_taxa_order( 'LGT', table, buglist, info[0] )
            lgt_status = False
            reduced_taxastring = taxastring.replace( 'C', '')
            if len( set( reduced_taxastring ) ) == 2:
                lgt_status = True
            if info[1] == twobugscore and info[2] >= twobugavg and lgt_status == True:
                newname = calc_donorrecip( info[0], element_list, taxastring, shift )
                twobuglist.append( newname )
                twobugsyn.append( taxastring )
            elif info[1] == twobugscore and info[2] >= twobugavg and lgt_status == False:
                # actually onebug explains 2 bug score
                twobugscore, twobuglist, twobugsyn = 0, ['-'], ['-']
                break
            elif info[1] < twobugscore or (info[1] == twobugscore and info[2] < twobugavg):
                break
    return twobugscore, twobuglist, twobugsyn


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

 
def det_multiple( bugscore, buglist, bugsyn, num ):
    """ Rename taxa if there are multiple """
    newscore, newbug, newsyn = bugscore, buglist[0], bugsyn[0]
    if len( buglist ) > 1:
        if num == 1:
            newbug = 'multiple'
            newsyn = bugsyn[0] #'-'
        else:
            conservedset, allset, conservedsyn = set([]), set([]), set([])
            for i in range( len( buglist ) ):
                pair = buglist[i]
                synteny = bugsyn[i]
                bug1, bug2 = re.split('[><\?-]', pair )
                if i == 0:
                    conservedset |= set( [bug1, bug2] )
                else:
                    conservedset &= set([bug1, bug2])
                allset |= set( [bug1, bug2] )
                conservedsyn |= set( [synteny] )
            if len( conservedset ) == 1 and len( conservedsyn ) <= 2:
                conservedbug = list( conservedset )[0] 
                otherset = allset - conservedset
                for i in reversed( range( len( conservedbug.split('|') ) ) ):
                    levelset = set([])
                    for bug in otherset:
                        levelset.add( bug.split('|')[i] )
                    if len( levelset ) == 1:
                        conservedlevel = list( levelset )[0]
                        conservedother = list(otherset)[0].split('|')[:i + 1]
                        if conservedbug.split('|')[i] == conservedlevel:
                            multiplebug = conservedother + [ wu.c__list_taxa[i+1] + '__multiple' ]
                            multiplebugname = '|'.join( multiplebug )
                        else:
                            multiplebugname = '|'.join( conservedother )
                        newbuglist = sorted( [conservedbug] + [multiplebugname] )
                        bug1, bug2 = re.split('[><\?-]', pair )
                        sign = re.search( '[><\?-]', pair ).group()
                        if bug1 == conservedbug:
                            if bug1 == newbuglist[0]:
                                newsyn = bugsyn[0]
                            else:
                                sign, newsyn = flip_annot( sign, bugsyn[0] )
                        else: 
                            if bug2 == newbuglist[0]:
                                sign, newsyn = flip_annot( sign, bugsyn[0] )
                            else: 
                                newsyn = bugsyn[0]
                        newbug = newbuglist[0] + sign + newbuglist[1]
                        break
                    if i == 0:
                        newbug = 'multiple'
                        newsyn = '-'
            else:
                newbug = 'multiple'
                newsyn = '-'
    return newscore, newbug, newsyn


def det_status( matrix_u, taxaorder_u, onebug_thresh, twobug_thresh ):
    """ Find whether the status is LGT or not """
    numbugs, numgenes = np.shape( matrix_u )
    onebugscore, onebuglist, onebugsyn = calc_onebug( matrix_u, taxaorder_u )
    found, status = False, 'Ambiguous'
    onebug, onesyn, twobugscore, twobug, twosyn = '', '', 0, '', ''
    if numgenes == 1 or numbugs == 1: #only 1 gene or 1 bug
        onebugscore, onebug, onesyn = det_multiple( onebugscore, onebuglist, onebugsyn, 1 )
        twobugscore, twobug, twosyn = 0, '-', '-'
        if onebug != 'multiple':
            found, status = True, 'NoLGT'
    else:
        twobugscore, twobuglist, twobugsyn = calc_twobug( matrix_u, taxaorder_u, twobug_thresh )
        onebugscore, onebug, onesyn = det_multiple( onebugscore, onebuglist, onebugsyn, 1 )
        twobugscore, twobug, twosyn = det_multiple( twobugscore, twobuglist, twobugsyn, 2 )
        if onebugscore >= onebug_thresh:
            if onebug != 'multiple':
                found, status = True, 'NoLGT'
        else:
            if twobugscore >= twobug_thresh:
                if not re.search( 'multiple', twobug ):
                    found, status = True, 'LGT'
    return onebugscore, onebug, onesyn, twobugscore, twobug, twosyn, found, status


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
    for i in range( lowest_known ):
        level = taxa1list[i].split('__')[0]
        if taxa1list[i] != taxa2list[i]:
            break
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
            
            #Scoring onebug/twobug scores at the species level.
            #Do not stop iteration until
            #1: Clear annotation (LGT/NoLGT)
            #2: One taxon/taxa-pair          
            for taxalevel in reversed( wu.c__list_taxa ): 
                #Generate table for contig
                matrix, taxaorder = generate_matrix( dict_genearray, taxalevel )
                matrix_u, taxaorder_u = spike_unknown( matrix, taxaorder, args.unknown )
                print( contig )
                onebugscore, onebug, onebugsyn, twobugscore, twobug, twobugsyn, found, status = det_status( matrix_u, taxaorder_u, args.onebug, args.twobug )
                if found == True:
                    break

            #Print results
            #Format uniref
            uniref90, uniref50 = find_top_uniref( uniref90list), find_top_uniref( uniref50list )
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


