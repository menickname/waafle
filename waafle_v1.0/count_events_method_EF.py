#!/usr/bin/python

"""
The goal of this script is to:
1) count the lgt
2) normalize by gene number from both taxa
3) output table
"""

from __future__ import print_function # python 2.7+ required
import os, sys, csv, re, argparse
import waafle_utils as wu
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
        "-i", "--input",
        required=True,
        help="output from waafle_lgtcaller"
        )
    parser.add_argument(
        "-d", "--dir",
        action='store_true',
        help="include contigs where directionality is known"
        )
    parser.add_argument(
        "-f", "--frac",
        action="store_true",
        help="print ints rather than fractions"
        )
    parser.add_argument(
        "-n", "--name",
        help="name of sample"
        )
    args = parser.parse_args()
    return args

def count_each_lvl( bugname, count, dict_genes ):
    for i in range( len( bugname.split('|') ) ):
        buglvl = '|'.join( bugname.split('|')[:i+1] )
        if buglvl in dict_genes.keys():
            dict_genes[ buglvl ] = dict_genes[ buglvl ] + count
        else:
            dict_genes[ buglvl ] = count
    return dict_genes

def add_counts( bug1name, bug2name, count, dict_lgtgenes ):
    if bug1name in dict_lgtgenes.keys(): #if bug1 is in lgt dict
        bug1dict = dict_lgtgenes[ bug1name ]
        if bug2name in bug1dict.keys():
            bug1dict[ bug2name ] = count + bug1dict[ bug2name ]
        else:
            bug1dict[ bug2name ] = count
        dict_lgtgenes[ bug1name ] = bug1dict
    else:
        dict_lgtgenes[ bug1name ] = { bug2name: count }
    return dict_lgtgenes
    

def assign_count( dict_count, bug, count ):
    if bug in dict_count.keys():
        dict_count[bug] = dict_count[bug] + count
    else:
        dict_count[bug] = count
    return dict_count
    

def get_name( index, bug ):
    name, namelist = '', []
    if len( bug.split('|') ) <= index: #don't know taxa
        for i in range( index + 1 ):
            level = wu.c__list_taxa[i]
            if i >= len( bug.split('|') ):
                name = level + '__Unk'
                namelist.append( name )
            else:
                name = bug.split('|')[i]
                namelist.append( name )
        name = '|'.join( namelist )
    else: #know taxa
        name = '|'.join( bug.split('|')[:index+1] )
    return name

# ---------------------------------------------------------------
# main
# ---------------------------------------------------------------

def main():
    """
    This script does the following:
    1) 
    """
    args = get_args()
    
    dict_lgt, dict_genes, dict_levelbugs = {}, {}, {}
    for astrline in open( args.input ):
        aastrline = astrline.strip().split('\t')
        if aastrline[0] != 'contig': #skip header
            result = wu.Result( aastrline )
            if re.search( 'LGT:', result.status ): #if there is lgt
                bug1, bug2 = re.split( '[><\?]', result.twobug )
                direction = re.search( '[><\?]', result.twobug ).group()
                num_bug1 = result.twosyn.count('A') + float( result.twosyn.count('C') )*0.5
                num_bug2 = result.twosyn.count('B') + float( result.twosyn.count('C') )*0.5
                lvl = result.status.split(':')[1]
                
                # count lgt genes
                most_known = min( len( bug1.split('|') ), len( bug2.split('|') ) )
                for i in range( len( wu.c__list_taxa ) ):
                    
                    # fix names if needed
                    bug1_abbr = get_name( i, bug1 )
                    bug2_abbr = get_name( i, bug2 )
                    
                    if i >= wu.c__list_taxa.index( lvl ): #lgt

                        #assign counts to lgt dict
                        if args.dir == True: #dir doesn't matter
                            dict_lgt = add_counts( bug1_abbr, bug2_abbr, num_bug2, dict_lgt )
                            dict_lgt = add_counts( bug2_abbr, bug1_abbr, num_bug1, dict_lgt )
                            dict_lgt = add_counts( bug1_abbr, bug1_abbr, num_bug1, dict_lgt )
                            dict_lgt = add_counts( bug2_abbr, bug2_abbr, num_bug2, dict_lgt )
                            dict_levelbugs.setdefault( wu.c__list_taxa[i], set([]) ).add( bug1_abbr )
                            dict_levelbugs.setdefault( wu.c__list_taxa[i], set([]) ).add( bug2_abbr )
                            dict_genes.setdefault( wu.c__list_taxa[i], [] ).append( num_bug1 + num_bug2)
                        else:
                            if direction == '?': #only add those with unknown dir
                                dict_lgt = add_counts( bug1_abbr, bug2_abbr, num_bug2, dict_lgt )
                                dict_lgt = add_counts( bug2_abbr, bug1_abbr, num_bug1, dict_lgt )
                                dict_lgt = add_counts( bug1_abbr, bug1_abbr, num_bug1, dict_lgt )
                                dict_lgt = add_counts( bug2_abbr, bug2_abbr, num_bug2, dict_lgt )
                                dict_levelbugs.setdefault( wu.c__list_taxa[i], set([]) ).add( bug1_abbr )
                                dict_levelbugs.setdefault( wu.c__list_taxa[i], set([]) ).add( bug2_abbr )
                                dict_genes.setdefault( wu.c__list_taxa[i], [] ).append( num_bug1 + num_bug2)

                    else: #no lgt, at this point bug1 and bug2 are identical
                        # assign counts to nolgt dict
                        dict_lgt = add_counts( bug1_abbr, bug2_abbr, num_bug1 + num_bug2, dict_lgt )
                        dict_levelbugs.setdefault( wu.c__list_taxa[i], set([]) ).add( bug1_abbr )
                        dict_genes.setdefault( wu.c__list_taxa[i], [] ).append( num_bug1 + num_bug2)

            elif re.search( 'NoLGT', result.status ): # no lgt
                bug = result.onebug
                num_bug = result.onesyn.count('A')
                most_known = len( bug.split('|') )
                for i in range( len( wu.c__list_taxa ) ):
                    bug_abbr = get_name( i, bug ) 
                    dict_lgt = add_counts( bug_abbr, bug_abbr, num_bug, dict_lgt )
                    dict_levelbugs.setdefault( wu.c__list_taxa[i], set([]) ).add( bug_abbr )
                    dict_genes.setdefault( wu.c__list_taxa[i], [] ).append( num_bug )
                    
            else: #ignore 'ambiguous'
                pass

    #get matrix of values, output table
    fh_table = open( args.name + '_table.txt', 'w' ) #get table (all phyl levels)
    dict_levelmatrix = {}
    for level in dict_levelbugs.keys():
        matrix, totalmatrix, headers = [], [], ['bugs']
        fh = open( level + '_matrix.txt', 'w' ) #get matrix
        
        buglist = list( dict_levelbugs[ level ] )
        for i in range( len( buglist ) ): #get matrix rows 
            bug1 = buglist[i]
            total_genes = sum( dict_genes[ level ] )
            headers.append( bug1 )
            line = [ bug1 ]
            for j in range( len( buglist ) ): #get matrix columns
                bug2 = buglist[j]
                value = 0

                if bug2 in dict_lgt[bug1]:
                    value = dict_lgt[bug1][bug2]
                if args.frac == True:
                    line.append( str(value) + '/' + str(total_genes) ) 
                    fh_table.write( '\t'.join( [bug1 + '<' + bug2, str(value) + '/' + str(total_genes)] ) + '\n' )
                else:
                    line.append( float(value)/total_genes )
                    fh_table.write( '\t'.join( str(x) for x in [bug1 + '<' + bug2, float(value)/total_genes] ) + '\n' )
            totalmatrix.append( line )
            matrix.append( line[1:] )
        dict_levelmatrix[level] = [matrix, buglist]

        # write matrix
        fh.write( '\t'.join( headers ) + '\n' )
        for astrline in totalmatrix:
            fh.write( '\t'.join( str(x) for x in astrline ) + '\n' )
            #print( '\t'.join( str(x) for x in astrline ) )
        fh.close()
    fh_table.close()

    #get donorness scores
    fh_donor = open( args.name + '_donorness.txt', 'w' )
    for level in dict_levelmatrix.keys():
        matrix, buglist = dict_levelmatrix[level]
        for i in range( len( buglist ) ):
            bug = buglist[i]
            diag = np.array( matrix )[i, i]
            row = np.sum( np.array( matrix )[i, :] ) - diag
            col = np.sum( np.array( matrix )[:, i] ) - diag
            if row == col == 0:
                continue
            elif col == 0:
                fh_donor.write( '\t'.join( str(x) for x in [bug, 'inf'] ) + '\n' )
            else:
                donorness = float(row)/col
                fh_donor.write( '\t'.join( str(x) for x in [bug, donorness] ) + '\n' )
    
if __name__ == "__main__":
    main()
