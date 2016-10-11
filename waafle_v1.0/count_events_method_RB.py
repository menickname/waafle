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
    

def add_lgtcounts( bug1name, bug2name, twosyn, direction, dict_lgt ):
    num_bug1 = twosyn.count('A') + 0.5*float( twosyn.count('C') )
    num_bug2 = twosyn.count('B') + 0.5*float( twosyn.count('C') )
    if direction == '>': #bug1 (A) is donor, bug2 (B) is recip
        dict_lgt = add_counts( bug2name, bug1name, num_bug1, dict_lgt )
        dict_lgt = add_counts( bug2name, bug2name, num_bug2, dict_lgt )
    elif direction == '<': #A is recip, B is donor
        dict_lgt = add_counts( bug1name, bug2name, num_bug2, dict_lgt )
        dict_lgt = add_counts( bug1name, bug1name, num_bug1, dict_lgt )
    else: #directionality not known
        like_bug1 = num_bug1/( float( num_bug1 ) + float( num_bug2 ) )
        like_bug2 = num_bug2/( float( num_bug1 ) + float( num_bug2 ) )
        dict_lgt = add_counts( bug2name, bug1name, num_bug1*like_bug1, dict_lgt )
        dict_lgt = add_counts( bug1name, bug2name, num_bug2*like_bug2, dict_lgt )
        dict_lgt = add_counts( bug2name, bug2name, num_bug2*like_bug2, dict_lgt )
        dict_lgt = add_counts( bug1name, bug1name, num_bug1*like_bug1, dict_lgt )
    return dict_lgt
    
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
                    bug1_abbr = get_name( i, bug1 )
                    bug2_abbr = get_name( i, bug2 )
                    
                    # fix names if needed
                    if i >= wu.c__list_taxa.index( lvl ): #lgt

                        #assign counts to lgt dict
                        dict_lgt = add_lgtcounts( bug1_abbr, bug2_abbr, result.twosyn, direction, dict_lgt )
                        dict_levelbugs.setdefault( wu.c__list_taxa[i], set([]) ).add( bug1_abbr )
                        dict_levelbugs.setdefault( wu.c__list_taxa[i], set([]) ).add( bug2_abbr )

                    else: #no lgt, at this point bug1 and bug2 are identical
                        # assign counts to nolgt dict
                        dict_lgt = add_counts( bug1_abbr, bug2_abbr, num_bug1 + num_bug2, dict_lgt )
                        dict_levelbugs.setdefault( wu.c__list_taxa[i], set([]) ).add( bug1_abbr )

                    #sum all genes            
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

    #get matrix of values
    for level in dict_levelbugs.keys():
        totalmatrix, headers = [], ['bugs']
        fh = open( level + '_table.txt', 'w' )
        
        buglist = list( dict_levelbugs[ level ] )
        for i in range( len( buglist ) ): #get rows 
            bug1 = buglist[i]
            total_genes = sum( dict_genes[ level ] )
            headers.append( bug1 )
            line = [ bug1 ]

            for j in range( len( buglist ) ): #get columns
                bug2 = buglist[j]
                value = 0

                if bug1 not in dict_lgt.keys(): #was only ever a donor
                    #line.append( str(value) + '/' + str(total_genes ) )
                    line.append( float(value)/total_genes )
                    print( '\t'.join( [bug1 + '<' + bug2, str(float(value)/total_genes)] ) )
                else:
                    if bug2 in dict_lgt[bug1]:
                        value = dict_lgt[bug1][bug2]
                    #line.append( str(value) + '/' + str(total_genes) ) 
                    line.append( float(value)/total_genes )
                    print( '\t'.join( [bug1 + '<' + bug2, str(float(value)/total_genes)] ) )
            totalmatrix.append( line )

        fh.write( '\t'.join( headers ) + '\n' )
        #print( '\t'.join( headers ) )

        for astrline in totalmatrix:
            fh.write( '\t'.join( str(x) for x in astrline ) + '\n' )
            #print( '\t'.join( str(x) for x in astrline ) )
        fh.close()
 
if __name__ == "__main__":
    main()
