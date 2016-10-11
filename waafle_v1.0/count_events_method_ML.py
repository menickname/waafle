#!/usr/bin/python

"""
The goal of this script is to:
1) count the lgt
2) normalize by gene number from both taxa
3) output results for merging table

This "count" is Method #1 from brainstorming.
Morgan Langille decided we should count # of LGT events per gene.

This makes it possible to take into account broken up contigs.
Consider the scenario below:

ABA
AAAAAAAAAAAAAAAAAAA

# LGT/ contig = 1/2
# LGT/ genes = 1/22

vs.

ABAAAAAAAAAAAAAAAAAAAA

# LGT/ contig = 1/1
# LGT/ genes = 1/22

vs.

AB
AAAAAAAAAAAAAAAAAAAA

# LGT/ contig = 1/1
# LGT/ genes = 1/22

If we chose to count # of events per contig, these 3 samples would
have different counts. 
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
        "-u", "--unknown",
        action="store_true",
        help="Count lgt at levels where taxa is unknown",
        )
    parser.add_argument(
        "-o", "--output",
        default="count.txt",
        help="output file for counts",
        )
    parser.add_argument(
        "-n", "--name",
        help="name of file",
        )
    args = parser.parse_args()
    return args

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
    1) Read in results file
    2) Count event if lgt, don't if no lgt
    3) Divide # of events by # of genes
    4) Output 
    """
    args = get_args()
    
    dict_count, dict_genes, dict_dir = {}, {}, {}
    
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
                
                # count lgt events
                if args.unknown == True:
                    most_known = len( wu.c__list_taxa )
                else:
                    most_known = min( len( bug1.split('|') ), len( bug2.split('|') ) )
                for i in range( most_known ):
                    bug1_abbr = get_name( i, bug1 )
                    bug2_abbr = get_name( i, bug2 ) 
                    if i >= wu.c__list_taxa.index( lvl ): #lgt

                        #assign event counts and genes
                        dict_count.setdefault( i, [] ).append( 1 )
                        dict_genes.setdefault( i, [] ).append( num_bug1 + num_bug2 )
                        if direction == '?':
                            dict_dir.setdefault( i, [] ).append( 0 )
                        else:
                            dict_dir.setdefault( i, [] ).append( 1 )

                    else: #no lgt at this level
                        # assign event counts and genes 
                        dict_count.setdefault( i, [] ).append( 0 )
                        dict_genes.setdefault( i, [] ).append( num_bug1 + num_bug2 )
                        dict_dir.setdefault( i, [] ).append( 0 )
            
            elif re.search( 'NoLGT', result.status ): # no lgt
                bug = result.onebug
                num_bug = result.onesyn.count('A')
                if args.unknown == True:
                    most_known = len( wu.c__list_taxa )
                else:
                    most_known = len( bug.split('|') )
                for i in range( most_known ):
                    bug_abbr = get_name( i, bug )
                    dict_count.setdefault( i, [] ).append( 0 )
                    dict_genes.setdefault( i, [] ).append( num_bug )
                    dict_dir.setdefault( i, [] ).append( 0 )
    
            else: #ignore 'ambiguous'
                pass
    
    # get all taxa at each taxonomic level along with total number of genes
    if args.name:
        name = args.name + '_' + args.output 
    else:
        name = args.output
    fh = wu.try_open( name, "w" )
    writer = csv.writer( fh, dialect="excel-tab" )

    for i in range( len( wu.c__list_taxa ) ):
        level = wu.c__list_taxa[i]

        # get lgt count and dir count
        count = 0
        dir = 0
        dir_rate = 0
        if i in dict_count.keys(): # false if no blast results
            count = sum( dict_count[i] )
            dir = sum( dict_dir[i] )
            if count != 0: # 0 if no lgt at level
                dir_rate = dir/float(count)

        # get genes
        genes = 0
        rate = 0
        contigs = 0
        if i in dict_genes.keys(): #false if no blast results
            genes = sum( dict_genes[i] )
            rate = sum( dict_count[i] )/float( sum( dict_genes[i] ) )
            contigs = len( dict_genes[i] )

        writer.writerow( [level + '__count', count] )
        writer.writerow( [level + '__genes', genes] )
        writer.writerow( [level + '__dir', dir] )
        writer.writerow( [level + '__rate', rate] )
        writer.writerow( [level + '__dirrate', dir_rate] )
        writer.writerow( [level + '__contigs', contigs] )
        
if __name__ == "__main__":
    main()
