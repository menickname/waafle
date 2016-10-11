#!/usr/bin/python

"""
Get a list of contigs to print
"""

import sys
import waafle_utils as wu
import argparse

# ---------------------------------------------------------------
# functions
# ---------------------------------------------------------------

def get_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument(
        "-results", "--results",
        required=True,
        help="output from waafle aggregator"
        )
    parser.add_argument(
        "-status", "--status",
        help="get the status you want",
        default="prefcall"
        )
    parser.add_argument(
        "-calltype", "--calltype",
        help="pref or unambig",
        default="LGT"
        )
    parser.add_argument(
        "-listname", "--listname",
        help="name of output file for list of contigs",
        default="contiglist.txt"
        )
    parser.add_argument(
        "-resultgrep", "--resultgrep",
        help="name of output file for contigs chosen",
        default="contigresults.txt"
        )
    args = parser.parse_args()
    return args


args = get_args()
contiglist = open( args.listname, 'w' )
resultgrep = open( args.resultgrep, 'w')

for astrline in open( args.results ):
    aastrline = astrline.strip().split('\t')
    if aastrline[0] != 'contig':
        info = wu.Result( aastrline )
        status, onebug, twobug, taxa, synteny = wu.split_info( info.pref_call )
        if args.status == 'unambig':
            status, onebug, twobug, taxa, synteny = wu.split_info( info.unambig_call )
        if status.split(':')[0] == args.calltype:
            contiglist.write( info.contig + '\n' )
            resultgrep.write( astrline )

contiglist.close()
resultgrep.close()
        
    
