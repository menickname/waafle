#!/usr/bin/env python

from __future__ import print_function # python 2.7+ required
import os, sys, csv, argparse
import waafle_utils as wu
from operator import itemgetter, attrgetter, methodcaller
from collections import Counter

# ---------------------------------------------------------------
# constants
# ---------------------------------------------------------------

phylevels = ['k', 'p', 'c', 'o', 'f', 'g', 's', 't']

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
        "-i", "--scoredcontigs",
        required=True,
        help="output from waafle_lgtscorer"
        )
    parser.add_argument(
        "-ans", "--answerkey",
        required=True,
        help="answerkey for synthetic contigs"
        )
    args = parser.parse_args()
    return args

def call_status( waafle_status, ans_status ):
    if waafle_status == ans_status:
        if waafle_status == 'LGT':
            call_status = 'TP'
        elif waafle_status == 'NoLGT':
            call_status = 'TN'
    else:
        if waafle_status == 'LGT':
            call_status = 'FP'
        elif waafle_status == 'NoLGT':
            call_status = 'FN'
        else:
            if ans_status == 'LGT':
                call_status = 'FN-ambig'
            else:
                call_status = 'TN-ambig'
    return call_status
# ---------------------------------------------------------------
# main
# ---------------------------------------------------------------

def main():
    args = get_args()
    dict_answerkey = {}
    for astrline in open( args.answerkey ):
        aastrline = astrline.strip().split('\t')
        taxalevel, contig, recip, donor, = aastrline[0], aastrline[1], aastrline[2], aastrline[3]
        dict_answerkey[contig] = [taxalevel, recip, donor]

    for bstrline in open( args.scoredcontigs ):
        bbstrline = bstrline.strip().split('\t')
        contig, waafle_status = bbstrline[0], bbstrline[1]
        LGTcall = ''
    
        #call status  
        ans_taxalevel = phylevels.index( dict_answerkey[contig][0] )
        waafle_taxalevel = phylevels.index( args.scoredcontigs[args.scoredcontigs.rindex('_')+1] )
        if ans_taxalevel <= waafle_taxalevel:
            ans_status = 'LGT'
        else:
            ans_status = 'NoLGT'
        LGTcall = call_status( waafle_status, ans_status )
        print( ans_taxalevel, waafle_taxalevel, ans_status, waafle_status, LGTcall )
        
        #call orgs
        if waafle_status == 'LGT' or waafle_status == 'ambiguous-LGT': 
            orgpairset = set( bbstrline[2].split('-') )

            # Check donor/recip
            if bbstrline[3] != 'NA' or bbstrline[4] != 'NA':
                recipient = bbstrline[3]
                donor = bbstrline[4]
            
        else:
            orgpairset = set( bbstrline[2] )

if __name__ == "__main__":
    main()

