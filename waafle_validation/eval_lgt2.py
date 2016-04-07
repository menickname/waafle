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

def call_orgs( ans_status, ans_reciplevel, ans_donorlevel, waafle_status, waafle_orgs ):
    # call orgs for ans
    ans_orgpairs = set([ans_reciplevel + '-' + ans_donorlevel, ans_donorlevel + '-' + ans_reciplevel])
    ans_pairs = set([ans_reciplevel + '-' + ans_donorlevel])
    ans_orgs = set([ans_reciplevel, ans_donorlevel])
    
    # call orgs for waafle
    if waafle_status == 'LGT' or waafle_status == 'ambiguous-LGT':
        waafle_orgpairs = set( waafle_orgs.split(';') )
        waafle_orgs = set([])
        for x in waafle_orgpairs:
            org1, org2 = x.split('-')
            waafle_orgs.add( org1 )
            waafle_orgs.add( org2 )
        if ans_status == 'LGT':
            if len( ans_orgpairs & waafle_orgpairs ) >= 1:
                status = ['Match', ans_pairs, waafle_orgpairs]
            else:
                status = ['No_Match', ans_pairs, waafle_orgpairs]
        else:
            if len( ans_orgs & waafle_orgs ) >= 1:
                status = ['Partial_Match', ans_orgs, waafle_orgpairs]
            else:
                status = ['No_Match', ans_orgs, waafle_orgpairs]
    else:
        waafle_orgs = set( waafle_orgs.split(';') )
        if ans_status == 'NoLGT':
            if len( ans_orgs & waafle_orgs ) >= 1: 
                status = ['Match', ans_orgs, waafle_orgs]
            else:
                status = ['No_Match', ans_orgs, waafle_orgs]
        else:
            if len( ans_orgs & waafle_orgs ) >= 1:
                status = ['Partial_Match', ans_pairs, waafle_orgs]
            else:
                status = ['No_Match', ans_pairs, waafle_orgs]
    return status

def call_dr( ans_recip, ans_donor, waafle_recip, waafle_donor ):
    waafle_pair = waafle_recip + '-' + waafle_donor
    if ans_recip == waafle_recip and ans_donor == waafle_donor:
        RD_call = 'Match_RD'
    elif ans_recip == waafle_recip and ans_donor != waafle_donor:
        RD_call = 'Match_R'
    elif ans_recip != waafle_recip and ans_donor == waafle_donor:
        RD_call = 'Match_D'
    else:
        RD_call = 'No_Match'
    return RD_call, waafle_pair

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

    set_waafle_contigs = set( [] )
    for bstrline in open( args.scoredcontigs ):
        bbstrline = bstrline.strip().split('\t')
        if bbstrline[1] == 'contiglen':
            continue
        contig, waafle_status = bbstrline[0], bbstrline[2]
        set_waafle_contigs.add( contig )
        LGTcall = ''
    
        #call status  
        ans_taxalevel = phylevels.index( dict_answerkey[contig][0] )
        waafle_taxalevel = phylevels.index( args.scoredcontigs[args.scoredcontigs.rindex('_')+1] )
        if ans_taxalevel <= waafle_taxalevel:
            ans_status = 'LGT'
        else:
            ans_status = 'NoLGT'
        LGTcall = call_status( waafle_status, ans_status )

        #call orgs
        ans_recip, ans_donor = dict_answerkey[contig][1], dict_answerkey[contig][2]
        ans_reciplevel = ans_recip.split('|')[waafle_taxalevel]
        ans_donorlevel = ans_donor.split('|')[waafle_taxalevel]
        callstatus = call_orgs( ans_status, ans_reciplevel, ans_donorlevel, waafle_status, bbstrline[4] )
        orgcall, ansorgs, waafleorgs = callstatus[0], ';'.join( callstatus[1] ), ';'.join( callstatus[2] ) 
        
        #call donor-recip
        waafle_recip, waafle_donor = bbstrline[5], bbstrline[6]
        if waafle_recip == 'NA' or waafle_donor == 'NA':
            RD_call = 'No_RDCall'
        else:
            RD_call, waafleorgs = call_dr( ans_reciplevel, ans_donorlevel, waafle_recip, waafle_donor )
        print('\t'.join( [contig, LGTcall, ans_status, waafle_status, orgcall, RD_call, ansorgs, waafleorgs] ) )
        
    for contig in dict_answerkey.keys():
        if contig not in set_waafle_contigs:
            ans_taxalevel = phylevels.index( dict_answerkey[contig][0] )
            waafle_taxalevel = phylevels.index( args.scoredcontigs[args.scoredcontigs.rindex('_')+1] )
            ans_status, LGTcall = "", ""
            if ans_taxalevel <= waafle_taxalevel:
                ans_status = 'LGT'
                LGTcall = 'FN-nogenes'
            else:
                ans_status = 'NoLGT'
                LGTcall = 'TN-nogenes'
            
            print( '\t'.join( [contig, LGTcall, ans_status, "No-genes-called", "NA", "NA", "NA", "NA"] ) )

if __name__ == "__main__":
    main()

