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
    sortedans = sorted( [ans_reciplevel, ans_donorlevel ] )
    ans_orgpairs = set([sortedans[0] + '-' + sortedans[1]])
    ans_pairs = set([ans_reciplevel + '-' + ans_donorlevel])
    ans_orgs = set([ans_reciplevel, ans_donorlevel])
    
    # call orgs for waafle
    if waafle_status == 'LGT' or waafle_status == 'ambiguous-LGT':
        waafle_orgpairs = set( [waafle_orgs] )
        waafle_orgpairs_mult = set( waafle_orgs.split(';') )
        waafle_singleorgs = set([])

        if len( waafle_orgpairs_mult ) > 1:
            #multiple orgs were called
            for x in waafle_orgpairs_mult:
                if len(x.split('-')) == 2:
                    org1, org2 = x.split('-')
                    waafle_singleorgs.add( org1 )
                    waafle_singleorgs.add( org2 )
                else:
                    continue
        else:
            org1, org2 = list(waafle_orgpairs)[0].split('-')
            waafle_singleorgs.add( org1 )
            waafle_singleorgs.add( org2 )
        
        if ans_status == 'LGT':
            # if WAAFLE is LGT/ambig LGT and ans is LGT:
            # get both organisms correct
            if len( ans_orgpairs & waafle_orgpairs ) >= 1:
                status = ['Match', ans_pairs, waafle_orgpairs]
            # get one out of two org pairs correct
            elif len( ans_orgpairs & waafle_orgpairs_mult ) >= 1:
                status = ['Partial-Multiple', ans_pairs, waafle_orgpairs_mult]
            # get no orgs correct
            else:
                status = ['No_Match', ans_pairs, waafle_orgpairs]
        else:
            # if WAAFLE is LGT/ambig LGT and ans is NoLGT:
            # get one out of two orgs correct
            if len( ans_orgs & waafle_singleorgs ) >= 1:
                status = ['Partial-WrongCall', ans_orgs, waafle_orgpairs]
            # get no orgs correct
            else:
                status = ['No_Match-WrongCall', ans_orgs, waafle_orgpairs]
    else:
        waafle_orgs_single = set( [waafle_orgs] )
        waafle_orgs_mult = set( waafle_orgs.split(';') )
        if ans_status == 'NoLGT':
            # if WAAFLE is NoLGT/ambig-NoLGT and ans is NoLGT:
            # get org correct
            if len( ans_orgs & waafle_orgs_single ) >= 1: 
                status = ['Match', ans_orgs, waafle_orgs_single]
            # get org partially right due to multiple
            elif len( ans_orgs & waafle_orgs_mult ) >=1:
                status = ['Partial-Multiple', ans_orgs, waafle_orgs_mult]
            # get org wrong
            else:
                status = ['No_Match', ans_orgs, waafle_orgs_single]
        else:
            # if WAAFLE is NoLGT/ambig-NoLGT and ans is LGT:
            if len( ans_orgs & waafle_orgs_single ) >= 1:
                status = ['Partial-WrongCall', ans_pairs, waafle_orgs_single]
            else:
                status = ['No_Match-WrongCall', ans_pairs, waafle_orgs_single]
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

    print('\t'.join( ['contig', 'call', 'answer_status', 'waafle_status', 'org_status', 'recipient/donor_status', 'answer_recipient/donor_orgs', 'top_waafle_orgs'] ) )

    for astrline in open( args.answerkey ):
        aastrline = astrline.strip().split('\t')
        taxalevel, contig, recip, donor, = aastrline[0], aastrline[1], aastrline[2], aastrline[3]
        dict_answerkey[contig] = [taxalevel, recip, donor]

    set_waafle_contigs = set( [] )
    for bstrline in open( args.scoredcontigs ):
        bbstrline = bstrline.strip().split('\t')
        
        #skip the first header
        if bbstrline[1] == 'contiglen':
            continue
        contig, waafle_status = bbstrline[0], bbstrline[2]
        set_waafle_contigs.add( contig )
        LGTcall = ''
    
        #call status  
        ans_taxalevel = phylevels.index( dict_answerkey[contig][0] )
        waafle_taxalevel = phylevels.index( bbstrline[4].split('_')[0] )
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
        RD_results = []
        RD_orgs = []
        RD_call = ''
        for orgpair in bbstrline[5].split(';'):
            if orgpair == 'NA':
                RD_call = 'No_RDCall'
            else:
                waafle_recip, waafle_donor = orgpair.split('-')[0], orgpair.split('-')[1]
                RD_call, RDorgs = call_dr( ans_reciplevel, ans_donorlevel, waafle_recip, waafle_donor )
                RD_results.append( RD_call )
                RD_orgs.append( RDorgs )
        if RD_call != 'No_RDCall':
            RD_final, waafleorgs  = ';'.join(RD_results), ';'.join(RD_orgs)
        else:
            RD_final = 'No_RDCall'
        print('\t'.join( [contig, LGTcall, ans_status, waafle_status, orgcall, ansorgs, waafleorgs, RD_final] ) )
        
    for contig in dict_answerkey.keys():
        if contig not in set_waafle_contigs:
            ans_taxalevel = phylevels.index( dict_answerkey[contig][0] )
            waafle_taxalevel = phylevels.index( dict_answerkey[contig][1].split('_')[0] )
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

