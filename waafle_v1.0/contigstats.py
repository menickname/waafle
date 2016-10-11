#!/usr/bin/python


from __future__ import print_function # python 2.7+ required
import os, sys, csv, argparse, re
import waafle_utils as wu
from operator import itemgetter, attrgetter, methodcaller
from collections import Counter

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
        "-c", "--contiglist",
        help="list of contigs to get stats from"
        )
    parser.add_argument(
        "-i", "--input",
        help="Results from WAAFLE",
        )
    parser.add_argument(
        "-gff", "--gff",
        help="GFF file for genes",
        )
    args = parser.parse_args()
    return args

def cov_by_bug( status, genelen, synteny, index ):
    onebug, twobug = 0, 0
    if re.search( 'LGT:', status ):
        syntenyletter = synteny[index]
        if syntenyletter == 'A':
            onebug = genelen
        else:
            twobug = genelen
    return onebug, twobug
    

# ---------------------------------------------------------------
# main
# ---------------------------------------------------------------
def main():
    args = get_args()

    # store results
    dict_presults = {}
    dict_uresults = {}
    for astrline in open( args.input ):
        aastrline = astrline.strip().split('\t')
        contig, length = aastrline[0], aastrline[1]
        pref_call, unambig_call = aastrline[4], aastrline[5]

        if contig != 'contig':
            p_status, p_onebug, p_twobug, p_taxa, p_synteny = pref_call.split(';')[0], pref_call.split(';')[1], pref_call.split(';')[2], pref_call.split(';')[3], pref_call.split(';')[4]
            u_status, u_onebug, u_twobug, u_taxa, u_synteny = unambig_call.split(';')[0], unambig_call.split(';')[1], unambig_call.split(';')[2], unambig_call.split(';')[3], unambig_call.split(';')[4]
            p_scorediff = float(p_twobug) - float(p_onebug)
            u_scorediff = float(u_twobug) - float(u_onebug)
            dict_presults[contig] = [length, p_status, p_onebug, p_twobug, p_scorediff, p_taxa, p_synteny ]
            dict_uresults[contig] = [length, u_status, u_onebug, u_twobug, u_scorediff, u_taxa, u_synteny ]

    # calculate contig coverage
    for contig, genelist in wu.iter_contig_genes( args.gff ):
        totalgenelen = 0
        numgenes = 0
        pbugone, pbugtwo = 0, 0
        ubugone, ubugtwo = 0, 0
        psmallfrac, usmallfrac = 0, 0

        for i in range( len( genelist ) ):
            gene = genelist[i]
            numgenes += 1
            genelen = gene.end - gene.start + 1
            totalgenelen += genelen

            # calculate coverage per bug
            p_one, p_two = cov_by_bug( dict_presults[contig][1], genelen, dict_presults[contig][6], i)
            u_one, u_two = cov_by_bug( dict_uresults[contig][1], genelen, dict_uresults[contig][6], i)
            pbugone += p_one
            pbugtwo += p_two
            ubugone += u_one
            ubugtwo += u_two

        if re.search( 'LGT:', dict_presults[contig][1] ):
            psmallfrac = min( pbugone, pbugtwo )/float( totalgenelen )

        if re.search( 'LGT:', dict_uresults[contig][1] ):
            usmallfrac = min( ubugone, ubugtwo )/float( totalgenelen )
                    
        genecov = totalgenelen / float( dict_uresults[contig][0] )
        pnewinfo = [genecov, pbugone, pbugtwo, psmallfrac]
        unewinfo = [genecov, ubugone, ubugtwo, usmallfrac]
        dict_presults[contig] = dict_presults[contig] + pnewinfo
        dict_uresults[contig] = dict_uresults[contig] + unewinfo
        print( '\t'.join( str(x) for x in [contig] + dict_presults[contig] ) )

if __name__ == "__main__":
    main()
