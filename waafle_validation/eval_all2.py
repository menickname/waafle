#!/usr/bin/env/python

from __future__ import print_function # python 2.7+ required
import os, sys, csv, argparse, subprocess, shlex

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
        "-w_genes", "--waafleeval",
        required=True,
        help="evaluated waafle genes"
        )
    parser.add_argument(
        "-ref_genes", "--refeval",
        required=True,
        help="evaluated ref genes"
        )
    parser.add_argument(
        "-contig_eval", "--contigseval",
        required=True,
        help="evaluated contig genes"
        )
    parser.add_argument(
        "-condition", "--currwd",
        required=True,
        help="get conditions under which this was generated"
        )
    args = parser.parse_args()
    return args

def count_grep_lines( myfile, grep_term ):
    p1 = subprocess.Popen( ["grep", grep_term, myfile], stdout=subprocess.PIPE )
    p2 = subprocess.Popen( ["wc", "-l"], stdin=p1.stdout, stdout=subprocess.PIPE )
    p1.stdout.close()
    output = p2.communicate()[0].strip()
    return output

# ---------------------------------------------------------------
# main
# ---------------------------------------------------------------

def main():
    args = get_args()
    
    called_wgenes = count_grep_lines( args.waafleeval, "called" )
    partial_wgenes = count_grep_lines( args.waafleeval, "partial" )
    misseds_wgenes = count_grep_lines( args.waafleeval, "missed_strand" )
    missedo_wgenes = count_grep_lines( args.waafleeval, "missed_overlap" )
    missedc_wgenes = count_grep_lines( args.waafleeval, "missed_nogene" )
    called_rgenes = count_grep_lines( args.refeval, "called" )
    partial_rgenes = count_grep_lines( args.refeval, "partial" )
    misseds_rgenes = count_grep_lines( args.refeval, "missed_strand" )
    missedo_rgenes = count_grep_lines( args.refeval, "missed_overlap" )
    missedc_rgenes = count_grep_lines( args.refeval, "missed_nogene" )

    total_wgenes = int( called_wgenes ) + int( partial_wgenes ) + int( misseds_wgenes ) + int( missedo_wgenes ) + int( missedc_wgenes)
    total_rgenes = int( called_rgenes ) + int( partial_rgenes ) + int( misseds_rgenes ) + int( missedo_rgenes ) + int( missedc_rgenes )
  
    answerheader, answerline = [ args.currwd ], [ args.currwd ] 
    TP_genes = called_rgenes
    FN_genes = int( misseds_rgenes ) + int( missedo_rgenes ) + int( missedc_rgenes )
    FP_genes = int( misseds_wgenes ) + int( missedo_wgenes )
    
    answerheader += [ 'called_wgenes', 'partial_wgenes', 'misseds_wgenes', 'missedo_wgenes', 'missedc_wgenes', 'total_wgenes', 'called_rgenes', 'partial_rgenes', 'misseds_rgenes', 'missedo_rgenes', 'missedc_wgenes', 'total_rgenes', 'TP_genes', 'FN_genes', 'FP_genes' ]
    answerline += [ called_wgenes, partial_wgenes, misseds_wgenes, missedo_wgenes, missedc_wgenes, str( total_wgenes ), called_rgenes, partial_rgenes, misseds_rgenes, missedo_rgenes, missedc_rgenes, str( total_rgenes ), str( TP_genes ), str( FN_genes ), str( FP_genes ) ]

    taxalist = ['k', 'p', 'c', 'o', 'f', 'g', 's']
    for taxa in taxalist:
        filename = args.contigseval + taxa + '.tsv'
        TP = count_grep_lines( filename, "TP" )
        TN = count_grep_lines( filename, "TN" )
        FP = count_grep_lines( filename, "FP" )
        FN = count_grep_lines( filename, "FN" )
        answerheader += ['TP_' + taxa, 'TN_' + taxa, 'FP_' + taxa, 'FN_' + taxa]
        answerline += [TP, TN, FP, FN]
    
    print( '\t'.join( answerheader ) )
    print( '\t'.join( answerline ) )
    
if __name__ == "__main__":
    main()
