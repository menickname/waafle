#!/usr/bin/env python

"""
WAAFLE MODULE: waafle_genecaller.py

Authors:
Tiffany Hsu
Eric Franzosa

This script interprets waafle blast output to make gene calls.
Run with the "-h" flag for usage help.

Details of the GTF/GFF file format:

Fields must be tab-separated. 
Empty columns should be denoted with a '.'.

0: seqname - name of the chromosome or scaffold
1: source - name of the program that generated this feature
2: feature - feature type name, e.g. Gene, Variation, Similarity
3: start - Start position of the feature, with sequence numbering starting at 1.
4: end - End position of the feature, with sequence numbering starting at 1.
5: score - A floating point value.
6: strand - defined as + (forward) or - (reverse).
7: frame - One of '0', '1' or '2'.
8: attribute - A semicolon-separated list of tag-value pairs.
"""

from __future__ import print_function # python 2.7+ required
import os, sys, csv, argparse
import waafle_utils as wu

# ---------------------------------------------------------------
# constants
# ---------------------------------------------------------------

c_dummy_gff_row = ["." for i in range(9)]
c_dummy_gff_row[1] = "WAAFLE"
c_dummy_gff_row[2] = "Gene"

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
        help="output from waafle_search"
        )
    parser.add_argument( 
        "-o", "--out",
        default="waafle_genes.gff",
        help="waafle gene calls",
        )
    args = parser.parse_args()
    return args

def hits2genes( blastoutfile ):
    """
    Silly placeholder code to illustrate hit object
    Will be replaced by Tiffany's ACTUAL code for merging blast
    hits into genes.
    """
    genes = []
    for contig, sorted_hits in wu.iter_contig_hits( blastoutfile ):
        """ file me in """
    return genes

def write_gff( genes, outfile ):
    """
    Write genes to a file
    """
    with wu.try_open( outfile, "w" ) as fh:
        writer = csv.writer( fh, dialect="excel-tab" )
        for contig, start, end, strand in genes:
            # important to copy [:]
            gene = c_dummy_gff_row[:] 
            gene[0] = contig
            gene[3] = start
            gene[4] = end
            gene[6] = strand
            writer.writerow( gene )

# ---------------------------------------------------------------
# main
# ---------------------------------------------------------------

def main():
    args = get_args()
    write_gff( hits2genes( args.input ), args.out )

if __name__ == "__main__":
    main()
