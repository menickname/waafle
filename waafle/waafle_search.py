#!/usr/bin/env python

"""
WAAFLE MODULE: waafle_search.py

Authors:
Tiffany Hsu
Eric Franzosa

This script calls a custom blast search of a set of contigs against 
the waafle database (a modified version of the chocophlan database). 
Run with the "-h" flag for usage help.
"""

from __future__ import print_function # python 2.7+ required
import os, sys, csv, argparse
import waafle_utils as wu

def get_args():
    """
    Get arguments passed to script
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument( 
        "-q", "--query",
        required=True,
        help="contigs file (fasta format)",
        )
    parser.add_argument( 
        "-d", "--db",
        required=True,
        help="path to chocophlan database",
        )
    parser.add_argument( 
        "-b", "--blastpath",
        default="blastn",
        help="path to blast",
        )
    parser.add_argument( 
        "-n", "--num_threads",
        default="1",
        help="number of cpu cores to use in search",
        )
    parser.add_argument( 
        "-o", "--out",
        default="waafle-blastout.tsv",
        help="name of blast output file",
        )
    parser.add_argument( 
        "-e", "--execute",
        default="0",
        help="Run the search",
        )
    args = parser.parse_args()
    return args

# ---------------------------------------------------------------
# main
# ---------------------------------------------------------------

def main():
    args = get_args()
    command = "{0} -query {1} -db {2} -out {3} -max_target_seqs {4} -num_threads {5} -outfmt \'{6}\'".format( 
        args.blastpath,
        args.query,
        args.db,
        args.out,
        wu.c_max_target_seqs,
        args.num_threads,
        wu.c_blast_format_string,
        )
    print( "Executing BLASTN command:\n{}".format( command ), file=sys.stderr )
    if int( args.execute ):
        os.system( command )

if __name__ == "__main__":
    main()
