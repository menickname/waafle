#!/usr/bin/env python

"""
This module is a part of:
WAAFLE, a [W]orkflow to [A]nnotate [A]ssemblies and [F]ind [L]GT [E]vents

Copyright (c) 2019 Harvard T.H. Chan School of Public Health

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
"""

from __future__ import print_function # Python 2.7+ required
import os
import sys
import csv
import argparse
import re
from collections import Counter

from waafle import dev_utils as wu

# ---------------------------------------------------------------
# description
# ---------------------------------------------------------------

description = wu.describe( """
{SCRIPT}: Applies junction results to QC WAAFLE calls

To be completed.
""" )

# ---------------------------------------------------------------
# output formats
# ---------------------------------------------------------------

c_formats = {}
for name, items in c_formats.items( ):
    c_formats[name] = [k for k in items.split( "\n" ) if k != ""]

# ---------------------------------------------------------------
# cli
# ---------------------------------------------------------------

def get_args( ):

    parser = argparse.ArgumentParser(
        description=description, 
        formatter_class=argparse.RawTextHelpFormatter,
        )

    g = parser.add_argument_group( "required inputs" )
    g.add_argument( 
        "contigs",
        help="contigs file (fasta format)",
        )
    g.add_argument( 
        "gff",
        help="GFF file for provided contigs",
        )

    g = parser.add_argument_group( "filtering parameters" )
    g.add_argument( 
        "--min-junction-hits",
        type=int,
        default=2,
        metavar="<int>",
        help="minimum read-hits to 'ok' a junction\n[default: 2]",
        )
    g.add_argument( 
        "--min-junction-ratio",
        type=int,
        default=0.5,
        metavar="<float>",
        help="minimum coverage (relative to flanking genes) to 'ok' a junction\n[default: 0.5]",
        )
    g.add_argument( 
        "--min-contig-genes",
        type=int,
        default=2,
        metavar="<int>",
        help="minimum gene count for the contig\n[default: 2]",
        )
    g.add_argument( 
        "--min-contig-length",
        type=int,
        default=500,
        metavar="<int>",
        help="minimum length for the contig\n[default: 500]",
        )
    g.add_argument( 
        "--max-contig-failures",
        type=float,
        default=0,
        metavar="<float>",
        help="allowed fraction of failing junctions\n[default: 0]",
        )

    args = parser.parse_args( )
    return args

# ---------------------------------------------------------------
# main
# ---------------------------------------------------------------

def main( ):

    args = get_args( )

    # write contig report
    wu.say( "Writing contig report." )
    with wu.try_open( p_contig_qc, "w" ) as fh:
        wu.write_rowdict( 
            format=c_formats["contig_qc"], 
            file=fh, )
        for c in sorted( contig_lengths ):
            ok = True
            rowdict = {}
            rowdict["contig"] = c
            rowdict["length"] = my_len = contig_lengths[c]
            if my_len < args.min_contig_length:
                ok = False
            rowdict["loci"] = my_loci = len( contig_loci.get( c, [] ) )
            if my_loci < args.min_contig_genes:
                ok = False
            rowdict["failing_junctions"] = my_fails = failures.get( c, 0 )
            # note pseudocount to avoid /0 on 1-gene contigs
            my_rate = my_fails / (my_loci - 1 + 1e-6)
            # define rate to be 0 on contigs without a junction
            rowdict["failure_rate"] = myrate = 0.0 if my_loci < 2 else my_rate
            if my_rate > args.max_contig_failures:
                ok = False
            rowdict["summary"] = "OK" if ok else "FAILED"
            # write
            wu.write_rowdict( 
                rowdict=rowdict,
                format=c_formats["contig_qc"],
                file=fh,
                precision=2,
                )

    # wrap-up
    wu.say( "Finished successfully." )           
                
if __name__ == "__main__":
    main( )
