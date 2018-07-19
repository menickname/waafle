#!/usr/bin/env python

"""
This module is a part of:
WAAFLE, the [W]orkflow to [A]nnotate [A]ssemblies and [F]ind [L]GT [E]vents

Copyright (c) 2018 Harvard T.H. Chan School of Public Health

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

import numpy as np

from waafle import utils2 as wu

# ---------------------------------------------------------------
# description
# ---------------------------------------------------------------

description = wu.describe( """
{SCRIPT}: Filter out low-quality contigs

This script maps reads to contigs (or analyzes an existing mapping)
to identify contigs whose junctions aren't supported by reads. A short
junction is supported if N+ read-pairs engulf the junction. A long junction
(i.e. too long to be engulfed by a read-pair) is supported if its coverage
is "reasonably similar" to the flanking genes' coverages (i.e. at least half
their mean). This script will also filter out contigs that are too short
or which contain too few genes.
""" )

# ---------------------------------------------------------------
# output formats
# ---------------------------------------------------------------

c_formats = {}

c_formats["site_hits"] = """
contig
mean
stdev
depths
"""

c_formats["gene_hits"] = """
contig
gene1
gene2
hits
"""

c_formats["junction_report"] = """
contig
gene1
gene2
hits_junction
hits_ok
coverage_gene1
coverage_gene2
coverage_junction
coverage_ratio
coverage_ok
acceptable
"""

c_formats["contig_report"] = """
contig
summary
length
loci
failing_junctions
failure_rate
"""

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
    parser.add_argument( 
        "contigs",
        help="contigs file (fasta format)",
        )
    parser.add_argument( 
        "gff",
        help="GFF file for provided contigs",
        )
    parser.add_argument( 
        "--reads1",
        help="1st-end sequencing reads",
        )
    parser.add_argument( 
        "--reads2",
        help="2nd-end sequencing reads",
        )
    parser.add_argument( 
        "--sam",
        help="sam file (if alignment already exists)",
        )
    parser.add_argument( 
        "--bowtie2-build",
        default="bowtie2-build",
        metavar="<path>",
        help="path to bowtie2-build\n[default: $PATH]",
        )
    parser.add_argument( 
        "--bowtie2",
        default="bowtie2",
        metavar="<path>",
        help="path to bowtie2\n[default: $PATH]",
        )
    parser.add_argument( 
        "--threads",
        default="1",
        metavar="<int>",
        help="number of cpu cores to use\n[default: 1]",
        )
    parser.add_argument( 
        "--tmpdir",
        default=".",
        metavar="<path>",
        help="where to place temp outputs\n[default: <.>]",
        )
    parser.add_argument( 
        "--outdir",
        default=".",
        metavar="<path>",
        help="where to place main outputs\n[default: <.>]",
        )
    parser.add_argument( 
        "--basename",
        metavar="<str>",
        help="basename for output files\n[default: <derived from input>]",
        )
    parser.add_argument( 
        "--read-gene-overlap",
        type=int,
        default=40,
        metavar="<1-N>",
        help="critical nucleotides for counting a read-gene overlap\n[default: <40>]",
        )
    parser.add_argument( 
        "--junction-fragments",
        type=int,
        default=2,
        metavar="<1-N>",
        help="critical fragments to 'ok' a junction\n[default: <2>]",
        )
    parser.add_argument( 
        "--junction-coverage-ratio",
        type=int,
        default=0.5,
        metavar="<float>",
        help="critical relative coverage to 'ok' a junction\n[default: <0.5>]",
        )
    parser.add_argument( 
        "--min-genes",
        type=int,
        default=2,
        metavar="<int>",
        help="minimum gene count for the contig\n[default: <2>]",
        )
    parser.add_argument( 
        "--min-length",
        type=int,
        default=500,
        metavar="<int>",
        help="minimum length for the contig\n[default: <500>]",
        )
    parser.add_argument( 
        "--allowed-failure-rate",
        type=float,
        default=0,
        metavar="<float>",
        help="allowed fraction of failing junctions\n[default: <0>]",
        )
    parser.add_argument( 
        "--write-detailed-output",
        action="store_true",
        help="write out the numbers that go into junction evaluation\n[default: off]",
        )
    args = parser.parse_args( )
    return args

# ---------------------------------------------------------------
# utils for running bowtie2
# ---------------------------------------------------------------

def bowtie2_build( p_bowtie2_build=None, p_contigs=None, p_index=None ):
    fields = {
        "PROG":    p_bowtie2_build,
        "CONTIGS": p_contigs, 
        "INDEX":   p_index, 
        }
    if os.path.exists( p_index + ".1.bt2" ):
        wu.say( "The index <{INDEX}> already exists.".format( **fields ) )
        wu.say( "(Delete to force rebuild.)" )
    else:
        wu.say( "Indexing <{CONTIGS}> to <{INDEX}>.".format( **fields ) )
        command = [
            "{PROG}",
            "{CONTIGS}",
            "{INDEX}",
            ]
        command = " ".join( command )
        command = command.format( **fields )
        os.system( "{PROG} {CONTIGS} {INDEX}".format( **fields ) )
        wu.say( "Build complete." )
    return None

def bowtie2_align( p_bowtie2=None, p_reads1=None, p_reads2=None, 
                   p_index=None, p_sam=None, threads=None ):
    fields = {
        "PROG":    p_bowtie2,
        "READS1":  p_reads1,
        "READS2":  p_reads2,
        "INDEX":   p_index,
        "SAM":     p_sam,
        "THREADS": threads,
        }
    wu.say( "Performing bowtie2 alignment." )
    command = [
        "{PROG}",
        "-x {INDEX}",
        "-1 {READS1}",
        "-2 {READS2}",
        "-S {SAM}",
        "--threads {THREADS}",
        "--no-mixed",
        "--no-discordant",
        ]
    command = " ".join( command )
    command = command.format( **fields )
    os.system( command )
    wu.say( "Alignment complete." )
    return None

# ---------------------------------------------------------------
# utils for parsing sam output
# ---------------------------------------------------------------

"""
@HD  VN:1.0                SO:unsorted
@SQ  SN:SRS011061_k119_3   LN:961
@SQ  SN:SRS011061_k119_5   LN:837
@SQ  SN:SRS011061_k119_11  LN:502
...
61NLYAAXX100508:5:100:10002:9010   83   SRS011061_k119_37319  25069  42  101M      =  25052  -118
61NLYAAXX100508:5:100:10002:9010   163  SRS011061_k119_37319  25052  42  46M1I43M  =  25069  118
61NLYAAXX100508:5:100:10002:9075   83   SRS011061_k119_14610  17113  42  101M      =  16942  -272
61NLYAAXX100508:5:100:10002:9075   163  SRS011061_k119_14610  16942  42  46M1I43M  =  17113  272
61NLYAAXX100508:5:100:10003:17250  77   *                     0      0   *         *  0      0
61NLYAAXX100508:5:100:10003:17250  141  *                     0      0   *         *  0      0
61NLYAAXX100508:5:100:10003:3146   99   SRS011061_k119_83764  304    42  101M      =  366    152
61NLYAAXX100508:5:100:10003:3146   147  SRS011061_k119_83764  366    42  90M       =  304    -152
"""

def cigar_length( cigar ):
    counts = [int( c ) for c in re.split( "[A-Z]+", cigar ) if c != ""]
    sigils = [s        for s in re.split( "[0-9]+", cigar ) if s != ""]
    # ignore read-only bands
    return sum( [c for c, s in zip( counts, sigils ) if s in "DHMNSX="] )

def good_pair( read1, read2 ):
    ret = True
    # mate pair
    if read1[0] != read2[0]:
        ret = False
    # aligned
    if read1[2] == "*":
        ret = False
    # concordantly
    if read1[2] != read2[2]:
        ret = False
    return ret

def process_pair( read1, read2 ):
    ret = None
    if good_pair( read1, read2 ):
        # read1 span
        start1 = int( read1[3] )
        end1   = start1 + cigar_length( read1[5] ) - 1
        span1  = [start1, end1]
        # read2 span
        start2 = int( read2[3] )
        end2   = start2 + cigar_length( read2[5] ) - 1
        span2  = [start2, end2]
        # format: read name, target, left span, right span
        span1, span2 = sorted( [span1, span2] )
        ret = [read1[0], read1[2], span1, span2]
    return ret

def iter_sam( p_sam ):
    counter = 0
    with wu.try_open( p_sam ) as fh:
        reader = csv.reader( fh, dialect="excel-tab" )
        for row in reader:
            if row[0][0] == "@":
                counter += 1
            else:
                counter += 2
                read1 = row
                read2 = reader.next( )
                hit = process_pair( read1, read2 )
                if hit is not None:
                    yield hit
            if counter % 10000 == 0:
                wu.say( "  {:,} lines processed".format( counter ) )

# ---------------------------------------------------------------
# utils for comparing SAM and GFF
# ---------------------------------------------------------------

def find_overlaps( span1=None, span2=None, loci=None, crit=1 ):
    hits = set( )
    for L in loci:
        a1, a2 = L.start, L.end
        for span in [span1, span2]:
            b1, b2 = span
            overlap = wu.calc_overlap( a1, a2, b1, b2, normalize=False )
            if overlap >= crit:
                hits.add( L.code )
    return hits

# ---------------------------------------------------------------
# utils for evaluating a contig
# ---------------------------------------------------------------

def code_start_stop( code ):
    items = code.split( ":" )
    start = int( items[0] )
    stop = int( items[1] )
    return [start, stop]

def evaluate_contig( loci=None, coverage=None, gene_hits=None, args=None ):
    rowdicts = []
    # determine the loci codes to consider
    codes = [L.code for L in sorted( loci, key=lambda x: x.start )]
    for i in range( len( codes ) - 1 ):
        rowdict = {}
        rowdict["gene1"]  = codes[i]
        rowdict["gene2"]  = codes[i+1]
        # check hits
        rowdict["hits_junction"]     = gh = gene_hits.get( (codes[i], codes[i+1]), 0 )
        rowdict["hits_ok"]           = hits_ok = gh >= args.junction_fragments
        # check coverage
        start1, stop1 = code_start_stop( codes[i] )
        start2, stop2 = code_start_stop( codes[i+1] )
        rowdict["coverage_gene1"]    = c1 = np.mean( coverage[start1-1:stop1] )
        rowdict["coverage_gene2"]    = c2 = np.mean( coverage[start2-1:stop2] )
        # define junction coverage as 0 if there's no gap
        cj = 0.0 if start2 <= stop1 else np.mean( coverage[stop1-1:start2] )
        rowdict["coverage_junction"] = cj
        # note: pseudocount to avoid /0
        rowdict["coverage_ratio"]    = cr = cj / ( np.mean( [c1, c2] ) + 1e-6 )
        rowdict["coverage_ok"]       = coverage_ok = cr >= args.junction_coverage_ratio
        # finalize
        rowdict["acceptable"]        = rowdict["hits_ok"] or rowdict["coverage_ok"]
        rowdicts.append( rowdict )
    return rowdicts

# ---------------------------------------------------------------
# utils for output formatting
# ---------------------------------------------------------------

def write_detailed_output( basename=None, tmpdir=None,
                           contig_coverage=None, contig_hits=None, ):

    p_site_hits = wu.name2path( basename, p_tmpdir, ".site_hits.tsv.gz" )
    p_gene_hits = wu.name2path( basename, p_tmpdir, ".gene_hits.tsv.gz" )

    # write: site_hits
    wu.say( "Writing site hits." )
    with wu.try_open( p_site_hits, "w" ) as fh:
        wu.write_rowdict( 
            format=c_formats["site_hits"], 
            file=fh, )
        for c in sorted( contig_coverage ):
            depths = contig_coverage[c]
            rowdict = {
                "contig" : c,
                "mean"   : np.mean( depths ),
                "stdev"  : np.std( depths ),
                "depths" : " ".join( ["{:.0f}".format( k ) for k in depths] ),
                }
            wu.write_rowdict( 
                rowdict=rowdict, 
                format=c_formats["site_hits"], 
                file=fh, )

    # write: gene_hits
    wu.say( "Writing gene-pair hits." )
    with wu.try_open( p_gene_hits, "w" ) as fh:
        wu.write_rowdict( 
            format=c_formats["gene_hits"], 
            file=fh, )
        for c in sorted( contig_hits ):
            for code1, code2 in sorted( contig_hits[c] ):
                value = contig_hits[c][(code1, code2)]
                rowdict = {
                    "contig" : c,
                    "gene1"  : code1,
                    "gene2"  : code2,
                    "hits"   : value,
                    }
                wu.write_rowdict( 
                    rowdict=rowdict, 
                    format=c_formats["gene_hits"], 
                    file=fh, )

# ---------------------------------------------------------------
# main
# ---------------------------------------------------------------

def main( ):

    args = get_args( )
    p_contigs = args.contigs
    p_gff = args.gff

    # define files
    p_outdir = args.outdir
    p_tmpdir = args.tmpdir
    basename = args.basename
    if basename is None:
        basename = wu.path2name( p_contigs )
    p_index           = wu.name2path( basename, p_tmpdir, ".index" )
    p_sam             = wu.name2path( basename, p_tmpdir, ".sam" )
    p_junction_report = wu.name2path( basename, p_outdir, ".junction_report.tsv" )
    p_contig_report   = wu.name2path( basename, p_outdir, ".contig_report.tsv" )

    # alignment workflow
    if args.sam is not None:
        p_sam = args.sam
        wu.say( "Using specified SAM file:", p_sam )
    elif args.reads1 is not None and args.reads2 is not None:
        # build process
        bowtie2_build( 
            p_bowtie2_build = args.bowtie2_build,
            p_contigs = args.contigs,
            p_index = p_index,
            )
        # alignment process
        bowtie2_align(
            p_bowtie2=args.bowtie2,
            p_reads1=args.reads1, 
            p_reads2=args.reads2,
            p_index=p_index,
            p_sam=p_sam,
            threads=args.threads,
            )
    else:
        wu.die( "Must provide reads for alignment or SAM file." )

    # load contig data
    contig_lengths = wu.read_contig_lengths( p_contigs )
    contig_coverage = {}
    contig_hits = {}
    wu.say( "Loading contig lengths." )
    for name, length in contig_lengths.items( ):
        contig_coverage[name] = np.zeros( length )
    contig_loci = {}
    wu.say( "Loading contig gene coordinates." )
    for name, loci in wu.iter_contig_loci( p_gff ):
        contig_loci[name] = loci

    # post-processing workflow
    wu.say( "Processing SAM file." )
    for read, contig, span1, span2 in iter_sam( p_sam ):
        inner = contig_hits.setdefault( contig, Counter( ) )
        # update pers-site coverage
        L = min( span1 + span2 ) - 1
        R = max( span1 + span2 ) - 1
        contig_coverage[contig][L:R+1] += 1
        # find hit loci
        hits = find_overlaps( 
            span1=span1, 
            span2=span2, 
            loci=contig_loci.get( contig, [] ),
            crit=args.read_gene_overlap,
            )
        # update self counts
        for code in hits:
            inner[(code, code)] += 1
        # update pair counts
        for code1 in hits:
            for code2 in hits:
                if code1 < code2:
                    inner[(code1, code2)] += 1

    # detailed output?
    if args.write_detailed_output:
        write_detailed_output(
            basename=basename,
            tmpdir=p_tmpdir,
            contig_coverage=contig_coverage,
            contig_hits=contig_hits,
            )

    # write junction report
    wu.say( "Writing junction report." )
    failures = {}
    with wu.try_open( p_junction_report, "w" ) as fh:
        wu.write_rowdict( 
            format=c_formats["junction_report"], 
            file=fh, )
        for c in sorted( contig_lengths ):
            rowdicts = evaluate_contig( 
                loci=contig_loci.get( c, [] ),
                coverage=contig_coverage[c], 
                gene_hits=contig_hits.get( c, {} ),
                args=args,
                )
            failures[c] = sum( [1 for r in rowdicts if not r["acceptable"]] )
            for rowdict in rowdicts:
                rowdict["contig"] = c
                wu.write_rowdict( 
                    rowdict=rowdict, 
                    format=c_formats["junction_report"], 
                    file=fh, )

    # write contig report
    wu.say( "Writing contig report." )
    with wu.try_open( p_contig_report, "w" ) as fh:
        wu.write_rowdict( 
            format=c_formats["contig_report"], 
            file=fh, )
        for c in sorted( contig_lengths ):
            ok = True
            rowdict = {}
            rowdict["contig"] = c
            rowdict["length"] = my_len = contig_lengths[c]
            if my_len < args.min_length:
                ok = False
            rowdict["loci"] = my_loci = len( contig_loci.get( c, [] ) )
            if my_loci < args.min_genes:
                ok = False
            rowdict["failing_junctions"] = my_fails = failures.get( c, 0 )
            # note pseudocount to avoid /0 on 1-gene contigs
            rowdict["failure_rate"] = my_rate = my_fails / (my_loci - 1 + 1e-6)
            # define rate to be 0 on contigs without a junction
            my_rate = 0 if my_loci < 2 else my_rate
            if my_rate > args.allowed_failure_rate:
                ok = False
            rowdict["summary"] = "OK" if ok else "FAILED"
            # write
            wu.write_rowdict( 
                rowdict=rowdict,
                format=c_formats["contig_report"],
                file=fh,
                )

    # wrap-up
    wu.say( "Finished successfully." )           
                
if __name__ == "__main__":
    main( )
