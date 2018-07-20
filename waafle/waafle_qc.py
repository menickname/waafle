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

from waafle import dev_utils as wu

# ---------------------------------------------------------------
# description
# ---------------------------------------------------------------

description = wu.describe( """
{SCRIPT}: Performs coverage-level quality control on contigs

This script maps reads to contigs (or analyzes an existing mapping)
to identify contigs whose junctions aren't supported by reads. A short
junction is supported if mate-pairs engulf the junction. A long junction
(i.e. too long to be engulfed by a mate-pair) is supported if its coverage
is "reasonably similar" to the flanking genes' coverages (e.g. at least half
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

c_formats["junctions"] = """
contig
gene1
gene2
gap
hits_junction
hits_ok
coverage_gene1
coverage_gene2
coverage_junction
coverage_ratio
coverage_ok
acceptable
"""

c_formats["contig_qc"] = """
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

    g = parser.add_argument_group( "required inputs" )
    g.add_argument( 
        "contigs",
        help="contigs file (fasta format)",
        )
    g.add_argument( 
        "gff",
        help="GFF file for provided contigs",
        )

    g = parser.add_argument_group( "provide reads or a sam file" )
    g.add_argument( 
        "--reads1",
        metavar="<path>",
        help="1st-end sequencing reads",
        )
    g.add_argument( 
        "--reads2",
        metavar="<path>",
        help="2nd-end sequencing reads",
        )
    g.add_argument( 
        "--sam",
        metavar="<path>",
        help="sam file (if alignment already exists)",
        )

    g = parser.add_argument_group( "output options" )
    g.add_argument( 
        "--tmpdir",
        default=".",
        metavar="<path>",
        help="where to place temp outputs\n[default: ./]",
        )
    g.add_argument( 
        "--outdir",
        default=".",
        metavar="<path>",
        help="where to place main outputs\n[default: ./]",
        )
    g.add_argument( 
        "--basename",
        metavar="<str>",
        help="basename for output files\n[default: <derived from input>]",
        )
    g.add_argument( 
        "--write-detailed-output",
        action="store_true",
        help="write out the numbers that go into junction evaluation\n[default: off]",
        )

    g = parser.add_argument_group( "filtering parameters" )
    g.add_argument( 
        "--min-overlap-sites",
        type=int,
        default=25,
        metavar="<int>",
        help="minimum nucleotide overlap for counting a read-gene hit\n[default: 25]",
        )
    g.add_argument( 
        "--min-junction-hits",
        type=int,
        default=2,
        metavar="<int>",
        help="minimum hits to 'ok' a junction\n[default: 2]",
        )
    g.add_argument( 
        "--min-junction-ratio",
        type=int,
        default=0.5,
        metavar="<float>",
        help="minimum relative coverage to 'ok' a junction\n[default: 0.5]",
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

    g = parser.add_argument_group( "bowtie2 options" )
    g.add_argument( 
        "--bowtie2-build",
        default="bowtie2-build",
        metavar="<path>",
        help="path to bowtie2-build\n[default: $PATH]",
        )
    g.add_argument( 
        "--bowtie2",
        default="bowtie2",
        metavar="<path>",
        help="path to bowtie2\n[default: $PATH]",
        )
    g.add_argument( 
        "--threads",
        default="1",
        metavar="<int>",
        help="number of cpu cores to use in bowtie2 alignment\n[default: 1]",
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
        wu.say( "(Move/delete to force rebuild.)" )
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
# utils for parsing SAM/GFF comparison
# ---------------------------------------------------------------

def concordant_hits( p_sam=None ):
    counter = 0
    mate1 = None
    mate2 = None
    for hit in wu.iter_sam_hits( p_sam ):
        # progress
        counter += 1
        if counter % int( 1e5 ) == 0:
            wu.say( "  SAM alignments processed: {:.1f}M".format( counter / 1e6 ) )
        # weave
        mate1 = mate2
        mate2 = hit
        # edge case
        if mate1 is None:
            continue
        # not a mate pair
        elif mate1.qseqid != mate2.qseqid:
            continue
        # not concordant
        elif mate1.sseqid != mate2.sseqid:
            continue
        # good pair
        else:
            yield [mate1, mate2]

def find_hit_loci( mate1=None, mate2=None, loci=None, args=None ):
    hits = set( )
    for L in loci:
        a1, a2 = L.start, L.end
        for read in [mate1, mate2]:
            b1, b2 = read.sstart, read.send
            overlap = wu.calc_overlap( a1, a2, b1, b2, normalize=False )
            #print( L.code, read.qseqid, a1, a2, b1, b2, overlap )
            if overlap >= args.min_overlap_sites:
                hits.add( L.code )
    return hits

# ---------------------------------------------------------------
# utils for evaluating a contig
# ---------------------------------------------------------------

def evaluate_contig( loci=None, coverage=None, gene_hits=None, args=None ):
    rowdicts = []
    loci = sorted( loci, key=lambda x: x.start )
    for i in range( len( loci ) - 1 ):
        rowdict = {}
        L1 = loci[i]
        L2 = loci[i+1]
        rowdict["gene1"] = code1 = L1.code
        rowdict["gene2"] = code2 = L2.code
        rowdict["gap"]   = L2.start - L1.end - 1
        # this ensures that lookup matches non-redundant storage
        pair_key = tuple( sorted( [code1, code2] ) )
        # check hits
        rowdict["hits_junction"] = my_hits = gene_hits.get( pair_key, 0 )
        rowdict["hits_ok"] = hits_ok = my_hits >= args.min_junction_hits
        # check coverage (note: base-0 start and pythonic end)
        rowdict["coverage_gene1"] = my_cov1 = np.mean( coverage[L1.start-1:L1.end] )
        rowdict["coverage_gene2"] = my_cov2 = np.mean( coverage[L2.start-1:L2.end] )
        # define junction coverage as 0 if there's no gap
        my_covj = 0.0 if L2.start <= L1.end else np.mean( coverage[L1.end-1:L2.start] )
        rowdict["coverage_junction"] = my_covj
        # note: pseudocount to avoid /0
        rowdict["coverage_ratio"] = my_ratio = my_covj / ( np.mean( [my_cov1, my_cov2] ) + 1e-6 )
        rowdict["coverage_ok"] = coverage_ok = my_ratio >= args.min_junction_ratio
        # finalize
        rowdict["acceptable"] = rowdict["hits_ok"] or rowdict["coverage_ok"]
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
    p_index     = wu.name2path( basename, p_tmpdir, ".index" )
    p_sam       = wu.name2path( basename, p_tmpdir, ".sam" )
    p_junctions = wu.name2path( basename, p_outdir, ".junctions.tsv" )
    p_contig_qc = wu.name2path( basename, p_outdir, ".contig_qc.tsv" )

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
        wu.die( "Must provide READS1/2 or SAM file." )

    # load contig data
    wu.say( "Loading contig lengths." )
    contig_lengths = wu.read_contig_lengths( p_contigs )
    contig_coverage = {}
    for name, length in contig_lengths.items( ):
        contig_coverage[name] = np.zeros( length )
    wu.say( "Loading contig gene coordinates." )
    contig_loci = {}
    for name, loci in wu.iter_contig_loci( p_gff ):
        contig_loci[name] = loci
    contig_hits = {}

    # post-processing workflow
    wu.say( "Processing SAM file." )
    for mate1, mate2 in concordant_hits( p_sam ):
        contig = mate1.sseqid
        inner = contig_hits.setdefault( contig, Counter( ) )
        # update pers-site coverage (note: base-0 start and pythonic end)
        coords = [mate1.sstart, mate1.send, mate2.sstart, mate2.send]
        L = min( coords ) - 1
        R = max( coords ) - 1
        contig_coverage[contig][L:R+1] += 1
        # find hit loci
        hits = find_hit_loci( 
            mate1=mate1, 
            mate2=mate2, 
            loci=contig_loci.get( contig, [] ),
            args=args,
            )
        # update self counts
        for code in hits:
            inner[(code, code)] += 1
        # update pair counts
        for code1 in hits:
            for code2 in hits:
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
    with wu.try_open( p_junctions, "w" ) as fh:
        wu.write_rowdict( 
            format=c_formats["junctions"], 
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
                    format=c_formats["junctions"], 
                    file=fh, )

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
