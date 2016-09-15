#!/usr/bin/env python

"""
WAAFLE MODULE: waafle_utils.py

Shared components of the waafle system.

Authors:
Tiffany Hsu
Eric Franzosa
"""

from __future__ import print_function # python 2.7+ required
import os, sys, csv, argparse, re
from operator import itemgetter, attrgetter, methodcaller
import numpy as np
from collections import Counter

"""
@codereview 9/2/2016

General comments: 

* Delete stuff here that we aren't using

* If a function or object is ONLY used in one script, it's better
to have it there

* The key reason for util is sharing object / functions across scripts
or forcing consistency as in the blast output

"""

# ---------------------------------------------------------------
# constants
# ---------------------------------------------------------------

c__list_taxa = ["k", "p", "c", "o", "f", "g", "s"]

c_blastfields = [
    ["qseqid", str],
    ["sseqid", str],
    ["qlen", int],
    ["slen", int],
    ["length", int],
    ["qstart", int],
    ["qend", int],
    ["sstart", int],
    ["send", int],
    ["pident", float],
    ["positive", int],
    ["gaps", int],
    ["evalue", float],
    ["bitscore", float],
    ["sstrand", str],
]

c_blast_format_string = " ".join( ["6"] + [fname for [fname, ftype] in c_blastfields] )
c_choco_header_delim = "|"
c_tax_delim = "."
# blast default is 500, which is sometimes too small for long contigs
c_max_target_seqs = 10000

c_gfffields = [
    ["seqname", str],
    ["source", str],
    ["feature", str],
    ["start", int],
    ["end", int],
    ["score", float],
    ["strand", str],
    ["frame", str],
    ["attribute", str],
]

c_taxafields = [
    ["contig", str],
    ["length", int],
    ["gene", int],
    ["strand", str],
    ["genestart", int],
    ["geneend", int],
    ["taxa", str],
    ["taxastart", str],
    ["taxaend", str],
    ["score", float],
    ["uniref50", str],
    ["uniref90", str],
    ["hits", int],
]

# ---------------------------------------------------------------
# classes for working with hits and gff (here to force equivalenence with output)
# ---------------------------------------------------------------

class Hit( ):
    """
    Processes the information from a single blast row;
    Row is provided already split by the csv reader.
    Hit MUST be compatible with the blast search defined above.

    Additional information is pulled into the hit object by parsing the
    hit chocophlan centroid's header. Chocophlan headers have the following
    fields delimited with a pipe ("|"):
    
    0: non-informative
    1: gi number
    2: non-informative
    3: refseq number
    4: coordinates
    5: ncbi tax id
    6: taxonomy ("."-delimited)
    7: UniRef90 hit
    8: UniRef50 hit
    """
    def __init__( self, blastrow ):
        assert len( blastrow ) == len( c_blastfields ), \
            "inconsistent blast row {}".format( str( blastrow ) )
        # pull values from blast line and coerce to appropriate types
        for [fname, ftype], value in zip( c_blastfields, blastrow ):
            setattr( self, fname, ftype( value ) )
        # derived coverage stats
        self.scov = ( abs( self.send - self.sstart ) + 1 ) / float( self.slen )
        self.qcov = ( abs( self.qend - self.qstart ) + 1 ) / float( self.qlen )
        # special scoverage that won't penalize hanging off contig end
        if self.sstrand == 'minus':
            sstart, send = self.slen - self.sstart + 1, self.slen - self.send + 1
        else:
            sstart, send = self.sstart, self.send
        self.ltrim = max( 0, sstart - self.qstart )
        self.rtrim = max( 0, self.slen - sstart - self.qlen + self.qstart )
        self.scov_modified = ( send - sstart + 1 ) / float( self.slen - self.ltrim - self.rtrim )
        # values extracted from chocophlan header
        chocoitems = self.sseqid.split( c_choco_header_delim )
        self.taxid = chocoitems[5]
        self.taxonomy = chocoitems[6].split( c_tax_delim )
        self.uniref90 = chocoitems[7].split('_')[1]
        self.uniref50 = chocoitems[8].split('_')[1]


class GFF( ):
    """
    Processes the information from a single gff line;
    Row is provided already split by the csv reader.

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
    def __init__( self, gffrow ):
        for [fname, ftype], value in zip( c_gfffields, gffrow ):
       	    setattr( self, fname, ftype( value ) )
        attritems = self.attribute.split(';')
        dict_attr = {}
        for items in attritems:
            if len( items.split('=') ) == 2:
                label, descriptor = items.split('=')
                dict_attr[label] = descriptor
        self.attributes = dict_attr

class Taxa( ):
    """
    Details of the scored orgs file:
    Fields should be tab-separated.
    
    0: contig - name of the contig
    1: length - length of the contig
    2: gene - gene number
    3: strand - defined as + (forward) or - (reverse).
    4: genestart - start position for the gene in the contig.
    5: geneend - end position for the end in the contig.
    6: taxa - taxa name at level chosen in waafle_orgscorer.py
    7: taxastart - start position for the taxa in the contig, may not equal gene start.
    8: taxaend - end position for the taxa in the contig, may not equal gene end.
    9: score -  ranges from 0-1. Average score across all positions for aggregated hits.
    10: uniref50 - uniref50 aggregated from all hits in gene.
    11: uniref90 - uniref90 aggregated from all hits in gene.
    12: hits - number of hits to taxa in gene.
    """
    def __init__( self, taxainfo ):
        for [fname, ftype], value in zip( c_taxafields, taxainfo ):
            setattr( self, fname, ftype( value ) )


class INode:
    """interval node: represents an interval + some network properties""" 
    def __init__( self, start, stop, strand="+" ):
        self.start, self.stop = sorted( [start, stop] )
        self.strand = strand
        self.neighbors = set()
        self.visited = False
    def __len__( self ):
        return self.stop - self.start + 1
    def attach( self, node ):
        self.neighbors.update( [node] )
    def get_connected_component( self ):
        # modified to use breadth-first search
        cc, front = {self}, {self}
        while any( [not inode.visited for inode in front] ):
            new_front = set()
            for inode in front:
                if not inode.visited:
                    inode.visited = True
                    new_front.update( inode.neighbors )
            cc.update( new_front )
            front = new_front
        return list( cc )
    def to_list( self ):
        return [self.start, self.stop, self.strand] 

# ---------------------------------------------------------------
# functions
# ---------------------------------------------------------------

def try_open( path, *args ):
    """
    Open a file; fail gracefully
    """
    try:
        fh = open( path, *args )
    except:
        sys.exit( "Can't open blast out file: {}".format( path ) )
    return fh

def iter_hits( blastoutfile ):
    """
    Iterate through the hits in a blast file
    """
    with try_open( blastoutfile ) as fh:                 
        for row in csv.reader( fh, dialect="excel-tab" ):
            yield Hit( row )

def iter_contig_hits( blastoutfile ):
    """
    Iterate through hits by contig (assumes file is sorted by query)
    """
    contig, hits = None, []
    with try_open( blastoutfile ) as fh:                 
        for row in csv.reader( fh, dialect="excel-tab" ):
            hit = Hit( row )
            if contig is not None and hit.qseqid != contig:
                yield contig, hits
                # reset
                hits = []
            contig = hit.qseqid
            hits.append( hit )
        # last case cleanup
        yield contig, hits

def iter_contig_genes( gfffile ):
    """
    Iterate through genes by contig (in gff)
    """
    contig, genes = None, []
    with try_open( gfffile ) as fh:
        for row in csv.reader( fh, dialect="excel-tab" ):
            if not re.search( r'^#', row[0] ):
                gene = GFF( row )
                if contig is not None and gene.seqname != contig:
                    yield contig, genes
                    # reset
                    genes = []
                contig = gene.seqname
                genes.append( gene )
	    # last case cleanup
        yield contig, genes

def iter_contig_taxa( orgfile ):
    """
    Iterate through taxa by contig (in org-scores file)
    """
    contig, taxa = None, []
    with try_open( orgfile ) as fh:
        for row in csv.reader( fh, dialect="excel-tab" ):
            if row[0] != "contig":
                taxon = Taxa( row )
                if contig is not None and taxon.contig != contig:
                    yield contig, taxa
                    # reset
                    taxa = []
                contig = taxon.contig
                taxa.append( taxon )
        #last case cleanup
        yield contig, taxa

def convert_strand( strandw ):
    strand = '-'
    if strandw == 'plus':
        strand = '+'
    return strand

def calc_overlap( a1, b1, a2, b2 ):
    """ compute overlap between two intervals """
    if b1 < a2 or b2 < a1:
        return 0
    else:
        outleft, inleft, inright, outright = sorted( [a1, b1, a2, b2] )
        denom = min( ( b1 - a1 + 1 ), ( b2 - a2 + 1 ) )
        return ( inright - inleft + 1 ) / float( denom )

def order_taxa( taxa ):
    """ Format the taxa class into an ordered list for printing. """
    orderedlist = [ str( taxa.contig ),
                    str( taxa.length ),
                    str( taxa.gene ),
                    str( taxa.strand ),
                    str( taxa.genestart ),
                    str( taxa.geneend ),
                    str( taxa.taxa ),
                    str( taxa.taxastart ),
                    str( taxa.taxaend ),
                    str( taxa.score ),
                    str( taxa.uniref50 ),
                    str( taxa.uniref90 ),
                    str( taxa.hits ),
                    ]
    return orderedlist

# ---------------------------------------------------------------
# tests
# ---------------------------------------------------------------

if __name__ == "__main__":
    counter = 0
    for hit in iter_hits( sys.argv[1] ):
        print( hit.qseqid, hit.length, hit.taxonomy[2:5], hit.evalue )
        counter += 1
        if counter > 10:
            break
    counter = 0
    for contig, hits in iter_contig_hits( sys.argv[1] ):
        print( contig, len( hits ) )
        counter += 1
        if counter > 10:
            break
