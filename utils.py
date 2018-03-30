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
import re
import gzip
import bz2
from collections import OrderedDict

# ---------------------------------------------------------------
# ---------------------------------------------------------------
# GENERIC HELPER FUNCTIONS
# ---------------------------------------------------------------
# ---------------------------------------------------------------

def say( *args ):
    print( " ".join( map( str, args ) ), file=sys.stderr )

def die( *args ):
    args = ["LETHAL ERROR:"] + list( args )
    say( *args )
    sys.exit( "EXITING." )

def path2name( p ):
    return os.path.split( p )[1].split( "." )[0]

def name2path( name, root=".", ext="" ):
    return os.path.join( root, name+ext )

def try_open( path, *args ):
    """ Open a file; fail gracefully """
    fh = None
    try:
        if path.endswith( ".gz" ):
            fh = gzip.GzipFile( path, *args )
        elif path.endswith( ".bz2" ):
            fh = bz2.BZ2File( path, *args )
        else:
            fh = open( path, *args )
    except:
        sys.exit( "Can't open file: {}".format( path ) )
    return fh

def describe( text, width=80, margin=2 ):
    # remove flanking whitespace
    text = text.strip( )
    lines = text.split( "\n" )
    margin = " " * margin
    # title
    rule = "=" * width
    newlines = [rule, margin + lines[0], rule, "\n"]
    newline = margin
    # description lines
    for line in lines[2:]:
        line = line.strip( )
        if line == "":
            newlines.append( newline )
            newlines.append( "\n" )
            newline = margin
            continue
        words = line.split( )
        for word in words:
            if len( word ) > width:
                newlines.append( newline )
                newlines.append( word )
                newline = margin
            elif len( newline + " " + word ) > width:
                newlines.append( newline )
                newline = margin + word
            else:
                newline += (" " if newline != margin else "") + word
    if len( newline ) > 0:
        newlines.append( newline )
    newlines += ["\n", rule]
    return "\n".join( [k for k in newlines] )

def read_contig_lengths( fasta ):
    data = OrderedDict( )
    header = None
    with try_open( fasta ) as fh:
        for line in fh:
            line = line.strip( )
            if line[0] == ">":
                header = line[1:].split( )[0]
                data[header] = 0
            else:
                data[header] += len( line )
    return data
 
# ---------------------------------------------------------------
# ---------------------------------------------------------------
# WORKING WITH BLAST
# ---------------------------------------------------------------
# ---------------------------------------------------------------

# ---------------------------------------------------------------
# blast constants
# ---------------------------------------------------------------

c_blast_fields = [
    ["qseqid"   , str],
    ["sseqid"   , str],
    ["qlen"     , int],
    ["slen"     , int],
    ["length"   , int],
    ["qstart"   , int],
    ["qend"     , int],
    ["sstart"   , int],
    ["send"     , int],
    ["pident"   , float],
    ["positive" , int],
    ["gaps"     , int],
    ["evalue"   , float],
    ["bitscore" , float],
    ["sstrand"  , str],
]
c_blast_format_string = " ".join( ["6"] + [fname for [fname, ftype] in c_blast_fields] )
# blast default is 500, which is sometimes too small for long contigs
c_max_target_seqs = 10000

# ---------------------------------------------------------------
# blast hit object
# ---------------------------------------------------------------

class Hit( ):

    """
    Processes the information from a single BLAST row.
    Row is provided already split by the csv reader.
    Hit MUST be compatible with the blast format defined above.

    Additional information is pulled into the hit object by parsing the
    subject sequence identifier. The expected format is:
    
    0 : unique gene id
    1 : taxon name
    2+: one or more functional annotations in "system=id" format
    """

    def __init__( self, blast_row ):
        if len( blast_row ) != len( c_blast_fields ):
            die( "inconsistent blast row: {}".format( str( blast_row ) ) )
        # pull values from blast line and coerce to appropriate types
        for [fname, ftype], value in zip( c_blast_fields, blast_row ):
            setattr( self, fname, ftype( value ) )
        # overrides
        self.sstrand = "-" if self.sstrand == "minus" else "+"
        # derived coverage stats
        self.scov = ( abs( self.send - self.sstart ) + 1 ) / float( self.slen )
        self.qcov = ( abs( self.qend - self.qstart ) + 1 ) / float( self.qlen )
        # special scoverage that doesn't penalize hanging off contig end
        if self.sstrand == "-":
            sstart = self.slen - self.sstart + 1
            send   = self.slen - self.send + 1
        else:
            sstart = self.sstart
            send   = self.send
        self.ltrim = max( 0, sstart - self.qstart )
        self.rtrim = max( 0, self.slen - sstart - self.qlen + self.qstart )
        self.scov_modified = (send - sstart + 1) / float( self.slen - self.ltrim - self.rtrim )
        # waafle score
        self.waafle_score = self.scov_modified * self.pident / 100.0 
        # values extracted from subject header
        items = self.sseqid.split( "|" )
        if len( items ) < 2:
            die( "bad subject id header:", self.sseqid )
        self.geneid = items[0]
        self.taxon = items[1]
        # process possible annotations
        self.annotations = {}
        if len( items ) > 2:
            for k in items[2:]:
                system, name = k.split( "=" )
                self.annotations[system] = name

# ---------------------------------------------------------------
# blast utils
# ---------------------------------------------------------------

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

# ---------------------------------------------------------------
# ---------------------------------------------------------------
# WORKING WITH GFFS
# ---------------------------------------------------------------
# ---------------------------------------------------------------

# ---------------------------------------------------------------
# constants
# ---------------------------------------------------------------

c_gff_fields = [
    ["seqname"   , str],
    ["source"    , str],
    ["feature"   , str],
    ["start"     , int],
    ["end"       , int],
    ["score"     , float],
    ["strand"    , str],
    ["frame"     , str],
    ["attribute" , str],
]

# ---------------------------------------------------------------
# gff line object (locus)
# ---------------------------------------------------------------

class Locus( ):

    def __init__( self, gff_row, attach_annotations=True ):
        # no name by default
        self.name = None
        # gff fields
        if len( gff_row ) != len( c_gff_fields ):
            die( "Bad GFF row:", gff_row )
        for [fname, ftype], value in zip( c_gff_fields, gff_row ):
       	    setattr( self, fname, ftype( value ) if value != "." else value )
        # annotations
        self.annotations = {}
        self.annotation_scores = {}
        for item in self.attribute.split( "; " ):
            match = re.search( "^(.*?) \"(.*)\"$", item )
            if attach_annotations and match:
                system, value = match.groups( )
                self.annotations[system] = value
                self.annotation_scores[system] = None
        # ignore status
        self.ignore = False

    def __len__( self ):
        return abs( self.end - self.start ) + 1

# ---------------------------------------------------------------
# gff utils
# ---------------------------------------------------------------

"""
Note: Currently we'll ignore any annotations in a user-supplied
GFF and focus on those assigned by WAAFLE.
Could update this later.
"""

def iter_loci( gff_file, attach_annotations=True ):
    with try_open( gff_file ) as fh:
        for row in csv.reader( fh, csv.excel_tab ):
            yield Locus( row, attach_annotations=attatch_annotations )

def iter_contig_loci( gff_file, attach_annotations=True ):
    contig, loci = None, []
    with try_open( gff_file ) as fh:
        for row in csv.reader( fh, csv.excel_tab ):
            l = Locus( row, attach_annotations=attach_annotations )
            if contig is not None and l.seqname != contig:
                yield contig, loci
                # reset
                loci = []
            contig = l.seqname
            loci.append( l )
        # last case cleanup
        yield contig, loci

# ---------------------------------------------------------------
# ---------------------------------------------------------------
# WORKING WITH TAXONOMY
# ---------------------------------------------------------------
# ---------------------------------------------------------------

# ---------------------------------------------------------------
# constants
# ---------------------------------------------------------------

c_unknown = "Unknown"
c_root    = "r__Root"

# ---------------------------------------------------------------
# taxonomy object
# ---------------------------------------------------------------

class Taxonomy( ):

    def __init__( self, path ):
        self.parents = {}
        self.children = {}
        with try_open( path ) as fh:
            for clade, parent in csv.reader( fh, csv.excel_tab ):
                self.parents[clade] = parent
                self.children.setdefault( parent, set( ) ).add( clade )
        # memoizer for get_leaf_count( )
        self.known_leaf_count = {}

    def get_parent( self, clade ):
        return self.parents.get( clade, c_root )

    def get_children( self, clade ):
        return self.children.get( clade, set( ) )

    def get_lineage( self, clade ):
        l = [clade]
        while l[-1] != c_root:
            parent = self.get_parent( clade )
            l.append( parent )
            clade = parent
        l.reverse( )
        return l

    def get_lca( self, *clades ):
        lca = c_root
        lineages = [self.get_lineage( c ) for c in clades]
        min_depth = min( [len( l ) for l in lineages] )
        for i in range( 0, min_depth ):
            level = {l[i] for l in lineages}
            if len( level ) == 1:
                lca = list( level )[0]
            else:
                break
        return lca

    def get_tails( self, clades, lca ):
        tails = []
        for c in clades:
            l = self.get_lineage( c )
            l.reverse( )
            t = []
            for c2 in l:
                if c2 == lca:
                    break
                else:
                    t.append( c2 )
            t.reverse( )
            tails.append( t )
        return tails

    def get_sisters( self, clade ):
        ret = set( )
        parent = self.get_parent( clade )
        for child in self.get_children( parent ):
            if child != clade:
                ret.add( child )
        return ret

    def get_leaf_count( self, clade ):
        ret = None
        if clade in self.known_leaf_count:
            ret = self.known_leaf_count[clade]
        elif clade not in self.children:
            ret = 1
        else:
            ret = 0
            for c in self.children[clade]:
                ret += self.get_leaf_count( c )
        self.known_leaf_count[clade] = ret
        return ret

# ---------------------------------------------------------------
# ---------------------------------------------------------------
# WORKING WITH INTERVALS
# ---------------------------------------------------------------
# ---------------------------------------------------------------

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

def calc_overlap( a1, a2, b1, b2 ):
    """ compute overlap between two intervals """
    a1, a2 = sorted( [a1, a2] )
    b1, b2 = sorted( [b1, b2] )
    if b1 > a2 or a1 > b2:
        return 0
    else:
        outleft, inleft, inright, outright = sorted( [a1, a2, b1, b2] )
        denom = min( ( a2 - a1 + 1 ), ( b2 - b1 + 1 ) )
        return ( inright - inleft + 1 ) / float( denom )

# ---------------------------------------------------------------
# ---------------------------------------------------------------
# TESTING
# ---------------------------------------------------------------
# ---------------------------------------------------------------

if __name__ == "__main__":
    pass
