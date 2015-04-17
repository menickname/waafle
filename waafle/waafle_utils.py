#!/usr/bin/env python

"""
WAAFLE MODULE: waafle_utils.py

Shared components of the waafle system.

Authors:
Tiffany Hsu
Eric Franzosa
"""

from __future__ import print_function # python 2.7+ required
import os, sys, csv, argparse

# ---------------------------------------------------------------
# constants
# ---------------------------------------------------------------

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
    ["positives", int],
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

# ---------------------------------------------------------------
# classes for working with hits (here to force equivalenence with output)
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
        self.scov = ( self.send - self.sstart + 1 ) / float( self.slen )
        self.qcov = ( self.qend - self.qstart + 1 ) / float( self.qlen )
        # special scoverage that won't penalize hanging off contig end
        self.ltrim = max( 0, self.sstart - self.qstart )
        self.rtrim = max( 0, self.slen - self.sstart - self.qlen + self.qstart )
        self.scov_modified = ( self.send - self.sstart + 1 ) \
            / float( self.slen - self.ltrim - self.rtrim )
        # values extracted from chocophlan header
        chocoitems = self.sseqid.split( c_choco_header_delim )
        self.taxid = chocoitems[5]
        self.taxonomy = chocoitems[6].split( c_tax_delim )
        self.uniref90 = chocoitems[7]
        self.uniref50 = chocoitems[8]

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

def hit_sort( hits ):
    return sorted( hits, key=lambda k: k.length, reverse=True )

def iter_contig_hits( blastoutfile ):
    """
    Iterate through hits by contig (assumes file is sorted by query)
    """
    contig, hits = None, []
    with try_open( blastoutfile ) as fh:                 
        for row in csv.reader( fh, dialect="excel-tab" ):
            hit = Hit( row )
            if contig is not None and hit.qseqid != contig:
                yield contig, hit_sort( hits )
                # reset
                hits = []
            contig = hit.qseqid
            hits.append( hit )
        # last case cleanup
        yield contig, hit_sort( hits )

# ---------------------------------------------------------------
# tests
# ---------------------------------------------------------------

if __name__ == "__main__":
    for contig, hits in iter_contig_hits( sys.argv[1] ):
        print( contig )
        for hit in hits:
            print( "\t", hit.length )
