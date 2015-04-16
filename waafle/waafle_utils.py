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

c_choco_header_delim = "|"
c_tax_delim = "."

# ---------------------------------------------------------------
# blast output configuration
# ---------------------------------------------------------------

blastfields = [
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
    ["positives", float],
    ["gaps", int],
    ["evalue", float],
    ["bitscore", float],
]
blast_format_string = " ".join( ["6"] + [fname for [fname, ftype] in blastfields] )

# ---------------------------------------------------------------
# classes for working with hits (here to force equivalenence with output)
# ---------------------------------------------------------------

class Hit( ):
    """
    Processes the information from a single blast row;
    Row is provided already split by the csv reader.
    Hit MUST be compatible with the blast search defined above.

    Addition information is pulled into the hit object by parsing the
    hit chocophlan centroid's header. Chocophlan headers have the following
    fields delimited with a pipe ("|"):
    
    0: non-informative
    1: gi number
    2: non-informative
    3: refseq number
    4: coordinates
    5: ncbi tax id
    6: taxonomy ("."-delimitted)
    7: UniRef90 hit
    8: UniRef50 hit
    """
    def __init__( self, blastrow ):
        assert len( blastrow ) == len( blastfields ), \
            "inconsistent blast row {}".format( str( blastrow ) )
        # pull values from blast line and coerce to appropriate types
        for [fname, ftype], value in zip( blastfields, blastrow ):
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

def iterhits( blastoutfile ):
    """
    Iterate through the hits in a blast file
    """
    try:
        fh = open( blastoutfile )
    except:
        sys.exit( "Can't open blast out file: {}".format( blastoutfile ) ) 
    for row in csv.reader( fh, dialect="excel-tab" ):
        yield Hit( row )

# ---------------------------------------------------------------
# tests
# ---------------------------------------------------------------

if __name__ == "__main__":
    for hit in iterhits( sys.argv[1] ):
        print( hit.qseqid, hit.length, hit.taxonomy[2:5], hit.evalue )
