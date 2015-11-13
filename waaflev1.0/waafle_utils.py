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
    ["gene", int],
    ["strand", str],
    ["start", int],
    ["end", int],
    ["taxa", str],
    ["score", float],
    ["percid", float],
    ["genecov", float],
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
        
    def genenum ( self, attribute ):
        attritems = attribute.split( ';' )
	for items in attritems:
            if len( items.split( '=' ) ) == 2:
                label, descriptor = items.split( '=' )[0], items.split( '=' )[1]
                if label == 'ID':
                    lastindex = descriptor.rindex( '_' )
                    self.genenum = int( descriptor[lastindex + 1:] )


class Taxa( ):
    """
    Processes the information from a single line from waafle_orgscorer.py script;
    Row is provided already split by the csv reader.

    Details of the scored orgs file:
    Fields should be tab-separated.
    
    0: contig - name of the contig
    1: gene - number of the gene annotated in the contig
    2: strand - defined as + (forward) or - (reverse).
    3: start - start position for the taxa in the contig. May not start of gene.
    4: end - end position for the taxa in the contig. May not end of gene.
    5: taxa - taxa name at level chosen in waafle_orgscorer.py
    6, 7, 8 - values are generated by aggregating BLAST hits that correspond to the taxa listed.
    Hits are scored by multiplying percid and gene coverage (for that hit).
    The highest scoring hits, along with their percid and gene coverage, are assigned to the coordinates they cover.
    At the end, the scores, percid, and gene coverages are averaged across all coordinates.
    6: score - ranges from 0-1. Average score across all positions for aggregated hits.
    7: percid - ranges from 0-1. Average percent ID across all positions for aggregrated hits.
    8: genecov - ranges from 0-1. Average gene coverage across all positions for aggregated hits.
    """
    def __init__( self, taxainfo ):
        for [fname, ftype], value in zip( c_taxafields, taxainfo ):
            setattr( self, fname, ftype( value ) )

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
                gene.genenum( gene.attribute )
		if contig is not None and gene.seqname != contig:
			yield contig, genes
			genes = []
	        contig = gene.seqname
		genes.append( gene )
	# last case cleanup
	yield contig, genes

def iter_contig_taxa( orgfile ):
    """
    Iterate thrugh taxa by contig (in org-scores file)
    """
    contig, taxa = None, []
    with try_open( orgfile ) as fh:
        for row in csv.reader( fh, dialect="excel-tab" ):
            taxon = Taxa( row )
            if contig is not None and taxon.contig != contig:
                yield contig, taxa
                #reset
                taxa = []
            contig = taxon.contig
            taxa.append( taxon )
        # last case cleanup
        yield contig, taxa
		
def calc_overlap( onestart, oneend, twostart, twoend ):
    """
    Calculate overlap between two hits or genes.
    """
    if onestart - twoend > 0 or twostart - oneend > 0:
        return 0
    else:
        coord_sorted = sorted( [onestart, oneend, twostart, twoend] )
        divisor = float( min( oneend - onestart + 1, twoend - twostart + 1 ) )
        overlap = float( ( coord_sorted[2] - coord_sorted[1] )/divisor )
        return overlap


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
