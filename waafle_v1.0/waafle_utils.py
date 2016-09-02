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

c_scorerfields = [
    ["contig", str],
    ["length", int],
    ["status", str],
    ["onescore", float],
    ["twoscore", float],
    ["taxa", str],
    ["synteny", str],
    ["uniref50", str],
    ["uniref90", str],
]

c_resultfields = [
    ["contig", str],
    ["length", int],
    ["uniref50", str],
    ["uniref90", str],
    ["pref_call", str],
    ["unambig_call", str],
    ["k", str],
    ["p", str],
    ["c", str],
    ["o", str],
    ["f", str],
    ["g", str],
    ["s", str],
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

"""
@codereview 9/2/2016

Split up the GFF attribute field into a dict of attributes
Interact with the dict via get, e.g. 

gff.attributes.get( "UniRef50", None )
"""

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
    
    #attributes are complex, format as dictionary
    def sepattr( self, attribute ):
        attritems = attribute.split( ';' )
        dict_attr = {}
        for items in attritems:
            if len(items.split('=')) == 2:
                label, descriptor = items.split('=')[0], items.split('=')[1]
                dict_attr[label] = descriptor
                if label == 'ID':
                    lastindex = descriptor.rindex( '_' )
                    self.genenum = int( descriptor[lastindex+1:] )
                if label == 'Uniref50':
                    self.uniref50 = descriptor
                if label == 'Uniref90':
                    self.uniref90 = descriptor
        self.attr = dict_attr

"""
@codereview 9/2/2016

def __repr__( self ):
	...
	return STRING

This is what will be called if you print an object

print( GFF )
"""

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

class Scores( ):
    """
    Processes the information from a single line after waafle_lgtscorer.py script;
    """
    def __init__( self, scorerinfo ):
        for [fname, ftype], value in zip( c_scorerfields, scorerinfo ):
            setattr( self, fname, ftype( value ) )

class Result( ):
    """
    Processes the information from a single line after waafle_aggregator.py script.
    """
    def __init__( self, resultinfo ):
        for [fname, ftype], value in zip( c_resultfields, resultinfo ):
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

"""
@codereview 9/2/2016
See if some of the overlap code can use the INode and genes functions
just to avoid duplication of code
"""
		
def calc_overlap( onestart, oneend, twostart, twoend ):
    """
    Calculate overlap between two hits or genes.
    """
    if onestart - twoend > 0 or twostart - oneend > 0:
        return 0
    else:
        divisor = float( min( oneend - onestart + 1, twoend - twostart + 1 ) )
        coord_sorted = sorted( [onestart, oneend, twostart, twoend] )
        overlap = float( ( coord_sorted[2] - coord_sorted[1] + 1 )/divisor )
        return overlap

def convert_strand( strandw ):
    strand = '-'
    if strandw == 'plus':
        strand = '+'
    return strand

def find_ends( indexscore_sort, status ):
    """
    Find start and stop sites for hits that correspond to a single taxon in a gene.
    There may be multiple starts and stops. The output corresponds to each other.
    For example, if a taxa covers a gene from bp 1 to 100, and 200 to 300, the output would be 1,200 and 100,300.
    """
    startlist = [indexscore_sort[0]]
    endlist = []
    for i in range( len( indexscore_sort ) - 1 ):
        findex = indexscore_sort[i]
        sindex = indexscore_sort[i + 1]
        if sindex - findex > 1:
            endlist.append( findex )
            if i < len( indexscore_sort ) - 2:
                startlist.append( sindex )
        if i == len( indexscore_sort ) - 2:
            endlist.append( sindex )
    if status == 'gene':
        return str(min(startlist)), str(max(endlist))
    else:
        return ','.join( [str(x) for x in startlist] ), ','.join( [str(y) for y in endlist] )

"""
@codereview 9/2/2016
I think we discussed replacing this with pure numpy operations for speed?
Looks like loops over coords are still in effect
"""

def score_hits( hitlist, genestart, geneend, status ):
    """
    Loop through set of hits that correspond to something (gene or taxon).
    Calculate a score for each hit. Score = Percent identity * Coverage.
    Assign the highest score (from any hit) to each base position covered by all hits.
    Calculate a final score by averaging these scores.
    """
    dict_indexscore = {}
    groupcov = 0
    finalscore, finalpercid, finalcov, calcstart, calcend = 0, 0, 0, 'NA', 'NA'
    genelen = geneend - genestart + 1
    if len( hitlist ) > 0:
        for hit in hitlist:
            cov = hit.scov_modified #scov
            #cov = min( abs(hit.qend - hit.qstart) + 1, genelen )/float( genelen ) #gcov
            percid = hit.pident/float( 100 )
            score = float( cov ) * percid
            info = [score, percid, cov]
            start = max( hit.qstart, genestart )
            end = min( hit.qend, geneend )
            for coordinate in range( start, end + 1 ):
                if dict_indexscore.get( coordinate, [0, 0, 0] ) != [0, 0, 0]:
                    oldscore = dict_indexscore[coordinate][0]
                    if oldscore < score:
                        dict_indexscore[coordinate] = info
                else:
                    dict_indexscore[coordinate] = info

        finalpercid = np.sum( [x[1][1] for x in dict_indexscore.items()] )/float( genelen )
        finalcov = np.sum( [x[1][2] for x in dict_indexscore.items()] )/float( genelen ) #subject coverage
        #finalcov = len( dict_indexscore.keys() )/float( genelen ) #gcov
        finalscore = np.sum( [x[1][0] for x in dict_indexscore.items()] )/float( genelen ) 
        index_sort = sorted( dict_indexscore.keys() )
        calcstart, calcend = find_ends( index_sort, status )
    return finalscore, finalpercid, finalcov, calcstart, calcend

def hits2genes( genestart, geneend, genesign, hitlist, overlap_thresh, scov): #Add in assign Unirefs
    """
    Count how many hits are associated with each gene with filters.
    Annotate genes with uniref50 and uniref90 annotations, as well as number of hits included.
    Return the gene coordinates (start and end), strand, and count as a list.
    """
    if genesign == '+':
        wsign = 'plus'
    elif genesign == '-':
        wsign = 'minus'
    else:
        wsign = '.'
    genehits = []
    counter, uniref50list, uniref90list, info = 0, [], [], []
    for hit in hitlist:
        overlap = calc_overlap( genestart, geneend, hit.qstart, hit.qend)
        if genesign == '.':
            wsign = hit.sstrand
        if overlap >= overlap_thresh and hit.sstrand == wsign and hit.scov_modified >= scov:
            genehits.append( hit )
            counter += 1
            uniref50, uniref90 = hit.uniref50[ hit.uniref50.rindex('_')+1 : ], hit.uniref90[ hit.uniref90.rindex('_')+1 : ]
            uniref50list.append( uniref50 )
            uniref90list.append( uniref90 )
    uniref50_format = ','.join( [ str(x) + ':' + str(y) for x, y in Counter( uniref50list ).most_common( 3 )] )
    uniref90_format = ','.join( [ str(x) + ':' + str(y) for x, y in Counter( uniref90list).most_common( 3 )] )
    info = [counter, uniref50_format, uniref90_format]
    return genehits, info

def write_gff( contig, genelist ):
    """
    For each contig, sort all genes by start coordinate regardless of strand.
    Output each gene in gff format to a new file.
    """
    gfflist = []
    counter = 1
    for gene in genelist:
        gffrow = GFF( [] )
        gffrow.seqname = contig
        gffrow.source = 'WAAFLE'
        gffrow.feature = 'CDS'
        gffrow.start = gene[0]
        gffrow.end = gene[1]
        gffrow.score = gene[3]
        gffrow.strand = gene[2]
        gffrow.frame = '0'
        ID = 'ID=' + contig + '_' + str(counter)
        Uniref50 = 'Uniref50=' + gene[4]
        Uniref90 = 'Uniref90=' + gene[5]
        gffrow.attribute = ';'.join( [ID, Uniref50, Uniref90] )
        counter += 1
        gffline = print_gff( gffrow )
        gfflist.append( gffline )
    return gfflist

def print_gff( gffrow ):
    """
    Format the gff class into an ordered list for printing.
    """
    gfflist = [ gffrow.seqname,
                gffrow.source,
                gffrow.feature,
                gffrow.start,
                gffrow.end,
                gffrow.score,
                gffrow.strand,
                gffrow.frame,
                gffrow.attribute
                ]
    return gfflist

def filter_genes( genelist, hitlist, lap_h, length, scov_genes, scov_hits ):
    """
    First, group hits into genes.
    Second, filter genes based on length, scoverage, or # of BLAST hits.
    """
    genelist_filtered = []
    for start, end, strand in genelist:
        gene_hits, gene_info = hits2genes( start, end, strand, hitlist, lap_h, scov_hits )
        genelen = end - start + 1
        if len(gene_hits) == 0:
            gene_score, gene_percid, gene_cov, gene_starts, gene_ends = 0, 0, 0, 0, 0
        else:
            gene_score, gene_percid, gene_cov, gene_starts, gene_ends = score_hits( gene_hits, start, end, 'gene' )
        if gene_cov >= scov_genes and genelen >= length:
            genelist_filtered.append( [start, end, strand, gene_score, gene_info[1], gene_info[2]] )
    return genelist_filtered

def coords2groups( coordslist, overlap_thresh ):
    """
    Bin a list of hits into groups by start/end coordinates.
    """
    groups = []
    for start, end, strand in coordslist:
        group_add = False
        if len( groups )== 0:
            groups.append( [start, end, strand] )
        else:
            for i in range( len( groups ) ):
                overlap = calc_overlap( groups[i][0], groups[i][1], start, end )
                if overlap >= overlap_thresh:
                   coord_sorted = sorted( [groups[i][0], groups[i][1], start, end] )
                   groups[i] = [coord_sorted[0], coord_sorted[3], strand]
                   group_add = True
                if i == len( groups ) - 1 and group_add == False:
                   groups.append( [start, end, strand] )
    return groups

def groups2genes( groups, overlap_thresh ):
    """
    Sort groups by start/end coordinates.
    Bin a list of groups into genes by start/end coordinates (of groups).
    """
    genes = []
    #Sort groups by end and start coordinates
    for gene in groups:
        length = int(gene[1]) - int(gene[0])
        gene.append( length )
    groups.sort( key=itemgetter( 3 ), reverse=True )
    group_sorted = [x[0:3] for x in groups]
    if len(groups) == 1:
        genes = group_sorted
    else:
        genes = coords2groups( group_sorted, overlap_thresh )
    return genes

def print_taxa( taxa ):
    """
    Format the taxa class into an ordered list for printing.
    """
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



def split_info( info ):
    """
    Split the info from 1 column in the final output.
    """
    status, onebug, twobug, taxa, synteny = info.split(';')
    if re.search( 'NoLGT', status ):
        if re.search( ',', taxa ):
            multipletaxa = taxa.split('|')
            taxa = multipletaxa
        else:
            taxa = [taxa]
    else:
        taxalist = []
        multipletaxa = re.split( '[><\?]', taxa )
        for i in range(0, len( multipletaxa )-1, 2):
            taxaone, taxatwo = multipletaxa[i].split('|'), multipletaxa[i+1].split('|') 
            taxalist += taxaone
            taxalist += taxatwo
        taxa = taxalist
    return status, onebug, twobug, taxa, synteny

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
