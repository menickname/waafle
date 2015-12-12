#!/usr/bin/python

"""
WAAFLE VALIDATION: fakemake.py

Authors:
Tiffany Hsu
Eric Franzosa

This script generates the synthetic contigs from a list of GCFs.
Taxa pairs are drawn from /n/huttenhower_lab/data/2014_repophlan/repophlan_microbes.txt
Run with the "-h" flag for usage help.
"""

from __future__ import print_function # python 2.7+ required
import argparse, random, re, math, pickle
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import waafle_utils as wu

# ---------------------------------------------------------------
# constants
# ---------------------------------------------------------------
# testing
#random.seed( 1 )
mingenelength=100
mincontiglen=300
mediancontiglen=905


# ---------------------------------------------------------------
# functions
# ---------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument( '--recipient', help='This is the recipient GCF number.' )
parser.add_argument( '--donor', help='This is the donor GCF number.' )
parser.add_argument( '--reciptaxa', help='This is the complete recipient taxonomy.' )
parser.add_argument( '--donortaxa', help='This is the complete donor taxonomy.' )
parser.add_argument( '--taxadiff', help='This is the phylogenetic difference between the 2.')
parser.add_argument( '--ngenes', type=int, default=1, help='This is the number of contigs we would like to generate for this donor-recipient pair.')
#I have not implemented this, but it could be a solution for generating the same donor/recip with multiple genes to see if the donor/recip pair is the problem, or not
args = parser.parse_args()

# ---------------------------------------------------------------
# functions
# ---------------------------------------------------------------
def load_scaffolds ( GCF ):
    """ returns dict with key=scaffold GI, item = SeqRecord """
    fileloc = '/n/huttenhower_lab/data/2014_repophlan/repophlan_microbes/fna/'
    genomefileloc = fileloc  + GCF + '.fna'
    dictScaffolds = SeqIO.index(genomefileloc, 'fasta')
    return dictScaffolds 

def load_bog ( GCF ):
    """ returns dict with key=gene GI, item=SeqRecord """
    fileloc = '/n/huttenhower_lab/data/2014_repophlan/repophlan_microbes/ffn/'
    bogfileloc = fileloc + GCF + '.ffn'
    dictBog = SeqIO.index(bogfileloc, 'fasta')
    return dictBog

def getcoords ( GI ):
    """ return coordinates in sorted order """
    # gi like ">bug|:1-1213,1213-154242 annotations" 
    # first string annotation
    # then isolate coordinates
    # then, in rare case of ,-separate coords, take the first pair
    coordinates = GI.split( " " )[0].split(':')[1].split( "," )[0]
    if re.search('c', coordinates):
        end, start = int(coordinates.split('-')[0][1:]), int(coordinates.split('-')[1])
    else:
        start, end = int(coordinates.split('-')[0]), int(coordinates.split('-')[1])
    return start, end

def bog2coords ( dictBog ):
    """ return dict with key=scaffold name, value=[ start, end, geneGI ] """
    coords = {}
    for geneGI in dictBog:
        start, end = getcoords(geneGI)
        scaffoldname = geneGI.split(':')[0]
        coords.setdefault(scaffoldname, []).append([start, end, geneGI])
    return coords

def calc_overlap( onestart, oneend, twostart, twoend ):
    """ Calculate overlap between two hits or genes. """
    if onestart - twoend > 0 or twostart - oneend > 0:
        return 0
    else:
        if onestart < twostart and twoend < oneend or onestart > twostart and oneend < twoend:
            divisor = float( min( oneend - onestart, twoend - twostart ) )
        else:
            divisor = float( min( oneend - onestart + 1, twoend - twostart + 1 ) )
        coord_sorted = sorted( [onestart, oneend, twostart, twoend] )
        overlap = float( ( coord_sorted[2] - coord_sorted[1] )/divisor )
        return overlap

def find_flanking_genes( rgenename, rcoords, rGI ):
    rgenes_before, rgenes_after, status = "", "", False 
    sorted_rgenes = list( enumerate( sorted( rcoords[rGI] ) ) )
    for i, info in sorted_rgenes:
        genename = info[2]
        if rgenename == genename:
            if i == 0:
                rgenes_before = None
                rgenes_after = sorted_rgenes[i+1: i+2]
            elif i == len( sorted_rgenes ) - 1:
                rgenes_before = sorted_rgenes[i-1: i]
                rgenes_after = None
            else:
                rgenes_before = sorted_rgenes[i-1: i]
                rgenes_after = sorted_rgenes[i+1: i+2]
                status = True
    return status, rgenes_before, rgenes_after

def make_contig( donor, recipient, dbog, rbog, rcoords, scaffolds ):
    rgenename, dGI, rGI = "", "", ""
    dgenename = random.choice(dbog.keys()) #choose a donor gene
    status = False
    while status == False: #repick recipient gene if it does not have flanking genes
        rgenename = random.choice(rbog.keys()) #choose a recipient gene
        dGI, rGI = dgenename.split(':')[0], rgenename.split(':')[0]
        rstart, rend = getcoords(rgenename) #get coords of recipient gene
        status, rgene_before, rgene_after = find_flanking_genes( rgenename, rcoords, rGI )
    contigstart, contigend = rgene_before[0][1][0], rgene_after[0][1][1]
    newseq = scaffolds[rGI].seq[ contigstart - 1: rstart - 1] + dbog[dgenename].seq + scaffolds[rGI].seq[rend: contigend]
    new_SeqRec = SeqRecord( newseq )
    new_SeqRec.description = 'donor:' + donor + '==' + dgenename + '||recipient:' + recipient + '==' + rgenename
    genelist = rgene_before + rgene_after
    return new_SeqRec, genelist

def get_coords( coords ):
    if re.search( '^:c', coords[0]):
        start = int( coords[1] )
        end = int( coords[0][2:] )
        frame = '-'
    else:
        start = int( coords[0][1:] )
        end = int( coords[1] )
        frame = '+'
    return [start, end, frame]

def printgff( gffrow ):
    """
    Format the gff class into a ordered list for printing.
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

def format_gff( dGI, rGI, genelist, contig ):
    dcoords = get_coords( dGI.split('|')[4].split('-') )
    rcoords = get_coords( rGI.split('|')[4].split('-') )
    first_gene = genelist[0]
    third_gene = genelist[1]
    first_name, third_name = first_gene[0], third_gene[0]
    first_coords, third_coords = get_coords( first_gene[1][2].split('|')[4].split('-') ), get_coords( third_gene[1][2].split('|')[4].split('-') )
    first_GI, third_GI = first_gene[1][2], third_gene[1][2]
    
    first_start = first_coords[0]-first_coords[0]+1
    first_end = first_coords[1]-first_coords[0]+1
    d_start = first_end+(rcoords[0]-first_coords[1])
    d_end = d_start+(dcoords[1]-dcoords[0])
    third_start = d_end+(third_coords[0]-rcoords[1])
    third_end = third_start+(third_coords[1]-third_coords[0])
    geneinfo = [[first_name, first_start, first_end, first_coords[2], first_GI],
                ['donor_' + str( int( first_name )+1 ), d_start, d_end, '+', dGI + '_' + rGI],
                [third_name, third_start, third_end, third_coords[2], third_GI],
                ] 
    genelist = []
    for item in geneinfo:
        gene = wu.GFF( [] )
        gene.seqname = contig
        gene.source = 'WAAFLE-VALIDATION'
        gene.feature = 'CDS'
        gene.start = item[1]
        gene.end = item[2]
        gene.score = 1
        gene.strand = item[3]
        gene.frame = 0 
        ID = 'ID=' + contig + '_' + str(item[0])
        desc = 'original_gene=' + item[4]
        gene.attribute = ';'.join([ID, desc])
        genelist.append( printgff( gene ) )
    return genelist

# ---------------------------------------------------------------
# main
# ---------------------------------------------------------------

def main():
    outname= "-".join( [args.recipient, args.donor] )

    fh_fnt = open( outname+".fnt", "w" )
    fh_ans = open( outname+"-ans.txt", "w" )
    fh_gff = open( outname+"-gff.txt", "w" )

    for i in range( args.ngenes ):

        # load the scaffolds and the genes
        scaffolds = load_scaffolds( args.recipient )
        rbog = load_bog( args.recipient )
        dbog = load_bog( args.donor )
        rcoords = bog2coords( rbog )
    
        # pick a donor and recipient genes and make a contig that contains surrounding genes
        mySeqRec, genelist = make_contig( args.donor, args.recipient, dbog, rbog, rcoords, scaffolds )
        donorinfo, recipinfo = mySeqRec.description.split('||')[0], mySeqRec.description.split('||')[1]
        dGI, rGI = donorinfo.split('==')[1], recipinfo.split('==')[1]

        # generate fasta
        mySeqRec.id = "%s|contig%05d" % ( outname, i )
        fh_fnt.write( ">" + mySeqRec.id + '\n')
        fh_fnt.write( str( mySeqRec.seq ) + '\n' )
        
        #generate answerkey
        fh_ans.write( '\t'.join( [str(x) for x in [args.taxadiff, mySeqRec.id, args.reciptaxa, args.donortaxa, len(mySeqRec.seq), dGI, rGI, genelist ] ] ) + '\n' )
        formatted_genelist = format_gff( dGI, rGI, genelist, mySeqRec.id )
        for line in formatted_genelist:
            fh_gff.write( '\t'.join( [str(x) for x in line] ) + '\n' )

    fh_fnt.close()
    fh_ans.close()
    fh_gff.close()

if __name__ == "__main__":
    main()
