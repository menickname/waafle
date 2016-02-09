#!/usr/bin/env python

"""
WAAFLE MODULE: waafle_orgscorer.py

Authors:
Tiffany Hsu
Eric Franzosa

This script combines gff output with BLAST hits and annotates genes with microbial taxa.
Run with the "-h" flag for usage help.

"""

from __future__ import print_function # python 2.7+ required
import os, sys, csv, argparse
import waafle_utils as wu
import numpy as np
from collections import Counter

# ---------------------------------------------------------------
# constants
# ---------------------------------------------------------------

c__dict_taxa = {
    "k": 0,
    "p": 1,
    "c": 2,
    "o": 3,
    "f": 4,
    "g": 5,
    "s": 6,
}

# ---------------------------------------------------------------
# functions
# ---------------------------------------------------------------
def get_args():
    """
    Get arguments passed to script
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument(
        "-g", "--gff",
        required=True,
        help="output from waafle_genecaller or user supplied gff",
        )
    parser.add_argument(
	    "-b", "--blast",
        required=True,
        help="output from waafle_search",
        )
    parser.add_argument(
        "-t", "--taxa",
        help="level of taxa to score",
        )
    parser.add_argument(
        "-o", "--out",
        default="waafle-scoredorgs.tsv",
        help="output for scored taxa",
        )
    parser.add_argument(
        "-scov_h", "--scov_hits",
        help="cutoff for gene coverage or subject coverage when grouping hits",
        default=0,
        type=float,
        )
    parser.add_argument(
        "-scov_g", "--scov_genes",
        help="cutoff for gene coverage or subject coverage for genes",
        default=0,
        type=float,
        )
    parser.add_argument(
        "-l", "--length",
        help="cutoff for gene length",
        default=0,
        type=float,
        )
    parser.add_argument(
        "-lap", "--overlap",
        help="amount of overlap to include hit in a gene",
        default=0.5,
        type=float,
        )
    parser.add_argument(
        "-merge-by", "--genemerge",
        help="options to merge genes for more conservative gene calls. Choose options 'none', 'strand', 'all'",
        default='none',
        type=str,
        )
    args = parser.parse_args()
    return args

def combine_genes( genelist, option, lap ):
    """
    Merge genes based on the option specified. 
    """
    if option == 'none':
        genecoords = []
        for gene in genelist:
            genecoords.append( [gene.start, gene.end, gene.strand] )
        return genecoords
    elif option == 'strand':
        posgenes, neggenes = [], []
        for gene in genelist:
            if gene.strand == '+':
                posgenes.append( [gene.start, gene.end, gene.strand] )
            else:
                neggenes.append( [gene.start, gene.end, gene.strand] )
        posgenes_group = wu.groups2genes( posgenes, lap )
        neggenes_group = wu.groups2genes( neggenes, lap )
        genecoords = sorted( posgenes_group + neggenes_group )
        return genecoords
    else:
        allgenes = []
        for gene in genelist:
            allgenes.append( [gene.start, gene.end, '.'] )
        allgenes_group = wu.groups2genes( allgenes, lap )
        genecoords = sorted( allgenes_group )
        return genecoords

def hits2taxa( contig, gene, genehits, taxalevel, uniref50, uniref90 ):
    """
    Group hits that correspond to specific taxa.
    If the gene did not have hits to begin with, output an "unknown" taxa with score 0.
    If the gene has hits that corresponds to orgs, score the orgs and output the taxa annotations.
    """
    dict_orghits = {}
    gene.sepattr( gene.attribute )
    for hit in genehits:
        taxa = hit.taxonomy[taxalevel]
        dict_orghits.setdefault( taxa, [] ).append( hit )
    taxalist = []
    if len( dict_orghits.keys() ) == 0:
        taxa = wu.Taxa( [] )
        taxa.contig, taxa.length = contig, 0
        taxa.gene, taxa.strand, taxa.genestart, taxa.geneend = gene.genenum, gene.strand, gene.start, gene.end
        taxa.taxa = "Unknown"
        taxa.score, taxa.percid, taxa.genecov, taxa.taxastart, taxa.taxaend = 0, 0, 0, 0, 0
        taxa.uniref50, taxa.uniref90 = uniref50, uniref90
        taxa.hits = 0
        taxalist.append( wu.print_taxa( taxa ) )
    else:
        for org in dict_orghits:
            orghits = dict_orghits[org]
            taxa = wu.Taxa( [] )
            taxa.contig, taxa.length = contig, hit.qlen
            taxa.gene, taxa.strand, taxa.genestart, taxa.geneend = gene.genenum, gene.strand, gene.start, gene.end
            taxa.taxa = org
            taxa.score, taxa.percid, taxa.genecov, taxa.taxastart, taxa.taxaend = wu.score_hits( orghits, int(gene.start), int(gene.end), 'taxa' )
            taxa.uniref50, taxa.uniref90 = uniref50, uniref90
            taxa.hits = len( dict_orghits[org] )
            taxalist.append( wu.print_taxa( taxa ) )
    return taxalist

# ---------------------------------------------------------------
# main
# ---------------------------------------------------------------

def main():
    """
    This script does:
    1) Groups hits into genes for contigs that have genes called.
    2) Groups hits into orgs per gene.
    3) Score each org.
    4) Prints out the taxa into a new file.
    """
    args = get_args()
    taxalevel = c__dict_taxa[args.taxa]

    fh = wu.try_open( args.out, "w" )
    writertaxa = csv.writer( fh, dialect="excel-tab" )
    fh2 = wu.try_open( 'waafle-genes_' + args.genemerge + '.gff', "w" )
    writergff = csv.writer( fh2, dialect="excel-tab" ) 
    writertaxa.writerow( ["contig", "contiglen", "genenum", "strand", "genestart", "geneend", "taxa", "taxastart", "taxaend", "score", "percid", "coverage", "uniref50", "uniref90", "orghitnum"] )
       
    #Merge genes based on options
    dict_contiggenes = {}
    for contig, genes in wu.iter_contig_genes( args.gff ):
        newgenecoords = combine_genes( genes, args.genemerge, args.overlap )
        dict_contiggenes[contig] = newgenecoords
        
    #Group taxa based on new genes
    for contig, contighits in wu.iter_contig_hits( args.blast ):
        if contig in dict_contiggenes:
            contiggenes = dict_contiggenes[contig]
            genelist = wu.filter_genes( contiggenes, contighits, args.overlap, args.length, args.scov_genes, args.scov_hits )
            gfflist = wu.write_gff( contig, genelist )
            for gff in gfflist:
                writergff.writerow( gff )
                newgene = wu.GFF( gff )
                start, end, strand = newgene.start, newgene.end, newgene.strand
                genehits, info = wu.hits2genes( start, end, strand, contighits, args.overlap, args.scov_hits )
                
                #score orgs 
                taxalist = hits2taxa( contig, newgene, genehits, taxalevel, info[1], info[2] )
                for taxa in taxalist:
                    writertaxa.writerow( taxa )
    fh.close()
    fh2.close()

if __name__ == "__main__":
    main()


