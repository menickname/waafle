#!/usr/bin/env/python

from __future__ import print_function
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import argparse
import waafle_utils as wu
import numpy as np
from matplotlib.pyplot import cm 
import random
import re

# ---------------------------------------------------------------
# functions
# ---------------------------------------------------------------

def get_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument(
        "-i", "--waafle",
        required=True,
        help="output from waafle"
        )
    parser.add_argument(
        "-refgenes", "--ref_genes",
        help="output for reference genes from Prodigal or NCBI in gff format, if available"
        )
    parser.add_argument(
        "-reforgs", "--ref_orgs",
        help="output for reference orgs from NCBI for synthetic contigs only, if available"
        )
    parser.add_argument(
        "-contigs", "--contiglist",
        help="list of contigs to output plots for"
        )
    parser.add_argument(
        "-taxa", "--taxalevel",
        help="taxa level"
        )
    args = parser.parse_args()
    return args

#----------------------------------------------------------------
# functions
#----------------------------------------------------------------



# ---------------------------------------------------------------
# main
# ---------------------------------------------------------------

def main():
    args = get_args()

    taxalevels = {'k':0, 'p':1, 'c':2, 'o':3, 'f':4, 'g':5, 's':6, 't':7}
    # Parse taxa data
    dict_mytaxa = {}
    for contig, taxalist in wu.iter_contig_taxa( args.waafle ):
        dict_mytaxa[contig] = taxalist

    # Parse gff for from ref (either Prodigal or AnswerKey)
    if args.ref_genes:
        dict_refgenes = {}
        for contig, gfflist in wu.iter_contig_genes( args.ref_genes ):
            dict_refgenes[contig] = gfflist

    # Parse taxa data for AnswerKey IF provided
    if args.ref_orgs:
        dict_reforgs = {}
        for astrline in open( args.ref_orgs ): 
            aastrline = astrline.strip().split('\t')
            index = taxalevels[args.taxalevel]
            donor = aastrline[3].split('|')[index]
            recip = aastrline[2].split('|')[index]
            dict_reforgs[ aastrline[1] ] = [recip, donor]

    # Loop through list of contigs
    for astrline in open( args.contiglist ):
        contig = astrline.strip()
        taxalist = dict_mytaxa[contig]
  
        # Arrange waafle genes in data structure
        orglabels = []
        dict_waaflegenes = {}
        for taxa in taxalist:
            orglabels.append( taxa.taxa )
            if taxa.strand == '+':
                dict_waaflegenes[taxa.gene] = [taxa.genestart, taxa.geneend, taxa.uniref50, taxa.uniref90]
            else:
                dict_waaflegenes[taxa.gene] = [taxa.geneend, taxa.genestart, taxa.uniref50, taxa.uniref90]

        # Make colors for taxa
        orgset = set(orglabels)
        colors = cm.Accent(np.linspace(0, 1, len(orgset)))
        dict_taxacolor = {}
        for i in range( len( list(orgset) ) ):
            dict_taxacolor[list(orgset)[i]] = colors[i]

        # Arrange taxa in data structure
        dict_waafletaxa = {}
        contiglen = 0
        for current_org in dict_taxacolor.keys():
            scorelist, startlist, endlist = [], [], []
            for taxa in taxalist:
                contiglen = float( taxa.length )
                if taxa.taxa == current_org:
                    newstart = taxa.taxastart.split(',')
                    newend = taxa.taxaend.split(',')
                    if taxa.strand == '+':
                        startlist += newstart
                        endlist += newend
                    else:
                        startlist += newend
                        endlist += newstart
                    scorelist += [taxa.score]*len(newstart) 
            dict_waafletaxa[current_org] = [scorelist, startlist, endlist]        

        # Initiate the plot      
        taxarows = len(orgset)/6 + 1
        generows = len( dict_waaflegenes.keys() )/6 + 1
        rows = 4 + taxarows + generows
        print( rows, 'rows' )
        fig = plt.figure(figsize=(12, rows*3))
        h_waafle = plt.subplot2grid( (rows,1), (0,0), rowspan=2, colspan=1 )
        h_waaflegenes = plt.subplot2grid( (rows,1), (2,0), rowspan=1, colspan=1 )        

        # Plot taxa and legend
        h_legend = plt.subplot2grid( (rows,1), (4,0), rowspan=taxarows, colspan=1)
        h_legend.set_xlim( 0, contiglen )
        h_legend.set_ylim( 0, len(dict_taxacolor)*2 )
        ylab = len(dict_taxacolor)*2
        for taxa in dict_waafletaxa.keys(): 
            scores, starts, ends = dict_waafletaxa[taxa][0], [int(x) for x in dict_waafletaxa[taxa][1]], [int(y) for y in dict_waafletaxa[taxa][2]]
            headlen = contiglen/50
            diffs = list(np.array(ends) - np.array(starts) - headlen )
            print( contig, taxa, starts, scores, diffs, headlen )
            for i in range( len( scores ) ):
                h_waafle.arrow( starts[i], scores[i], diffs[i], 0, width=0.02, color=dict_taxacolor[taxa], head_width=0.02, head_length=headlen, alpha=0.4 )
            ylab -= 1
            p = patches.Rectangle( (0, ylab), contiglen/20, 1, color=dict_taxacolor[taxa], alpha=0.4 )
            h_legend.add_patch(p)
            h_legend.text( contiglen/20, ylab + 0.5, str(taxa), fontsize=contiglen/200 )
            ylab -= 1
        h_waafle.set_ylim( -0.1, 1.1 )
        h_waafle.set_xlim( 0, contiglen ) 
        h_waafle.set_ylabel( 'scores' )
        h_waafle.get_xaxis().set_ticks([])
        h_legend.set_axis_off()
                   
        # Plot waafle genes
        score = 1
        ylab = len( dict_waaflegenes.keys() )*2
        h_geneleg = plt.subplot2grid( (rows,1), (4+taxarows,0), rowspan=generows, colspan=1 )
        h_geneleg.set_xlim( 0, contiglen )
        h_geneleg.set_ylim( 0, len( dict_waaflegenes.keys() )*2 )
        for gene in dict_waaflegenes.keys():
            if score <= 0:
                score = 1
            else:
                score -= 0.1
            diff = dict_waaflegenes[gene][1] - dict_waaflegenes[gene][0]
            h_waaflegenes.arrow( dict_waaflegenes[gene][0], score, diff, 0, width=0.04, color="blue", head_width=0.04, head_length=headlen, alpha=0.4 )
            ylab -= 1
            h_geneleg.text(0, ylab, 'Gene' + str(gene) + ' Uniref50:' + dict_waaflegenes[gene][2], fontsize=8 )
            ylab -= 0.5
            h_geneleg.text(0, ylab, 'Gene' + str(gene) + ' Uniref90:' + dict_waaflegenes[gene][3], fontsize=8 )
        
        h_waaflegenes.set_ylim( 0, 1 )
        h_waaflegenes.set_xlim( 0, contiglen )
        h_waaflegenes.get_xaxis().set_ticks([])
        h_waaflegenes.get_yaxis().set_ticks([0, 0.5, 1]) 
        h_waaflegenes.set_ylabel( 'WAAFLE-genes' )
        h_geneleg.set_axis_off()
        
        # Plot reference genes (Prodigal or AnswerKey) and taxa (if AnswerKey) 
        DR = False
        if args.ref_orgs:
            DR = True
    
        if args.ref_genes:
            # Initiate the plot
            h_refgenes = plt.subplot2grid( (rows,2), (3,0), rowspan=1, colspan=2 )
            gfflist = dict_refgenes[contig]
	    score = 1
            for gene in gfflist:
                ID = gene.attribute.split(';')[0]
                if re.search( 'donor', ID ) and DR == True:
			colors = 'red'
		elif not re.search('donor', ID ) and DR == True:
                        colors = 'blue'
                else:
                        colors = 'black'
                
                # Set y-axis
                if score <= 0:
                    score = 1
                else:
                    score -= 0.1		
                if gene.strand == '+':
                    h_refgenes.arrow( gene.start, score, gene.end-gene.start, 0, width=0.04, color=colors, head_width=0.04, head_length=headlen, alpha=0.4 )
                    if re.search('donor', ID ) and DR == True:
                        h_refgenes.text( (gene.start+gene.end)/2, score, 'donor:' + dict_reforgs[contig][1] )
                    elif not re.search('donor', ID ) and DR == True:
                        h_refgenes.text( (gene.start+gene.end)/2, score, 'recip:' + dict_reforgs[contig][0] )
                else:
                    h_refgenes.arrow( gene.end, score, gene.start-gene.end, 0, width=0.04, color=colors, head_width=0.04, head_length=headlen, alpha=0.4)
                    if re.search('donor', ID ) and DR == True:
                        h_refgenes.text( (gene.start+gene.end)/2, score, 'donor:' + dict_reforgs[contig][1] )
                    elif not re.search('donor', ID ) and DR == True:
                        h_refgenes.text( (gene.start+gene.end)/2, score, 'recip:' + dict_reforgs[contig][0] )
            
            h_refgenes.set_ylim( 0, 1 )
            h_refgenes.set_xlim( 0, contiglen )
            h_refgenes.get_yaxis().set_ticks([0, 0.5, 1])
            h_refgenes.set_xlabel( 'coordinates (bp)' )
            h_refgenes.set_ylabel( 'answers' )

        filename = contig + '.png'
        plt.savefig( filename )
        
if __name__ == "__main__":
    main()
