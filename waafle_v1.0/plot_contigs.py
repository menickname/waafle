#!/usr/bin/env/python

from __future__ import print_function
import matplotlib
matplotlib.use('Agg')

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
        "-scores", "--orgscorer",
        required=True,
        help="output from waafle orgscorer"
        )
    parser.add_argument(
        "-results", "--results",
        required=True,
        help="output from waafle aggregator"
        )
    parser.add_argument(
        "-reforgs", "--reforgs",
        default=False,
        help="output for reference orgs from NCBI for synthetic contigs only, if available"
        )
    parser.add_argument(
        "-refgenes", "--refgenes",
        default=False,
        help="output for reference genes from NCBI for synthetic contigs only, if available"
        )
    parser.add_argument(
        "-contigs", "--contiglist",
        help="list of contigs to output plots for"
        )
    parser.add_argument(
        "-taxa", "--taxalevel",
        help="taxa level"
        )
    parser.add_argument(
        "-o", "--type",
        help="Format of output (pdf, png)",
        default="pdf"
        )
    args = parser.parse_args()
    return args

#----------------------------------------------------------------
# constants
#----------------------------------------------------------------
taxalevels = {'k':0, 'p':1, 'c':2, 'o':3, 'f':4, 'g':5, 's':6, 't':7}

#----------------------------------------------------------------
# functions
#----------------------------------------------------------------

def arrange_taxagenes( taxalist ):
    dict_waaflegenes = {}
    for taxa in taxalist:
        if taxa.strand == '+':
            dict_waaflegenes[taxa.gene] = [taxa.genestart, taxa.geneend, taxa.uniref50, taxa.uniref90]
        else:
            dict_waaflegenes[taxa.gene] = [taxa.geneend, taxa.genestart, taxa.uniref50, taxa.uniref90]
    return dict_waaflegenes

def arrange_gffgenes( gfflist ):
    dict_genes = {}
    counter = 1
    for gene in gfflist:
        dict_genes[counter] = [gene.start, gene.end ]
        counter += 1
    return dict_genes

# ---------------------------------------------------------------
# main
# ---------------------------------------------------------------

def main():
    args = get_args()

    # Parse orgscorer data
    dict_mytaxa = {}
    for contig, taxalist in wu.iter_contig_taxa( args.orgscorer ):
        dict_mytaxa[contig] = taxalist

    # Parse answer data IF provided
    dict_reforgs = {}
    if args.reforgs:
        for astrline in open( args.reforgs ): 
            aastrline = astrline.strip().split('\t')
            index = taxalevels[args.taxalevel]
            donor = aastrline[3].split('|')[index]
            recip = aastrline[2].split('|')[index]
            dict_reforgs[ aastrline[1] ] = [recip, donor]
    dict_refgenes = {}
    if args.refgenes:
        for contig, gfflist in wu.iter_contig_genes( args.refgenes ): 
            dict_refgenes[contig] = gfflist

    # Parse final results into dict
    dict_results = {}
    for contiginfo in open( args.results ):
        line = contiginfo.strip().split('\t')
        if line[0] == 'contig': #remove header
            continue
        else:
            result = wu.Result( line )
            dict_results[ result.contig ] = result.__dict__[args.taxalevel]

    # Loop through list of contigs
    for astrline in open( args.contiglist ):
        contig = astrline.strip()
        taxalist = dict_mytaxa[contig]
        results = dict_results[contig]

        dict_waaflegenes = arrange_taxagenes( taxalist )

        # Make colors for taxa and arrange taxa in data structure
        orgset = set([])
        coloredorgs = wu.split_info( results )[3]
        bcolors = cm.RdYlBu(np.linspace(0, 1, len( coloredorgs )))
        gcolors = cm.gray( np.linspace(0, 1, 1 ) )
        dict_taxacolor, dict_waafletaxa, counter = {}, {}, 0
        contiglen = 0
        for info in taxalist:
            contiglen = float( info.length )
            abbrtaxa = info.taxa.split('|')[taxalevels[args.taxalevel]]
            if abbrtaxa in coloredorgs and abbrtaxa not in dict_waafletaxa.keys():
                dict_taxacolor[abbrtaxa] = bcolors[counter]
                counter += 1
            elif abbrtaxa not in coloredorgs: 
                dict_taxacolor[abbrtaxa] = gcolors[0] #colors[i]
            if abbrtaxa in dict_waafletaxa.keys():
                scorelist, startlist, endlist = [], [], []
                oldscores, oldstarts, oldends = dict_waafletaxa[ abbrtaxa ]
                if info.strand == '+':
                    startlist = info.taxastart.split(',')
                    endlist = info.taxaend.split(',')
                else:
                    startlist = info.taxaend.split(',')
                    endlist = info.taxastart.split(',')
                scorelist = [info.score]*len(startlist)
                dict_waafletaxa[ abbrtaxa ] = [ oldscores + scorelist, oldstarts + startlist, oldends + endlist ]
            else:
                scorelist, startlist, endlist = [], [], []
                if info.strand == '+':
                    startlist = info.taxastart.split(',')
                    endlist = info.taxaend.split(',')
                else:
                    startlist = info.taxaend.split(',')
                    endlist = info.taxastart.split(',')
                scorelist = [info.score]*len(startlist)
                dict_waafletaxa[abbrtaxa] = [scorelist, startlist, endlist]
        print( dict_waafletaxa )
 
        # Initiate the plots      
        taxarows = len( dict_waafletaxa.keys() )
        generows = len( dict_waaflegenes.keys() )
        ansrows = 0
        if not args.reforgs:
            rows = (taxarows + generows)*5/2 + 4 #make taxa and gene rows 2/5 of the plot, plus 4 for spaces in between plots
            plotrows = rows - taxarows - generows - 4
        else:
            rows = (taxarows + generows)*6/2 + 5 #make taxa and gene rows 2/6 of plot, plus 6 for spaces in between plots 
            remainrows = rows - taxarows - generows - 5
            plotrows = remainrows*3/4
            ansrows = remainrows*1/4
        fig = plt.figure(figsize=(12, 12), dpi=300)
        
        # Plot for orgscorer hits
        h_waafle = plt.subplot2grid( (rows,1), (0,0), rowspan=plotrows*2/3, colspan=1 ) #first 2/3rds of the plot
        h_waafle.set_ylim( -0.1, 1.1 )
        h_waafle.set_xlim( 0, contiglen )
        h_waafle.set_ylabel( 'scores' )
        #h_waafle.get_xaxis().set_ticks( np.arange(0, contiglen, contiglen/10) )

        # Plot for genes
        h_genes = plt.subplot2grid( (rows,1), (1+plotrows*2/3,0), rowspan=plotrows*1/3, colspan=1 )
        h_genes.set_ylim( 0, 1 )
        h_genes.set_xlim( 0, contiglen )
        #h_genes.get_xaxis().set_ticks(np.arange( 0, contiglen, contiglen/10) )
        h_genes.get_yaxis().set_ticks([0, 0.5, 1])
        h_genes.set_ylabel( 'genes' )

        if args.reforgs:
            next_coord = 3+plotrows + ansrows 
        else:
            next_coord = 2+plotrows

        # Plot for taxa legend
        h_taxaleg = plt.subplot2grid( (rows,1), (next_coord,0),rowspan=taxarows, colspan=1 ) #for taxa annotations
        h_taxaleg.set_xlim( 0, contiglen )
        h_taxaleg.set_ylim( 0, len(dict_taxacolor.keys()) )

        # Plot for Uniref annotations
        h_geneleg = plt.subplot2grid( (rows,1), (1 + next_coord + taxarows, 0), rowspan=generows, colspan=1 ) #for Uniref annotations
        h_geneleg.set_xlim( 0, contiglen )
        h_geneleg.set_ylim( 0, len(dict_waaflegenes.keys()) )

        # Plot taxa and taxa legend
        xlab = 0
        ylab = len(dict_taxacolor)
        ylab -= 2
        for taxa in dict_waafletaxa.keys():
            #Plot all taxa
            scores, starts, ends = dict_waafletaxa[taxa]
            diffs = list( np.array( [int(x) for x in ends] ) - np.array( [int(y) for y in starts] ) )
            adjdiffs = []
            for i in range( len( diffs ) ):
                value = diffs[i]
                arrowhead = abs( 0.05 * value )
                if float(value) < 0:
                    newvalue = value + arrowhead
                else:
                    newvalue = value - arrowhead
                if taxa in coloredorgs:
                    h_waafle.arrow( float( starts[i] ), float( scores[i] ), newvalue, 0, width=0.02, color=dict_taxacolor[taxa], head_width=0.02, head_length=arrowhead, alpha=0.8, zorder=2 )
                else:
                    h_waafle.arrow( float( starts[i] ), float( scores[i] ), newvalue, 0, width=0.02, color=dict_taxacolor[taxa], head_width=0.02, head_length=arrowhead, alpha=0.2 )
            #Plot taxa legend
            if taxa in coloredorgs:
                p = patches.Rectangle( (xlab, ylab), contiglen/20, 1, color=dict_taxacolor[taxa], alpha=0.8)
            else:
                p = patches.Rectangle( (xlab, ylab), contiglen/20, 1, color=dict_taxacolor[taxa], alpha=0.2)
            h_taxaleg.add_patch(p)
            h_taxaleg.text( xlab + contiglen/20, ylab + 0.5, str(taxa), fontsize=5 ) #contiglen/200 )
            xlab += contiglen/4
            if xlab >= contiglen:
                xlab = 0
                ylab -= 2
        h_taxaleg.set_axis_off()
       
        # Plot waafle genes
        score = 1
        ylab = len( dict_waaflegenes.keys() )
        for gene in dict_waaflegenes.keys():
            start, end, uniref50, uniref90 = dict_waaflegenes[ gene ]
            if score <= 0:
                score = 1
            else:
                score -= 0.1
            diff = float( end ) - float( start )
            arrowhead = abs( 0.05*diff )
            if float(diff) < 0:
                newdiff = diff + arrowhead
            else:
                newdiff = diff - arrowhead
            h_genes.arrow( dict_waaflegenes[gene][0], score, newdiff, 0, width=0.04, color="blue", head_width=0.04, head_length=arrowhead, alpha=0.9 )
            #Plot unirefs
            h_geneleg.text(0, ylab, 'Gene' + str(gene) + ' Uniref50:' + uniref50, fontsize=5 )
            h_geneleg.text(contiglen/4, ylab, 'Gene' + str(gene) + ' Uniref90:' + uniref90, fontsize=5 )
            ylab -= 2
        h_geneleg.set_axis_off()
         
        # Plot answers: In this case answers are genes! 
        if args.reforgs and args.refgenes:
            h_ref = plt.subplot2grid( (rows,1), ( 2 + plotrows,0), rowspan=ansrows, colspan=1 ) #initiate plot
            h_ref.set_ylim( 0, 1.1 )
            h_ref.set_xlim( 0, contiglen )
            h_ref.get_yaxis().set_ticks([0, 0.5, 1])
            h_ref.set_xlabel( 'coordinates (bp)' )
            h_ref.set_ylabel( 'answers' )

            recipient, donor = dict_reforgs[contig]
            genelist = dict_refgenes[contig]
            score = 1
            for gene in genelist:
                genelen = 0
                if gene.strand == '+':
                    genelen = gene.end - gene.start + 1
                    start = gene.start
                else:
                    genelen = (gene.end - gene.start + 1)*(-1)
                    start = gene.end
                arrowhead = abs(genelen*0.05)
                if genelen < 0:
                    newdiff = genelen + arrowhead
                else:
                    newdiff = genelen - arrowhead
                ID = gene.attribute.split(';')[0]
                if re.search( 'donor', ID ):
                    mycolor = 'red'
                    h_ref.text( (gene.start+gene.end)/2, score, 'donor:' + donor )
                else:
                    mycolor = 'blue'
                    h_ref.text( (gene.start+gene.end)/2, score, 'recip:' + recipient )
                h_ref.arrow( start, score, newdiff, 0, width=0.04, color=mycolor, head_width=0.04, head_length=arrowhead, alpha=0.4)
                # Set y-axis
                if score <= 0:
                    score = 1
                else:
                    score -= 0.1		
            
        filename = contig + "." + args.type
        plt.savefig( filename )
        plt.close(fig)
        

if __name__ == "__main__":
    main()
