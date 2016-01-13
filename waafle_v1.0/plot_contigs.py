#!/usr/bin/env/python

from __future__ import print_function
import matplotlib.pyplot as plt
import argparse
import waafle_utils as wu
import numpy as np
from matplotlib.pyplot import cm 
import pandas as pd

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
        "-ref", "--reference",
        help="output from reference, if available"
        )
    parser.add_argument(
        "-contigs", "--contiglist",
        help="list of contigs to output plots for"
        )
    args = parser.parse_args()
    return args

# ---------------------------------------------------------------
# main
# ---------------------------------------------------------------

def main():
    args = get_args()
    dict_mytaxa = {}
    for contig, taxalist in wu.iter_contig_taxa( args.waafle ):
        dict_mytaxa[contig] = taxalist
    
    dict_reftaxa = {}
    if args.reference:
        for contig, taxalist in wu.iter_contig_taxa( args.reference ):
            dict_reftaxa[contig] = taxalist

    for astrline in open( args.contiglist ):
        contig = astrline.strip()
        taxalist = dict_mytaxa[contig]
        
        fig = plt.figure(figsize=(12, 6))
        hax = fig.add_subplot(122)
 
        # Make colors
        orglabels = []
        for taxa in taxalist:
            orglabels.append( taxa.taxa )
        orgset = set(orglabels)
        colors = cm.Accent(np.linspace(0, 1, len(orgset)))

        dict_color = {}
        for i in range( len( list(orgset) ) ):
            dict_color[ list(orgset)[i]] = colors[i]

        dict_taxa = {}
        for current_org in dict_color.keys():
            ycoord, xstart, xend = [], [], []
            for taxa in taxalist:
                if taxa.taxa == current_org:
                    ycoord.append( float( taxa.score ) )
                    xstart.append( int( taxa.taxastart ) )
                    xend.append( int( taxa.taxaend ) )
            dict_taxa[ current_org ] = [ycoord, xstart, xend]
        
        # Plot
        for taxa in dict_taxa.keys(): 
            score, start, end = dict_taxa[ taxa ][0], dict_taxa[taxa][1], dict_taxa[taxa][2]
            hax.hlines( score, start, end, lw=10, colors = dict_color[taxa], label = taxa  )
        hax.set_ylim( 0, 1 ) 
        hax.set_xlabel( 'coordinates (bp) ')
        hax.set_ylabel( 'scores' )
        hax.legend( bbox_to_anchor=(-1.05, 1), loc=2, borderaxespad=0. )           
        
        filename = contig + '.png'
        plt.savefig( filename )        
        

if __name__ == "__main__":
    main()
