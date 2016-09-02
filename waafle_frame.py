#!/usr/bin/env python

import os
import sys
import re
import argparse

class WaafleFrame( ):

    """
    Track the per-nucleotide and per-gene scores
    for all taxa contributing to a single contig
    """

    def __init__( self, contig_hits, gene_coords ):

        # dict of 2d numpy arrays (one per taxon)
        self.nt_scores = {}
        
        # dict of 2d numpy arrays (one per taxon)
        self.gene_scores = {}

        for hit in contig_hits:
            # determine bug
            # determine gene coords
            # augment nucleotide-level values

        self.gene_scores = make_gene_scores( )           
        
    def tax_collapse( self, taxlevel ):

        # reduce nt_scores to a specified level of the taxonony
        # rebuild the gene_scores structure

    def make_gene_scores( self, include_unknown=False ):
        # collapse the per-nucleotide scores to per-gene scores
