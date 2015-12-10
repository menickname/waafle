#!/usr/bin/env python

"""
WAAFLE VALIDATION: mkcontiglist.py

Authors:
Tiffany Hsu
Eric Franzosa

This script generates a list taxa-pairs with phylogenic diversity.
Taxa pairs are drawn from /n/huttenhower_lab/data/2014_repophlan/repophlan_microbes.txt
Run with the "-h" flag for usage help.
"""

from __future__ import print_function # python 2.7+ required
import random, subprocess

# ---------------------------------------------------------------
# constants
# ---------------------------------------------------------------

phylevels = ['k', 'p', 'c', 'o', 'f', 'g', 's', 't']
microbeGCF = '/n/huttenhower_lab/data/2014_repophlan/repophlan_microbes.txt'

# ---------------------------------------------------------------
# functions
# ---------------------------------------------------------------

def find_uniqtaxa( microbeGCF ):
    dict_GCFtaxa = {}
    dict_taxaGCF = {}
    for line in open( microbeGCF ):
        GCF, microbelist = line.split('\t')[0], line.strip().split('\t')[1].split('|')
        dict_GCFtaxa[GCF] = microbelist
        for taxa in microbelist:
            taxalevel = taxa[0]
            dict_taxaGCF.setdefault( taxa, [] ).append( GCF )
    dict_k, dict_p, dict_c, dict_o, dict_f, dict_g, dict_s = {}, {}, {}, {}, {}, {}, {}
    dict_list = [ dict_k, dict_p, dict_c, dict_o, dict_f, dict_g, dict_s ]
    for line in open( microbeGCF ):
        GCF, microbelist = line.split('\t')[0], line.strip().split('\t')[1].split('|')
        for i in range( len( microbelist ) - 1 ):
            highertaxa = microbelist[i]
            lowertaxa = microbelist[i+1]
            dict_list[i].setdefault( highertaxa, set() ).add( lowertaxa )
    return dict_GCFtaxa, dict_taxaGCF, dict_list



# ---------------------------------------------------------------
# main
# ---------------------------------------------------------------

def main():
    dict_GCFtaxa, dict_taxaGCF, dict_list = find_uniqtaxa( microbeGCF )
    for i in range( len( phylevels ) ):
        if i == 0:
            dict_of_interest = dict_list[i]
            donor, recipient = random.sample( dict_of_interest.keys(), 2 )
            dGCF, rGCF = random.choice( dict_taxaGCF[ donor ] ), random.choice( dict_taxaGCF[ recipient ] )
        if i > 0:
            dict_of_interest = dict_list[i-1]
            level_len = 0
            while level_len < 2:
                highlevel = random.choice( dict_of_interest.keys() )
                level_len = len( dict_of_interest[highlevel] )
            donor, recipient = random.sample( dict_of_interest[ highlevel ], 2 )
            dGCF, rGCF = random.choice( dict_taxaGCF[ donor ] ), random.choice( dict_taxaGCF[ recipient ] )
        print( '\t'.join( [phylevels[i], dGCF, rGCF, '|'.join( dict_GCFtaxa[ dGCF ] ), '|'.join( dict_GCFtaxa[ rGCF ] )] ) )



if __name__ == "__main__":
    main()
