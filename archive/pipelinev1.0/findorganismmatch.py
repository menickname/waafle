#!/usr/bin/python

'''
This script is part of the HGT pipeline.
After BLAST has been run on the contigs, this script should be run to annotate the results with the organism names. The organism names will be inserted into the 3rd column of the BLAST results. This is done by:

1) Creating a dictionary of {ref names: taxa names}.
2) Matching the ref names in BLAST to the dictionary in Step 1.
'''

#Import
import sys
import subprocess
import json
import argparse
import re

#Define arguments
parser = argparse.ArgumentParser()
parser.add_argument('blast_results', help = 'Location and file of blast results') 
args = parser.parse_args()

#Set up dictionary of taxa genus/species names and their gi numbers
#Set the below up as a json file
#Taxa = {}
#for entry in open('/n/huttenhower_lab_nobackup/data/hgt/repophlan_pangenomes/info_repophlan_gene_lib.txt'): #Create a dictionary to match our entries to.
	#nentry = entry.strip().split("\t")
	#Taxa.update({nentry[0]:nentry[1]})
with open('/n/home05/thsu/hgt/data/dictGI_genusspecies.json') as infile1:
        dictGI_taxa = json.load(infile1)
	

#Since the above list only contains species/genera names, we will add in the rest using the below
#full_taxonomy = open('/n/huttenhower_lab/data/repophlan_chocophlan_pangenomes/single_celled_orgs.tax.txt')
for hit in open(args.blast_results):
	nhit = hit.strip().split('\t')
	species = dictGI_taxa[nhit[1]] #This returns the genus/species names
	#nspecies = re.sub(r'.centroids.ffn', '', species)
	#fulltaxonomy = subprocess.check_output('grep --max-count=1 \"'+nspecies+'\" /n/huttenhower_lab/data/2014_repophlan/repophlan_microbes.txt', shell = True)
	#nfulltaxonomy = fulltaxonomy.strip().split('\t')[1]
	nhit.insert(2, species)
	print '\t'.join(nhit)
	


	
	
