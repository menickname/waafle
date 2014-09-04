#!/usr/bin/python

'''
This script will:
1. Take in 2 species names.
2. 
'''

#Import
import argparse
import random
import re
import subprocess
import math
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

#Define arguments
parser = argparse.ArgumentParser()
parser.add_argument('fakelist', help = "Text file delineating donor and recipient species. This should be generated from 'generatelist.py'.")
args = parser.parse_args()

#Functions
def getRandomGene(GCF):
	microbes_genes = '/n/huttenhower_lab/data/2014_repophlan/repophlan_microbes/ffn/'
	geneloc = microbes_genes + GCF + '.ffn'
	genelist = list(SeqIO.parse(geneloc, 'fasta'))
	randomGENE = random.choice(genelist) 
	return(randomGENE)

def getGenome(GCF2):
	microbes_genomes = '/n/huttenhower_lab/data/2014_repophlan/repophlan_microbes/fna/'
	genomeloc = microbes_genomes + GCF2 + '.fna'
	genome = SeqIO.index(genomeloc, 'fasta')
	return(genome)

#Determine number of contigs that will be: 1) HGT (P) 2) Truncated HGT 3) No HGT (N)
count = 0
mynewrecords = []
contig_info_list = []

for astrline in open(args.fakelist):
	count += 1
	recipientGCF = astrline.strip().split('\t')[0]
	recipienttaxa = astrline.strip().split('\t')[1]
	phyldistance = astrline.strip().split('\t')[2]
	donorGCF = astrline.strip().split('\t')[3]
	donortaxa = astrline.strip().split('\t')[4]
	recipientGENE, donorGENE = getRandomGene(recipientGCF), getRandomGene(donorGCF)
	recipientGI, donorGI = str(recipientGENE.id).split(':')[0], str(donorGENE.id).split(':')[0]
	
	#Determine gene coordinates
	recipient_coord, donor_coord = str(recipientGENE.id).split(':')[1].split('-'), str(donorGENE.id).split(':')[1].split('-')
	if re.search('c', recipient_coord[0]):
		recipient_start, recipient_end = recipient_coord[1], str(recipient_coord[0])[1:]
	else:
		recipient_start, recipient_end = recipient_coord[0], recipient_coord[1]
	if re.search('c', donor_coord[0]):
		donor_start, donor_end = donor_coord[1], str(donor_coord[0])[1:]
        else:
                donor_start, donor_end = donor_coord[0], donor_coord[1]		
	
	#Determine contig size
	contig_len = max(int(random.expovariate(1/float(905))), 300) #905 was the average for all contigs in the tongue dorsum set
	contig_status = random.choice(['hgt', 'nohgt'])
	whichend = random.choice(['start', 'end'])
	#print recipientGCF, recipientGI, recipient_start, recipient_end
        #print donorGCF, donorGI, donor_start, donor_end
	
	#Get recipient genome and contig sections
	recipientGENOME = getGenome(recipientGCF)
	section1 = recipientGENOME[recipientGI].seq[0:int(recipient_start)]
	section3 = recipientGENOME[recipientGI].seq[int(recipient_end): len(recipientGENOME[recipientGI].seq)]
	
	#Slice the contig by contig_status and insert the correct section
	if contig_status == 'hgt':
		section2 = donorGENE.seq
	else:
		newrecipientGENE = getRandomGene(recipientGCF)
                section2 = newrecipientGENE.seq

	if whichend == 'end':
		donorstart = int(recipient_start)
	else:
		donorstart = int(recipient_end)
	
	newseq = section1 + section2 + section3
	
	contigstart = donorstart - contig_len/2
	contigend = donorstart + contig_len/2
	if contigstart < 0:
		contigstart = 0
	if contigend > len(recipientGENOME[recipientGI].seq):
		contigend = len(recipientGENOME[recipientGI].seq) - 1
	
	#print recipientGCF, recipientGI, recipient_start, recipient_end
	#print donorGCF, donorGI, donor_start, donor_end, len(donorGENE.seq)
	#print contig_len, contig_status, whichend
	#print recipient_start, donorstart, recipient_end
	#print contigstart, contigend
	
	#Create a new seq record to output a fasta file
	newseq_sliced = newseq[contigstart:contigend]
	new_SeqRec = SeqRecord(newseq_sliced)
	new_SeqRec.id = 'contig' + str(count)
	new_SeqRec.description = 'donor:' + donorGCF + '|recipient:' + recipientGCF
	new_SeqRec.annotations['source'] = recipienttaxa.split('|')
	new_SeqRec.annotations['taxonomy'] = donortaxa.split('|')
	mynewrecords.append(new_SeqRec) 
	
	#Generate an answer key
	contig_info = [new_SeqRec.id, new_SeqRec.description, donortaxa, recipienttaxa, phyldistance, contig_status, whichend, contig_len]
	contig_info_list.append(contig_info) 
	#print contig_info

file_table = open('contig_info.txt', 'w')
for line in contig_info_list:
        file_table.write('\t'.join(str(i) for i in line) + '\n')
file_table.close()

SeqIO.write(mynewrecords, "fakecontigs.fasta", "fasta")

