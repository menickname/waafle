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

def getGeneList(GCF, GI, recipientstart, recipientend):
        microbes_genes = '/n/huttenhower_lab/data/2014_repophlan/repophlan_microbes/ffn/'
        geneloc = microbes_genes + GCF + '.ffn'
        genelist = list(SeqIO.parse(geneloc, 'fasta'))
	expectedgenes = []
	for i in range(len(genelist)):
		if GI in genelist[i].id:
			genecoord = str(genelist[i].id).split(':')[1].split('-')
			if re.search('c', genecoord[0]):
				start, end = int(genecoord[1]), int(str(genecoord[0])[1:])
			else:
				start, end = int(genecoord[0]), int(genecoord[1])
			#print start, end, recipientstart, recipientend
			if start >= recipientstart and start <= recipientend:
				if end <= recipientend:
					genelen = end - start
					#print start, end, recipientstart, recipientend, genelen, 'shorterlen'
				elif end > recipientend:
					genelen = recipientend - start
					#print start, end, recipientstart, recipientend, genelen, 'greaterlen'
				if genelen != 0:
					expectedgenes.append([GI, genelen, start, recipientend])
			if start > recipientend:
				#print start, end, recipientstart, recipientend, 'itneedstoend'
				break
	return(expectedgenes)

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
	section1 = recipientGENOME[recipientGI].seq[0:int(recipient_start)-1] #Python is 0-based, but sequence coordinates are not.
	section3 = recipientGENOME[recipientGI].seq[int(recipient_end): len(recipientGENOME[recipientGI].seq)]
	
	#Slice the contig by contig_status and insert the correct section
	if contig_status == 'hgt':
		section2 = donorGENE.seq
	else:
		donortaxa = recipienttaxa
		donorGCF = recipientGCF
		donorGENE = getRandomGene(recipientGCF)
                section2 = donorGENE.seq

	if whichend == 'end':
		center = int(recipient_start) - 1
	else:
		center = int(recipient_start) - 2 + len(section2)

	newseq = section1 + section2 + section3
	
	contigstart = center - contig_len/2
	contigend = center + contig_len/2
	if contigstart < 0:
		contigstart = 0
	if contigend > len(newseq):
		contigend = len(newseq)
	
	#Determine if this will pass the threshhold
	#Get the length of the donor gene within the contig
	if contig_status == 'hgt':
		if whichend == 'end':
			donorgenelen = contigend - center + 1
			#print contigstart, recipient_start, 'end'
			recipientgenes = getGeneList(recipientGCF, recipientGI, int(contigstart), int(recipient_start))
			print recipientgenes
		else:
			donorgenelen = center - contigstart + 1
			#print recipient_end, str(int(recipient_end) + contig_len/2), 'start'
			recipientgenes = getGeneList(recipientGCF, recipientGI, int(recipient_end), int(recipient_end) + contig_len/2)
			print recipientgenes
		#Determine whether status should be changed
		rgenelen_list = []
                for i in range(len(recipientgenes)):
                	recipientgenelen = int(recipientgenes[i][1])
                        if recipientgenelen > 100:
                        	rgenelen_list.append(i)
		if donorgenelen < 100 or len(rgenelen_list) == 0:
                     	contig_status = 'nohgt_changed'
		

	#Get the length of the recipient genes within the contig
	#print recipientGCF, getGeneList(recipientGCF, recipientGI, recipient_)


	#print recipientGCF, recipientGI, recipient_start, recipient_end
	#print donorGCF, donorGI, donor_start, donor_end, len(donorGENE.seq)
	#print contig_len, contig_status, whichend
	#print recipient_start, center, recipient_end
	#print contigstart, contigend
	
	#Create a new seq record to output a fasta file
	newseq_sliced = newseq[contigstart:contigend]
	new_SeqRec = SeqRecord(newseq_sliced)
	new_SeqRec.id = 'contig' + str(count)
	new_SeqRec.description = 'donor:' + donorGCF + '|' + donorGENE.id + '||recipient:' + recipientGCF + '|' + recipientGENE.id
	new_SeqRec.annotations['source'] = recipienttaxa.split('|')
	new_SeqRec.annotations['taxonomy'] = donortaxa.split('|')
	mynewrecords.append(new_SeqRec) 
	
	
	#Generate an answer key
	contig_info = [new_SeqRec.id, new_SeqRec.description, donortaxa, recipienttaxa, phyldistance, contig_status, whichend, contig_len, contigend - contigstart + 1, center, contigstart, contigend]
	contig_info_list.append(contig_info)
	print contig_info, '\n'

#file_table = open('contig_info.txt', 'w')
#for line in contig_info_list:
#        file_table.write('\t'.join(str(i) for i in line) + '\n')
#file_table.close()

#SeqIO.write(mynewrecords, "fakecontigs.fasta", "fasta")

