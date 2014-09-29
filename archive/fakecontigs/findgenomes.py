#!/usr/bin/python

'''
This script will:
1. Pick two random microbial species.
'''

#Import
import argparse
import random
import re
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from ete2 import Tree
import json

#Define arguments
parser = argparse.ArgumentParser()
parser.add_argument('numbercontigs', type =int, help = 'Number of fake contigs you want to create.')
args = parser.parse_args()

#Read in Newick Tree
tree = Tree("/n/home05/thsu/hgt/data/bs_tree.reroot.nwk.txt")

#Read in IMG to NCBI conversion
#dict_table = {}
#img_ncbi_table = open('/n/home05/thsu/hgt/data/taxontable2884_29-jul-2014.xls')
#for astrline in img_ncbi_table:
#	aastrline = astrline.strip().split('\t')
#	ncbi = aastrline[7]
#	img = aastrline[0]
#	dict_table.setdefault(ncbi, []).append(img)
#with open("IMG_taxon_table.json", "w") as outfile:
#	json.dump(dict_table, outfile)

#Aggregate contigs and results
mynewrecords = []
contig_info_list = []
for i in range(args.numbercontigs + 1):

	#Pick two random genomes
	microbes_list = '/n/huttenhower_lab/data/2014_repophlan/repophlan_microbes.txt'
	donor = ''
	recipient = ''
	while donor == recipient: #To avoid picking the same microbes (which should rarely happen)
		donor = random.choice(open(microbes_list).readlines())
		recipient = random.choice(open(microbes_list).readlines())
	donor_GCF = donor.strip().split('\t')[0]
	recipient_GCF = recipient.strip().split('\t')[0]
	donor_name = donor.strip().split('\t')[1]
	recipient_name = recipient.strip().split('\t')[1]
	
	#For both donor and recipient, pick the random genes to be swapped
	microbes_genes = '/n/huttenhower_lab/data/2014_repophlan/repophlan_microbes/ffn/'
	donor_loc = microbes_genes + donor_GCF + '.ffn'
	recipient_loc = microbes_genes + recipient_GCF + '.ffn'
	donor_genes = list(SeqIO.parse(donor_loc, "fasta"))
	recipient_genes = list(SeqIO.parse(recipient_loc, "fasta"))

	#Get the coordinates to swap 
	donorGENE = random.choice(donor_genes)
	recipientGENE = random.choice(recipient_genes)
	donor_gi = str(donorGENE.id).split(':')[0]
	donor_coord = str(donorGENE.id).split(':')[1].split('-')
	if re.search('c', donor_coord[0]):
		donor_start, donor_end = donor_coord[1], str(donor_coord[0])[1:]
	else:
		donor_start, donor_end = donor_coord[0], donor_coord[1]
	recipient_gi = str(recipientGENE.id).split(':')[0]
	recipient_coord = str(recipientGENE.id).split(':')[1].split('-')
	if re.search('c', recipient_coord[0]):
		recipient_start, recipient_end = recipient_coord[1], str(recipient_coord[0])[1:]
	else:
		recipient_start, recipient_end = recipient_coord[0], recipient_coord[1]

	#Go back to whole genomes to get chunk of "donor"
	microbes_genomes = '/n/huttenhower_lab/data/2014_repophlan/repophlan_microbes/fna/'
	recipiome_loc = microbes_genomes + recipient_GCF + '.fna'
	recipient_genome = SeqIO.index(recipiome_loc, "fasta")
	donor_gene_len = int(donor_end) - int(donor_start) + 1
	#print donor_gene_len, donor_start, donor_end, len(donor_genome[donor_gi].seq)
	rcontig_len = random.randint(donor_gene_len, 200000) #The 200,000 bp limit for upper contig is arbitrary. Here, we are basing the recipient contig size off the size of the donor gene, but choosing the start/end sites through the start/end sites of the recipient gene.
	fromstart = int(recipient_start) - rcontig_len + donor_gene_len
	if fromstart > 0:
		rstart = random.randint(fromstart, int(recipient_start))
		rend = rcontig_len - donor_gene_len - (int(recipient_start) - rstart) + int(recipient_end)
	else:
		rstart = random.randint(0, int(recipient_start))
		rend = rcontig_len - donor_gene_len - (int(recipient_start) - rstart) + int(recipient_end)
	if rend > len(recipient_genome[recipient_gi].seq):
		rend = len(recipient_genome[recipient_gi].seq)
		rcontig_len = rend - rstart + 1 
	#print rstart, recipient_start, recipient_end, rend, rcontig_len, donor_gene_len

	#Insert the donor into the recipient for a new contig
	donor_gene = donorGENE.seq
	recipient_cont1 = recipient_genome[recipient_gi].seq[rstart: int(recipient_start)]
	recipient_cont2 = recipient_genome[recipient_gi].seq[int(recipient_end): rend]
	newseq = recipient_cont1 + donor_gene + recipient_cont2
	new_SeqRec = SeqRecord(newseq)
	new_SeqRec.id = 'contig' + str(i)
	new_SeqRec.description = 'donor:' + donor_GCF + '|recipient:' + recipient_GCF
	new_SeqRec.annotations['source'] = recipient_name.split('|')
	new_SeqRec.annotations['taxonomy'] = donor_name.split('|') 
	mynewrecords.append(new_SeqRec)

	#Calculate parameter X, which is the distance between the 2 genes
	SeqIO.write(donorGENE, "donorgene.faa", "fasta")
	SeqIO.write(recipientGENE, "recipientgene.faa", "fasta")
	#from Bio.Blast.Applications import NcbitblastxCommandline
	#cline =  NcbitblastxCommandline(query='donorgene.faa', subject='recipientgene.faa', outfmt=6, out="blastresults.txt") 
	from Bio.Blast.Applications import NcbiblastnCommandline
	cline =  NcbiblastnCommandline(query='donorgene.faa', subject='recipientgene.faa', outfmt=6, out="blastresults.txt")
	stdout, stderr = cline()
	#print cline
	blast_output = open("blastresults.txt")
	blastout_list = blast_output.readline().strip().split('\t')
	if len(blastout_list) > 1:
		pident, length, evalue, bitscore = blastout_list[2], blastout_list[3], blastout_list[10], blastout_list[11]
	else:
		pident, length, evalue, bitscore = 'NA', 'NA', 'NA', 'NA'
	#Most of the lengths suck and are very short, not sure we can measure this...


	#Calculate parameter Y, which is the phylogenetic distance between the 2 genomes
	recipient_taxonomy = new_SeqRec.annotations['source']
	donor_taxonomy = new_SeqRec.annotations['taxonomy']
	for i in range(min(len(recipient_taxonomy), len(donor_taxonomy))):
		if recipient_taxonomy[i] == donor_taxonomy[i]:
			continue
		else:
			difference_let = re.search(r'^[A-z]', recipient_taxonomy[i]).group()
			#print i, difference_let
			break
	#print recipient_taxonomy
	#print donor_taxonomy
	

	#Calculate parameter Y in a second way
	with open('/n/home05/thsu/hgt/results/fakecontigs/IMG_taxon_table.json') as infile1:
        	dictTable = json.load(infile1)
	recipient_ncbi = recipient_gi.split('|')[1]
	donor_ncbi = donor_gi.split('|')[1]
	
	#from cogent.parse.ncbi_taxonomy import NcbiTaxonomyFromFiles
	import subprocess
	subprocess.call('rm taxid.txt', shell = True)
	subprocess.call('grep --max-count=1 \"'+recipient_ncbi+'\" /n/home05/thsu/hgt/data/gi_taxid_nucl.dmp | tee -a taxid.txt', shell=True)
	subprocess.call('grep --max-count=1 \"'+donor_ncbi+'\" /n/home05/thsu/hgt/data/gi_taxid_nucl.dmp | tee -a taxid.txt', shell = True)
	recipient_taxid = open('taxid.txt').readlines()[0].split('\t')[1]
	donor_taxid = open('taxid.txt').readlines()[1].split('\t')[1]
	recipient_IMG = dictTable[recipient_taxid.strip()]
	donor_IMG = dictTable[donor_taxid.strip()]
	nrecipient_IMG = 't' + str(recipient_IMG[0]) #donor_IMG
	ndonor_IMG = 't' + str(donor_IMG[0])
	print recipient_IMG, donor_IMG
	R = tree&nrecipient_IMG	
	D = tree&ndonor_IMG
	print R, D, tree.get_distance(R, D)
	
	#print recipient_GCF
	#if recipient_ncbi in dictTable:
	#	print recipient_ncbi, dictTable[recipient_ncbi]
	#else:
	#	print recipient_ncbi, "Nothing"

	'''
	#calculate parameter X in a 2nd way
	#from Bio import pairwise2
	#from Bio.pairwise2 import format_alignment
	#alignments = pairwise2.align.globalms(str(donorGENE.seq), str(recipientGENE.seq), 2, -3, -5, -2) #These are in accordance w/BLASTN defaults
	#for i in range(len(alignments)):
	#pairwise_topscore = alignments[0][2]


	#Create table of contig information
	contig_info = [new_SeqRec.id, new_SeqRec.description, donor_name, recipient_name, pident, length, evalue, bitscore, difference_let, len(new_SeqRec.seq)]
	contig_info_list.append(contig_info)


file_table = open('contig_info.txt', 'w')
for line in contig_info_list:
	file_table.write('\t'.join(line) + '\n')

file_table.close()
SeqIO.write(mynewrecords, "fakecontigs.fasta", "fasta")
'''
