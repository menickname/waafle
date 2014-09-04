#!/usr/bin/python

'''
This script will pull out sections within whole genomes to create fake contigs. This should be done by:
1) Use what was generated from the shell script pull sequences for the donor organism.
2) Use 
'''


#Import
import argparse
import sys
import re
import collections


#Define arguments
parser = argparse.ArgumentParser()
parser.add_argument('GCFlist_recipient', help = 'Location and file of GCF list for recipient')
parser.add_argument('GCFlist_donor', help = 'Location and file of GCF list for donor')
parser.add_argument('fastalist_recipient', help = 'Location and file of FASTA to take')
parser.add_argument('fastalist_donor', help = 'Location and file of FASTA to take for donor')
args = parser.parse_args()


#Get the GCF list for the recipient into a list format
GCFlist_recipient = []
for GCF in open(args.GCFlist_recipient):
	rGCF = GCF.strip()
	rrGCF = '/n/huttenhower_lab/data/2014_repophlan/repophlan_microbes/fna/' + rGCF + '.fna'
	GCFlist_recipient.append(rrGCF)


#Get the GCF list for the donor into a list format
GCFlist_donor = []
for GCF in open(args.GCFlist_donor):
	dGCF = GCF.strip()
	ddGCF = '/n/huttenhower_lab/data/2014_repophlan/repophlan_microbes/fna/' + dGCF + '.fna'
	GCFlist_donor.append(ddGCF)

#Get the fasta headers from the shell script into a list format for the recipient
headerlistr = []
coordlistr = []
for fastaheader in open(args.fastalist_recipient):
	nfastaheader = fastaheader.strip().split(':')
	coordlistr.append(nfastaheader[1])
	headerlistr.append(nfastaheader[0])


#Get the fasta headers from the shell script into a list format for the donor
headerlistd = []
coordlistd = []
for fastaheader in open(args.fastalist_donor):
	nfastaheader = fastaheader.strip().split(':')
	coordlistd.append(nfastaheader[1])
	headerlistd.append(nfastaheader[0]) 

#print headerlistd
#print headerlistr
#print coordlistd
#print coordlistr
#print GCFlist_recipient
#print GCFlist_donor

#Get the whole genome sequences from the fasta headers
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


#Get the sequences for the recipient
seqListR = []
for i in GCFlist_recipient:
	for seq_record in SeqIO.parse(i, 'fasta'):
		for j in range(len(headerlistr)):
			#print seq_record.id, headerlistr[j]
			if seq_record.id == headerlistr[j]:
				coord = coordlistr[j].split('-')
				if re.search(r'c', coord[0]):
					endcoord = int(str(coord[0]).replace('c','')) + 1
					startcoord = int(str(coord[1]).replace('c','')) - 1
					complement = 'c'
				else:
					startcoord = int(coordlistr[j].split('-')[0]) - 1
					endcoord = int(coordlistr[j].split('-')[1]) + 1
					complement = ''
				#print seq_record.id # startcoord, endcoord
				seqlen = len(seq_record.seq[startcoord:endcoord])
				real_startcoord = startcoord - 50
				real_endcoord = real_startcoord + 500 #Currently, this is an arbitrary contig length of 600 bp (at 500 b/c we will insert 100 bp)
				newseq = seq_record.seq[real_startcoord: real_endcoord]
				newseq_r = SeqRecord(newseq)
				newseq_r.id = headerlistr[j] + complement + str(real_startcoord) + ':' + str(startcoord) + '-' +  str(endcoord) + ':' + str(real_endcoord)
				#print newseq_r.id
				#print newseq_r.seq	
				seqListR.append(newseq_r)

output_handle = open("recipient.fasta", "w")
SeqIO.write(seqListR, output_handle, "fasta")
output_handle.close()

#Get the sequences for the donor
seqListD = []
for i in GCFlist_donor:
	for seq_record in SeqIO.parse(i, 'fasta'):
		for j in range(len(headerlistd)):
			#print seq_record.id, headerlistd[j]
			if seq_record.id == headerlistd[j]:
				coordD = coordlistd[j].split('-')
				if re.search(r'c', coordD[0]):
					endcoordD = int(str(coordD[0]).replace('c','')) + 1
					startcoordD = int(str(coordD[1]).replace('c','')) - 1
					complement = 'c'
				else:
					startcoordD = int(coordlistd[j].split('-')[0]) - 1
					endcoordD = int(coordlistd[j].split('-')[1]) + 1
					complement = ''
				seqlenD = len(seq_record.seq[startcoordD:endcoordD])
				real_startcoordD = startcoordD
				real_endcoordD = real_startcoordD + 100
				newseqD = seq_record.seq[real_startcoordD:real_endcoordD]
				newseqD_r = SeqRecord(newseqD)
				newseqD_r.id = headerlistd[j] + complement + str(real_startcoordD) + ':' + str(startcoordD) + '-' + str(endcoordD) + ':' + str(real_endcoordD)
				#print newseqD_r.id
				#print newseqD_r.seq
				seqListD.append(newseqD_r)
output_handleD = open("donor.fasta", "w")
SeqIO.write(seqListD, output_handleD, "fasta")
output_handleD.close()

#Insert the donor into the recipient
newHGT = []
recipient_r = list(SeqIO.parse('recipient.fasta', 'fasta'))
donor_r = list(SeqIO.parse('donor.fasta', 'fasta'))
for i in range(len(recipient_r)):
	r_seq = recipient_r[i].seq
	d_seq = donor_r[i].seq
	new_seq = r_seq[:401] + d_seq + r_seq[401:]
	new_seq_record = SeqRecord(new_seq)
	new_seq_record.id = recipient_r[i].id + '===' + donor_r[i].id
	newHGT.append(new_seq_record)
	#print new_seq, len(new_seq) 
output_handleN = open("fakeHGT.fasta", "w")
SeqIO.write(newHGT, output_handleN, "fasta")
output_handleN.close()
