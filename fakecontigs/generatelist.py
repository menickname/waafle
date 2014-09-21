#!/usr/bin/python

'''
Update: 9/20/14
Need a way to terminate this if it takes too long trying to grep the lines.

This script generates a list of donor and recipient GCFs and taxonomies at 5 different levels. 
'''

# Import
import argparse
import random
import re
import subprocess
import collections

# Define arguments
parser = argparse.ArgumentParser()
args = parser.parse_args()

# Go into microbes list and pick a "recipient" species.
microbeslist = open('/n/huttenhower_lab/data/2014_repophlan/repophlan_microbes.txt')
recipline = subprocess.check_output(['shuf', '-n 1'], stdin=microbeslist)
recipGCF, reciptaxa = recipline.strip().split('\t')[0], recipline.strip().split('\t')[1]
cnt = collections.Counter()

# Get differing levels of phylogenetic distance for this particular microbe
leveldiff = ['t', 's', 'g', 'f', 'o', 'c', 'p']
donorGCFlist = []
count = 0

#Do not stop the loop until we have at least 5 different levels of phylogenetic distance.
while len(cnt.keys()) < 5:
	for level in leveldiff:
		count += 1
		searchstring = '.+' + level + '__'
		if re.search(searchstring, reciptaxa) != None:
			matchstring = re.search(searchstring, reciptaxa).group() # This increases the chance of finding something at multiple levels
		else:
			break
	
		# Get the donor gene.
		donorlines = subprocess.Popen(['grep', matchstring], stdin=open('/n/huttenhower_lab/data/2014_repophlan/repophlan_microbes.txt'), stdout=subprocess.PIPE)
		donorline = subprocess.check_output(['shuf', '-n 1'], stdin=donorlines.stdout) 
		donorGCF, donortaxa = donorline.strip().split('\t')[0], donorline.strip().split('\t')[1]
		reciptaxalist, donortaxalist = reciptaxa.split('|'), donortaxa.split('|')
		
		# Get the taxonomic difference between the donor and recip species.
		if donorGCF == recipGCF:
			taxadiff = 't'
		else:
			for j in range(min(len(reciptaxalist),len(donortaxalist))):
				if reciptaxalist[j] == donortaxalist[j]:
					continue
				else:
					if j <= 6:
						taxadiff = leveldiff[len(leveldiff)-j]
					else:
						taxadiff = 't'
					break
		
		# Output the importnat information
		newlist = [taxadiff, donorGCF, recipGCF, donortaxa, reciptaxa]
		if cnt[taxadiff] < 1 and donorGCF not in donorGCFlist:
			cnt[taxadiff] += 1
			donorGCFlist.append(donorGCF)
			newline = '\t'.join(newlist)
			print newline
