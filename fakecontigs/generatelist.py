#!/usr/bin/python

#Import
import argparse
import random
import re
import subprocess
import collections

#Define arguments
parser = argparse.ArgumentParser()
args = parser.parse_args()

#Go into microbes list and pick a "seed" species.
#To randomly pick:
microbeslist = open('/n/huttenhower_lab/data/2014_repophlan/repophlan_microbes.txt')
donorline = subprocess.check_output(['shuf', '-n 1'], stdin=microbeslist)
donorGCF, donortaxa = donorline.strip().split('\t')[0], donorline.strip().split('\t')[1]
cnt = collections.Counter()

#Get differing levels of phylogenetic distance for this particular microbe
leveldiff = ['s', 'g', 'f', 'o', 'c', 'p']

#Do not stop the loop until we have at least 4 different levels
while len(cnt.keys()) <= 4:
	for level in leveldiff:
		searchstring = '.+' + level + '__'
		matchstring = re.search(searchstring, donortaxa).group() #Everything up until the level diff should match
		
		#Get diff
		recipientlines = subprocess.Popen(['grep', matchstring], stdin=open('/n/huttenhower_lab/data/2014_repophlan/repophlan_microbes.txt'), stdout=subprocess.PIPE)
		recipientline = subprocess.check_output(['shuf', '-n 1'], stdin=recipientlines.stdout)  #microbeslist)
		recipientGCF, recipienttaxa = recipientline.strip().split('\t')[0], recipientline.strip().split('\t')[1]
		donortaxonomy, recipienttaxonomy = donortaxa.split('|'), recipienttaxa.split('|')
		for j in range(min(len(donortaxonomy),len(recipienttaxonomy))):
			if donortaxonomy[j] == recipienttaxonomy[j]:
				continue
			else:
				if j <= 6:
					taxonomydiff = leveldiff[len(leveldiff)-j]
				else:
					taxonomydiff = 'Below species'
				break
		if recipientGCF != donorGCF and taxonomydiff != 'Below species':
			newlist = [donorGCF, donortaxa, taxonomydiff, recipientGCF, recipienttaxa]
			if cnt[taxonomydiff] <= 4:
				cnt[taxonomydiff] += 1
				newline = '\t'.join(newlist)
				print newline

