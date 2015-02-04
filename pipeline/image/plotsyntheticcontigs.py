#!/usr/bin/python

"""
This script is a wrapper script that calls an R script. The R script then plots all the graphs.
"""

# import
import argparse
import subprocess

# arguments
parser = argparse.ArgumentParser()
parser.add_argument( '--contigs', help='List of contigs to generate plots for.' )
parser.add_argument( '--answerkey', help='Location of answerkey.' )
parser.add_argument( '--genetable', help='Location of genetable.' )
args = parser.parse_args()

# R script location
scriptloc = '/n/home05/thsu/bitbucket/hgt_project/pipeline/image/runfunction.R'

# loop through list of contigs
for astrline in open(args.contigs):
	contigname = astrline.strip()
	filename = '_'.join(contigname.split('|'))
	subprocess.call(["Rscript", scriptloc, contigname, args.answerkey, args.genetable, filename])
