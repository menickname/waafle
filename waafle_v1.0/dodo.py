from __future__ import print_function
import collections
import re
import sys
from doit.tools import config_changed

#-------------------------------------
#Set variables
#-------------------------------------

waafle_loc = '/n/home05/thsu/bitbucket/waafle/waafle_v1.0/'
CONTIGS = '/n/home05/thsu/bitbucket/waafle/waafle_validation/test-fasta.ffn'
LAP = str(0.5)
LENGTH = str(0)
SCOV_H = str(0)
SCOV_G = str(0)
TAXA = ['k', 'p', 'c', 'o', 'f', 'g', 's']
ONEBUG = str(0.8)
TWOBUG = str(0.8)

#-------------------------------------
#Tasks
#-------------------------------------

def task_blast_contigs( ):
    # BLAST the contigs against Repophlan
	repophlan = '/n/huttenhower_lab_nobackup/data/hgt/blast/blast_db/repophlan_31122013_speciescentroids_uniprot.db'
	script = 'waafle_search.py'
	scriptloc = waafle_loc + script
	blastpath = '/usr/local/bin/blastn'
	return {
        'actions': [ 'python ' + scriptloc + ' -q ' + CONTIGS + ' -d ' + repophlan + ' -b ' + blastpath + ' -e 1 ' + '-o ' + '%(targets)s' ],
        'targets': ['waafle-blastout.tsv'],
        }

def task_call_genes( ):
	# Take BLAST hits and call genes
	script = "waafle_genecaller.py"
	scriptloc = waafle_loc + script
	return {
		'actions': ['python ' + scriptloc + ' -i ' + 'waafle-blastout.tsv' + ' -o ' + '%(targets)s ' + ' -lap ' + LAP + ' -l ' + LENGTH + ' -scov_h ' + SCOV_H + ' -scov_g ' + SCOV_G ],
		'targets': ['waafle-genes.gff'],
		'file_dep': ['waafle-blastout.tsv']
		}

def task_score_orgs( ):
	# Take gff and call orgs	
	script = "waafle_orgscorer.py"
	scriptloc = waafle_loc + script
	for taxa in TAXA:
		yield {
			'name': 'scoreorg_iteration_%s' %taxa,
			'actions': ['python ' + scriptloc + ' -g ' + 'waafle-genes.gff' + ' -b ' + 'waafle-blastout.tsv' + ' -t ' + taxa + ' -o ' + '%(targets)s' + ' -cov ' + SCOV_H + ' -lap ' + LAP],
			'targets': ['waafle-scoredorgs_%s.tsv' %taxa],
			'file_dep': ['waafle-genes.gff', 'waafle-blastout.tsv']
			}

def task_find_lgt( ):
	# Call which contigs have LGT
	script = "waafle_lgtscorer.py"
	scriptloc = waafle_loc + script
	unknown = str(False)
	for taxa in TAXA:
		input_file = 'waafle-scoredorgs_%s.tsv' %taxa
		yield {
			'name': 'lgt_iteration_%s' %taxa,
			'actions': ['python ' + scriptloc + ' -i ' + input_file  + ' -o ' + '%(targets)s' + ' -s1 ' + ONEBUG + ' -s2 ' + TWOBUG + ' -u ' + unknown ],
			'targets': ['waafle-scoredcontigs_%s.tsv' %taxa],
			'file_dep': ['waafle-scoredorgs_%s.tsv' %taxa],
			}
