from __future__ import print_function
import collections
import re
import sys

#-------------------------------------
#Set variables
#-------------------------------------

num_unique_taxa = 63
num_unique_contigs = 2

#-------------------------------------
#Tasks
#-------------------------------------

def task_generate_random_taxa( ):
    # Call script to generate a tab-separated file of GCFs at different phylogenetic levels.
    generate_random_taxa = '/n/home05/thsu/bitbucket/waafle/waafle_validation/mkcontiglist.py'
    for i in range( num_unique_taxa ):
        yield {
            'name': 'iteration%d' %i,
            'actions': ['python ' + generate_random_taxa + '> %(targets)s'],
            'targets': ['file_%s.tsv' %i]
            }

    
def task_concatenate_taxa( ):
    # Concatenate similar files
    files_sim = 'file_*'
    return {
        'actions': ['cat ' + files_sim +  ' > %(targets)s'],
        'targets': ['concatenated_taxa.tsv'],
        'file_dep': ['file_0.tsv']
        }


def task_generate_contigs( ):
    # Generate the contigs
    generate_contigs = '/n/home05/thsu/bitbucket/waafle/waafle_validation/run.py'
    return {
        'actions': [ 'python ' + generate_contigs + ' concatenated_taxa.tsv ' + num_unique_contigs ],
        'file_dep': ['concatenated_taxa.tsv'],
        }

def task_concatenate_contigs( ):
	# Concatenate similar files
	files_fnt = '*.fnt'
	return {
		'actions': ['cat ' + files_fnt + '> %(targets)s'], 	
		'targets': ['concatenated_fasta.ffn']
		}	

def task_concatenate_answer( ):
	# Concatenate similar files
    files_ans = '*-ans*'
    return {
        'actions': ['cat ' + files_ans + '> %(targets)s'],
        'targets': ['concatenated_answers.tsv']
        }

def task_rm_contigs( ):
	return {
		'actions': ['rm ' + 'GCF*'],
		'file_dep':['concatenated_fasta.ffn', 'concatenated_taxa.tsv']
	}

def task_rm_files( ):
	return {
		'actions': ['rm ' + 'file_*'],
		'file_dep':['concatenated_fasta.ffn', 'concatenated_taxa.tsv']
	}
