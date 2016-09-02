from __future__ import print_function
import collections
import re
import os
import errno
from doit.tools import config_changed
from doit import get_var

#-------------------------------------
#Comments
#-------------------------------------
"""
08/19/16
This will be used to run all the real samples.
"""

#-------------------------------------
#Set variables
#-------------------------------------

waafleloc = '/n/home05/thsu/bitbucket/waafle/waafle_v1.0/'
workingdir = os.getcwd() + '/'
#blastdir = '/n/hutlab11_nobackup/data/waafle/blastn_output_completed'
resultdir = '/n/regal/huttenhower_lab/thsu/waafle/results/throat/'
taxalist = ['k', 'p', 'c', 'o', 'f', 'g', 's']

SCOVH = str(0.75) 
LAP = str(0.1)
GENELEN = str(200) 
ONEBUG = str(0.5)
TWOBUG = str(0.8)

blastfile = get_var( 'blastfile' )
name = re.search( 'SRS[0-9]+', blastfile ).group()


#-------------------------------------
#Functions
#-------------------------------------

def make_sure_path_exists(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

#-------------------------------------
#Tasks
#-------------------------------------

def task_call_genes( ):
    # Take BLAST hits and call genes
    script = "waafle_genecaller.py"
    scriptloc = waafleloc + script
    blastloc = blastfile
    resultloc = resultdir + name
    make_sure_path_exists( resultloc )
    yield {
        'name': name,
        'actions': ['python ' + scriptloc + ' -i ' + blastloc + ' -o ' + '%(targets)s' + ' -lap ' + LAP + ' -l ' + GENELEN + ' -scov ' + SCOVH ],
        'file_dep': [blastloc],
        'targets': [resultloc + '/waafle-genes.gff']
		}   

def task_score_orgs( ):
    # Take gff and call orgs
    script = "waafle_orgscorer.py"
    scriptloc = waafleloc + script
    blastloc = blastfile
    resultloc = resultdir + name
    for taxa in taxalist:
        yield {
            'name': name + '_' + taxa,
            'actions': ['python ' + scriptloc + ' -g ' + resultloc + '/waafle-genes.gff' + ' -b ' + blastloc + ' -t ' + taxa + ' -o ' + '%(targets)s' + ' -lap ' + LAP + ' -scov ' + '0' ],
            'targets': [resultloc + '/waafle-scoredorgs_%s' %taxa + '.tsv'],
            'file_dep': [resultloc + '/waafle-genes.gff', blastloc],
            }

def task_find_lgt( ):
	# Call which contigs have LGT
    script = "waafle_lgtscorer.py"
    scriptloc = waafleloc + script
    resultloc = resultdir + name
    for taxa in taxalist:
        orgfile = resultloc + '/waafle-scoredorgs_%s' %taxa + '.tsv'
        yield {
            'name': name + '_' + taxa,
            'actions': ['python ' + scriptloc + ' -i ' + orgfile + ' -o ' + '%(targets)s' + ' -s1 ' + ONEBUG + ' -s2 ' + TWOBUG],
            'targets': [resultloc + '/waafle-scoredcontigs_%s' %taxa + '.tsv'],
            'file_dep': [resultloc + '/waafle-scoredorgs_%s' %taxa + '.tsv']
            }

def task_concatenate( ):
    # Concatenate all results
    filelist = []
    resultloc = resultdir + name
    for taxa in taxalist:
        files = resultloc + '/waafle-scoredcontigs_%s' %taxa + '.tsv'
        filelist.append( files )
    concatfile = resultloc + '/concatfile.txt'
    yield {
        'name': name,
        'actions': ['cat ' + ' '.join( filelist ) + ' > ' + concatfile],
        'targets': [concatfile],
        'file_dep': filelist
        }

def task_aggregate( ):
    # Aggregate contigs
    script = "waafle_aggregator.py"
    scriptloc = waafleloc + script
    resultloc = resultdir + name
    concatfile = resultloc + '/concatfile.txt'
    yield {
        'name': name,
        'actions': ['python ' + scriptloc + ' -i ' + concatfile + ' >%(targets)s '],
        'targets': [resultloc + '/' + name + '-scoredcontigs.tsv'],
        'file_dep': [concatfile]
        }
