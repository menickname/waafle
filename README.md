WAAFLE
====================
*Workflow to Annotate Assemblies and Find LGT Events*

Last updated on January 13, 2016.

Authors: Tiffany Hsu and Eric A. Franzosa

Huttenhower Lab, Harvard School of Public Health,
Boston, MA

-----------------------------------------------------
Overview
=======
WAAFLE is a tool that will take metagenomic shotgun reads, assemble them into
contigs, and then characterize contigs with potential lateral gene transfer events. 

Installation
=======
TBD.

Usage
=======
Currently, there are 4 scripts for the pipeline, and _ scripts
for the validation pipeline.

###WAAFLE Pipeline
The WAAFLE pipeline is located in the 'waafle_v1.0' folder. It consists of
several steps:

1. waafle_search.py
```
$ DB=/n/huttenhower_lab_nobackup/data/hgt/blast/blast_db_updated/repophlan_31122013_speciescentroids.db
```

2. Sort the BLAST hits by name, length (in decreasing order), and bitscore (in decreasing order).
```
$ sort -k1,1 -k4,4nr -k12,12nr blastn.out > blastn.sorted
```

3. Identify gene groups from the BLAST results by grouping BLAST hits that overlap, merging groups that overlap, and then removing groups < some length.
```
$ python blast2groups2.py --blastoutput blastn.sorted --hitoverlap 0.5 --groupoverlap 0.5 --length 100 > groupedcombhits.out
```

4. Score each taxa within each group by coverage_over_group x percID across all BLAST hits that cover that taxon at the specified phylogenetic level (species, genus, class, etc). Output those that are A) Explained by 1 organism only B) Explained by at least 1 organism across all genes on the contig at high confidence.
```
$ python scoreorgs.py --fasta contigs.fasta --blastoutput groupedcombhits.out --taxa g --delta 0.75
```

5. Detect which organisms are likely to have LGT by picking contigs with >=2 organisms with high to low drops in scores. This is currently being revised.
```
$ python detectlgt.py --dict1 dddictCOGS_g.json --dict2 dddictCGOS_g.json --delta 0.75 --epsilon 0.25
```

The above have now been linked into one python script runpipeline.py.
```
usage: runpipeline.py [-h] [--delta DELTA] [--epsilon EPSILON]
                      [--length LENGTH] [--taxa TAXA]
                      [--hitoverlap HITOVERLAP] [--groupoverlap GROUPOVERLAP]

optional arguments:
  -h, --help            show this help message and exit
  --delta DELTA         Upper threshold for calling high confidence BLAST hits
  --epsilon EPSILON     Lower threshold for calling low confidence BLAST hits
  --length LENGTH       Length of groups to exclude.
  --taxa TAXA           Taxa level to detect LGT at.
  --hitoverlap HITOVERLAP
                        Amount hits should overlap to join a group.
  --groupoverlap GROUPOVERLAP
                        Amount groups should overlap to merge.
```

###Validation Pipeline

The validation pipeline is located in the 'fakecontigs' folder. It consists of several steps:

1. Generate a list of donor-recipient pairs at 5 different phylogenetic
levels (in which all the recipients are identical). This list includes GCF
numbers, full taxonomies, and phylogenetic distance for each donor-recipient
pair. 
```
$ python generatelist.py > donorreciplist.txt
```

2. Make n number of contigs corresponding to a donor-recipient taxa pair.
```
$ python fakemake.py --recipient 'GCF_000415265' --donor 'GCF_000347175' --reciptaxa 'k__Bacteria|p__Firmicutes|c__Bacilli|o__Lactobacillales|f__Enterococcaceae|g__Enterococcus|s__Enterococcus_faecium|t__Enterococcus_faecium_SD3B_2' --donortaxa 'k__Bacteria|p__Spirochaetes|c__Spirochaetia|o__Spirochaetales|f__Leptospiraceae|g__Leptospira|s__Leptospira_alstoni|t__Leptospira_alstoni_serovar_Sichuan_str_79601' --ngenes 5
```

3. Wrapper script which takes list from Step 1, inputs donor-recipient pairs
into fakemake.py, and outputs a fasta file and
answer key.  
```
$ python run.py donorreciplist.txt
```
