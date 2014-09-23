WAAFLE
===================== 
*Widget to Annotate Assemblies and Find LGT Events*

Last updated on Sept 23 2014.

Authors: Tiffany Hsu and Eric Franzosa
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
Currently, there are 2 sets of scripts for the actual pipeline, and the
validation pipeline.

The actual pipeline is located in the 'pipeline' folder. It consists of
several steps:
1. BLAST
1. 

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