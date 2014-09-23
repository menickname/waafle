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

###Current Pipeline
The actual pipeline is located in the 'pipeline' folder. It consists of
several steps:

1. Run BLAST against the contigs.
```
$ DB=/n/huttenhower_lab_nobackup/data/hgt/blast/blast_db_updated/repophlan_31122013_speciescentroids.db
$ INPUT=contigs.fasta
$ OUTPUT=blastn.out 
$ blastn -db $DB -query $INPUT -out $OUTPUT -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp'
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
$ python aggregatebylen_matchsets.py dddictOrgGroupScore.json dddictGroupOrgScore.json 0.75 0.25 100
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