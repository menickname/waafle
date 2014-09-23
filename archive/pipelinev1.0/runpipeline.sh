#!/bin/bash
#This script should join all the pipeline scripts together to run the pipeline smoothly.


#Run BLASTN on the contigs you would like to search for LGT.
#These contigs should be in FASTA format.
module load bio/ncbi-blast-2.2.28+
INPUT=$1 #FASTA file for BLAST
OUTPUT='blastn.out'
DB='/n/huttenhower_lab_nobackup/data/hgt/blast/blast_db/repophlan_gene_clusters.db'
#blastn -db $DB -query $INPUT -out $OUTPUT -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp'


SCRIPT_LOC='/n/home05/thsu/hgt/bsub_and_scripts/pipeline'


#Step 1: Annotate the BLAST results with their organism names
python $SCRIPT_LOC/findorganismmatch.py $OUTPUT > blastn_annot.out


#Step 2: Sort the output by contig name, length, and bitscore.
sort -k1,1 -k5,5nr -k13,13n blastn_annot.out > blastn_annot.sorted


#Step 3: Group the BLASTN hits
OVERLAP=$2
python $SCRIPT_LOC/aggregatebylen.py blastn_annot.sorted $OVERLAP > grouped$OVERLAP.out
python $SCRIPT_LOC/aggregatebylen_comboverlaps.py grouped$OVERLAP.out $OVERLAP > grouped$OVERLAP.comb$OVERLAP.out 


#Step 4: Score the organism hits
TAXA_LEVEL=$3
python $SCRIPT_LOC/aggregatebylen_scoreorgs.py grouped$OVERLAP.comb$OVERLAP.out $TAXA_LEVEL


#Step 5: Find potential LGT/HGT hits
DELTA=$4
EPSILON=$5
LENCUT=$6

#Determine contig lengths
python /n/home05/thsu/hgt/bsub_and_scripts/tools/findcontiglen.py $INPUT 
python $SCRIPT_LOC/aggregatebylen_matchsets.py dddictOrgGroupScore_$TAXA_LEVEL.json dddictGroupOrgScore_$TAXA_LEVEL.json $DELTA $EPSILON $LENCUT > results.txt 
