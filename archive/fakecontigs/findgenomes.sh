#!/bin/bash
#This script should 1)Get the the gi and accession numbers for the genes we would like to use for combining contigs.

SPECIES=$1 #This allows us to specify some species name, such as: 'Prevotella_oris'
MICROBES_REF='/n/huttenhower_lab/data/2014_repophlan/repophlan_microbes.txt'
GCF_PATH='/n/huttenhower_lab/data/2014_repophlan/repophlan_microbes/ffn/'
CENTROID_PATH='/n/huttenhower_lab/data/repophlan_chocophlan_pangenomes/repophlan_31122013_speciescentroids_ffn/g__Prevotella.s__'$SPECIES'.centroids.ffn'
COUNT=0
#ALLFASTA=''
ALLFASTA_HEADS=''

#Pick two random genomes in Repophlan.
shuf -n2 $MICROBES_REF > myspecies

#Concatenate all GCF files for a list of all genes that belong to this species.
cat $MICROBES_REF | grep $SPECIES | cut -f1 > $SPECIES'gcf_list'
while read GCF; do
	COUNT=$((COUNT+1))
	#cat $GCF_PATH$GCF'.ffn' > $SPECIES'_GCFALL'$COUNT	
	cat $GCF_PATH$GCF'.ffn' | grep '>' | sed 's/ [[:graph:]]*//g' > $SPECIES'_GCF'$COUNT
	#ALLFASTA=$ALLFASTA' '$SPECIES'_GCFALL'$COUNT
	ALLFASTA_HEADS=$ALLFASTA_HEADS' '$SPECIES'_GCF'$COUNT 
done < $SPECIES'gcf_list'
#cat $ALLFASTA > $SPECIES'_GCF_ALL.fasta'
cat $ALLFASTA_HEADS > $SPECIES'_GCF_ALL'
sort $SPECIES'_GCF_ALL' > $SPECIES'_GCF_ALL.sorted'
ALLGENES=$SPECIES'_GCF_ALL.sorted'
rm $ALLFASTA_HEADS $SPECIES'_GCF_ALL' $ALLFASTA 


#Get a list of all centroids that belong to this species.
cat $CENTROID_PATH | grep '>' | sort > $SPECIES'_centroids.sorted' 
CENTROIDS=$SPECIES'_centroids.sorted'


#Determine which centroids are shared between all genes vs. the centroids
comm -12 $ALLGENES $CENTROIDS > $SPECIES'.common' #Are common between centroids and genomes. Use for best BLAST results.
comm -13 $ALLGENES $CENTROIDS > $SPECIES'.centroidsonly' #Are unique to the centroids. This should be investigated, since it should == 0.
comm -23 $ALLGENES $CENTROIDS > $SPECIES'.genesonly' #Are unique to the genes of all species w/$SPECIES name in Repophlan. Use for testing on BLAST.
COMMON=$SPECIES'.common'
CENTROIDSONLY=$SPECIES'.centroidsonly'
GENESONLY=$SPECIES'.genesonly'


#List out the number of genes that come from the genome, the centroids, and their union/intersections.
wc -l $ALLGENES
wc -l $CENTROIDS
wc -l $COMMON
wc -l $GENESONLY
wc -l $CENTROIDSONLY
rm $ALLGENES $CENTROIDS


#Determine which genes to use for pasting genes together.
awk 'rand() < 0.01' $COMMON | sed 's/>//g' | head -n20 > $SPECIES'.randcommon'
awk 'rand() < 0.01' $GENESONLY | sed's/>//g' | head -n20 > $SPECIES'.randgenesonly'
#head -n5  $SPECIES'.randcommon'

