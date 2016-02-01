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

### WAAFLE Pipeline
The WAAFLE pipeline is located in the 'waafle_v1.0' folder. It consists of
several steps:

1. __waafle_search.py__
    * This script BLASTs contigs against the Repophlan database.
    * Usage:  
    ``$ repophlan_db= /n/huttenhower_lab_nobackup/data/hgt/blast/blast_db/repophlan_31122013_speciescentroids.db  ``
    ``$ blastpath= '/usr/local/bin/blastn'``  
    ``$ python waafle_search.py -q CONTIGFILE.fasta -d repophlan_db -b blastpath -e 1 -o waafle-blastout.tsv``  

2. __waafle_genecaller.py__
    * This script takes the BLAST output, calls genes, and outputs the gene calls into [gff3 format](http://www.sequenceontology.org/gff3.shtml).
    * The parameters you must pass include overlap of hits (lap_h), overlap of genes (lap_g), subject coverage of hits (scov_h), and subject coverage of genes (scov_g). We are still optimizing which default parameter values to set.
    * Usage:  
    ```$ python waafle_genecaller.py -i waafle-blastout.tsv -o waafle-genes.gff -lap_h 0.5 -lap_g 0.5 -l 200 -scov_h 0.25 -scov_g 0.25```  

3. __waafle_orgscorer.py__
    * This script takes any gff3 file (can be from external gene callers such as Prodigal) and the BLAST output and outputs a tab-delimited file containing scores for each taxa within each gene.
    * The parameters that must be passed include overlap of hits and genes (lap), as well as subject coverage of hits (scov_h). These are used to group BLAST hits into genes before hits are grouped by taxa. The last parameter is the phylogenetic level to score taxa at, these are specified by 'k', 'p', 'c', 'o', 'f', 'g', or 's', which stand for kingdom, phylum, class, order, family, genus, or species.
    * Usage:  
    ```$ python waafle_orgscorer.py -g waafle-genes.gff -b waafle-blastout.tsv -o waafle-scoredorgs_g.tsv -t g -scov_h 0.25 -lap 0.5```  

4. __waafle_lgtscorer.py__
    * This script takes the taxa scores, scores contigs for LGT, determines potential donor/recipient taxa, and outputs each contig and its predicted status. 
    * The parameters that must be passed include the score threshold for 1 taxa (s1), score threshold for 2 taxa (s2), and whether or not we want to call "unknown" taxa (u). "True" turns it on, "False" turns it off.
    * Usage:  
    ```$ python waafle_lgtscorer.py -i waafle-scoredorgs_g.tsv -o waafle-scoredcontigs_g.tsv -s1 0.8 -s2 0.8 -u False```  


The above have now been linked into one doit script. So far, this doit script runs only 1 sample. It requires you specify the contig file for each run, though this could potentially be bypassed by including a task that loops through all files, and sends each file as a variable to the next set of tasks.


### Validation Pipeline

The validation pipeline is located in the 'waafle_validation' folder. There are two main steps, generating synthetic contigs and evaluating them.  

#### Generating Synthetic Contigs
1. __mkcontiglist.py__
    * This script generates a list of donor-recipient pairs at 8 different phylogenetic levels (in which all the donors and recipients are random). 
    * This list includes GCF numbers, full taxonomies, and phylogenetic distance for each donor-recipient pair. 
    * Usage:
    ```$ python mkcontiglist.py > donorrecip_list.txt```  

2. __mksyncontig.py__
    * This script takes the list of donor-recipient pairs generated in the first script and generates synthetic contigs. It outputs three files: 1) a gff format (to show where the genes of the contig are), 2) a fasta file containing the sequences 3) an answer key for the level the contig has LGT.
    * The parameters needed include the recipient taxa's GCF number (recipient), the donor taxa's GCF number (donor), the recipient taxa's taxonomy name across all levels (reciptaxa), the donor taxa's taxonomy name across all levels (donortaxa), the phylogenetic level at which the two taxa differ (taxadiff), and the number of contigs to make (ngenes). The latter is not a typo, we were just stupid in naming.
    * This script is never run alone. Instead, we run a wrapper called runmkcontig.py. 
    * Usage:
    ```$ python mksyncontig.py --recipient GCF_RECIP --donor GCF_DONOR --reciptaxa TAXA_RECIP_PHYLEVEL --donortaxa TAXA_DONOR_PHYLEVEL --taxadiff g --ngenes 1 ```  

3. __run_mksyncontig.py__
    * This is a wrapper script which takes list from step 1 and runs `mksyncontig.py` multiple times to generate contigs.
    * After this script is run, all gff, fasta, and answer key files should be concatenated into one.  
    * This script has two parameters. The first positional argument is the list of donor and recipient GCFs generated from step 1. The second positional argument is for 'ngenes' in step 2, or the number of contigs to generate per donor/recipient pair.
    * This script will only run as many donor/recipient pairs as are in the list of donor and recipient GCFs. It may be helpful to run step 1 multiple times to concatenate a longer list.
    * Usage:  
    ```$ python run_mksyncontig.py donorrecip_list.txt 1```

### Evaluation Pipeline
1. __eval_genes3.py__
    * This is a script that evaluates whether WAAFLE-called genes are similar relative to NCBI. It does so by comparing GFF overlaps.
    * The parameters needed include the WAAFLE genes gff and the synthetic contigs gff.
    * The `-o_gff` flag is the reference gff. Genes in the `-w_gff` flag will be compared to genes in the `-o_gff` flag, allowing us to find true positives, false positives, and false negative genes.
    * In the below, the first line of code compares WAAFLE-called genes to the synthetic genes, allowing us to call false positives. The second line of code compares reference genes to WAAFLE-called genes, allowing us to call true positives and false negative genes.
    * Usage:  
    ```$ python eval_genes3.py -w_gff waafle-genes.gff -o_gff syntheticcontig-genes.gff > waaflegenes_eval.tsv```
    ```$ python eval_genes3.py -w_gff syntheticcontig-genes.gff -o_gff waafle-genes.gff > syntheticgenes_eval.tsv```

2. __eval_lgt2.py__
    * This is a script that evaluates whether WAAFLE correctly called LGT or not. For each contig, it reports its status of true positive (TP), true negative (TN), false positive (FP), and false negative (FN). It also reports whether the taxa called match the answers, as well as whether the donor/recipients were called correctly.
    * The parameters needed include the scored contigs from the last step in the WAAFLE pipeline and the answer key generated when creating synthetic contigs.
    * Usage:  
    ```$ python eval_lgt2.py -i waafle-scoredcontigs_g.tsv -ans answerkey.tsv > contigs_evaluated.tsv```

3. __eval_all2.py__
    * This is a script that counts the number TP, FP, and FN genes that WAAFLE called, as well as TP, TN, FP, and FN LGTs WAAFLE called at each phylogenetic level. 
    * The output is 2 lines, the first is the header, the second is the numbers.
    * This allows us to concatenate large number of files to create a table of these values if WAAFLE is tested under multiple conditions.
    * The parameters needed include the evaluated WAAFLE and reference genes from step 1, as well as the evaluated contigs from step 2.
    * Usage:  
    ```$ python eval_all2.py -w_genes waaflegenes_eval.tsv -ref_genes syntheticgenes_eval.tsv -contig_eval contigs_evaluated.tsv -condition WAAFLE_conditions_as_str > results.txt```