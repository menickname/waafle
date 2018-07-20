# WAAFLE

**WAAFLE** (the **W**orkflow to **A**nnotate **A**ssemblies and **F**ind **L**GT Events) is a method for finding novel LGT (Lateral Gene Transfer) events in assembled metagenomes, including those from human microbiomes. WAAFLE was developed in the [Huttenhower Lab](http://huttenhower.sph.harvard.edu) at the Harvard T.H. Chan School of Public Health by **Tiffany Hsu** and **Eric A. Franzosa**. Please direct questions to the [WAAFLE user group](waafle-users@googlegroups.com).

## Citation

(*A publication describing WAAFLE and its applications is currently in-prep. In the meantime, if you use WAAFLE in your work, please cite the WAAFLE website:* http://huttenhower.sph.harvard.edu/waafle)

## Installation

WAAFLE is a Python package and will *eventually* be installable via pip. In the meantime, if you have mercurial (`hg`) available, you can clone the WAAFLE repository to get started:

```
$ hg clone https://bitbucket.org/biobakery/waafle
```

Alternatively, you can download and extract the WAAFLE project directly:

```
$ wget https://bitbucket.org/biobakery/waafle/get/default.zip
$ unzip default.zip
```
To use WAAFLE, you'll need to add the src folder `waafle/waafle` to your `$PATH` and `$PYTHONPATH`.

## Software requirements

* Python 2.7+
* Python `numpy` (*version TBD*)
* NCBI `blastn` (*version TBD*)
* `bowtie2` (for `waafle_qc`; *version TBD*) 

## Database requirements

WAAFLE requires two input databases: 1) a **nucleotide sequence database** and 2) a corresponding **taxonomy file**.

### Formatting a sequence database for WAAFLE

WAAFLE performs a nucleotide-level search of metagenomic contigs against a collection of taxonomically annotated protein-coding genes (*not* complete genomes). A common way to build such a database is to combine a collection of microbial pangenomes of interest. The protein-coding genes should be organized in FASTA format and then indexed for use with `blastn`. For example, the FASTA file `waafledb.fnt` would be indexed as:

```
$ makeblastdb -in waafledb.fnt -dbtype nucl
```

WAAFLE expects the input FASTA sequence headers to follow a specific format. At a minimum, the headers must contain a unique sequence ID (`GENE_123` below) followed by `|` (pipe) followed by a taxon name or taxonomic ID (`SPECIES_456` below):

```
>GENE_123|SPECIES_456
```

Headers can contain additional `|`-delimited fields corresponding to functional annotations of the gene. These fields have the format `SYSTEM=IDENTIFIER` and multiple such fields can be included per header, as in:

```
>GENE_123|SPECIES_456|PFAM=P00789|EC=1.2.3.4
```
Headers are allowed to contain different sets of functional annotation fields. WAAFLE currently expects at most one annotation per annotation system per gene; this will be improved in future versions. (We currently recommend annotating genes against the [UniRef90 and UniRef50](http://www.uniprot.org/help/uniref) databases to enable link-outs to more detailed functional annotations in downstream analyses.)

WAAFLE assumes that the taxa listed in the sequence database file are all at the same taxonomic level (for example, all genera or all species or all strains).

### Formatting a WAAFLE taxonomy file

WAAFLE requires a taxonomy file to understand the taxonomic relationships among the taxa whose genes are included in the sequence database. The taxonomy file is a tab delimited list of child-parent relationships, for example:

```
k__Animalia      r__Root
p__Chordata      k__Animalia
c__Mammalia      p__Chordata
o__Primates      c__Mammalia
f__Hominidae     o__Primates
g__Homo          f__hominidae
s__Homo sapiens  g__Homo
```

While the format of this file is relatively simple, it has a number of critical structural constraints that must be obeyed:

* All taxon names/IDs used in the sequence database must appear is the taxonomy file.
 
* The file must contain a root taxon from which all other taxa descend (the root taxon should be named `r__Root`, as above).

* All taxon names/IDs used in the sequence database must be the same distance from the root.

The following properties of the taxonomy file are optional:

* Additional taxa *below* the level of the taxa in the sequence file can be included in the taxonomy file. For example, a species-level sequence database could contain isolate genomes as children of the species-level clades in the taxonomy file. (WAAFLE can use such entries as "leaf support" for LGT events.)

* We recommend prefixing taxonomic clades to indicate their level. For example, `g__Homo` identifies *Homo* as a genus above.

### Getting databases

(*WAAFLE-formatted databases compatible with the microbial pangenomes from http://metaref.org will be made available when the WAAFLE publication is submitted. In the meantime, the [WAAFLE demo](https://bitbucket.org/biobakery/waafle/src/default/docs/demo.md) provides small databases that are illustrative of the above formats and suitable for experimenting with the method.*)

## Inputs

An individual WAAFLE run requires two inputs (above and beyond the fixed sequence database and taxonomy files described above): 1) a file containing **metagenomic contigs** and (optionally) 2) a GFF file describing **gene coordinates** along those contigs.

Contigs should be provided as nucleotide sequences in FASTA format. Contigs are expected to have unique, BLAST-compatible headers. WAAFLE is optimized for working with fragmentary contigs from partially assembled metagenomes (spanning 2-20 genes, or roughly 1-20 kb). WAAFLE is not optimized to work with extremely long contigs (100s of kbs), scaffolds, or closed genomes. The WAAFLE developers recommend [MEGAHIT](https://github.com/voutcn/megahit) as a general-purpose metagenomic assembler.

The optional GFF file, if provided, should conform to the [GFF format]([https://useast.ensembl.org/info/website/upload/gff.html).

## Performing a WAAFLE analysis

Analyzing a set of contigs with WAAFLE consistent of three steps, one of which is optional. These three steps are carried out by three independent scripts: `waafle_search`, `waafle_genecaller` (optional), and `waafle_orgscorer`.

### Step 1: running `waafle_search`

`waafle_search` is a light wrapper around `blastn` to help guide the nucleotide-level search of your metagenomic contigs against a WAAFLE-formatted database (for example, it ensures that all of the non-default BLAST output fields required for downstream processing are generated).

A sample call to `waafle_search` with input contigs `contigs.fna` and a blast database located in `waafledb` would be:

```
$ waafle_search.py contigs.fna waafledb
```

By default, this produces an output file `contigs.blastout` in the same location as the input contigs. See the `--help` menu for additional configuration options:

```
usage: waafle_search.py [-h] [--blastn <path>] [--threads <int>]
                        [--out <path>]
                        query db

================================================================================
  ./hg/waafle/waafle/waafle_search.py: Step 1 in the WAAFLE pipeline
================================================================================

  This script executes a custom BLAST search of a set of contigs against a
  WAAFLE-formatted database.

================================================================================

positional arguments:
  query            contigs file (fasta format)
  db               path to WAAFLE BLAST database

optional arguments:
  -h, --help       show this help message and exit
  --blastn <path>  path to blastn binary
                   [default: $PATH]
  --threads <int>  number of CPU cores to use in blastn search
                   [default: 1]
  --out <path>     path for blast output file
                   [default: <derived from input>]
```

### Step 1.5 (optional): running `waafle_genecaller`

If the user chooses not to provide a GFF file along with their contigs, WAAFLE can identify gene coordinates of interest directly from the BLAST output produced in the previous step:

`$ waafle_search.py contigs.blastout`

This produces a file in GFF format called  `contigs.gff` for use in the next and last WAAFLE step. See the `--help` menu for additional configuration options:

```
usage: waafle_genecaller.py [-h] [--gff <path>] [--min-overlap <float>]
                            [--min-gene-length <int>] [--min-scov <float>]
                            [--stranded]
                            blastout

================================================================================
  ./hg/waafle/waafle/waafle_genecaller.py: (Optional) Step 1.5 in the WAAFLE pipeline
================================================================================

  Use the results of waafle_search to identify candidate gene loci in a set of
  contigs and output them as a GFF file for use in the next step. Users can
  optionally supply their own (independently-generated) GFF file.

================================================================================

positional arguments:
  blastout              (custom) blast output from waafle_search

optional arguments:
  -h, --help            show this help message and exit
  --gff <path>          path for (output) waafle gene calls (.gff)
                        [default: <derived from input>]
  --min-overlap <float>
                        if a large hit covers this fraction of a smaller hit, consider them part of the same gene group
                        [default: 0.1]
  --min-gene-length <int>
                        minimum allowed gene length
                        [default: 200]
  --min-scov <float>    (modified) scoverage filter for hits to gene catalog
                        [default: 0.75]
  --stranded            only merge hits into hits/genes of the same strandedness
                        [default: off]

```

### Step 2 : running `waafle_orgscorer`

The final and most critical step of a WAAFLE analysis is combining the BLAST output generated in Step 1 with a GFF file to 1) taxonomically score genes along the length of the input contigs and then 2) identify those contigs as derived from a single clade or a pair of clades (i.e. putative LGT). Assuming you have run steps 1 and 1.5 as described above, a sample call to `waafle_orgscorer` would be:

```
$ waafle_orgscorer.py \
  contigs.fna \
  contigs.blastout \
  contigs.gff \
  taxonomy.tsv
```

This will produce three output files which divide and describe your contigs as putative LGT contigs, single-clade (no-LGT) contigs, and unclassified contigs (e.g. those containing no genes):

* `contigs.lgt.tsv`
*  `contigs.no_lgt.tsv`
* `contigs.unclassified.tsv`

These files and their formats are described in more detailed below (see "Interpretting WAAFLE outputs").

`waafle_orgscorer` offers many options for fine-tuning your analysis. The various analysis parameters have been pre-optimized for maximum specificity on both short contigs (containing as little as two partial genes) and longer contigs (10s of genes). These options are detailed in the `--help` menu:

```
usage: waafle_orgscorer.py [-h] [--outdir <path>] [--basename <str>]
                           [--write-details] [--quiet] [-k1 <0.0-1.0>]
                           [-k2 <0.0-1.0>]
                           [--disambiguate-one <report-best/meld>]
                           [--disambiguate-two <report-best/jump/meld>]
                           [--range <float>] [--jump-taxonomy <1-N>]
                           [--allow-lca] [--ambiguous-fraction <0.0-1.0>]
                           [--sister-penalty <0.0-1.0>] [--clade-genes <1-N>]
                           [--clade-leaves <1-N>]
                           [--weak-loci <ignore/penalize/assign-unknown>]
                           [--transfer-annotations <lenient/strict/very-strict>]
                           [--min-overlap <0.0-1.0>] [--min-gene-length <int>]
                           [--min-scov <float>] [--stranded]
                           contigs blastout gff taxonomy

================================================================================
  ./hg/waafle/waafle/waafle_orgscorer.py: Step 2 in the WAAFLE pipeline
================================================================================

  Merges blast hits into genes on contigs-of-interest. Uses corresponding
  taxonomy file, and the WAAFLE algorithm, to identify contigs that are best
  explained by a single clade vs. a pair of clades. The latter events correpond
  to putative LGTs.

================================================================================

optional arguments:
  -h, --help            show this help message and exit

required inputs:
  contigs               contigs file (.fasta format)
  blastout              output of waafle_search for one set of contigs (.blastout)
  gff                   gene calls (from waafle_genecaller or user-supplied) for <contigs> (.gff)
  taxonomy              taxonomy file for the blast database used to make <blastout>

output formatting:
  --outdir <path>       directory for writing output files
                        [default: .]
  --basename <str>      basename for output files
                        [default: derived from input]
  --write-details       make an additional output file with per-gene clade scores
                        [default: off]
  --quiet               don't show running progress
                        [default: off]

main parameters:
  -k1 <0.0-1.0>, --one-clade-threshold <0.0-1.0>
                        minimum per-gene score for explaining a contig with a single clade
                        [default: 0.5]
  -k2 <0.0-1.0>, --two-clade-threshold <0.0-1.0>
                        minimum per-gene score for explaining a contig with a pair of clades (putative LGT)
                        [default: 0.8]
  --disambiguate-one <report-best/meld>
                        what to do when other one-clade explanations fall within <--range> of the best explanation
                        [default: meld]
  --disambiguate-two <report-best/jump/meld>
                        what to do when other two-clade explanations fall within <--range> of the best explanation
                        [default: meld]
  --range <float>       when disambiguating, consider explanations within <--range> of the best explanation
                        [default: 0.05]
  --jump-taxonomy <1-N>
                        before starting, perform 1+ 'jumps' up the taxonomy (e.g. species->genus)
                        [default: off]

post-detection LGT filters:
  --allow-lca           when melding LGT clades, allow the LGT LCA to occur as a melded clade
                        [default: off]
  --ambiguous-fraction <0.0-1.0>
                        allowed fraction of ambiguous (A OR B) gene length in a putative A+B contig
                        [default: 0.1]
  --sister-penalty <0.0-1.0>
                        allowed mean prevalence of missing genes in sisters of LGT clades (or just recipient if known)
                        [default: 0.0]
  --clade-genes <1-N>   required minimum genes assigned to each LGT clade
                        [default: off]
  --clade-leaves <1-N>  required minimum leaf count supporting each LGT clade (or just recipient if known)
                        [default: off]

gene-hit merge parameters:
  --weak-loci <ignore/penalize/assign-unknown>
                        method for handling loci that are never assigned to known clades
                        [default: ignore]
  --transfer-annotations <lenient/strict/very-strict>
                        stringency of gene annotation transfer to loci
                        [default: strict]
  --min-overlap <0.0-1.0>
                        only merge hits into genes if the longer of the two covers this portion of the shorter
                        [default: 0.1]
  --min-gene-length <int>
                        minimum allowed gene length
                        [default: 200]
  --min-scov <float>    (modified) scoverage filter for hits to gene catalog
                        [default: 0.75]
  --stranded            only merge hits into hits/genes of the same strandedness
                        [default: off]

```

## Interpreting WAAFLE outputs

The `contigs.lgt.tsv` output file lists the details of putative LGT contigs. Its fields are a superset of the types of fields included in the other output files. The following represents the first two lines/rows of a `contigs.lgt.tsv` file *transposed* such that first line (column headers) is shown as the first column and the details of the first LGT contig (second row) are shown as the second column:

```
CONTIG_NAME          12571                                                                                                                       
CALL                 lgt                                                                                                                         
CONTIG_LENGTH        9250                                                                                                                        
MIN_MAX_SCORE        0.856                                                                                                                       
AVG_MAX_SCORE        0.965                                                                                                                       
SYNTENY              AABAAAA                                                                                                                     
DIRECTION            B>A                                                                                                                         
CLADE_A              s__Ruminococcus_bromii                                                                                                      
CLADE_B              s__Faecalibacterium_prausnitzii                                                                                             
LCA                  f__Ruminococcaceae                                                                                                          
MELDED_A             --                                                                                                                          
MELDED_B             --                                                                                                                          
TAXONOMY_A           r__Root|k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Ruminococcaceae|g__Ruminococcus|s__Ruminococcus_bromii  
TAXONOMY_B           r__Root|k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Ruminococcaceae|g__Faecalibacterium|s__Faecalibacterium_prausnitzii
LOCI                 252:668:-|792:1367:-|1557:2360:-|2724:3596:-|4540:5592:+|5608:7977:+|8180:8425:+                                            
ANNOTATIONS:UNIPROT  R5E4K6|D4L7I2|D4JXM0|D4L7I1|D4L7I0|None|D4L7H8 
```

The fields in detail:

* **`CONTIG_NAME`**: the name of the contig from the input FASTA file.

* **`CALL`**: indicates that this was an LGT contig.

* **`CONTIG_LENGTH`**: the length of the contig in nucleotides.

* **`MIN_MAX_SCORE`**: the minimum score for the pair of clades explaining the contig along the length of the contig. (i.e. the score for identifying this as a putative LGT contig, with a default threshold of 0.8.)

* **`AVG_MAX_SCORE`**: the average score for the pair of clades explaining the contig (used for ranking multiple potential explanations of the contig).

* **`SYNTENY`**: the pattern of genes assigned to the A or B clades along the contig. `*` indicates that the gene could be contributed by either clade; `~` indicates an ignored gene; `!` indicates a problem (should not occur).

* **`DIRECTION`**: indicates this as a putative B-to-A transfer event, as determined from synteny (A genes flank the inserted B gene). `A?B` indicates that the direction was uncertain.

* **`CLADE_A`**: the name of clade A.

* **`CLADE_B`**: the name of clade B.

* **`LCA`**: the lowest common ancestor of clades A and B. A higher-level LCA indicates a more remote LGT event.
 
* **`MELDED_A`**: if using a meld reporting option, the individual melded clades will be listed here. For example, if a contig could be explained by a transfer from *Genus_1 species_1.1* to either *Genus_2 species_2.1* or *Genus_2 species_2.2*, this field would list `species_2.1; species 2.2` and *Genus 2* would be given as `CLADE_A`.

* **`MELDED_B`**: *see above*.

* **`TAXONOMY_A`**: the full taxonomy of `CLADE_A`.

* **`TAXONOMY_B`**: the full taxonomy of `CLADE_B`.

* **`LOCI`**: Ccordinates of the loci (genes) that were considered for this contig in format `START:STOP:STRAND`.

* **`ANNOTATIONS:UNIPROT`**: indicates that UniProt annotations were provided for the genes in the input sequence database (in format `UNIPROT=IDENTIFIER`). The best-scoring UniProt annotation for each gene is given here. (Additional annotations would appear as additional, similarly-formatted columns in the output.)

## Filtering questionable contigs

The WAAFLE workflow described above has been optimized to distinguish LGTs from other biological events (e.g. gene deletion). However, it cannot intrinsically identify spurious LGTs resulting from misassembly (e.g. chimerism). For this, we provide a separate method, `waafle_qc`.

`waafle_qc` maps reads to contigs to evaluate gene-gene junctions (involved in LGTs or otherwise). `waafle_qc` considers a contig to be "OK" if all of its junctions are supported. A junction is supported if one of the followed two conditions to be true:

1. Multiple sequencing fragments (paired reads) engulf the junction. This criterion is most applicable to junctions that are smaller than typical sequencing insert sizes (0-200 nts).

2. The junction is well covered relative to its flanking genes. This criterion is most applicable to longer junctions (>200 nts).

A sample call to `waafle_qc` might look like:

```
$ waafle_qc.py \
  contigs.fna \
  contigs.gff \
  --reads1 contigs_reads.1.fq \
  --reads2 contigs_reads.2.fq \
```

With this call, `waafle_qc` will use [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) to index the contigs and then align the input reads (pairwise) against the index to produce a SAM file. (`waafle_qc` can also interpret a mapping from an existing SAM file.) The alignment results are then interpreted to score junctions and contigs, producing an output file for each.

A sample report for an individual junction looks like:

```
CONTIG             contig_k119_11579 
GENE1              4975:6210:-          
GENE2              6279:6896:+          
GAP                68                   
HITS_JUNCTION      18                   
HITS_OK            True                 
COVERAGE_GENE1     32.4256              
COVERAGE_GENE2     33.3706              
COVERAGE_JUNCTION  43.0000              
COVERAGE_RATIO     1.3071               
COVERAGE_OK        True                 
ACCEPTABLE         True
```

This junction was supported by 18 sequencing fragments (`HITS_JUNCTION`) and was similar in coverage to its flanking genes (`COVERAGE_RATIO`), thus passing both of the `waafle_qc` filters outlined above.

A sample report for a pair of contigs (one passing, one failing) looks like:

```
CONTIG             contig_k119_10097  contig_k119_10099 
SUMMARY            OK                    FAILED               
LENGTH             2176                  4154                 
LOCI               3                     5                    
FAILING_JUNCTIONS  0                     1                    
FAILURE_RATE       0.00                  0.25
```

`waafle_qc` can be tuned to filter junctions more or less stringently, as well as to filter contigs based on other properties (e.g. length). Consult the `--help` menu for a full list of options:

```
usage: waafle_qc.py [-h] [--reads1 <path>] [--reads2 <path>] [--sam <path>]
                    [--tmpdir <path>] [--outdir <path>] [--basename <str>]
                    [--write-detailed-output] [--min-overlap-sites <int>]
                    [--min-junction-hits <int>] [--min-junction-ratio <float>]
                    [--min-contig-genes <int>] [--min-contig-length <int>]
                    [--max-contig-failures <float>] [--bowtie2-build <path>]
                    [--bowtie2 <path>] [--threads <int>]
                    contigs gff

================================================================================
  waafle_qc.py: Performs coverage-level quality control on contigs
================================================================================

  This script maps reads to contigs (or analyzes an existing mapping) to
  identify contigs whose junctions aren't supported by reads. A short junction
  is supported if mate-pairs engulf the junction. A long junction (i.e. too long
  to be engulfed by a mate-pair) is supported if its coverage is "reasonably
  similar" to the flanking genes' coverages (e.g. at least half their mean).
  This script will also filter out contigs that are too short or which contain
  too few genes.

================================================================================

optional arguments:
  -h, --help            show this help message and exit

required inputs:
  contigs               contigs file (fasta format)
  gff                   GFF file for provided contigs

provide paired reads or a .sam file:
  --reads1 <path>       sequencing reads (mate-1)
  --reads2 <path>       sequencing reads (mate-2)
  --sam <path>          sam file (from existing alignment)

output options:
  --tmpdir <path>       where to place temp outputs
                        [default: ./]
  --outdir <path>       where to place main outputs
                        [default: ./]
  --basename <str>      basename for output files
                        [default: <derived from input>]
  --write-detailed-output
                        write out coverage values for all sites and all junctions
                        [default: off]

filtering parameters:
  --min-overlap-sites <int>
                        minimum nucleotide overlap for counting a read-gene hit
                        [default: 25]
  --min-junction-hits <int>
                        minimum read-hits to 'ok' a junction
                        [default: 2]
  --min-junction-ratio <float>
                        minimum coverage (relative to flanking genes) to 'ok' a junction
                        [default: 0.5]
  --min-contig-genes <int>
                        minimum gene count for the contig
                        [default: 2]
  --min-contig-length <int>
                        minimum length for the contig
                        [default: 500]
  --max-contig-failures <float>
                        allowed fraction of failing junctions
                        [default: 0]

bowtie2 options:
  --bowtie2-build <path>
                        path to bowtie2-build
                        [default: $PATH]
  --bowtie2 <path>      path to bowtie2
                        [default: $PATH]
  --threads <int>       number of cpu cores to use in bowtie2 alignment
                        [default: 1]
``` 