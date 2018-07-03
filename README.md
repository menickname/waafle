# WAAFLE

WAAFLE (*W*orkflow to *A*nnotate *A*ssemblies and *F*ind *L*GT Events) is a method for finding novel LGT (*L*ateral *G*ene *T*ransfer) events in assembled metagenomes.

## Authors

WAAFLE was developed in the Huttenhower Lab at the Harvard T.H. Chan School of Public Health by **Tiffany Hsu** and **Eric A. Franzosa**.

## Citation

*A publication describing WAAFLE and its applications is currently in-prep. In the meantime, if you use WAAFLE in your work, please cite the WAAFLE website:* http://huttenhower.sph.harvard.edu/waafle.

## Installation

WAAFLE is a Python package and will eventually be installable via pip. In the meantime, if you have mercurial (`hg`) available, you can clone the WAAFLE repository to get started:

```
$ hg clone https://bitbucket.org/biobakery/waafle
```

Alternatively, you can download and extract the WAAFLE project directly:

```
$ wget https://bitbucket.org/biobakery/waafle/get/default.zip
$ unzip default.zip
```

## Software requirements

* Python 2.7+
* Python numpy (versions **TBD**)
* NCBI BLAST (versions **TBD**)

## Database requirements

WAAFLE requires two input databases: a **nucleotide sequence database** and a **taxonomy file**. 

### Formatting a WAAFLE sequence database

The nucleotide sequence database should be organized in FASTA format and then indexed for use with BLAST. For example, the FASTA file `waafle.fnt` would be indexed as:

```
$ makeblastdb -in waafle.fnt -dbtype nucl
```

WAAFLE expects the sequence headers to follow a specific format. At a minimum, the headers must contain a unique sequence ID followed by a taxon name or taxonomic ID:

```
>GENE_123|SPECIES_456
```

Headers can contain additional piped fields corresponding to functional annotations of the gene, as in:

```
>GENE_123|SPECIES_456|PFAM=P00789|EC=1.2.3.4
```

There is no limit to the number of functional annotation fields per header. Headers are allowed to contain different sets of functional annotation fields.

### Formatting a WAAFLE taxonomy file

The taxonomy file is a tab delimited list of child-parent relationships, for example:

```
Animalia      Root
Chordata      Animalia
Mammalia      Chordata
Primates      Mammalia
Hominidae     Primates
Homo          hominidae
Homo sapiens  Homo
```

A number of requirements must be satisfied by this file:

* The file must contain a root taxon from which all others descend (listed as `Root` above).
* All taxon names/IDs used in the sequence database must appear is the taxonomy file.
* All taxon names/IDs used in the sequence database must be the same distance from the root (e.g. all species or all genera).
* The taxonomy file can contain additional levels below the level used in the sequence file (for example, individual genome IDs within a parent species).

### Getting databases

WAAFLE-formatted databases compatible with the microbial pangenomes from http://metaref.org will be made available when the WAAFLE publication is submitted. In the meantime, the WAAFLE demo provides small databases that are illustrative of the above formats and suitable for experimenting with the method.

## Inputs

An individual WAAFLE run requires two inputs (above and beyond the fixed database and taxonomy files described above): a file containing metagenomic contigs and (optionally) a GFF file describing gene coordinates along those contigs.

The contigs should be provided as nucleotide sequences in FASTA format. Contigs are expected to have unique, BLAST-compatible headers. WAAFLE is optimized for working with fragmentary contigs from partially assembled metagenomes (spanning 2-20 genes, or roughly 1-20 kb). WAAFLE is not optimized to work with extremely long contigs (100s of kbs), scaffolds, or closed genomes.

The optional GFF file, if provided, should conform to the established [GFF format]([https://useast.ensembl.org/info/website/upload/gff.html).

## Workflow

### 1. WAAFLE search

...

### 1.5 (Optional) WAAFLE gene caller 

...

### 2. WAAFLE analysis

...

## Interpreting output

* lgt
* no_lgt
* ambiguous
* details