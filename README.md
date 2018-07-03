# WAAFLE

WAAFLE (*W*orkflow to *A*nnotate *A*ssemblies and *F*ind *L*GT Events) is a method for
finding novel LGT (*L*ateral *G*ene *T*ransfer) events in assembled metagenomes.

## Authors

WAAFLE was developed in the Huttenhower Lab at the Harvard T.H. Chan School of Public Health by **Tiffany Hsu** and **Eric A. Franzosa**.

## Citation

*A publication describing WAAFLE and its applications is currently in-prep. In the meantime, if you use WAAFLE in your work, please cite the WAAFLE website: * http://huttenhower.sph.harvard.edu/waafle.

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
* Numpy
* NCBI BLAST

## Database requirements

### Getting databases

* ChocoPhlAn -> WAAFLE db
* Taxonomy

### Making databases / database format

* Sequences
* Taxonomy

## Inputs

* Contigs
* Optional GFF

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