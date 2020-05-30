# torkel
Aligns long reads to splicing graphs and predicts FSM, ISM, NIC and NNC

Dependencies: gffutils


uLTRA
===========

uLTRA is a tool for aligning and annotating Long Transcriptomic Reads to a database of known transcripts. It is tailored for finding splice matches to annotations and provided genomic alignment coordinates to a sam file. The SAM file includes information of which splice sites are found and if the reads is a full splice match, incomplete splice match, Novel in catalogue, or novel not in catalogue, as defined in [cite Sqanti here]. 


uLTRA is distributed as a python package supported on Linux / OSX with python v>=3.4. [![Build Status](https://travis-ci.org/ksahlin/uLTRA.svg?branch=master)](https://travis-ci.org/ksahlin/uLTRA).

Table of Contents
=================

  * [INSTALLATION](#INSTALLATION)
    * [Using conda](#Using-conda)
    * [Using pip](#Using-pip)
    * [Downloading source from GitHub](#Downloading-source-from-github)
    * [Dependencies](#Dependencies)
  * [USAGE](#USAGE)
    * [constructing](#Constructing)
    * [aligning](#Aligning)
    * [Output](#Output)
    * [Parameters](#Parameters)
  * [CREDITS](#CREDITS)
  * [LICENCE](#LICENCE)



INSTALLATION
----------------

### Using conda
Conda is the preferred way to install uLTRA.

1. Create and activate a new environment called ultra

```
conda create -n ultra python=3 pip 
source activate ultra
```

2. Install uLTRA 

```
pip install uLTRA
```
3. You should now have 'uLTRA' installed; try it:
```
uLTRA --help
```

Upon start/login to your server/computer you need to activate the conda environment "ultra" to run uLTRA as:
```
source activate ultra
```

### Using pip 

To install uLTRA, run:
```
pip install uLTRA
```
`pip` will install the dependencies automatically for you. `pip` is pythons official package installer and is included in most python versions. If you do not have `pip`, it can be easily installed [from here](https://pip.pypa.io/en/stable/installing/) and upgraded with `pip install --upgrade pip`. 


### Downloading source from GitHub

#### Dependencies

Make sure the below listed dependencies are installed (installation links below). Versions in parenthesis are suggested as uLTRA has not been tested with earlier versions of these libraries. However, uLTRA may also work with earliear versions of these libaries.
* [parasail](https://github.com/jeffdaily/parasail-python)
* [pysam](http://pysam.readthedocs.io/en/latest/installation.html) (>= v0.11)
* dill


With these dependencies installed. Run

```sh
git clone https://github.com/ksahlin/uLTRA.git
cd uLTRA
./uLTRA
```


USAGE
-------

uLTRA can be used with either Iso-Seq or ONT reads. 
 


### Constructing

Fist we construct splicing structure and extract reference sequnces from a genome.

```
uLTRA prep_splicing  all_genes.gtf outfolder/  [parameters]
```

```
uLTRA prep_seqs  genome.fasta  outfolder/  [parameters]
```


### Aligning

```
uLTRA align  genome.fasta  reads.fast[a/q] outfolder/  [parameters]
```


### Pipeline

Perforns all thee steps in one

```
uLTRA pipeline  all_genes.gtf   genome.fasta  reads.fast[a/q] outfolder/  [parameters]
```

#### Output

uLTRA outputs a sam file with alignments to the genome. In addition, it outputs to extra tags: 



CREDITS
----------------

Please cite [1] when using uLTRA.

1. Kristoffer Sahlin, Veli Makinen (2019) "Long read transcript annotation with uLTRA using co-linear chaining".

Bib record: 


LICENCE
----------------

GPL v3.0, see [LICENSE.txt](https://github.com/ksahlin/uLTRA/blob/master/LICENCE.txt).


