uLTRA
===========

uLTRA is a tool for splice alignment of long transcriptomic reads to a genome, guided by a database of exon annotations. uLTRA takes reads in fast(a/q) and a genome annotation as input and outputs a SAM-file. The SAM-file includes information on which splice sites are found and if the read is a full splice match (and to which transcript), incomplete splice match, Novel in catalog, or novel not in the catalog, as defined in [SQANTI](https://github.com/ConesaLab/SQANTI). uLTRA is highly accurate when aligning to small exons [see some examples](https://github.com/ksahlin/ultra/tree/master/data/images). 

uLTRA is distributed as a python package supported on Linux / OSX with python v>=3.4. [![Build Status](https://travis-ci.org/ksahlin/uLTRA.svg?branch=master)](https://travis-ci.org/ksahlin/uLTRA).

### New since v0.0.2
Since v0.0.2, uLTRA can be used as an **end-to-end aligner for annotation and detection of novel genes or isoforms** (default mode). This is because uLTRA (>=v0.0.2) now incorporates [minimap2](https://github.com/lh3/minimap2). [minimap2](https://github.com/lh3/minimap2) is run upon start of uLTRA, and the results are used both for (i) not aligning reads with uLTRA which had a primary alignment to regions not indexed by uLTRA (e.g. genomic regions or unannotated genes) and (ii) to consult at end of program which aligner had a better fit (based on cigar) of the primary alignment and chose this alignment to be primary. uLTRA still uses its own alignment algorithm to align to and around all annotated gene regions. uLTRA can therefore, at worst, be seen as an advanced wrapper around minimap2 that refines alignments around annotated regions. See updated `CREDITS` when using this version. uLTRA can still be used as a stand alone aligner as presented in our [preprint](https://www.biorxiv.org/content/10.1101/2020.09.02.279208v1) by specifying `--disable_mm2`.

### New since v0.0.3

uLTRA is now used less than half of the memory used in previous versions and is about 20% faster.



Table of Contents
=================

  * [INSTALLATION](#INSTALLATION)
    * [Using conda](#Using-conda)
    * [Downloading source from GitHub](#Downloading-source-from-github)
  * [USAGE](#USAGE)
    * [Indexing](#Indexing)
    * [aligning](#Aligning)
    * [Output](#Output)
    * [Parameters](#Parameters)
  * [CREDITS](#CREDITS)
  * [LICENCE](#LICENCE)



INSTALLATION
----------------

## Using conda
Conda is the preferred way to install uLTRA. You can either clone this repository and 
run the script `INSTALL.sh` or you can perform step 1-6 below manually for more control.

### Installation with INSTALL.sh script

```
git clone https://github.com/ksahlin/uLTRA.git --depth 1
cd uLTRA
./INSTALL.sh [An install directory in your PATH]
```
The install script is tested in bash environment and will perform the steps 1-6 below automatically
for you. 

You need to activate the conda environment "ultra" to run uLTRA as:
```
conda activate ultra
```

### Manual installation


#### 1. Create and activate a new environment called ultra

```
conda create -n ultra python=3 pip 
conda activate ultra
```

#### 2. Install uLTRA 

```
pip install ultra-bioinformatics
```

#### 3. Install third party MEM finder [slaMEM](https://github.com/fjdf/slaMEM) and aligner [minimap2](https://github.com/lh3/minimap2)

```
git clone git@github.com:fjdf/slaMEM.git
cd slaMEM
make 
```
Place the generated binary `slaMEM` in your path. Minimap2 can be installed through conda with `conda install -c bioconda minimap2`, or [manually](https://github.com/lh3/minimap2). 

#### 4. You should now have 'uLTRA' installed; try it:

```
uLTRA --help
```

Upon start/login to your server/computer you need to activate the conda environment "ultra" to run uLTRA as:
```
conda activate ultra
```

#### 5. Test uLTRA

Download/use test data available in this repository [here](https://github.com/ksahlin/ultra/tree/master/test) and run: 

```
uLTRA pipeline [/your/local/directory/to/test]/SIRV_genes.fasta  \
               [/your/local/directory/to/test]/SIRV_genes_C_170612a.gtf  \
               [/your/local/directory/to/test]/reads.fa outfolder/  [optional parameters]
```
Specify full path to annotation, otherwise `gffutils` will complain.


#### 6. (Optional) Install of MUMmer 

While MUMmer is usually not used in uLTRA, if slaMEM [fails](https://github.com/fjdf/slaMEM/issues/3), uLTRA falls back on finding MEMs with MUMmer until the slaMEM bug has been fixed. In this corner case, uLTRA needs MUMmer avaialble in the path. MUMmer can be installed with

```
conda install --yes -c bioconda mummer
```


## Downloading source from GitHub


Make sure the below-listed dependencies are installed (installation links below). Versions in parenthesis are suggested as uLTRA has not been tested with earlier versions of these libraries. However, uLTRA may also work with earlier versions of these libraries. All below dependencies except `slaMEM` can be installed as `pip install X` or through conda.
* [parasail](https://github.com/jeffdaily/parasail-python)
* [edlib](https://github.com/Martinsos/edlib)
* [pysam](http://pysam.readthedocs.io/en/latest/installation.html) (>= v0.11)
* [dill](https://pypi.org/project/dill/)
* [intervaltree](https://github.com/chaimleib/intervaltree/tree/master/intervaltree)
* [gffutils](https://pythonhosted.org/gffutils/)
* [slaMEM](https://github.com/fjdf/slaMEM)


With these dependencies installed. Run

```sh
git clone https://github.com/ksahlin/uLTRA.git
cd uLTRA
./uLTRA
```


USAGE
-------

uLTRA can be used with either PacBio Iso-Seq or ONT cDNA/dRNA reads. 


### Indexing

First, we construct the data structures used in uLTRA using a genome annotation GTF file and a genome fasta file.
Make sure to specify full path to annotation, otherwise `gffutils` will complain.

```
uLTRA index genome.fasta  /full/path/to/annotation.gtf  outfolder/  [parameters]
```


### Aligning

For example

```
uLTRA align genome.fasta reads.[fa/fq] outfolder/  --ont --t 48   # ONT cDNA reads using 48 cores
uLTRA align genome.fasta reads.[fa/fq] outfolder/  --isoseq --t 48 # PacBio isoseq reads
uLTRA align genome.fasta reads.[fa/fq] outfolder/  --k 14  --t 48 # PacBio dRNA reads or reads with >10-12% error rate
```

uLTRA's index takes about 7Gb for human, and each instance needs a separate copy of the index (if parallelized, that is, `--t` greater than 1). So if you have a computer/cluster with 8Gb per core and n cores it is straightforward to set `--t n-1` (n-1 to leave some space for the main process). 

### Pipeline

Performs all the steps in one

```
uLTRA pipeline genome.fasta annotation.gtf reads.fa outfolder/  [parameters]
```

#### Output

uLTRA outputs a SAM-file with alignments to the genome. In addition, it outputs to extra tags describing whether all the splices sites are known and annotated (FSM), new splice combinations (NIC), etc. For details see the definitions of notations in the [Sqanti paper](https://genome.cshlp.org/content/28/7/1096).



CREDITS
----------------

Please cite [1] when using uLTRA. If you are using uLTRA v0.0.2 or later **please also cite** [minimap2](https://github.com/lh3/minimap2) as uLTRA incorporates minimap2 for alignment of some reads. For example "We aligned reads to the genome using uLTRA [1], which incorporates minimap2 [CIT].".

1. Kristoffer Sahlin, Veli Makinen. 2020. "Accurate spliced alignment of long RNA sequencing reads" [preprint available here](https://www.biorxiv.org/content/10.1101/2020.09.02.279208v1).

Bib record: 

@article {Sahlin2020.09.02.279208,
  author = {Sahlin, Kristoffer and Makinen, Veli},
  title = {Accurate spliced alignment of long RNA sequencing reads},
  elocation-id = {2020.09.02.279208},
  year = {2020},
  doi = {10.1101/2020.09.02.279208},
  publisher = {Cold Spring Harbor Laboratory},
  abstract = {Long-read RNA sequencing techniques are quickly establishing themselves as the primary sequencing technique to study the transcriptome landscape. Many such analyses are dependent upon splice alignment of reads to the genome. However, the error rate and sequencing length of long-read technologies create new challenges for accurately aligning these reads. We present an alignment method uLTRA that, on simulated and synthetic data, shows higher accuracy over state-of-the-art with substantially higher accuracy for small exons. We show several examples on biological data where uLTRA aligns to known and novel isoforms with exon structures that are not detected with other aligners. uLTRA is available at https://github.com/ksahlin/ultra.Competing Interest StatementThe authors have declared no competing interest.},
  URL = {https://www.biorxiv.org/content/early/2020/09/03/2020.09.02.279208},
  eprint = {https://www.biorxiv.org/content/early/2020/09/03/2020.09.02.279208.full.pdf},
  journal = {bioRxiv}
}


LICENCE
----------------

GPL v3.0, see [LICENSE.txt](https://github.com/ksahlin/uLTRA/blob/master/LICENCE.txt).


