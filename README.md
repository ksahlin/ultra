uLTRA
===========
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/ultra_bioinformatics/README.html) [![Build Status](https://travis-ci.org/ksahlin/uLTRA.svg?branch=master)](https://travis-ci.org/ksahlin/uLTRA)


uLTRA is a tool for splice alignment of long transcriptomic reads to a genome, guided by a database of exon annotations. uLTRA takes reads in fast(a/q) and a genome annotation as input and outputs a SAM-file. The SAM-file includes information on which splice sites are found and if the read is a full splice match (and to which transcript), incomplete splice match, Novel in catalog, or novel not in the catalog, as defined in [SQANTI](https://github.com/ConesaLab/SQANTI). uLTRA is particularly accurate when aligning to small exons [see some examples](https://github.com/ksahlin/ultra/tree/master/data/images). 

uLTRA is distributed as a python package supported on Linux / OSX with python v>=3.4. 


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
=================

There is a [bioconda recipe](https://bioconda.github.io/recipes/ultra_bioinformatics/README.html), [docker image](https://quay.io/repository/biocontainers/ultra_bioinformatics?tab=tags), and a [singularity container](https://depot.galaxyproject.org/singularity/ultra_bioinformatics%3A0.0.4--pyh5e36f6f_1) of uLTRA v0.0.4 created by [sguizard](https://github.com/sguizard). You can use, e.g., the bioconda recipe for an easy automated installation. 

If a newer version of uLTRA is not available through bioconda (or you simply want more control of/customize your installation), alternative ways of installations are provided below. Current version of uLTRA is 0.0.4 (see changelog at end of this readme).

## Using conda

Conda is the preferred way to install uLTRA. You can either clone this repository and 
run the script `INSTALL.sh` or you can perform step 1-6 below manually for more control.

### Installation with INSTALL.sh script

```
git clone https://github.com/ksahlin/uLTRA.git --depth 1
cd uLTRA
./INSTALL.sh <install_directory>
```
The install script is tested in bash environment and will perform the steps 1-6 below automatically
for you. You need to have the `<install_directory>` included in your shell path or, alternatively, move the binaries 
`minimap2`, `slaMEM`, `StrobeMap` installed to `<install_directory>` to a directory in your path so that uLTRA finds them.

To run uLTRA, you need to activate the conda environment "ultra":
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
uLTRA pipeline [/your/full/path/to/test]/SIRV_genes.fasta  \
               /your/full/path/to/test/SIRV_genes_C_170612a.gtf  \
               [/your/full/path/to/test]/reads.fa outfolder/  [optional parameters]
```

Specify the **absolute path** to the GTF-file on your system, otherwise `gffutils` will complain and giva a cryptic `ValueError: unknown url type:` error message. Outfile will be `outfolder/reads.sam`, unless you specify your custom prefix filename with `--prefix`. **Update: As from version 0.0.4.1 and upwards you don't need the absolute path.** 


#### 6.(Optional) Install of StrobeMap

Using NAM seeds is new since version 0.0.4. It can reduce runtime, disk usage and provide fixed memory usage to default MEM finding. See changlog at end of this README. [StrobeMap](https://github.com/ksahlin/strobemers) is installed on Linux with:

```
wget https://github.com/ksahlin/strobemers/raw/main/strobemers_cpp/binaries/Linux/StrobeMap-0.0.2
mv StrobeMap-0.0.2 StrobeMap
chmod +x StrobeMap
```
Place the generated binary `StrobeMap` in your path.

#### 7. (Optional) Install of MUMmer 

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
* [StrobeMap](https://github.com/ksahlin/strobemers)

With these dependencies installed. Run

```sh
git clone https://github.com/ksahlin/uLTRA.git
cd uLTRA
./uLTRA
```


USAGE
=================

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
uLTRA align genome.fasta reads.[fa/fq] outfolder/  --ont --t 8   # ONT cDNA reads using 8 cores
uLTRA align genome.fasta reads.[fa/fq] outfolder/  --isoseq --t 8 # PacBio isoseq reads
```

You can set a custom location of where to get the index from using `--index [PATH]`. Otherwise, uLTRA will try to read the index from the `outfolder/` by default. The aligned reads will be written to `outfolder/reads.sam` unless `--prefix` is set. For example, `--prefix sample_X` will output the reads in `outfolder/sample_X.sam`.

### Pipeline

Performs all the steps in one

```
uLTRA pipeline genome.fasta /full/path/to/annotation.gtf reads.fa outfolder/  [parameters]
```

#### Output

uLTRA outputs a SAM-file with alignments to the genome. In addition, it outputs to extra tags describing whether all the splices sites are known and annotated (FSM), new splice combinations (NIC), etc. For details see the definitions of notations in the [Sqanti paper](https://genome.cshlp.org/content/28/7/1096).



CREDITS
----------------

Please cite [1] when using uLTRA. **Please also cite** [minimap2](https://github.com/lh3/minimap2) as uLTRA incorporates minimap2 for alignment of some genomic reads outside indexed regions. For example "We aligned reads to the genome using uLTRA [1], which incorporates minimap2 [CIT].".

1. Kristoffer Sahlin, Veli Mäkinen, Accurate spliced alignment of long RNA sequencing reads, Bioinformatics, Volume 37, Issue 24, 15 December 2021, Pages 4643–4651, https://doi.org/10.1093/bioinformatics/btab540


LICENCE
----------------

GPL v3.0, see [LICENSE.txt](https://github.com/ksahlin/uLTRA/blob/master/LICENCE.txt).


VERSION INFO
---------------

### New since v0.0.4.1

Removed requirement to specify full path to GTF annotation, making implementation in [nf-core modules](https://github.com/nf-core/modules/pull/871#pullrequestreview-787413998) and pipelines easier.

### New since v0.0.4

An option `--use_NAM_seeds` is added. This parameter changes the seeding of MEMs to NAMs (with strobemers). If `--use_NAM_seeds` is specified, uLTRA calls `StrobeMap` (binary can be acquired [here](https://github.com/ksahlin/strobemers/tree/main/strobemers_cpp/binaries) and put `StrobeMap` in your path). NAM seeding makes uLTRA faster and produces much smaller intermediate files. The speed improvement is soemwhere between 15-70%, depending on the number of threads. The more threads the better speed improvement.

The memory usage with `--use_NAM_seeds` is "fixed" regardless of number of cores/threads (about ~80-90Gb for human genome). This is in contrast to the default version where memory grows with number of cores. This makes `--use_NAM_seeds` less memory consuming than the default option (MEMs) when uLTRA is given about 18 cores or more (with `--t`), and more memory consuming than the default version for `t < 18`. 

The alignment accuracy is largely the same. Using NAM seeds decreases the accuracy of about 0.01%-0.05% compared to MEMs (i.e., 1 alignment in every 2,000-10,000). Accuracy is measured as `(correct alignments)/(total reads)`. An alignment is correct if all splice sites are aligned to correctly and exactly without offset allowed. Everything else is classified as incorrect. This is the most stringent critera for correct and was evaluated on simulated data.

Due to the "fixed" memory usage, faster runtime, and smaller intermediate file size, I recommend `--use_NAM_seeds` option for large datasets (>5M reads) if running on nodes with more than 20 cores.

### New since v0.0.3

uLTRA now uses less than half of the memory used in previous versions and is about 20% faster.

### New since v0.0.2
Since v0.0.2, uLTRA can be used as an **end-to-end aligner for annotation and detection of novel genes or isoforms** (default mode). This is because uLTRA (>=v0.0.2) now incorporates [minimap2](https://github.com/lh3/minimap2). [minimap2](https://github.com/lh3/minimap2) is run upon start of uLTRA, and the results are used both for (i) not aligning reads with uLTRA which had a primary alignment to regions not indexed by uLTRA (e.g. genomic regions or unannotated genes) and (ii) to consult at end of program which aligner had a better fit (based on cigar) of the primary alignment and chose this alignment to be primary. uLTRA still uses its own alignment algorithm to align to and around all annotated gene regions. uLTRA can therefore, at worst, be seen as an advanced wrapper around minimap2 that refines alignments around annotated regions. See updated `CREDITS` when using this version. uLTRA can still be used as a stand alone aligner as presented in our [paper](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btab540/6327681) by specifying `--disable_mm2`.




