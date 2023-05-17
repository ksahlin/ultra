uLTRA
===========
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/ultra_bioinformatics/README.html) [![Build Status](https://travis-ci.org/ksahlin/uLTRA.svg?branch=master)](https://travis-ci.org/ksahlin/uLTRA)


uLTRA is a tool for splice alignment of long transcriptomic reads to a genome, guided by a database of exon annotations. uLTRA is particularly accurate when aligning to small exons [see some examples](https://github.com/ksahlin/ultra/tree/master/data/images). 

uLTRA is distributed as a python package supported on Linux / OSX with python (versions 3.4 or above). 


Table of Contents
=================

  * [INSTALLATION](#INSTALLATION)
  * [USAGE](#USAGE)
    * [Indexing](#Indexing)
    * [aligning](#Aligning)
    * [Output](#Output)
    * [Parameters](#Parameters)
  * [CREDITS](#CREDITS)
  * [LICENCE](#LICENCE)



INSTALLATION
=================

## Conda recipe

There is a [bioconda recipe](https://bioconda.github.io/recipes/ultra_bioinformatics/README.html), [docker image](https://quay.io/repository/biocontainers/ultra_bioinformatics?tab=tags), and a [singularity container](https://depot.galaxyproject.org/singularity/ultra_bioinformatics%3A0.0.4--pyh5e36f6f_1) of uLTRA v0.0.4 created by [sguizard](https://github.com/sguizard). You can use, e.g., the bioconda recipe for an easy automated installation. 

If a newer version of uLTRA is not available through bioconda (or you simply want more control of/customize your installation), alternative ways of installations are provided below. Current version of uLTRA is 0.1 (see changelog at end of this readme).

## Using the INSTALL.sh script

You can clone this repository and 
run the script `INSTALL.sh` as

```
git clone https://github.com/ksahlin/uLTRA.git --depth 1
cd uLTRA
./INSTALL.sh <install_directory>
```

The install script is tested in bash environment. 

To run uLTRA, you need to activate the conda environment "ultra":

```
conda activate ultra
```

## Without the INSTALL.sh script

You can also manually perform below steps for more control.

#### 1. Create conda environment

Create a conda environment called ultra and activate it

```
conda create -n ultra python=3 pip 
conda activate ultra
```

#### 2. Install uLTRA 

```
pip install ultra-bioinformatics
```

#### 3. Install third party tools 

Install [namfinder](https://github.com/ksahlin/namfinder) and [minimap2](https://github.com/lh3/minimap2) and
place the generated binaries `namfinder` and `minimap2` in your path. 

#### 4. Verify installation

You should now have 'uLTRA' installed; try it

```
uLTRA --help
```

Upon start/login to your server/computer you need to activate the conda environment "ultra" to run uLTRA as:
```
conda activate ultra
```

You can also download and use test data available in this repository [here](https://github.com/ksahlin/ultra/tree/master/test) and run: 

```
uLTRA pipeline [/your/full/path/to/test]/SIRV_genes.fasta  \
               /your/full/path/to/test/SIRV_genes_C_170612a.gtf  \
               [/your/full/path/to/test]/reads.fa outfolder/  [optional parameters]
```



## Entirly from source


Make sure the below-listed dependencies are installed (installation links below). All below dependencies except `namfinder` can be installed as `pip install X` or through conda.
* [parasail](https://github.com/jeffdaily/parasail-python)
* [edlib](https://github.com/Martinsos/edlib)
* [pysam](http://pysam.readthedocs.io/en/latest/installation.html) (>= v0.11)
* [dill](https://pypi.org/project/dill/)
* [intervaltree](https://github.com/chaimleib/intervaltree/tree/master/intervaltree)
* [gffutils](https://pythonhosted.org/gffutils/)
* [namfinder](https://github.com/ksahlin/namfinder)

With these dependencies installed. Run

```sh
git clone https://github.com/ksahlin/uLTRA.git
cd uLTRA
./uLTRA
```


USAGE
----------------

uLTRA can be used with either PacBio Iso-Seq or ONT cDNA/dRNA reads. 

### Indexing

```
uLTRA index genome.fasta  /full/path/to/annotation.gtf  outfolder/  [parameters]
```

Important parameters: 

1. `--disable_infer` can speed up the indexing considerably, but it only works if you have the `gene feature` and `transcript feature` in your GTF file.

### Aligning

For example

```
uLTRA align genome.fasta reads.[fa/fq] outfolder/  --ont --t 8   # ONT cDNA reads using 8 cores
uLTRA align genome.fasta reads.[fa/fq] outfolder/  --isoseq --t 8 # PacBio isoseq reads
```

Important parameters:

1. `--index [PATH]`: You can set a custom location of where to get the index from using, otherwise, uLTRA will try to read the index from the `outfolder/` by default. 
2. `--prefix [PREFIX OF FILE]`: The aligned reads will be written to `outfolder/reads.sam` unless `--prefix` is set. For example, `--prefix sample_X` will output the reads in `outfolder/sample_X.sam`.

### Pipeline

Perform all the steps in one

```
uLTRA pipeline genome.fasta /full/path/to/annotation.gtf reads.fa outfolder/  [parameters]
```

### Common errors

Not having a properly formatted GTF file. Before running uLTRA, notice that it reqires a _properly formatted GTF file_. If you have a GFF file or other annotation format, it is adviced to use [AGAT](https://github.com/NBISweden/AGAT) for file conversion to GTF as many other conversion tools do not respect GTF format. For example, you can run AGAT as:

```
agat_convert_sp_gff2gtf.pl --gff annot.gff3 --gtf annot.gtf
```



CREDITS
----------------

Please cite 

1. Kristoffer Sahlin, Veli Mäkinen, Accurate spliced alignment of long RNA sequencing reads, Bioinformatics, Volume 37, Issue 24, 15 December 2021, Pages 4643–4651, https://doi.org/10.1093/bioinformatics/btab540

when using uLTRA. **Please also cite** [minimap2](https://github.com/lh3/minimap2) as uLTRA incorporates minimap2 for alignment of some genomic reads outside indexed regions. For example "We aligned reads to the genome using uLTRA [1], which incorporates minimap2 [CIT].".




LICENCE
----------------

GPL v3.0, see [LICENSE.txt](https://github.com/ksahlin/uLTRA/blob/master/LICENCE.txt).




