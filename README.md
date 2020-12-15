uLTRA
===========

uLTRA is a tool for splice alignment of long transcriptomic reads to a genome, guided by a database of exon annotations. uLTRA takes reads in fast(a/q) and a genome annotation as input and outputs a SAM-file. The SAM-file includes information of which splice sites are found and if the read is a full splice match (and to which transcript), incomplete splice match, Novel in catalogue, or novel not in catalogue, as defined in [SQANTI](https://github.com/ConesaLab/SQANTI). uLTRA is highly accurate when aligning to small exons [see some examples](https://github.com/ksahlin/ultra/data/images/).

uLTRA is distributed as a python package supported on Linux / OSX with python v>=3.4. [![Build Status](https://travis-ci.org/ksahlin/uLTRA.svg?branch=master)](https://travis-ci.org/ksahlin/uLTRA).

Table of Contents
=================

  * [INSTALLATION](#INSTALLATION)
    * [Using conda](#Using-conda)
    * [Downloading source from GitHub](#Downloading-source-from-github)
    * [Dependencies](#Dependencies)
  * [USAGE](#USAGE)
    * [Indexing](#Indexing)
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
conda activate ultra
```

2. Install uLTRA 

```
pip install uLTRA
```

3. You should now have 'uLTRA' installed; try it:
```
uLTRA --help
```

4. Install [slaMEM](https://github.com/fjdf/slaMEM)

```
git clone git@github.com:fjdf/slaMEM.git
cd slaMEM
make 
```
And either place the generated binary `slaMEM`in your path or run `export PATH=$PATH:$PWD/` if you are in the slaMEM folder).


Upon start/login to your server/computer you need to activate the conda environment "ultra" to run uLTRA as:
```
source activate ultra
```

### Downloading source from GitHub

#### Dependencies

Make sure the below listed dependencies are installed (installation links below). Versions in parenthesis are suggested as uLTRA has not been tested with earlier versions of these libraries. However, uLTRA may also work with earliear versions of these libaries.
* [parasail](https://github.com/jeffdaily/parasail-python)
* [pysam](http://pysam.readthedocs.io/en/latest/installation.html) (>= v0.11)
* dill
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

uLTRA can be used with either Iso-Seq or ONT reads. 


### Indexing

Fist we construct the datastructures used in uLTRA using a genome annotation GTF file and a genome fasta file.

```
uLTRA prep_splicing  all_genes.gtf outfolder/  [parameters]
```

```
uLTRA prep_seqs  genome.fasta  outfolder/  [parameters]
```


### Aligning

For example, to align ONT cDNA reads using 48 cores, run

```
uLTRA align  genome.fasta  reads.[fa/fq] outfolder/  --ont --t 48   # ONT cDNA reads using 48 cores
uLTRA align  genome.fasta  reads.[fa/fq] outfolder/  --isoseq --t 48 # PacBio isoseq reads
uLTRA align  genome.fasta  reads.[fa/fq] outfolder/  --k 14  --t 48 # PacBio dRNA reads or reads with >10-12% error rate
```


### Pipeline

Perforns all the steps in one

```
uLTRA pipeline test/SIRV_genes_C_170612a.gtf  test/SIRV_genes.fasta  test/reads.fa outfolder/  [parameters]
```

#### Output

uLTRA outputs a SAM-file with alignments to the genome. In addition, it outputs to extra tags describing whether all the splices sites are known and annotated (FSM), new splice combinations (NIC), etc. For details see the defitions of notations in the [Sqanti paper](https://genome.cshlp.org/content/28/7/1096).



CREDITS
----------------

Please cite [1] when using uLTRA.

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


