Installing the Python implementation
====================================

uLTRA is now written in Rust — see [README.md](README.md) for that, and
[RUST-PORT.md](RUST-PORT.md) for what changed. The Python implementation is kept as the reference
the Rust version is verified against, and these are its original installation instructions.

It needs parasail-python, pysam, dill, gffutils, intervaltree, edlib, namfinder and minimap2.
Note that **namfinder has no osx-arm64 conda build**, so this route does not work on Apple
Silicon; the Rust version compiles namfinder in and does.

---


## Conda recipe

There is a [bioconda recipe](https://bioconda.github.io/recipes/ultra_bioinformatics/README.html), [docker image](https://quay.io/repository/biocontainers/ultra_bioinformatics?tab=tags), and a [singularity container](https://depot.galaxyproject.org/singularity/ultra_bioinformatics%3A0.0.4--pyh5e36f6f_1) of uLTRA created by [sguizard](https://github.com/sguizard). You can use, e.g., the bioconda recipe for an easy automated installation. 

Alternative ways of installations are provided below.

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


