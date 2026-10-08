uLTRA
===========
[![build](https://github.com/ksahlin/ultra/actions/workflows/build.yml/badge.svg)](https://github.com/ksahlin/ultra/actions/workflows/build.yml) [![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/ultra_bioinformatics/README.html)

uLTRA is a tool for splice alignment of long transcriptomic reads to a genome, guided by a database
of exon annotations. It is particularly accurate when aligning to small exons
[(examples)](data/images). See the
[paper](https://doi.org/10.1093/bioinformatics/btab540), or this
[YouTube video](https://www.youtube.com/watch?v=M7cK80kXXMU).

## v0.3.0 — uLTRA is now written in Rust

Same command line and output files, so existing pipelines do not need editing. A minimap2 binary is
the only dependency. **1.4–2.9x faster** depending on the annotation, and peak memory no longer
grows with `--t` — 3.8 GB against the Python implementation's 10.4 GB on 200 000 Drosophila reads
at `--t 8`. Alignment accuracy is slightly improved. Every input read now appears in
`reads.sam` — the Python implementation silently discards the ones minimap2 could not place.

The Python implementation is still here ([INSTALL-python.md](INSTALL-python.md)), kept as the
reference the Rust version is verified against. More detailed stats in
**[RUST-PORT.md](RUST-PORT.md)**.

Install
-------

A prebuilt binary — needs only glibc 2.17 (RHEL/CentOS 7, 2012).
[Releases](https://github.com/ksahlin/ultra/releases) also has `linux-aarch64`, `macos-arm64` and
`macos-x86_64`:

```bash
curl -L https://github.com/ksahlin/ultra/releases/download/v0.3.0/uLTRA-0.3.0-linux-x86_64.tar.gz | tar xz
```

Or from source, with a Rust toolchain (1.74+):

```bash
git clone https://github.com/ksahlin/ultra.git
cd ultra
cargo build --release --manifest-path rust/Cargo.toml
```

That builds `rust/target/release/uLTRA`. Copy it onto your `PATH`.

[minimap2](https://github.com/lh3/minimap2) must also be on your `PATH`; it is used to detect reads
aligning outside the annotated regions. `--disable_mm2` skips that step and needs nothing at all.

> The bioconda package still installs the Python implementation. The update to this version is
> [in review](https://github.com/bioconda/bioconda-recipes/pull/70047).

#### Test the installation

```bash
uLTRA pipeline test/SIRV_genes.fasta test/SIRV_genes_C_170612a.gtf test/reads.fa /tmp/ultra_test
```

This writes four aligned reads to `/tmp/ultra_test/reads.sam`.

Usage
-----

uLTRA works with PacBio Iso-Seq and ONT cDNA/dRNA reads.

```bash
# everything in one command
uLTRA pipeline genome.fasta /full/path/to/annotation.gtf reads.fa outfolder/

# or index once and align many times
uLTRA index genome.fasta /full/path/to/annotation.gtf outfolder/
uLTRA align genome.fasta reads.fq outfolder/ --ont --t 8      # ONT cDNA
uLTRA align genome.fasta reads.fq outfolder/ --isoseq --t 8   # PacBio Iso-Seq
```

Alignments are written to `outfolder/reads.sam`.

| parameter | |
| --- | --- |
| `--t` | threads (default 3) |
| `--index PATH` | read the index from somewhere other than `outfolder/` |
| `--prefix NAME` | write `outfolder/NAME.sam` instead of `reads.sam` |
| `--disable_infer` | much faster indexing, if your GTF has `gene` and `transcript` features |
| `--disable_mm2` | skip the minimap2 step entirely |

`uLTRA --help` lists the rest.

#### A properly formatted GTF file

This is the most common cause of failure. If you have GFF or another format, convert it with
[AGAT](https://github.com/NBISweden/AGAT) — many other converters do not respect the GTF format:

```bash
agat_convert_sp_gff2gtf.pl --gff annot.gff3 --gtf annot.gtf
```

Credits
-------

Please cite:

1. Kristoffer Sahlin, Veli Mäkinen, Accurate spliced alignment of long RNA sequencing reads,
   *Bioinformatics*, Volume 37, Issue 24, 15 December 2021, Pages 4643–4651,
   https://doi.org/10.1093/bioinformatics/btab540

**Please also cite** [minimap2](https://github.com/lh3/minimap2), which uLTRA uses to align genomic
reads outside the indexed regions. For example: "We aligned reads to the genome using uLTRA [1],
which incorporates minimap2 [CIT]."

Licence
-------

GPL v3.0, see [LICENSE.txt](LICENSE.txt).
