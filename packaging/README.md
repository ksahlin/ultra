# Distribution

Four prebuilt binaries and a conda recipe. `build-release.sh` produces the
binaries; nothing here tags, pushes or publishes.

```bash
packaging/build-release.sh          # -> packaging/dist/*.tar.gz + SHA256SUMS
```

| target | tarball | glibc floor |
| --- | --- | --- |
| linux x86_64 | `uLTRA-<v>-linux-x86_64.tar.gz` | 2.17 |
| linux aarch64 | `uLTRA-<v>-linux-aarch64.tar.gz` | 2.17 |
| macOS arm64 | `uLTRA-<v>-macos-arm64.tar.gz` | — |
| macOS x86_64 | `uLTRA-<v>-macos-x86_64.tar.gz` | — |

Each is a single file of 650–870 KB. namfinder, edlib and zlib are compiled in
and linked statically; the Linux binaries need only `libc`, `libpthread` and
`libdl`, and the macOS ones only `libSystem` and `libc++`. **minimap2 is the
only runtime dependency**, and only for the default path — `--disable_mm2`
needs nothing at all.

## Installing

```bash
# conda, the recommended route -- conda does not set the quarantine attribute
conda install -c bioconda -c conda-forge ultra_bioinformatics

# or a binary, piped straight out of curl so nothing is ever quarantined
curl -L https://github.com/ksahlin/ultra/releases/download/vX.Y.Z/uLTRA-X.Y.Z-macos-arm64.tar.gz | tar xz
```

### macOS will kill a downloaded binary

This is not a hypothetical. See PORTING.md *Finding 42*: a binary that carries
`com.apple.quarantine` is **SIGKILLed on sight** — exit 137, no message, no
dialog. Measured:

| how it arrives | quarantined? | runs? |
| --- | --- | --- |
| `curl … \| tar xz` | no | **yes** |
| `.tar.gz`, extracted with command-line `tar` | no | **yes** |
| `.zip`, double-clicked (Archive Utility) | **yes** | **no — killed** |
| the bare binary, downloaded in a browser | **yes** | **no — killed** |

Hence: **releases ship `.tar.gz`, never `.zip`, and never a bare binary.**
Anyone who ends up with a quarantined copy anyway can clear it:

```bash
xattr -d com.apple.quarantine ./uLTRA
```

The real fix is Apple notarisation, which needs a paid Apple Developer
account and signing credentials. Until then conda is the route that works
without the user having to know any of this.

## What is verified, and what is not

Verified on this machine:

- all four binaries build, and the Linux pair's highest required glibc symbol
  is `GLIBC_2.17` — checked with `objdump -T`, not assumed from the triple
- both macOS binaries run from an **empty environment** (`env -i`, `PATH`
  reduced to `/usr/bin:/bin`), x86_64 through Rosetta
- `reads.sam` is **byte-identical across architectures** — macOS arm64 vs
  x86_64, on the smoke corpus and on 10 000 real SIRV reads, which is what
  makes a single expected hash a valid cross-platform assertion
- the quarantine behaviour in the table above, each row measured

**Not** verified here, because this is a macOS machine with no Linux and no
container runtime: that the Linux binaries *execute*. They are cross-compiled
and statically checked only. `.github/workflows/build.yml` is what closes
that gap — it runs each Linux binary inside a `manylinux2014` image, whose
glibc is exactly the 2.17 floor being claimed, and asserts the same hash. That
workflow has not run yet either; it lands with the first push.

This distinction is PORTING.md's rule 1 and the reason it exists: a wheel
tested under an interpreter that already had the dependencies, and a release
binary fetched with `gh release download`, which does not set the quarantine
attribute, both passed and both were broken for real users.

## Before the recipe can be submitted

`packaging/conda/meta.yaml` is complete except for two things that are the
author's to decide:

1. **There is no licence file in the repository.** `README.md` says "GPL
   v3.0, see LICENCE.txt" and links to it on `master`; that link is a 404 and
   no such file has ever been committed. bioconda requires `license_file`, so
   the GPL-3.0 text has to be added before submission.
2. **Whether this replaces `ultra_bioinformatics` or becomes a new package.**
   The existing one is the Python tool, with parasail-python, pysam, dill,
   gffutils and intervaltree at runtime and no osx-arm64 build at all. This
   recipe's only runtime dependency is minimap2 and it builds on all four
   subdirs.
