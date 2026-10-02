#!/usr/bin/env bash
# A real alignment, not just --version: a binary that starts but cannot align
# passes --version happily. Uses the repository's own SIRV fixtures, pulled in
# by test.source_files.
set -euxo pipefail
uLTRA pipeline --disable_mm2 \
    test/SIRV_genes.fasta test/SIRV_genes_C_170612a.gtf test/reads.fa out
test -s out/reads.sam
test "$(grep -vc '^@' out/reads.sam)" = "4"
