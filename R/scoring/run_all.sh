#!/usr/bin/env bash
# Rebuild every result from scratch. Expects R with `survival` and `data.table`.
# Network access is required: expression data is pulled from NCBI GEO and the NCI GDC.
set -euo pipefail
cd "$(dirname "$0")"

echo "== unit tests =="
Rscript -e 'testthat::test_dir("tests/testthat")'

echo "== corpus composition (uses cached headers under bronze/headers) =="
Rscript R/scan_corpus.R

echo "== score the corpus (resumable; skips studies already done) =="
Rscript R/score_corpus.R

echo "== permutation test, 95 contrasts =="
Rscript R/null_corpus_test.R 1000

echo "== figures =="
Rscript R/figure_corpus_ranking.R
Rscript R/figure_null_tests.R

echo "done. See results/"
