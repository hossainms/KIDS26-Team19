---
date: 2026-09-16
time: "13:41"
time_zone: America/Chicago (CDT)
person: Ian Maurer
coding_agent: Claude Code
model: Claude Opus 5
machine: beelink
branch: ian
commits: [1d52d2e]
---

# Bronze download of 588 studies

## What ran

`Rscript R/run_bronze.R` from the repo root, twice. The second run retried one failure and skipped files already on disk.

## Result

- 777 manifest rows for 588 studies. 774 files downloaded, 1.1 GB, in `bronze/geo/`.
- GSE126933 failed on the first run and succeeded on the retry.
- GDS4288, GDS4290, and GDS6083 have no files. These are GEO curated dataset IDs, not study IDs, so they have no matrix folder. Their parent study IDs still need to be looked up.

## Check against Nobel's inventory

- `GSM_Expression_Inventory (1).xlsx` has 1,332 rows but only 774 unique files. 558 rows are duplicates.
- Those 774 files match the 774 downloaded here exactly.
- The summary sheet's sample total of 51,795 likely counts the duplicates too.

## Next

Silver step: parse each file into study, sample, and expression tables, and flag files with no expression values.
