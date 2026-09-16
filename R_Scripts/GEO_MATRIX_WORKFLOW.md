# One-script GEO workflow

Use **05_Download_Metadata_Inventory.R**. It contains keyword search, matrix-file
discovery, downloads, GPL resolution, and reporting. It replaces scripts 03 and 04
and does not source other project scripts.

Packages: `xml2`; `readxl` for Excel; `GEOquery` only when headers cannot resolve
the GPL (`BiocManager::install("GEOquery")`).

Run from the repository root:

```sh
Rscript R_Scripts/05_Download_Metadata_Inventory.R
# Optional input workbook/output directory:
Rscript R_Scripts/05_Download_Metadata_Inventory.R GSE_Metadata_Inventory.xlsx downloads/geo_metadata_inventory
# Or keyword search followed by downloading:
Rscript R_Scripts/05_Download_Metadata_Inventory.R --query '(AML OR "acute myeloid leukemia") AND (DMSO OR vehicle)' --keywords 'AML,acute myeloid leukemia,DMSO,vehicle' --out downloads/geo_search
```

Workbook input uses the Metadata sheet and SeriesAccession column. Search mode
supports `--max-studies 10` for a bounded trial; truncation is recorded.

**downloaded_matrices.csv** contains one row per validated matrix: study accession,
GPL IDs in Platform, filename, URL, local path, size, title, organism, sample count,
platform source, and errors. Local headers supply actual sample counts and study
metadata where available. Workbook annotations and associated search terms are
retained. PubMed IDs come from workbook annotations when available. Missing
metadata is left blank, never invented.

GPL resolution keeps known IDs, otherwise reads local sample/series headers,
then uses GEOquery to match matrix GSMs to platforms. Multi-GPL studies retain
separate files. Ambiguous mappings stay unresolved with PlatformError.
A valid downloaded matrix does not necessarily contain expression values.

Other outputs:

- all_results.csv: combined report including errors and prior results.
- failures.csv: discovery/download failures and invalid GSE inputs.
- platform_failures.csv: downloaded files with unresolved GPL IDs.
- summary.csv: file counts, platform completeness, bytes, and current-run progress.
- matrices/: downloaded compressed files.
- search_*.csv: search provenance when using keyword mode.

Reports are checkpointed after every file and replaced atomically. Reruns merge
existing rows by study/filename, preserving unprocessed results and reusing valid
local files. Use a new output directory for a separate report. Rows already lost
through older report overwrites are recovered as their studies are processed
again. Missing local files are flagged and excluded from the successful report.

Sourcing the script defines functions without starting downloads:

```r
source("R_Scripts/05_Download_Metadata_Inventory.R")
path <- "downloads/geo_metadata_inventory/downloaded_matrices.csv"
report <- fill_geo_platforms(read.csv(path, stringsAsFactors = FALSE))
write.csv(report, path, row.names = FALSE, na = "")
```

Offline tests:

```sh
Rscript R_Scripts/tests/test_geo_matrix_download.R
Rscript R_Scripts/tests/test_geo_keyword_discovery.R
Rscript R_Scripts/tests/test_geo_platform_resolution.R
Rscript R_Scripts/tests/test_geo_pipeline.R
```

The opt-in test_geo_matrix_live.R test makes network requests and downloads files.
Reference: [NCBI programmatic access](https://www.ncbi.nlm.nih.gov/geo/info/geo_paccess.html#FTP).
