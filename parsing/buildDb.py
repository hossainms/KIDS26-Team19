"""Read all GEO series matrix files in Toy-Datasets and load them into an in-memory DuckDB.

For each *_series_matrix.txt.gz file, this extracts a small set of fields into a
`samples` table: Series_geo_accession, Series_platform_id, Sample_geo_accession,
Sample_organism_ch1, and Sample_data_row_count.
"""

import csv
import glob
import gzip
import os

import duckdb

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "Toy-Datasets")
TABLE_BEGIN_MARKER = "!series_matrix_table_begin"

WANTED_KEYS = (
    "Series_geo_accession",
    "Series_platform_id",
    "Sample_geo_accession",
    "Sample_organism_ch1",
    "Sample_data_row_count",
)


def parse_tsv_line(line):
    return next(csv.reader([line], delimiter="\t", quotechar='"'))


def parse_series_matrix(path):
    fields = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith(TABLE_BEGIN_MARKER):
                break  # stop before the (large) expression matrix
            if not line.strip().startswith("!"):
                continue
            key, *values = parse_tsv_line(line)
            key = key.lstrip("!")
            if key in WANTED_KEYS:
                fields[key] = values

    series_accession = fields["Series_geo_accession"][0]
    platform_id = fields["Series_platform_id"][0]
    sample_ids = fields["Sample_geo_accession"]
    organisms = fields["Sample_organism_ch1"]
    row_counts = fields["Sample_data_row_count"]

    return [
        (series_accession, platform_id, sample_id, organism, int(row_count))
        for sample_id, organism, row_count in zip(sample_ids, organisms, row_counts)
    ]


def build_db():
    con = duckdb.connect()
    con.execute(
        "CREATE TABLE samples ("
        "series_accession VARCHAR, series_platform_id VARCHAR, "
        "sample_geo_accession VARCHAR, sample_organism_ch1 VARCHAR, "
        "sample_data_row_count INTEGER)"
    )

    paths = sorted(glob.glob(os.path.join(DATA_DIR, "*_series_matrix.txt.gz")))
    for path in paths:
        rows = parse_series_matrix(path)
        con.executemany("INSERT INTO samples VALUES (?, ?, ?, ?, ?)", rows)
        print(f"Loaded {rows[0][0]}: {len(rows)} samples")

    return con


if __name__ == "__main__":
    con = build_db()
    print(con.execute("SELECT COUNT(*) FROM samples").fetchone())

