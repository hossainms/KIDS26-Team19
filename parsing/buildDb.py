"""Load GEO sample metadata into a DuckDB file for a diagnosis's downloaded matrices.

Intended to run against the `matrices/` directory produced by
R_Scripts/05_Download_Metadata_Inventory.R for a given diagnosis, e.g.:
    python parsing/buildDb.py downloads/geo_aml/matrices --db-path data/geo.duckdb

When downloaded_matrices.csv (or all_results.csv) exists next to data_dir, it is
loaded into a `studies` table (study-level inventory from run_geo_pipeline()).
The `samples` table is always built by parsing *_series_matrix.txt.gz files under
data_dir for per-sample fields (series_accession, sample_geo_accession, etc.).
"""

import argparse
import csv
import glob
import gzip
import os

import duckdb

REPO_ROOT = os.path.join(os.path.dirname(__file__), "..")
DEFAULT_DB_PATH = os.path.join(REPO_ROOT, "data", "geo.duckdb")
TABLE_BEGIN_MARKER = "!series_matrix_table_begin"
REPORT_FILENAMES = ("downloaded_matrices.csv", "all_results.csv")

WANTED_KEYS = (
    "Series_geo_accession",
    "Series_platform_id",
    "Sample_geo_accession",
    "Sample_organism_ch1",
    "Sample_data_row_count",
    "Sample_characteristics_ch1",
    "Sample_molecule_ch1",
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
    n_samples = len(sample_ids)

    def per_sample(key):
        values = fields.get(key)
        if not values:
            return [None] * n_samples
        if len(values) < n_samples:
            values = values + [None] * (n_samples - len(values))
        return values[:n_samples]

    organisms = per_sample("Sample_organism_ch1")
    row_counts = per_sample("Sample_data_row_count")
    characteristics = per_sample("Sample_characteristics_ch1")
    molecules = per_sample("Sample_molecule_ch1")

    rows = []
    for sample_id, organism, row_count, characteristic, molecule in zip(
        sample_ids, organisms, row_counts, characteristics, molecules
    ):
        parsed_count = None
        if row_count is not None and str(row_count).strip():
            try:
                parsed_count = int(row_count)
            except ValueError:
                parsed_count = None
        rows.append(
            (
                series_accession,
                platform_id,
                sample_id,
                organism,
                parsed_count,
                characteristic,
                molecule,
            )
        )
    return rows


def find_report_path(data_dir):
    """Find run_geo_pipeline()'s inventory report next to a matrices/ directory."""
    parent = os.path.dirname(os.path.normpath(data_dir))
    for name in REPORT_FILENAMES:
        candidate = os.path.join(parent, name)
        if os.path.exists(candidate):
            return candidate
    return None


def load_samples_from_matrices(con, data_dir):
    con.execute(
        "CREATE TABLE samples ("
        "series_accession VARCHAR, series_platform_id VARCHAR, "
        "sample_geo_accession VARCHAR, sample_organism_ch1 VARCHAR, "
        "sample_data_row_count INTEGER, sample_characteristics_ch1 VARCHAR, "
        "sample_molecule_ch1 VARCHAR)"
    )

    paths = sorted(glob.glob(os.path.join(data_dir, "*_series_matrix.txt.gz")))
    for path in paths:
        rows = parse_series_matrix(path)
        con.executemany("INSERT INTO samples VALUES (?, ?, ?, ?, ?, ?, ?)", rows)
        print(f"Loaded {rows[0][0]}: {len(rows)} samples")

    return len(paths)


def build_db(data_dir, db_path=None):
    """Build studies (optional) and samples tables for the matrices in data_dir.

    db_path=None keeps the database in memory; otherwise the file at db_path is
    deleted first so every build starts from a clean database.
    """
    if db_path:
        os.makedirs(os.path.dirname(db_path), exist_ok=True)
        if os.path.exists(db_path):
            os.remove(db_path)
        con = duckdb.connect(db_path)
    else:
        con = duckdb.connect()

    report_path = find_report_path(data_dir)
    if report_path:
        con.execute("CREATE TABLE studies AS SELECT * FROM read_csv_auto(?)", [report_path])
        count = con.execute("SELECT COUNT(*) FROM studies").fetchone()[0]
        print(f"Loaded studies from {report_path}: {count} rows")

    matrix_count = load_samples_from_matrices(con, data_dir)
    sample_count = con.execute("SELECT COUNT(*) FROM samples").fetchone()[0]
    print(f"Parsed {matrix_count} matrix file(s); {sample_count} sample row(s) in samples")

    return con


def main():
    parser = argparse.ArgumentParser(
        description="Build the GEO samples DuckDB from series matrix files."
    )
    parser.add_argument(
        "data_dir",
        help="Directory of *_series_matrix.txt.gz files, e.g. downloads/geo_aml/matrices.",
    )
    parser.add_argument(
        "--db-path", default=DEFAULT_DB_PATH,
        help="DuckDB file to write; deleted and rebuilt from scratch each run.",
    )
    args = parser.parse_args()

    con = build_db(args.data_dir, args.db_path)
    print(con.execute("SELECT COUNT(*) FROM samples").fetchone())
    con.close()


if __name__ == "__main__":
    main()

