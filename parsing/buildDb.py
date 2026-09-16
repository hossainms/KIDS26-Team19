"""Read all GEO series matrix files in Toy-Datasets and load them into an in-memory DuckDB.

For each *_series_matrix.txt.gz file, this creates:
  - rows in `series_metadata` (accession, key, value)
  - rows in `sample_metadata` (series_accession, sample_accession, key, value)
  - a wide table `expr_<ACCESSION>` holding the expression matrix (ID_REF + one column per sample)
"""

import csv
import glob
import gzip
import os

import duckdb

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "Toy-Datasets")
TABLE_BEGIN_MARKER = "!series_matrix_table_begin"
TABLE_END_MARKER = "!series_matrix_table_end"


def parse_tsv_line(line):
    return next(csv.reader([line], delimiter="\t", quotechar='"'))


def parse_series_matrix(path):
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        lines = f.readlines()

    begin_idx = next(i for i, line in enumerate(lines) if line.startswith(TABLE_BEGIN_MARKER))
    end_idx = next(i for i, line in enumerate(lines) if line.startswith(TABLE_END_MARKER))

    meta_fields = [
        parse_tsv_line(line)
        for line in lines[:begin_idx]
        if line.strip().startswith("!")
    ]

    accession = next(
        values[0]
        for key, *values in meta_fields
        if key.lstrip("!") == "Series_geo_accession"
    )
    sample_ids = next(
        values
        for key, *values in meta_fields
        if key.lstrip("!") == "Sample_geo_accession"
    )

    series_rows = [
        (accession, key.lstrip("!"), value)
        for key, *values in meta_fields
        if key.lstrip("!").startswith("Series_")
        for value in values
    ]
    sample_rows = [
        (accession, sample_id, key.lstrip("!"), value)
        for key, *values in meta_fields
        if key.lstrip("!").startswith("Sample_")
        for sample_id, value in zip(sample_ids, values)
    ]

    header = parse_tsv_line(lines[begin_idx + 1])
    sample_columns = header[1:]
    data_rows = []
    for line in lines[begin_idx + 2 : end_idx]:
        if not line.strip():
            continue
        fields = parse_tsv_line(line)
        probe_id, values = fields[0], fields[1:]
        data_rows.append(
            [probe_id] + [float(v) if v not in ("", "NA") else None for v in values]
        )

    return accession, series_rows, sample_rows, sample_columns, data_rows


def load_expression_table(con, accession, sample_columns, data_rows):
    table_name = f"expr_{accession}"
    columns_sql = ", ".join(
        ['"ID_REF" VARCHAR'] + [f'"{col}" DOUBLE' for col in sample_columns]
    )
    con.execute(f'CREATE TABLE "{table_name}" ({columns_sql})')
    placeholders = ", ".join(["?"] * (1 + len(sample_columns)))
    con.executemany(f'INSERT INTO "{table_name}" VALUES ({placeholders})', data_rows)
    return table_name


def build_db():
    con = duckdb.connect()
    con.execute("CREATE TABLE series_metadata (accession VARCHAR, key VARCHAR, value VARCHAR)")
    con.execute(
        "CREATE TABLE sample_metadata "
        "(series_accession VARCHAR, sample_accession VARCHAR, key VARCHAR, value VARCHAR)"
    )

    paths = sorted(glob.glob(os.path.join(DATA_DIR, "*_series_matrix.txt.gz")))
    for path in paths:
        accession, series_rows, sample_rows, sample_columns, data_rows = parse_series_matrix(path)
        con.executemany("INSERT INTO series_metadata VALUES (?, ?, ?)", series_rows)
        con.executemany("INSERT INTO sample_metadata VALUES (?, ?, ?, ?)", sample_rows)
        table_name = load_expression_table(con, accession, sample_columns, data_rows)
        print(f"Loaded {accession}: {len(data_rows)} probes x {len(sample_columns)} samples -> {table_name}")

    return con


if __name__ == "__main__":
    con = build_db()
    print(con.execute("SELECT COUNT(*) FROM series_metadata").fetchone())
    print(con.execute("SELECT COUNT(*) FROM sample_metadata").fetchone())

