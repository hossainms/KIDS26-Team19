"""Shared schema and streaming header parser for the GEO DuckDB.

Two scripts use this module:
  parsing/initDb.py    creates an empty database with the fixed tables
  parsing/updateDb.py  adds series matrix files to an existing database

Tables: diagnosis -> dataset (one row per series matrix file) -> sample (one row
per GEO sample). `sample` carries diagnosis_id as well so it can be filtered by
diagnosis without joining through dataset.
"""

from __future__ import annotations

import csv
import gzip
import io
import os
import re

TABLE_BEGIN_MARKER = "!series_matrix_table_begin"

# Persistent, repo-tracked database. Kept separate from data/geo.duckdb, which
# the legacy buildDb.py deletes and rebuilds on every run.
DEFAULT_DB_PATH = os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "geo.db")
)

WANTED_KEYS = (
    "Series_geo_accession",
    "Series_platform_id",
    "Sample_geo_accession",
    "Sample_organism_ch1",
    "Sample_data_row_count",
    "Sample_characteristics_ch1",
    "Sample_molecule_ch1",
    "Series_pubmed_id",
)

SERIES_KEYS = tuple(k for k in WANTED_KEYS if k.startswith("Series_"))
SAMPLE_KEYS = tuple(k for k in WANTED_KEYS if k.startswith("Sample_"))

# GEO omits !Series_pubmed_id until a study is published, so it is only
# required under --strict.
OPTIONAL_KEYS = ("Series_pubmed_id",)

KEY_COLUMNS = {key: key.lower() for key in WANTED_KEYS}
INTEGER_COLUMNS = frozenset({"sample_data_row_count"})

DATASET_COLUMNS = ("diagnosis_id", "source_file", *(KEY_COLUMNS[k] for k in SERIES_KEYS))
SAMPLE_COLUMNS = ("dataset_id", "diagnosis_id", *(KEY_COLUMNS[k] for k in SAMPLE_KEYS))


def _column_type(column: str) -> str:
    return "BIGINT" if column in INTEGER_COLUMNS else "VARCHAR"


_DDL = f"""
CREATE SEQUENCE IF NOT EXISTS diagnosis_id_seq START 1;
CREATE SEQUENCE IF NOT EXISTS dataset_id_seq START 1;
CREATE SEQUENCE IF NOT EXISTS sample_id_seq START 1;

CREATE TABLE IF NOT EXISTS diagnosis (
    diagnosis_id BIGINT PRIMARY KEY DEFAULT nextval('diagnosis_id_seq'),
    diagnosis_name VARCHAR NOT NULL UNIQUE
);

CREATE TABLE IF NOT EXISTS dataset (
    dataset_id BIGINT PRIMARY KEY DEFAULT nextval('dataset_id_seq'),
    diagnosis_id BIGINT NOT NULL REFERENCES diagnosis(diagnosis_id),
    source_file VARCHAR NOT NULL UNIQUE,
    {", ".join(
        f"{KEY_COLUMNS[key]} {_column_type(KEY_COLUMNS[key])}" for key in SERIES_KEYS
    )}
);

CREATE TABLE IF NOT EXISTS sample (
    sample_id BIGINT PRIMARY KEY DEFAULT nextval('sample_id_seq'),
    dataset_id BIGINT NOT NULL REFERENCES dataset(dataset_id),
    diagnosis_id BIGINT NOT NULL REFERENCES diagnosis(diagnosis_id),
    {", ".join(
        f"{KEY_COLUMNS[key]} {_column_type(KEY_COLUMNS[key])}" for key in SAMPLE_KEYS
    )}
);
"""

INSERT_DATASET_SQL = (
    "INSERT INTO dataset (dataset_id, {columns}) "
    "VALUES (nextval('dataset_id_seq'), {placeholders}) RETURNING dataset_id"
).format(
    columns=", ".join(DATASET_COLUMNS),
    placeholders=", ".join("?" for _ in DATASET_COLUMNS),
)

INSERT_SAMPLE_SQL = (
    "INSERT INTO sample (sample_id, {columns}) "
    "VALUES (nextval('sample_id_seq'), {placeholders})"
).format(
    columns=", ".join(SAMPLE_COLUMNS),
    placeholders=", ".join("?" for _ in SAMPLE_COLUMNS),
)


class MatrixError(Exception):
    """A series matrix file cannot be loaded."""


def create_schema(con) -> None:
    con.execute(_DDL)


def upsert_diagnosis(con, name: str) -> int:
    """Return the diagnosis_id for `name`, creating the row if needed."""
    normalized = name.strip().lower()
    if not normalized:
        raise ValueError("diagnosis name must not be empty")
    con.execute(
        "INSERT INTO diagnosis (diagnosis_id, diagnosis_name) "
        "VALUES (nextval('diagnosis_id_seq'), ?) ON CONFLICT DO NOTHING",
        [normalized],
    )
    return con.execute(
        "SELECT diagnosis_id FROM diagnosis WHERE diagnosis_name = ?", [normalized]
    ).fetchone()[0]


def open_matrix(path: str) -> io.TextIOBase:
    opener = gzip.open if path.endswith(".gz") else open
    return opener(path, "rt", encoding="utf-8", errors="replace")


def read_header(path: str) -> dict[str, list[str]]:
    """Read only the `!`-prefixed header, stopping before the expression matrix."""
    fields: dict[str, list[str]] = {}
    with open_matrix(path) as handle:
        for line in handle:
            if line.startswith(TABLE_BEGIN_MARKER):
                break
            if not line.startswith("!"):
                continue
            key, *values = next(csv.reader([line], delimiter="\t", quotechar='"'))
            key = key.lstrip("!")
            if key in KEY_COLUMNS:
                fields[key] = values
    return fields


def missing_keys(fields: dict[str, list[str]], strict: bool = False) -> list[str]:
    required = (
        WANTED_KEYS if strict else tuple(k for k in WANTED_KEYS if k not in OPTIONAL_KEYS)
    )
    return [key for key in required if not fields.get(key)]


def _scalar(fields: dict[str, list[str]], key: str) -> str | None:
    values = fields.get(key)
    if not values or not str(values[0]).strip():
        return None
    return values[0]


def _as_int(value) -> int | None:
    if value is None:
        return None
    match = re.search(r"-?\d+", str(value))
    return int(match.group()) if match else None


def _cast(column: str, value):
    return _as_int(value) if column in INTEGER_COLUMNS else value


def parse_matrix(path: str, diagnosis_id: int, strict: bool = False):
    """Return (dataset_row, sample_rows) for one series matrix file.

    Raises MatrixError when required keys are missing so the caller can report
    the file instead of writing a partial study.
    """
    fields = read_header(path)
    absent = missing_keys(fields, strict=strict)
    if absent:
        raise MatrixError("missing required keys: " + ", ".join(absent))

    # Platform-specific matrix filenames are authoritative when present. Some
    # GEO headers repeat the wrong or shared Series_platform_id across files in
    # a multi-GPL study, while the filename remains platform-specific.
    filename_platform = re.search(r"(?:^|-)GPL([0-9]+)(?:_series_matrix|$)", os.path.basename(path))
    series_platform = (
        f"GPL{filename_platform.group(1)}"
        if filename_platform
        else _scalar(fields, "Series_platform_id")
    )
    dataset_row = (
        diagnosis_id,
        os.path.basename(path),
        *(_scalar(fields, key) for key in SERIES_KEYS),
    )
    platform_index = SERIES_KEYS.index("Series_platform_id")
    dataset_row = dataset_row[:2 + platform_index] + (series_platform,) + dataset_row[3 + platform_index:]

    n_samples = len(fields["Sample_geo_accession"])

    def per_sample(key: str) -> list:
        values = fields.get(key) or []
        if len(values) < n_samples:
            values = values + [None] * (n_samples - len(values))
        return values[:n_samples]

    columns = {key: per_sample(key) for key in SAMPLE_KEYS}
    sample_rows = [
        tuple(_cast(KEY_COLUMNS[key], columns[key][i]) for key in SAMPLE_KEYS)
        for i in range(n_samples)
    ]
    return dataset_row, sample_rows


def series_accession(path: str) -> str | None:
    return _scalar(read_header(path), "Series_geo_accession")
