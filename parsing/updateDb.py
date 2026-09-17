"""Add GEO series matrix files to an existing DuckDB built by parsing/initDb.py.

    python parsing/updateDb.py downloads/geo_aml/matrices \
        --db-path data/geo.duckdb --diagnosis aml --report logs/update_aml.txt

Each file is validated for the required header keys and skipped if its filename
or series accession is already present. Skips and failures are written to the
report (stdout by default); the exit code is 1 when any file is rejected.
"""

from __future__ import annotations

import argparse
import contextlib
import os
import sys

import duckdb

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from schema import (  # noqa: E402
    INSERT_DATASET_SQL,
    INSERT_SAMPLE_SQL,
    MatrixError,
    parse_matrix,
    upsert_diagnosis,
)

MATRIX_SUFFIXES = ("_series_matrix.txt.gz", "_series_matrix.txt")


def iter_matrix_files(data_dir: str, recursive: bool = False):
    if not os.path.isdir(data_dir):
        raise NotADirectoryError(data_dir)
    if recursive:
        for root, _dirs, files in os.walk(data_dir):
            for name in sorted(files):
                if name.endswith(MATRIX_SUFFIXES):
                    yield os.path.join(root, name)
        return
    for name in sorted(os.listdir(data_dir)):
        path = os.path.join(data_dir, name)
        if name.endswith(MATRIX_SUFFIXES) and os.path.isfile(path):
            yield path


def _existing_dataset(con, source_file: str, accession: str | None):
    return con.execute(
        "SELECT source_file FROM dataset "
        "WHERE source_file = ? OR (? IS NOT NULL AND series_geo_accession = ?) LIMIT 1",
        [source_file, accession, accession],
    ).fetchone()


def load_file(con, path: str, diagnosis_id: int, strict: bool = False) -> int:
    """Insert one file's dataset row and sample rows; return the sample count."""
    dataset_row, sample_rows = parse_matrix(path, diagnosis_id, strict=strict)
    accession = dataset_row[2]

    duplicate = _existing_dataset(con, dataset_row[1], accession)
    if duplicate:
        raise MatrixError(f"already in the database (as {duplicate[0]})")
    if not sample_rows:
        raise MatrixError("no samples found in header")

    con.execute("BEGIN TRANSACTION")
    try:
        dataset_id = con.execute(INSERT_DATASET_SQL, list(dataset_row)).fetchone()[0]
        con.executemany(
            INSERT_SAMPLE_SQL,
            [(dataset_id, diagnosis_id, *row) for row in sample_rows],
        )
        con.execute("COMMIT")
    except Exception:
        con.execute("ROLLBACK")
        raise
    return len(sample_rows)


@contextlib.contextmanager
def open_report(path: str | None):
    if not path or path == "-":
        yield sys.stdout
        return
    parent = os.path.dirname(os.path.abspath(path))
    if parent:
        os.makedirs(parent, exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        yield handle


def update_db(
    data_dir: str,
    db_path: str,
    diagnosis: str,
    report=sys.stdout,
    strict: bool = False,
) -> tuple[int, int]:
    """Load every matrix file in data_dir; return (loaded, rejected) counts."""
    if not os.path.exists(db_path):
        raise FileNotFoundError(f"{db_path} not found; run parsing/initDb.py first")

    loaded = rejected = 0
    con = duckdb.connect(db_path)
    try:
        diagnosis_id = upsert_diagnosis(con, diagnosis)
        for path in iter_matrix_files(data_dir):
            try:
                n_samples = load_file(con, path, diagnosis_id, strict=strict)
            except (MatrixError, OSError, KeyError) as error:
                rejected += 1
                print(f"SKIP\t{os.path.basename(path)}\t{error}", file=report)
            else:
                loaded += 1
                print(f"OK\t{os.path.basename(path)}\t{n_samples} samples", file=report)
    finally:
        con.close()
    return loaded, rejected


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("data_dir", help="Directory of *_series_matrix.txt(.gz) files.")
    parser.add_argument("--db-path", required=True, help="DuckDB file to update.")
    parser.add_argument(
        "--diagnosis", required=True, help="Diagnosis these files belong to, e.g. aml."
    )
    parser.add_argument(
        "--report", default="-", help="Write per-file results here ('-' for stdout)."
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Also require Series_pubmed_id, which GEO omits until publication.",
    )
    args = parser.parse_args(argv)

    try:
        with open_report(args.report) as report:
            loaded, rejected = update_db(
                args.data_dir,
                args.db_path,
                args.diagnosis,
                report=report,
                strict=args.strict,
            )
            print(f"# loaded {loaded} file(s), rejected {rejected}", file=report)
    except (FileNotFoundError, NotADirectoryError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2

    print(f"Loaded {loaded} file(s) into {args.db_path}; rejected {rejected}.")
    return 1 if rejected else 0


if __name__ == "__main__":
    raise SystemExit(main())
