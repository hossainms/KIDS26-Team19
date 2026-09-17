"""Create an empty GEO DuckDB with the fixed diagnosis/dataset schema.

    python parsing/initDb.py --db-path data/geo.duckdb --diagnosis aml

Run once. Use parsing/updateDb.py to load series matrix files into it.
"""

from __future__ import annotations

import argparse
import os
import sys

import duckdb

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from schema import create_schema, upsert_diagnosis  # noqa: E402

DEFAULT_DB_PATH = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "data", "geo.duckdb"
)


def init_db(db_path: str, diagnoses=(), overwrite: bool = False) -> str:
    db_path = os.path.abspath(db_path)
    if os.path.exists(db_path):
        if not overwrite:
            raise FileExistsError(
                f"{db_path} already exists; pass --overwrite to recreate it"
            )
        os.remove(db_path)

    parent = os.path.dirname(db_path)
    if parent:
        os.makedirs(parent, exist_ok=True)

    con = duckdb.connect(db_path)
    try:
        create_schema(con)
        for name in diagnoses:
            upsert_diagnosis(con, name)
    finally:
        con.close()
    return db_path


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--db-path", default=DEFAULT_DB_PATH)
    parser.add_argument(
        "--diagnosis",
        action="append",
        default=[],
        metavar="NAME",
        help="Seed a diagnosis row; repeatable (e.g. --diagnosis aml).",
    )
    parser.add_argument(
        "--overwrite", action="store_true", help="Delete an existing database first."
    )
    args = parser.parse_args(argv)

    try:
        path = init_db(args.db_path, args.diagnosis, overwrite=args.overwrite)
    except FileExistsError as error:
        print(f"error: {error}", file=sys.stderr)
        return 1

    seeded = ", ".join(d.strip().lower() for d in args.diagnosis) or "none"
    print(f"Initialized {path} (tables: diagnosis, dataset; diagnoses seeded: {seeded})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
