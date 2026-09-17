"""Create the persistent GEO DuckDB with the fixed diagnosis/dataset/sample schema.

    python parsing/initDb.py --diagnosis aml

Run this once. The resulting data/geo.db is committed to the repo; every later
change goes through parsing/updateDb.py, which opens the same file in place.
"""

from __future__ import annotations

import argparse
import os
import sys

import duckdb

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from schema import DEFAULT_DB_PATH, create_schema, upsert_diagnosis  # noqa: E402


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
    parser.add_argument(
        "--db-path",
        default=DEFAULT_DB_PATH,
        help="Persistent DuckDB file to create (default: %(default)s).",
    )
    parser.add_argument(
        "--diagnosis",
        action="append",
        default=[],
        metavar="NAME",
        help="Seed a diagnosis row; repeatable (e.g. --diagnosis aml).",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Delete the existing database first; this discards all loaded data.",
    )
    args = parser.parse_args(argv)

    try:
        path = init_db(args.db_path, args.diagnosis, overwrite=args.overwrite)
    except FileExistsError as error:
        print(f"error: {error}", file=sys.stderr)
        return 1

    seeded = ", ".join(d.strip().lower() for d in args.diagnosis) or "none"
    print(
        f"Initialized {path} (tables: diagnosis, dataset, sample; "
        f"diagnoses seeded: {seeded})"
    )
    print("Commit this file so the whole team shares the same database.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
