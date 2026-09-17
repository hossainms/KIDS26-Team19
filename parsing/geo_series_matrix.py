"""Extract GEO series matrix header fields via DuckDB (direct .txt query)."""

from __future__ import annotations

from pathlib import Path

import duckdb
import pandas as pd

SERIES_HEADER_SQL = """
WITH lines AS (
  SELECT line
  FROM read_csv(
    ?,
    columns = {'line': 'VARCHAR'},
    delim = '',
    quote = '',
    escape = '',
    header = false
  )
),
series_lines AS (
  SELECT line
  FROM lines
  WHERE line LIKE '!Series_%'
),
parsed AS (
  SELECT
    split_part(line, chr(9), 1) AS field,
    trim(both '"' from split_part(line, chr(9), 2)) AS value
  FROM series_lines
  WHERE strpos(line, chr(9)) > 0
)
SELECT
  max(CASE WHEN field = '!Series_geo_accession' THEN value END) AS series_accession,
  max(CASE WHEN field = '!Series_title' THEN value END) AS series_title
FROM parsed
"""


def series_metadata_from_matrix(path: str | Path) -> pd.DataFrame:
    """Query a GEO series_matrix .txt file and return accession + title as a DataFrame."""
    resolved = str(Path(path).resolve())
    con = duckdb.connect()
    return con.execute(SERIES_HEADER_SQL, [resolved]).df()
