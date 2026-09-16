# Agent guide for KIDS26-Team19

Public repo for the St. Jude KIDS26 BioHackathon. The team searches GEO for AML cell-line studies with a drug and a DMSO or vehicle control. It scores samples with LSC6 and LSC17 and ranks drugs for lab follow-up. The team works in R.

## Folders

The `ian` branch uses a medallion layout. Each layer is rebuilt from the one before it.

- `bronze/`: raw files exactly as downloaded. Never edit them. `bronze/manifest.csv` is committed and lists every file with URL, size, MD5, and status. `bronze/geo/<GSE>/` holds the series matrix files and is ignored by git.
- `silver/`: parsed, tidy tables, one per concept: studies, samples, expression, QC.
- `gold/`: analysis-ready tables: sample scores and drug rankings.
- `R/`: functions, one file per layer (`bronze.R`), plus `run_*.R` scripts run from the repo root.
- `tests/testthat/`: tests. Run `Rscript -e 'testthat::test_dir("tests/testthat")'` from the repo root.
- `logs/`: one Markdown log per commit or group of commits.
- `R_Scripts/`, `*.xlsx`, `project-management/`, `docs/`: files from team members. Do not change them without asking their author.

## Logs

Every commit gets a log in `logs/`. Name it `YYYY-MM-DD-HHMM-short-topic.md` using Central time (America/Chicago). Start it with this frontmatter:

```
---
date: 2026-09-16
time: "13:41"
time_zone: America/Chicago (CDT)
person: <who asked for the work>
coding_agent: <tool, or "none">
model: <model, or "none">
machine: <host>
branch: <branch>
commits: [<short sha>]
---
```

The body says what changed, what ran, what was found, and what comes next. State counts from commands, not guesses. The tests check the frontmatter.

## Rules

- The repo is public. Commit only public data and never patient-level data.
- Keep raw data out of git. Commit the manifest and the code that rebuilds it.
- Write a failing test first, then the code.
- Ask before installing R packages.
- No attribution, co-author, or "generated with" lines in commits.
- Branching follows `docs/git-github-basics.md`.
