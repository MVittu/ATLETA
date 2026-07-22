# Paper

Structure:

- `latex/` — LaTeX source. `main.tex` is the entry point; individual sections live in `latex/sections/`; `references.bib` holds the bibliography.
- `figures/` — Plots and images referenced by the LaTeX source.
- `docs/` — Reference material, notes, and source papers used while writing.
- `data/` — Processed data snapshots used to generate figures/tables (not raw experiment data).

To build: compile `latex/main.tex` (e.g. `latexmk -pdf main.tex` from within `latex/`).
