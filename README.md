# ATLETA

Reproducible analysis for a 2026 cross-sectional study of carbon-plate insole use and self-reported injury prevalence among Italian competitive athletics athletes.

The main result is an adjusted association, not a causal effect: among competition-period carbon users, higher reported carbon share was associated with higher competition-period injury prevalence. Adoption itself and preparation-period dose showed no clear association. The retrospective design and model sensitivity make the competition finding suggestive rather than confirmatory.

## Main outputs

- [`paper/latex/main.pdf`](paper/latex/main.pdf): compiled manuscript.
- [`paper/latex/`](paper/latex/): reproducible manuscript source.
- [`experiments/multilevel/`](experiments/multilevel/): primary Bayesian multilevel analysis, population weighting, diagnostics, and results.
- [`table/table1/`](table/table1/): descriptive Table 1 pipeline.
- [`experiments/coda/`](experiments/coda/): exploratory compositional analysis retained for reproducibility; it is not part of the final paper's conclusions.

## Data privacy

Individual survey responses, record identifiers, expert-validation exports, and detailed recruitment-monitoring cells are deliberately excluded from Git. The repository contains only analysis code and disclosure-reviewed aggregate outputs. See [`DATA_AVAILABILITY.md`](DATA_AVAILABILITY.md).

The recall scripts write participant-level intermediates to ignored paths under `data/processed/recall/`. Do not change those outputs to a tracked directory.

## Use the fitted posterior without refitting

The expensive Stan sampling results are already committed as compact posterior draws. Load the primary fit directly with base R:

```r
draws <- readRDS("experiments/multilevel/runs/01_primary/association/models/posterior_draws.rds")
names(draws)
```

Saved draws for the baseline, sensitivity, population-weighted, and subgroup analyses are in the neighboring `runs/` directories. The ignored `model.rds` files are machine-specific compiled-model caches; they save compilation time, not sampling time, and are therefore not distributed.

## Reproduce

Run commands from the repository root with R 4.5 or newer. The primary models require RStan and the packages loaded by each script.

```text
Rscript scripts/recall/01_clean_recall_data.r
Rscript scripts/recall/02_estimate_population_and_recall.r
Rscript experiments/multilevel/analysis.R
Rscript experiments/multilevel/population_analysis.R
Rscript table/table1/create_table1.R
```

Raw inputs are not distributed, so reproducing from source requires authorized local copies in the documented `data/raw/` paths. Committed aggregate outputs allow inspection of the reported results without those files. Use `--quick` with the main multilevel script for a shorter pipeline check.

## Repository status

This snapshot accompanies the near-final manuscript. Code is licensed under the MIT License; the manuscript and figures remain copyright of their authors unless otherwise stated.
