# Operative guide: Bayesian multilevel carbon-use and injury analysis

## Purpose

This experiment estimates cross-sectional associations between reported carbon-insole use and current-season injury. It separates:

- **adoption:** a median-share carbon user versus a nonuser;
- **dose among users:** the 75th versus the 25th percentile of carbon share.

Preparation exposure is paired with preparation injury, and competition exposure with competition injury. The design improves period matching and accounts for the two records contributed by each athlete, but it does not establish causality or temporal order.

## Data and analysis population

The source is `data/raw/recall/DSCBATLETAIT-CompleteCase_DATA_FINAL.csv`. Inclusion requires consent, calculated eligibility, and a complete survey.

- 355 athletes are eligible.
- One is excluded because injury timing is incomplete.
- Two are excluded because a period has zero total eligible training-frequency score.
- The balanced analysis contains 352 athletes and 704 athlete-period records.

Period injury is derived from the reported injury timing. Codes 1 or 3 identify preparation injury; codes 2 or 3 identify competition injury. There are 115 preparation injuries and 149 competition injuries.

For each period, carbon share is

\[
\text{carbon share} = \frac{\text{sum of five carbon frequency items}}
{\text{sum of carbon and non-carbon frequency items}}.
\]

The five domains are long endurance, middle endurance, lactate, maximum velocity, and technique. Each item is a top-coded weekly frequency count (0, 1, 2, 3, 4, or 5+ sessions/week); the exported data encode 5+ as 5. The sums therefore approximate reported type-specific weekly frequency, not minutes, workload, or necessarily unique sessions. `Palestra` is kept as a separate non-carbon adjustment variable because it has no matching carbon item.

## Primary model

Each athlete contributes a preparation and a competition Bernoulli outcome. The log-odds model contains:

- period-specific intercepts;
- period-specific indicators for any carbon use;
- a period-specific three-degree-of-freedom natural spline for carbon share among users;
- a shared athlete random intercept;
- age, competition sex, height, weight, athletics experience, and prior injury;
- period-matched total frequency scores for the five training domains and `palestra`;
- hierarchically shrunk effects for the seven observed discipline indicators.

Continuous adjustment variables are standardized. The athlete and discipline effects use non-centered parameterizations. Log-odds priors are regularizing: exposure SD 0.7, adjustment SD 0.5, discipline scale half-normal(0, 0.5), and athlete scale half-normal(0, 1).

The baseline model removes both carbon-use terms while retaining the same outcome, adjustment set, discipline structure, and athlete random intercept.

## Estimands

Posterior predictions are standardized over the observed adjustment-variable and discipline distributions. The athlete random intercept is marginalized using 15-point Gauss-Hermite quadrature.

For each period, the reported estimands are:

- adjusted prevalence for nonusers and for low-, median-, and high-share users;
- adoption prevalence difference and ratio: median user versus nonuser;
- dose prevalence difference and ratio: 75th versus 25th percentile among users;
- posterior probability that each association is positive;
- posterior probability that a difference exceeds 5 percentage points or a ratio exceeds 1.10;
- posterior probability inside a practical-null region of +/-2 percentage points or ratio 0.95-1.05.

The user quartiles are 28.6%, 42.1%, and 57.1% in preparation, and 36.0%, 50.0%, and 71.4% in competition.

## Validation and sensitivity analyses

Predictive comparison uses PSIS-LOO with the two period records grouped by athlete. The athlete random intercept is integrated out before computing each held-out athlete's joint likelihood. Posterior predictive checks compare observed and replicated injury counts by period.

Two sensitivity models test the main modeling choices:

1. a linear carbon-share term instead of the spline;
2. wider exposure and adjustment priors (SD 1.2 and 1.0).

The primary association is compared only with the matched no-exposure baseline. Sensitivity models are not treated as competing substantive hypotheses.

## Population weighting and conditional estimates

The primary analysis is unweighted and targets the observed analytic sample. `population_analysis.R` adds a Bayesian pseudo-posterior sensitivity using post-stratification weights for six sex-by-age-category cells. The target frame comes from `tables/recall/population_history_by_age_sex.csv`, derived from `data/raw/recall/raccolta_dati_campione.xlsx`:

- Under 20 uses the available 2026 championship counts;
- Under 23 and Senior use the rounded mean of 2023-2025 counts because 2026 counts are incomplete;
- target size is 1,328 athletes: 743 men and 585 women;
- model weights are normalized to mean one and range from 0.72 to 1.26;
- Kish effective sample size is 340.5 of 352.

This weighting can improve representation of that championship frame if selection is ignorable within the six cells. It does not make the survey representative of all Italian athletics participants and cannot correct selection within cells.

Sex- and age-conditional estimands standardize the weighted model within men, women, Under 20, Under 23, and Senior athletes using common exposure coefficients. They are descriptive risk-scale estimates, not tests of exposure-effect interaction.

## Table 1

`table/table1/table1_by_sex.csv` follows Hayes-Larson et al., [Who is in this study, anyway?](https://pmc.ncbi.nlm.nih.gov/articles/PMC6773463/): it includes all analysis variables, unweighted and target-weighted distributions, sex strata because sex-conditional results are examined, signed SMDs rather than p-values, and missing counts. The SMD direction is women minus men.

## Folder system

```text
experiments/multilevel/
|-- analysis.R
|-- model.stan
|-- population_analysis.R
|-- toolchain.R
|-- weighted_model.stan
|-- visualize_results.R
|-- visualize_distributions.R
|-- model.rds
|-- docs/
|   |-- guide.md
|   `-- results.md
`-- runs/
    |-- 00_data_audit/tables/
    |-- 01_primary/
    |   |-- association/{models,plots,tables}/
    |   |-- baseline/{models,plots,tables}/
    |   `-- tables/
    |-- 02_sensitivity/
    |   |-- linear_share/{models,plots,tables}/
    |   |-- wide_prior/{models,plots,tables}/
    |   `-- tables/
    |-- 03_population_weighting/{models,tables}/
    |-- 04_subgroups/tables/
    `-- run_config.csv
```

`model.rds` is a local compiled Stan cache and is ignored by Git.

On Windows, `toolchain.R` discovers Rtools 4.5 in its standard installation directory. Set the `RTOOLS45_HOME` environment variable before running the pipeline only when Rtools is installed elsewhere.

## Run

From the repository root:

```powershell
Rscript experiments/multilevel/analysis.R
```

For a short pipeline check:

```powershell
Rscript experiments/multilevel/analysis.R --quick
```

Run the population-weighted and conditional sensitivity with:

```powershell
Rscript experiments/multilevel/population_analysis.R
Rscript table/table1/create_table1.R
```

The full run uses seed `20260719`, four chains, 2,000 iterations per chain, and 1,000 warmup iterations. It overwrites only files under `experiments/multilevel/runs/`.

## Interpretation rule

The defensible claim is an adjusted association, not an effect: among carbon users, greater reported carbon share during competition was associated with greater competition-period injury prevalence. The result persists after age-sex post-stratification to the defined championship frame. Adoption itself and preparation-period dose showed no clear association. The competition magnitude is model-sensitive and the exposure model did not clearly improve held-out prediction, so the result is suggestive rather than confirmatory.
