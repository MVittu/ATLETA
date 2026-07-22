# Operative guide: CoDa models of training allocation

> **Status:** Exploratory analysis retained for reproducibility. It uses the 17 July 2026 survey snapshot and is not part of the final multilevel paper or its conclusions.

## Purpose

This experiment describes how athletes allocate their reported insole-eligible training across four jointly constrained components:

1. preparation with carbon insoles;
2. preparation without carbon insoles;
3. competition with carbon insoles;
4. competition without carbon insoles.

Each component is the sum of the five matching 0–5 frequency items for `fondo_lungo`, `fondo_medio`, `lattacido`, `max_velocity`, and `tecnica`. `Palestra` is excluded because the survey has no matching carbon-insole item. The four scores are closed to sum to one for CoDa modeling.

This is an exploratory composition of ordinal frequency scores, not a measured proportion of training time. It must not be interpreted as a causal injury model or as exact training volume.

## Model scope and rationale

The survey contains no spatial or group identifiers, so the analysis does not include spatial or group random effects. The implemented model has these features:

- all four components are modeled jointly rather than separately;
- ALR coordinates map the four-part simplex to three real-valued outcomes;
- the three coordinates have an unrestricted residual covariance matrix;
- DIC and WAIC use each athlete's joint three-coordinate likelihood;
- LCPO uses exact leave-one-athlete-out posterior prediction, excluding all three coordinates together;
- component probabilities are obtained by inverse-ALR transformation.

This experiment uses a logistic-normal model rather than a Dirichlet model. Its unrestricted 3×3 covariance lets the data identify the dependence among the three ALR coordinates.

## Population, zeros, and ALR reference

The source is `data/raw/recall/DSCBATLETAIT-CompleteCase_DATA_2026-07-17_1419.csv`. Eligible rows require consent, calculated eligibility, and a complete survey.

- 337 athletes are eligible and complete.
- 240 have positive values in all four components and form the primary analysis.
- 97 contain at least one zero and enter three sensitivity analyses only.
- Zero sensitivities replace zeros with 0.25, 0.5, or 1 frequency-score unit, leave positive scores unchanged, and then close each row.

ALR requires positive components. In the primary data, `prep_noncarbon` has the lowest variance of log composition (0.153) and is therefore the reference. A separate check repeats exact LOO ranking with every possible reference; model ranking, not raw density values across different coordinate systems, is compared.

## Joint Bayesian model

For athlete \(i\), the response is

\[
z_i = \left[
\log(y_{i,prep\_carbon}/y_{i,prep\_noncarbon}),
\log(y_{i,competition\_carbon}/y_{i,prep\_noncarbon}),
\log(y_{i,competition\_noncarbon}/y_{i,prep\_noncarbon})
\right].
\]

Each candidate fits

\[
z_i \sim \mathcal{N}_3(x_i B, \Sigma),
\]

with a conjugate matrix-normal/inverse-Wishart prior. Predictors are standardized where continuous. The prior is regularizing and permits exact posterior simulation and exact multivariate Student-t leave-one-out prediction without extra packages.

Current-season injury is not used as a predictor because it may follow or alter the recalled training allocation. Prior two-season injury is used only in candidate structures that explicitly include history.

## Candidate structures

| Type | Predictors |
|---|---|
| I | Intercept only |
| II | Age and competition sex |
| III | Age, competition sex, height, and weight |
| IV | Age, competition sex, athletics experience, and prior injury |
| V | Eleven discipline indicators |
| VI | All main effects above plus log total frequency score |
| VII | Type VI plus quadratic age and experience |
| VIII | Type VI plus age-by-sex and experience-by-prior-injury interactions |

LCPO is the primary selection criterion because it evaluates genuinely held-out athletes. WAIC is secondary; DIC is retained only because it was requested. Differences in expected log predictive density and their standard errors must be reported so small ranking differences are not overstated.

## Folder system

```text
experiments/coda/
├── analysis.R
├── docs/
│   └── guide.md
└── runs/
    ├── 00_data_audit/
    │   └── tables/
    ├── 01_primary_positive/
    │   ├── models/
    │   ├── plots/
    │   └── tables/
    ├── 02_zero_sensitivity/
    │   ├── pseudocount_0p25/
    │   ├── pseudocount_0p5/
    │   ├── pseudocount_1/
    │   └── tables/
    └── 03_reference_check/
        └── tables/
```

Component results stay in joint tables with a `component` or `coordinate` column. Separate GC1–GC4 directories are not created because they would imply four independent analyses and repeat files unnecessarily.

## Run

From the repository root:

```powershell
Rscript experiments/coda/analysis.R
```

For a short pipeline check using 200 rather than 1,000 posterior draws:

```powershell
Rscript experiments/coda/analysis.R --quick
```

The run is deterministic under seed `20260719` and overwrites only its own output files.

## Current full-run result

Type V (discipline profile) has the best athlete-level LCPO in the primary analysis. It also ranks first under all three zero replacements and all four ALR references. Its advantage over Type VIII is small: ELPD difference 1.54 with SE 9.38, so the result supports parsimony rather than a decisive predictive separation.

Among the 240 all-positive athletes, observed mean allocations are 21.5% preparation-carbon, 31.4% preparation-noncarbon, 23.1% competition-carbon, and 24.0% competition-noncarbon. The selected model's posterior population means are 20.6%, 32.4%, 23.2%, and 23.8%, respectively.

The discipline profile is useful for describing allocation differences, but the experiment does not establish that allocation causes or prevents injury. Structural zeros, ordinal frequency coding, retrospective measurement, and multi-selected discipline indicators remain the main limitations.
