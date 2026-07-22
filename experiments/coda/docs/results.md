# CoDa experiment results

> **Status:** Exploratory analysis retained for reproducibility. It uses the 17 July 2026 survey snapshot and is not part of the final multilevel paper or its conclusions.

## Main finding

Athletes' discipline profiles were the simplest useful predictors of how reported training was divided across:

- preparation with carbon insoles;
- preparation without carbon insoles;
- competition with carbon insoles;
- competition without carbon insoles.

The discipline-only model ranked first by athlete-level leave-one-out prediction. Its advantage over the more complex interaction model was small (ELPD difference 1.54, SE 9.38), so this supports choosing the simpler model, not claiming a clear predictive victory.

## Data used

- 337 eligible, complete surveys.
- 240 athletes had positive values in all four components and formed the primary analysis.
- The other 97 athletes were included in zero-replacement sensitivity checks.

## Average reported allocation

| Component | Observed | Model estimate (95% interval) |
|---|---:|---:|
| Preparation, carbon | 21.5% | 20.6% (19.4-22.0%) |
| Preparation, non-carbon | 31.4% | 32.4% (30.9-33.9%) |
| Competition, carbon | 23.1% | 23.2% (21.9-24.5%) |
| Competition, non-carbon | 24.0% | 23.8% (22.5-25.2%) |

The fitted averages closely reproduce the observed composition. This does not mean individual allocations are predicted precisely.

## Robustness

The discipline-only model remained first when:

- zeros were replaced with 0.25, 0.5, or 1 frequency-score unit;
- each of the four components was used as the ALR reference.

The broad conclusion is therefore insensitive to the tested zero handling and reference choice.

## Interpretation

Discipline appears more useful than demographics, body measurements, experience, prior injury, nonlinear terms, or the tested interactions for describing training-allocation patterns. The differences between the best complex models are uncertain, so individual coefficient estimates should remain exploratory.

## Limits

- The inputs are ordinal 0-5 frequency scores, not measured training time.
- Zero replacement is an assumption; structural non-use may have a different meaning.
- Discipline codes can be multi-selected and currently lack descriptive labels in the dataset.
- This analysis describes training composition. It does not estimate injury risk or a causal effect of carbon insoles.

Detailed outputs are under `experiments/coda/runs/`; the complete method and run command are in `experiments/coda/docs/guide.md`.
