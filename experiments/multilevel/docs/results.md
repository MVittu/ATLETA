# Bayesian multilevel experiment results

## Main finding

Among carbon users, moving from the 25th to the 75th percentile of competition carbon share was associated with **14.3 percentage points higher adjusted injury prevalence** (95% credible interval 2.7 to 25.8; 99.0% posterior probability of a positive association). The adjusted prevalence ratio was 1.41 (95% credible interval 1.06 to 1.86).

This is the only clear primary association. It is not a causal estimate.

## Primary adjusted associations

| Comparison | Prevalence difference (95% CrI) | Probability positive |
|---|---:|---:|
| Preparation: median user vs nonuser | +5.4 pp (-6.3 to +16.2) | 82.0% |
| Preparation: high vs low share among users | -0.8 pp (-10.6 to +9.1) | 43.8% |
| Competition: median user vs nonuser | +4.8 pp (-8.7 to +17.7) | 77.1% |
| Competition: high vs low share among users | +14.3 pp (+2.7 to +25.8) | 99.0% |

Therefore, adoption of carbon use was not clearly associated with injury in either period, and preparation-period dose was also unclear.

## Sensitivity

The competition dose association remained positive, but its size depended on specification:

| Model | Difference (95% CrI) | Probability positive |
|---|---:|---:|
| Primary nonlinear share | +14.3 pp (+2.7 to +25.8) | 99.0% |
| Linear share | +6.0 pp (-0.3 to +12.6) | 96.9% |
| Wider priors | +18.7 pp (+5.2 to +31.9) | 99.5% |

The consistent direction supports a possible positive association during competition. The changing magnitude and the linear interval crossing zero argue against presenting 13 percentage points as a stable effect size.

## Population weighting

The original model was not survey-weighted. Post-stratification to the estimated 2026 championship frame by sex and age category changed the competition dose estimate minimally:

| Analysis | Difference (95% CrI) | Probability positive |
|---|---:|---:|
| Unweighted | +14.3 pp (+2.7 to +25.8) | 99.0% |
| Age-sex weighted | +14.9 pp (+3.3 to +26.4) | 99.4% |

The normalized weights ranged from 0.72 to 1.26 and reduced the effective sample size from 352 to 340.5. The weighted result applies only to the championship frame represented by the recall workbook, subject to no residual selection within the six post-strata.

## Conditional estimates

For competition-period users, the weighted Q75-versus-Q25 prevalence differences were:

- men: +15.4 pp (95% CrI +3.3 to +26.9);
- women: +14.4 pp (+3.2 to +25.5);
- Under 20: +14.3 pp (+3.1 to +25.3);
- Under 23: +15.3 pp (+3.3 to +26.9);
- Senior: +15.5 pp (+3.4 to +27.1).

These similar conditional estimates use common exposure coefficients. They do not demonstrate absence of sex or age interaction; an interaction model would be required for that claim.

## Table 1

`table/table1/table1_by_sex.csv` reports the analytic sample overall and by sex, alongside target-weighted distributions and signed SMDs. The largest expected sex differences are height and weight; weighting reduces the age SMD from 0.198 to 0.137. No p-values are used.

## Model checks

- All four models had maximum R-hat at or below 1.004, no divergent transitions, no maximum-tree-depth hits, and acceptable E-BFMI.
- Posterior predictive means reproduced the observed counts: 115.5 versus 115 preparation injuries and 149.5 versus 149 competition injuries.
- Every athlete-level PSIS-LOO Pareto-k value was at most 0.5.
- The exposure model improved held-out ELPD by only 2.22 (SE 2.30) over the no-exposure model. This is not a clear predictive advantage.
- The weighted model had maximum R-hat 1.002, no divergences or tree-depth hits, and minimum E-BFMI 0.599.

## Conclusion and limits

The dataset supports one carefully worded result: **higher carbon share among competition-period users was associated with higher competition-period injury prevalence after adjustment**. It does not support claims that starting carbon use increases injury, that the same pattern exists in preparation, or that carbon use causes injury.

The main limits are retrospective cross-sectional reporting, uncertain temporal order, residual confounding, top-coded weekly frequency counts rather than measured use time or workload, only 44 competition nonusers, sensitivity to the dose-response form, and a target frame limited to championship athletes. Longitudinal data would be needed for a stronger causal claim.

Detailed tables and plots are under `experiments/multilevel/runs/`; the complete process is in `experiments/multilevel/docs/guide.md`.
