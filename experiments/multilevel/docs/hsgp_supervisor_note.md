# Note for supervisor: Gaussian-process dose-response model

## Short answer

The supervisor's question is whether our main result depends on forcing the carbon-share effect into a too-simple shape.

The answer is no: the current primary analysis already uses a flexible Bayesian smooth for carbon share. It is a Bayesian multilevel logistic model with:

- the same adjustment variables as before;
- an athlete-level random intercept;
- partially pooled discipline effects;
- a separate adoption term for any carbon use;
- a smooth dose-response curve for carbon share among users.

In practical terms, this is a Bayesian generalized additive mixed model. We implemented the smooth curve directly in Stan as a Hilbert-space approximate Gaussian process (HSGP), using 40 basis functions.

## What questions does the model answer?

The model separates the exposure question into three simpler questions.

1. **Adoption:** is a median-share carbon user different from a non-user?
2. **Dose among users:** among carbon users, does injury prevalence change as carbon share increases?
3. **Robustness:** does the competition-period association remain when the curve is allowed to be nonlinear, instead of being forced to be linear or fixed as a 3-df spline?

This distinction is important because non-users do not have a meaningful carbon-share dose. We therefore model non-use separately from dose among users.

## Model in one equation

For athlete $i$ in period $p$, the model is:

$$
\text{logit}\{P(y_{ip}=1)\}
= \alpha_p
+ \beta^{\mathrm{any}}_p A_{ip}
+ f_p(s_{ip})
+ X_{ip}\beta
+ Z_i\delta
+ u_i .
$$

Here:

- $A_{ip}$ is any carbon use;
- $s_{ip}$ is carbon share;
- $f_p(s_{ip})$ is the smooth carbon-share curve;
- $u_i$ is the athlete random intercept;
- $Z_i\delta$ is the partially pooled discipline contribution;
- $X_{ip}\beta$ is the adjustment set.

The GP curve is zero for non-users. Among users, carbon share is centered at the pooled median user share, 0.50. This makes $\beta^{\mathrm{any}}_p$ interpretable as non-user versus median-share user.

## What is Bayesian about the curve?

The model does not estimate one curve and then treat it as fixed.

Each posterior draw contains one possible dose-response curve, together with the other model parameters. Across all posterior draws, we get a posterior distribution over plausible curves.

For one posterior draw, the curve is:

$$
f_p(x) = \sum_{m=1}^{40} \phi_m(x)\,\beta^{\mathrm{share}}_{p,m}.
$$

The $\phi_m(x)$ terms are fixed basis functions. The $\beta^{\mathrm{share}}_{p,m}$ terms are estimated by Stan.

So the intuition is: the final uncertainty is an average over many plausible smooth curves, weighted by their posterior support. Operationally, this is done by MCMC sampling from the joint posterior, not by manually assigning weights to separate curves.

## What does the HSGP approximation do?

A full Gaussian process would require a covariance matrix across all observed share values. That is more expensive and less convenient inside this custom repeated-period multilevel model.

The HSGP approximation replaces the full GP with a finite basis expansion:

$$
f_p(x) \approx \sum_{m=1}^{40} \phi_m(x)b_{p,m}.
$$

This keeps the model interpretable as a GP smooth, but makes it computationally manageable in Stan.

We checked that 40 basis functions were enough. In the primary model, the highest-frequency basis term had negligible remaining spectral weight:

```text
preparation: max ratio = 0.0019
competition: max ratio = 0.0034
```

Both are below the 0.01 adequacy threshold used in the script.

## Why Matern 5/2?

The Matern 5/2 kernel allows smooth curves, but not unrealistically smooth ones.

That fits this setting. Carbon share is a bounded proportion, the sample size is modest, and we expect a continuous exposure-response relationship. A squared-exponential kernel would assume a very smooth curve. A rougher Matern kernel would allow more local wiggle. Matern 5/2 is a conservative middle choice.

It is also close in spirit to replacing the earlier natural spline with a more adaptive smooth: still smooth, but with flexibility learned from the data.

## Why write it directly in Stan?

We wrote it directly in Stan because the rest of the analysis is already a custom Stan model.

That lets us keep the exact same structure for:

- the two period records per athlete;
- the shared athlete random intercept;
- period-specific exposure effects;
- discipline partial pooling;
- posterior-standardized prevalence contrasts;
- athlete-grouped PSIS-LOO;
- the matched no-exposure baseline;
- the population-weighted sensitivity model.

Using `brms::gp(share)` would be reasonable for a simpler model. Here, direct Stan code makes the separation between adoption and dose explicit, keeps non-users at zero GP contribution, and lets all downstream estimands use the same posterior draws.

## What changed from the earlier spline model?

The earlier model used a 3-df natural spline for carbon share. That meant the amount of curvature was fixed before fitting.

The HSGP model keeps the same multilevel structure but replaces the fixed spline with a GP smooth. The data estimate both:

- how large the curve can be, through $\alpha^{\mathrm{gp}}_p$;
- how quickly it can bend, through $\rho^{\mathrm{gp}}_p$.

The main gain is not complexity for its own sake. The gain is that uncertainty about the shape of the dose-response curve is included in the final uncertainty intervals.

## Interpretation

The GP does not discover thresholds and does not create exposure groups.

It estimates one continuous curve among carbon users. The P25-versus-P75 contrast and the plotted 25%, 50%, and 75% anchors are summaries of that curve.

The conclusion should therefore be phrased as a flexible dose-response sensitivity of the multilevel model, not as evidence for a threshold.

## Code map

The implementation is split into two parts: R builds the basis; Stan estimates the posterior.

Main R file:

- `experiments/multilevel/analysis_hsgp.R:85`: sets Matern 5/2 with `matern_nu <- 2.5`.
- `experiments/multilevel/analysis_hsgp.R:93`: defines the Hilbert-space sine basis.
- `experiments/multilevel/analysis_hsgp.R:206`: centers carbon share at the pooled user median.
- `experiments/multilevel/analysis_hsgp.R:212-216`: sets the boundary, 40 basis functions, and basis frequencies.
- `experiments/multilevel/analysis_hsgp.R:221-225`: builds the GP basis only for users; non-users get zero GP contribution.
- `experiments/multilevel/analysis_hsgp.R:273-291`: passes the HSGP objects and priors to Stan.
- `experiments/multilevel/analysis_hsgp.R:418-426`: checks whether 40 basis functions are adequate.
- `experiments/multilevel/analysis_hsgp.R:475-476`: writes the approximation diagnostic table.

Main Stan file:

- `experiments/multilevel/model_hsgp.stan:5-11`: defines the Matern spectral density.
- `experiments/multilevel/model_hsgp.stan:23-33`: declares the HSGP basis and GP prior inputs.
- `experiments/multilevel/model_hsgp.stan:39-41`: declares GP basis coefficients and hyperparameters.
- `experiments/multilevel/model_hsgp.stan:53-57`: scales the basis coefficients using the Matern spectral density.
- `experiments/multilevel/model_hsgp.stan:64-66`: applies priors to the GP components.
- `experiments/multilevel/model_hsgp.stan:74-79`: inserts adoption plus GP dose into the multilevel logistic likelihood.

Weighted sensitivity:

- `experiments/multilevel/weighted_model_hsgp.stan:72-76`: repeats the same HSGP exposure term in the population-weighted pseudo-posterior model.
