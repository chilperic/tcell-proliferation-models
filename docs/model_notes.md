# Model notes

## Data format

The synthetic CFSE-like datasets have one row per sampled time point. Most files use 21 generations, indexed from 0 to 20. Each generation has two columns:

- `Nkd`
- `Nkp`

where `k` is the generation index. The refactored reader sums the pair into total generation count `Nk`.

Some files include a `Time` column and some do not. If no `Time` column is available, the reader assumes sampled rows correspond to 12, 24, ..., hours and prepends the initial condition at time zero.

## Smith-Martin model

The implementation follows the delayed differential-equation structure from the thesis. The state vector contains:

```text
NA0, NA1, ..., NAg, N0, N1, ..., Ng
```

where `NAk` is the A-phase count for generation `k` and `Nk` is the total count for generation `k`.

The parameter vector is:

```text
lambda_0, lambda, d_0, d, delta
```

The legacy name `deltat` is accepted and internally normalized to `delta`.

## Cyton model

Version 0.2 defaults to a stable exponential competing-risk ODE. The historical analytic formula is retained with method="legacy". See [migration notes](MIGRATION-0.2.md).

The parameter vector is:

```text
f0, f1, mu0d, mud, mu0p, mup
```

where `mu` parameters represent mean waiting times.

## Fitting

The fitting functions use box-constrained `optim(..., method = "L-BFGS-B")` with multiple random starts. The defaults are intentionally modest. For serious reproduction, use many more starts and record the random seed.
