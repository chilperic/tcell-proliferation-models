# Reproduce and inspect

The validated environment uses R 4.4 and packages recorded in
`environment-lock-linux-64.txt`. Recreate the exact conda package builds on
Linux x86-64 with micromamba:

```bash
micromamba create -p ./r-env --file environment-lock-linux-64.txt
micromamba run -p ./r-env R CMD INSTALL .
micromamba run -p ./r-env Rscript -e 'testthat::test_local()'
micromamba run -p ./r-env Rscript scripts/run_model.R examples/reproducible_run.json runs/example-001
```

The CSV and run.json record parameters, time coordinates, solver method/history,
runtime/package versions and output checksum. Preserve these together.
The explicit lock is Linux-specific. On other platforms install the package's
DESCRIPTION dependencies and rerun tests; this does not guarantee bitwise
identical results. Docker instructions preserve the earlier R 4.4.3 snapshot
workflow; it is separate from the tested conda environment.

Read [migration notes](MIGRATION-0.2.md) before reusing thesis fits. The Cyton
numerical backend has changed; the historical analytic formula remains
available with method="legacy". Smith-Martin's pre-zero history is a scientific
assumption, not merely a solver option. Select it explicitly. All original
synthetic observations are fingerprinted in data-manifest.json.
