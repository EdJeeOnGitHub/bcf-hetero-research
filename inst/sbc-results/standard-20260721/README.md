# Joint BCF SBC: standard profile, 2026-07-21

This directory records the compact outputs from the 250-replicate joint
prior-predictive simulation-based calibration run:

```sh
Rscript tools/sbc/joint-sbc.R \
  --profile standard \
  --workers 8 \
  --output sbc-output/standard-20260721
```

The run used the default seed `20260721`, 240 observations per replicate, two
chains, 500 burn-in iterations, and 200 retained iterations per chain. It
finished in 7.4 minutes. Each replicate fit the correctly specified
homoscedastic, shared-variance, and variance-ratio models. Homoscedastic fits to
the same shared and ratio draws provide paired misspecification benchmarks and
are not included in SBC ranks.

## Calibration result

Seventeen of 21 monitored estimands passed every pre-specified gate. All rank
uniformity checks passed, including every variance estimand. Four aggregate
gates failed:

- `shared/mean_tau`: 80% coverage was 219/250 (87.6%), three above the 99%
  Monte Carlo acceptance limit of 216/250. Its 50% and 90% coverage and rank
  checks passed.
- `ratio/point_log_var_ratio`: 50% coverage was 101/246 (41.1%), two below the
  99% acceptance limit of 103/246. Its 80% and 90% coverage and rank checks
  passed.
- `homoscedastic/point_mu`: 17/250 replicates (6.8%) were excluded by the
  `R-hat < 1.05` and `ESS >= 100` convergence gate, above the 5% limit.
- `ratio/point_mu`: 30/250 replicates (12.0%) were excluded by the same gate.

The two point-`mu` failures are convergence-exclusion failures rather than
rank or coverage failures. Without excluding those chains, homoscedastic
point-`mu` coverage was 51.2%, 81.2%, and 90.8%, and ratio point-`mu` coverage
was 52.0%, 80.8%, and 90.4% at the nominal 50%, 80%, and 90% levels.

The result supports joint mean and variance calibration, with longer chains
recommended for pointwise mean-function inference. The isolated marginal
coverage misses should be interpreted in the context of 21 estimands and
three coverage levels per estimand.

## Paired estimator comparison

Differences below are correctly specified heteroscedastic fit minus the
homoscedastic fit on the identical heteroscedastic draw. Lower RMSE and higher
predictive log score are better.

| Scenario | Mean-function RMSE | Treatment RMSE | Log-variance RMSE | Predictive log score |
|---|---:|---:|---:|---:|
| Ratio | -0.00783 | -0.00321 | -0.21425 | +0.04620 |
| Shared | -0.00390 | -0.00191 | -0.10382 | +0.01720 |

The paired summary CSV includes Monte Carlo standard errors. The 95% Monte
Carlo intervals exclude zero for both log-variance RMSE and predictive score
improvements, as well as both mean-function RMSE improvements. The shared
treatment-RMSE interval includes zero.

## Committed outputs

- `rank-summary.csv`: convergence, rank-uniformity, and coverage gates.
- `benchmark-summary.csv`: aggregate recovery and predictive metrics.
- `paired-benchmark-summary.csv`: heteroscedastic-minus-base differences and
  Monte Carlo standard errors.
- `rank-histograms.png` and `rank-ecdfs.png`: rank diagnostics.

Raw ranks, replicate-level benchmarks, sampler logs, and the complete RDS fit
are intentionally omitted. They can be regenerated with the command above.
