# Release preparation: 2.0.2.9011

Prepared on 2026-10-08 from the targeting-cash pinned 9011 source archive
(SHA256 `3094b5bdbf3b5508eb86a53c89603fb707e6e0524e7aedb23bdbc9d1021b6eec`).
Existing fork changes are preserved where they agree with this snapshot;
the snapshot supplies the final sampler and generated bindings. README and
NEWS are updated for release. Historical simulation outputs are retained.

## Checks

Run from the repository root with the installed release package:

```sh
R CMD build --no-build-vignettes --no-manual .
# Set BCF_TEST_LIBRARY to a fresh library containing declared dependencies.
export R_LIBS_USER="$BCF_TEST_LIBRARY"
R CMD INSTALL --library="$BCF_TEST_LIBRARY" bcf_2.0.2.9011.tar.gz
Rscript -e 'stopifnot(as.character(packageVersion("bcf")) == "2.0.2.9011")'
BCF_EXPECT_VERSION=2.0.2.9011 BCF_JOINT_MEAN_EVERY=5 Rscript tests/hetero-mean-conditionals.R
Rscript tests/joint-variance-reference.R .
Rscript tests/paired-variance-reference.R .
Rscript tests/variance-change-reference.R .
Rscript tests/shared-sampler-replay.R
```

Release-checkout validation on 2026-10-08 used R 4.6.1 and GCC 16.2.1,
with a newly installed BCF in a fresh library and dependencies selected from
an existing local R cache. Package version and all four disabled optional
sampler defaults were asserted explicitly.

- Source build and installation: passed. Regenerating Rcpp bindings reproduced
  both checked-in binding files byte for byte.
- `R CMD check --no-manual --no-build-vignettes`: **Status: OK**, with
  `BCF_EXPECT_VERSION=2.0.2.9011`, `BCF_JOINT_MEAN_EVERY=5` and
  `_R_CHECK_FORCE_SUGGESTS_=false`. All executable sampler tests and vignette
  code passed. There were no check errors, warnings or notes; informational
  messages reported missing suggestions, installed size and GNU make.
  `testthat`, `spelling`, `latex2exp`, `rpart.plot` and `partykit` were unavailable
  offline. The empty testthat runner had no tests; spelling was skipped.
  Full suggested-package coverage, PDF manuals and vignette rebuilding remain
  unverified. No real-data convergence claim follows from this result.
- Analytical mean-conditionals with `joint_mean_every = 5`: passed, including
  fixed Gaussian and half-Cauchy scales.
- Joint and paired variance numerical references: passed, including weighted
  likelihoods and saved-tree prediction replay.
- Variance split-change reference: passed for single/multiple covariates and
  an ancestor-constrained region.
- Combined sampler replay: passed for fixed and half-Cauchy mean scales.
- Heteroscedastic smoke/replay, weighted variance, mean-scale prior,
  mean split-change reference and joint SBC smoke: passed.

Reference scripts now resolve absolute source paths before Rcpp compilation
and locate the source automatically during package checks. They also clear
the checker's relative `R_TESTS` startup path before child compilation. The obsolete empty
`testthat.R` runner was removed; it had no active test invocation. Missing
joint/paired variance help entries were added. These changes do not alter
sampler code. All 22 nongenerated sampler files match the pinned archive;
regenerated Rcpp bindings remain the documented exception.

These are implementation and synthetic checks, not real-data convergence
certification. Historical SBC failures remain documented separately.

## Publication review

All 138 file blobs reachable from the locally available Git refs were rescanned
for private-key headers, common GitHub/AWS/Google API credentials, and quoted
credential assignments. No matches were found. Filename review identified
only vignette build metadata and aggregate synthetic calibration CSVs; their
provenance is documented in the SBC results README. No participant datasets
were found in that review. The 92-file prepared source was also scanned,
including binary bytes; no
credential or absolute machine-path patterns matched. The built source archive
was scanned separately with no matches. CSV outputs contain aggregate synthetic
calibration summaries; the RDS contains vignette build metadata. Upstream
authors, citation, GPL-3 declaration and Thomas Wiemann attribution are retained.
This is a bounded review of locally available history, not a guarantee about
remote refs that cannot be fetched. Do not change repository visibility until
remote-only history has also been reviewed.

Remote publication is pending. The current environment cannot resolve GitHub
or CRAN hostnames, its default SSH configuration fails a permissions check,
and GitHub CLI reports an invalid authentication token. Bypassing the local SSH
configuration still fails DNS resolution. Remote branches, history, protection
rules and tag existence therefore remain unverified.

Before publication: fetch all remote branches/tags, reconcile newer work,
review remote-only history, verify the release tag is unused, push through the
repository's normal review process, verify remote SHAs, and install the tagged
source with authenticated access in a fresh R library. No release tag is
created until the remote tag check succeeds. Visibility and any rename remain
separate owner actions; anonymous access/install must be checked afterwards.
