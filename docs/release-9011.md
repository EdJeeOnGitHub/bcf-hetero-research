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
Remote publication review on 2026-10-08 fetched all fork branches and tags,
plus PR heads. The default branch was `master` at `54aba577`; the release
branch matched prepared commit `ca23fd1`. No newer default-branch work existed.
Both open PRs were integrated with their original commits and authors: #4
for parallel prediction worker library paths and #5 for the macOS FLIBS
installation workaround. No sampler source was changed during publication.
All optional sampler defaults remain disabled. The targeting-cash campaign
configuration remains joint mean every 5 iterations, variance split changes,
joint/paired variance every iteration, and 4000 burn-in/2000 retained draws.

The refreshed reachable-history scan covered 172 unique file blobs before
publication documentation edits, with no credential or personal machine-path
matches. Content review covered the historical file inventory, branch diffs,
implementation notes, synthetic vignette generators, aggregate calibration CSV
schemas and provenance, vignette RDS metadata, and calibration graphics.
No participant datasets, compiled executables or private operational notes
were identified. Upstream authors, citation, GPL-3 and Thomas Wiemann's
bayesm.HART attribution are preserved.

GitHub review covered open and closed issues, issue comments, PR bodies,
inline comments, reviews, releases, Actions runs and artifacts. There were
no Actions runs, artifacts, releases, inline comments or PR reviews. Discussions
and Pages were disabled; the enabled wiki had no accessible Git repository.
Issue content contains technical reproduction details and campaign-relative
paths, not participant data or personal absolute machine paths. The repository
had no branch protection or repository rulesets. Publication uses a release
PR and merge without force pushes. The requested name was available and the
existing repository was renamed to `EdJeeOnGitHub/bcf-hetero-research`; the
`hetero` remote and README installation command were updated. `origin` and its
disabled push URL are retained.

Final validation results are recorded below before tagging. Authenticated
and anonymous tagged installations and remote SHA verification are performed
as publication gates; their outcome is reported with the public handoff.

### Final source validation

Publication validation on 2026-10-08 rebuilt and installed the reconciled
source in a fresh library using native R 4.6.1 and GCC 16.2.1. Declared
dependencies were available, including RcppArmadillo 15.6.0.1 and
RcppParallel 6.2.1. Version, selected installation path and all four disabled
optional defaults were asserted. All five required sampler checks passed:
mean conditionals with joint mean every 5, joint variance reference, paired
variance reference, variance split-change reference and combined replay.
Heteroscedastic smoke, weighted variance, mean-scale prior, mean split-change
and joint SBC smoke also passed. Parallel prediction with two workers passed
saved-tree replay for both shared and ratio models when the fresh library was
selected in-script and library environment overrides were removed.

`R CMD check --no-manual --no-build-vignettes` with
`_R_CHECK_FORCE_SUGGESTS_=false`, version 2.0.2.9011 and joint mean every 5
finished with **Status: OK**, including all executable tests. The five missing
suggestions remain testthat, spelling, latex2exp, rpart.plot and partykit.
Full suggestions, spelling, PDF manuals and vignette rebuilding remain
unverified; these synthetic checks do not certify real-data convergence.
An initial standalone reference compile encountered the sandbox's read-only
compiler cache; rerunning with `CCACHE_DISABLE=1` passed without source edits.
No sampler files differ from the prepared, archive-validated commit.
