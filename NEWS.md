# plantstressR 0.7.0

Two more trials to work on, so that every part of the package can be tried
without bringing your own data first.

## New data

* `wheat_salinity`: six cultivars grown at three sites under sodium chloride,
  in complete blocks within each site, 162 plots and 14 traits. This is the
  data set the multi-environment functions were written for. It carries a real
  genotype-by-environment structure -- `W2` holds a narrow band of ranks
  everywhere, `W4` runs from near the top on the upland site to near the bottom
  on the coastal one -- so `stress_stability()` has something to report. Its
  ion-relation traits (`Na`, `Cl`, `K`) also exercise the direction dictionary,
  which resolves rising sodium and falling potassium as damage without being
  told.
* `maize_heat`: 38 plants, four hybrids, one heat treatment, and every common
  defect on purpose -- a cell with two control plants, a trait abandoned
  partway through the campaign at 37 per cent missing, and a constant column
  that is numeric but is not a trait. Until now every bundled trial was clean,
  which left `validate_stress_data()` with nothing to report and users with no
  idea what its warnings look like.

Both are simulated and documented as such, and the scripts that build them are
in `data-raw/` in the package repository.

## Documentation

* The multi-environment example in the README and the vignette used to
  manufacture two sites by splitting the blocks of the drought trial, which
  confounded site with block. It now uses `wheat_salinity`, where the sites are
  sites.
* `validate_stress_data()` gains a second example, on a trial that actually
  fails its checks.
* The dashboard lets you pick which of the three trials to load.

## Bug fixes

* The manual claimed `maize_heat` had 11 columns. It has 12.

# plantstressR 0.6.0

Found by an adversarial re-read of the code added in 0.3.0-0.5.0, prompted by
two earlier defects of the same shape: check-green, plausible number, wrong
conclusion.

## Bug fixes

* The standard error of a standardized index assumed that the control and the
  stress group have the same variance. They frequently do not -- an inflated
  variance under stress is one of the things Glass's delta exists to cope with
  -- and the assumption contradicted the Welch test already used for
  `p_value`. The numerator variance is now taken from the data,
  `sqrt(s_c^2 / n_c + s_s^2 / n_s)`, which reduces **exactly** to the former
  `sqrt((n_c + n_s) / (n_c * n_s))` when the two spreads agree, so balanced
  homoscedastic trials are unaffected.
* That assumption also made `weights = "precision"` behave differently
  depending on whether `block` had been supplied: the blocked path noticed
  unequal dispersion and the unblocked path could not. On a two-trait example
  where the stress treatment inflated one variance sixfold, the unblocked path
  reported both traits as equally precise; it now gives the erratic one three
  per cent of the weight. `se` and the confidence limits widen accordingly for
  traits whose variance the treatment changed.

## Documentation

* The claim that `"precision"` coincides with `"equal"` in any balanced trial
  was too strong: that holds when the dispersion of the two groups is also
  comparable. Corrected in the manual and the vignette.
* `cv_isi` in `stress_stability()` is only interpretable where `mean_isi` is
  clearly positive, the integrated index being signed. Stated in the manual.
* The internal note claiming the block adjustment preserves treatment means
  "exactly" was true only for balanced layouts. Under an unbalanced one the
  means do move, and should: a treatment over-represented in the better blocks
  was being flattered by them. The treatment contrast is left at its
  least-squares value either way, which was verified against `lm()`.

# plantstressR 0.5.1

Fixes found by the Linux job of the continuous integration, which failed while
the Windows and macOS jobs passed.

## Bug fixes

* The test suite asserted that the bootstrap interval contains the point
  estimate. It is not required to: a percentile interval built from a skewed
  bootstrap distribution can legitimately exclude the observed value, and the
  assertion held on two platforms only by numerical luck. The invariant was
  wrong, not the interval, and the documentation now states the property.
* `stress_index_ci()` warns when a unit is laid out in fewer than four blocks.
  Blocked trials are resampled block by block, so it is the number of blocks and
  not the number of plants that decides how many distinct resamples exist; with
  two or three of them a large share of the draws repeat the same block and the
  interval is coarse.
* `plot()` on a bootstrap interval no longer calls the deprecated
  `ggplot2::geom_errorbarh()`.
* British spellings corrected in the documentation, which declares `en-US`, and
  the German title of the Wricke reference added to the word list.

# plantstressR 0.5.0

Multi-environment trials.

## New features

* `by` now accepts several columns. `by = c("genotype", "site")` analyses each
  genotype at each site against the control of that cell, which is how a
  multi-site or multi-year trial is described. The levels are pasted into one
  `unit` label so that every downstream table keeps a single key.
* `stress_stability()` summarizes such a trial: how severely each genotype was
  affected on average, and how consistently. It reports the spread of the index
  across environments (`sd_isi`, `cv_isi`), the spread of the genotype's
  position (`mean_rank`, `rank_min`, `rank_max`) and Wricke's `ecovalence`,
  the share of the genotype-by-environment interaction each genotype carries --
  computed on the integrated index rather than on yield. A genotype that ranks
  first at one site and last at the next is a different proposition from one
  that ranks second everywhere, and a mean hides the difference.
* `plot()` on the result places severity against consistency, so the material
  that is both tolerant and predictable sits in one corner.
* The dashboard gains an environment picker and an "Across environments" tab.

## Details worth knowing

* Ecovalence needs a complete genotype-by-environment table and at least two of
  each; genotypes missing from some environments are kept, with `n_env` saying
  how many they appeared in, but their ecovalence is `NA`.
* Ranks inside `stress_stability()` are recomputed within each environment, so
  they describe a genotype's standing among the material it was actually grown
  beside.
* `validate_stress_data()` checks every grouping column, and reports cells that
  have no control observations using the full combination of levels.

# plantstressR 0.4.0

An uncertainty statement for the ranking.

## New features

* `stress_index_ci()` resamples the trial and recomputes indices, weights,
  aggregation and ranking on every draw, returning a percentile interval for
  the integrated index together with the range of ranks each unit took and
  `p_best`, the share of resamples in which it came out first. A ranking table
  looks decisive; `p_best` says whether it is. The resampling respects the
  layout: whole blocks are drawn in a blocked trial, plants inside each
  unit-by-treatment cell otherwise.
* `plot()` on the result draws the index with its interval, shaded by `p_best`.
* The dashboard gains an "Is the ranking real?" tab wired to the same function.

## Details worth knowing

* `stress_index_ci(scale = )` decides whether the scaling standard deviation is
  resampled along with the response. It is an estimate too, and on few control
  plants it can collapse -- on the bundled trial the smallest control standard
  deviation is `0.0027` for `Fv/Fm` and reaches `0.0006` in the worst resample,
  a fivefold shrinkage of the ruler that inflates every index measured against
  it and produces upper limits three times the estimate. The default
  `"fixed"` holds each trait at its observed scaling, so the interval describes
  the plants' response rather than the yardstick; `"resampled"` recomputes
  everything and is both slower and heavier-tailed.
* `calculate_sri()` now keeps the columns it analyzed inside the returned
  object, which is what lets `stress_index_ci()` resample without being handed
  the table again.

# plantstressR 0.3.0

Designs with blocks, an interactive dashboard, and a correction to precision
weighting that changes published rankings.

## New features

* All of `calculate_sri()`, `stress_network()`, `stress_ordination()` and
  `validate_stress_data()` gain a `block` argument. Trait values are cleared of
  their additive block effect through a `trait ~ treatment + block` model, so a
  bench, strip or run that sat slightly better than the others no longer
  inflates the residual scale. The blocked path estimates one model per trait
  and unit, pools the residual variance over the whole trial, and reports the
  model's own t-test and degrees of freedom.
* `validate_stress_data()` checks the block layout: too few block levels, a
  treatment confined to a single block (which cannot be separated from it), and
  empty treatment-by-block cells.
* `run_plantstress_app()` returns, rewritten against the current API. The
  dashboard covers the whole workflow -- design checks, signature, ranking,
  network, ordination and tolerance indices -- with CSV upload and download.
  `shiny` is an optional dependency in `Suggests`; the app lives in
  `inst/shiny` and calls nothing but exported functions, and a test enforces
  that so it cannot drift out of step with the package again.
* `calculate_sri()` returns a new `se_sampling` column, the part of the
  standard error that comes from sampling the experimental units, with the
  scaling standard deviation treated as fixed.

## Bug fixes

* `integrated_stress_index(weights = "precision")` no longer penalizes the
  traits that responded to the stress. The full standard error of a
  standardized effect contains a `d^2` term, so weighting by `1/se^2` demoted
  precisely the traits carrying the signal: on the bundled example the weight
  was correlated `-0.92` with the size of the response, and the least affected
  trait took 42 per cent of the total weight. Weighting now uses
  `se_sampling`. **Rankings produced with `weights = "precision"` before this
  release should be recomputed.** Note that once effects are standardized,
  precision depends only on replication, so in a balanced trial with no missing
  data this scheme now coincides with `"equal"`.
* `simulate_brachiaria_stress()` no longer leaves a stray `names` attribute on
  29 of the trait columns, inherited from the genotype vector used to build the
  genotype effects.

## Data

* `brachiaria_stress` gains a `block` column with four complete blocks, and
  every trait carries an additive block shift, so the block machinery can be
  demonstrated and tested on the bundled data.

# plantstressR 0.2.0

Scope pivot. The package now answers a single question -- what is the
physiological stress signature of a plant, and how do its traits interact --
instead of covering multimodal integration, classification and reporting.

## New features

* `calculate_sri()` computes a signed, standardized Stress Response Index per
  trait, with four effect-size definitions (`"glass"`, `"cohen"`, `"hedges"`,
  `"relative"`), confidence intervals, Welch tests and multiplicity adjustment.
* `trait_directions()` resolves the stress direction of each trait so that a
  positive index always means impairment.
* `integrated_stress_index()` collapses the per-trait indices into one weighted
  index and ranks units within each stress level, with `"equal"`,
  `"precision"`, `"pca"` or manual weights. `stress_weights()` and
  `stress_contributions()` expose how a ranking was produced.
* `plot_stress_signature()` draws the signature as a diverging heat map with
  hierarchically ordered traits, or as a radar plot.
* `stress_network()` estimates a regularized partial-correlation network among
  traits and detects modules with Louvain, walktrap or fast-greedy clustering.
  `plot()` renders it as a `ggplot2` object.
* `stress_ordination()` provides the sample-level principal component view,
  optionally through `FactoMineR`.
* `stress_module_scores()` scores the network modules with the trait indices,
  turning a long trait signature into a handful of module-level statements.
* `stress_tolerance_index()` computes the classical productivity-based selection
  indices (`TOL`, `MP`, `GMP`, `HM`, `SSI`, `STI`, `YI`, `YSI`, `RDI`, `SSPI`)
  with the trial stress intensity and an average rank across them, so that the
  physiological ranking can be cross-checked against the agronomic one.
* `plot()` methods for the index objects: `plot()` on a `plantstress_sri` draws
  the signature, and `plot()` on a `plantstress_isi` draws either the ranking
  (`type = "ranking"`) or the trait decomposition behind it
  (`type = "contribution"`).
* `validate_stress_data()` runs the design checks used by every other function.

## Breaking changes

* Removed the classification, multimodal-integration, marker-selection,
  reporting and Shiny functions of 0.1.0: `classify_stress_state()`,
  `merge_modalities()`, `define_layers()`, `layer_contribution()`,
  `compare_genotypes()`, `select_markers()`, `explain_signature()`,
  `report_stress()`, `plot_multiomics_style()`, `preprocess_traits()`,
  `stress_signature()`, `trait_network()`, `run_plantstress_app()` and the
  `validate_*` pair.
* `trait_network()` is superseded by `stress_network()`, which conditions on the
  remaining traits instead of reporting marginal correlations.
* `stress_signature()` is superseded by `calculate_sri()` plus
  `stress_ordination()`.
* License changed from `MIT + file LICENSE` to `GPL (>= 3)`.
* Dropped the `caret`, `e1071`, `randomForest`, `nnet`, `shiny` and
  `visNetwork` dependencies.

## Data

* `brachiaria_stress` was regenerated. Traits now share latent physiological
  factors (photochemical, stomatal, water status, nutrition, growth), so the
  data set has a genuine module structure for the network examples.

# plantstressR 0.1.0

* First internal release.
