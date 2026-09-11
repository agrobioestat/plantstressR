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
* `calculate_sri()` now keeps the columns it analysed inside the returned
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

* `integrated_stress_index(weights = "precision")` no longer penalises the
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
