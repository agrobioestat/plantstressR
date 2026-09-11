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
