## Test environments

* local Windows 11, R 4.5.3
* GitHub Actions: ubuntu-latest, macOS-latest and windows-latest (R release)

## R CMD check results

0 errors | 0 warnings | 1 note

    * checking CRAN incoming feasibility ... NOTE
    Maintainer: 'Joao Everthon da Silva Ribeiro <everthonribeiro10@gmail.com>'
    New submission

This is the expected new-submission note.

## Notes for the reviewer

* This is a major scope change relative to 0.1.0. The package no longer provides
  classification or multimodal integration; it now covers only multivariate
  stress phenotyping. The removed functions are listed in NEWS.md.
* The license changed from `MIT + file LICENSE` to `GPL (>= 3)` to match the
  other packages of the maintainer. The maintainer is the sole copyright holder
  of all the code in the package.
* `FactoMineR` is the only heavy optional dependency. It sits in Suggests, its
  single use in `stress_ordination()` is guarded by `requireNamespace()`, and the
  function falls back to `stats::prcomp()` when it is absent.
* All examples, tests and the vignette run on the bundled synthetic data set and
  take a few seconds in total.
