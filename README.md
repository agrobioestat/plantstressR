---
output: github_document
---



# plantstressR

<!-- badges: start -->
[![R-CMD-check](https://github.com/agrobioestat/plantstressR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/agrobioestat/plantstressR/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

`plantstressR` answers one question:

> **What is the physiological stress signature of a plant, and how do its traits
> interact?**

You bring one tidy table of physiological traits measured under a control
treatment and one or more stress treatments. The package does the rest: it
puts every trait on a common signed scale, integrates them into a single
rankable index, draws the signature, and maps the network of traits that
respond together.

## Installation

```r
# install.packages("remotes")
remotes::install_github("agrobioestat/plantstressR")
```

## The workflow

| Function | What it gives you |
|---|---|
| `validate_stress_data()` | Design checks before anything is computed |
| `calculate_sri()` | A signed, standardized Stress Response Index per trait |
| `integrated_stress_index()` | One weighted index per genotype, plus a ranking |
| `plot_stress_signature()` | The signature as a heat map or a radar plot |
| `stress_network()` | A regularized partial-correlation network with trait modules |
| `stress_module_scores()` | The modules scored with the trait indices |
| `stress_tolerance_index()` | The classical selection indices (`STI`, `SSI`, `GMP`, `TOL`, ...) |
| `stress_ordination()` | The sample-level view of the trait space |

## Example


``` r
library(plantstressR)

traits <- c("Fv_Fm", "PIabs", "DI0_RC", "A", "gs", "Ci", "RWC", "SPAD", "shoot_biomass")

sri <- calculate_sri(
  brachiaria_stress,
  treatment = "drought_level",
  control = "control",
  traits = traits,
  by = "genotype",
  verbose = FALSE
)

summary(sri)
#> # A tibble: 8 × 8
#>   unit  group    n_traits n_significant mean_sri mean_abs_sri most_impaired
#>   <chr> <chr>       <int>         <int>    <dbl>        <dbl> <chr>        
#> 1 G1    moderate        9             8     3.19         3.20 RWC          
#> 2 G1    severe          9             8     6.23         6.23 RWC          
#> 3 G2    moderate        9             9     3.22         3.22 RWC          
#> 4 G2    severe          9             9     6.08         6.08 shoot_biomass
#> 5 G3    moderate        9             8     3.55         3.55 PIabs        
#> 6 G3    severe          9             9     7.05         7.05 PIabs        
#> 7 G4    moderate        9             8     2.87         2.87 A            
#> 8 G4    severe          9             9     5.98         5.98 A            
#> # ℹ 1 more variable: most_resilient <chr>
```

Positive values always mean *stress impaired this trait*, whichever direction
the trait moves. `Fv/Fm` falling and `DI0/RC` rising are both damage, and
`trait_directions()` keeps that bookkeeping straight.

Rank the genotypes:


``` r
integrated_stress_index(sri, weights = "precision")
#> <plantstress_isi>
#>   Weighting: precision
#>   Aggregate: mean
#>   Ranking:   tolerance (rank 1 = most tolerant)
#> # A tibble: 8 × 6
#>   unit  group      isi isi_scaled n_traits  rank
#>   <chr> <chr>    <dbl>      <dbl>    <int> <dbl>
#> 1 G1    moderate  1.84        0          9     1
#> 2 G4    moderate  1.98       33.4        9     2
#> 3 G3    moderate  2.23       91.7        9     3
#> 4 G2    moderate  2.27      100          9     4
#> 5 G1    severe    3.86        0          9     1
#> 6 G2    severe    4.18       37.4        9     2
#> 7 G4    severe    4.33       55.3        9     3
#> 8 G3    severe    4.72      100          9     4
```

Draw the signature:


``` r
plot_stress_signature(sri, units = "G1")
```

<div class="figure">
<img src="man/figures/README-signature-1.png" alt="plot of chunk signature" width="100%" />
<p class="caption">plot of chunk signature</p>
</div>

Map how the traits interact:


``` r
net <- stress_network(
  brachiaria_stress,
  traits = traits,
  treatment = "drought_level",
  level = c("moderate", "severe"),
  sri = sri
)

net
#> <plantstress_network>
#>   Method:      partial (lambda = 0)
#>   Traits:      9 | complete observations: 56
#>   Edges kept:  7 of 36 (|r| >= 0.1)
#>   Modules:     3 (modularity = 0.536)
#>    - module 1 [hub: PIabs] PIabs, Fv_Fm, DI0_RC
#>    - module 2 [hub: gs] gs, A, Ci
#>    - module 3 [hub: RWC] RWC, SPAD, shoot_biomass
```


``` r
plot(net, color_by = "sri")
```

<div class="figure">
<img src="man/figures/README-network-plot-1.png" alt="plot of chunk network-plot" width="100%" />
<p class="caption">plot of chunk network-plot</p>
</div>

And score those modules with the signature:


``` r
stress_module_scores(net, sri)
#> # A tibble: 24 × 8
#>    unit  group    module label n_traits hub   score traits                  
#>    <chr> <chr>     <int> <chr>    <int> <chr> <dbl> <chr>                   
#>  1 G1    moderate      3 M3           3 RWC    4.34 RWC, shoot_biomass, SPAD
#>  2 G1    moderate      1 M1           3 PIabs  3.56 DI0_RC, Fv_Fm, PIabs    
#>  3 G1    moderate      2 M2           3 gs     1.68 A, Ci, gs               
#>  4 G1    severe        3 M3           3 RWC    8.20 RWC, shoot_biomass, SPAD
#>  5 G1    severe        1 M1           3 PIabs  6.75 DI0_RC, Fv_Fm, PIabs    
#>  6 G1    severe        2 M2           3 gs     3.73 A, Ci, gs               
#>  7 G2    moderate      3 M3           3 RWC    3.71 RWC, shoot_biomass, SPAD
#>  8 G2    moderate      1 M1           3 PIabs  3.28 DI0_RC, Fv_Fm, PIabs    
#>  9 G2    moderate      2 M2           3 gs     2.68 A, Ci, gs               
#> 10 G2    severe        3 M3           3 RWC    6.89 RWC, shoot_biomass, SPAD
#> # ℹ 14 more rows
```

Cross-check against the classical selection indices:


``` r
stress_tolerance_index(
  brachiaria_stress,
  trait = "shoot_biomass",
  treatment = "drought_level",
  control = "control",
  by = "genotype",
  indices = c("STI", "GMP", "SSI")
)
#> <plantstress_sti>
#>   Trait:   shoot_biomass
#>   Units:   4 (genotype)
#>   Indices: STI, GMP, SSI
#>   Stress intensity [moderate]: 0.272
#>   Stress intensity [severe]: 0.565
#>   rank_overall 1 = most tolerant on the average of the indices
#> # A tibble: 8 × 10
#>   unit  group       Yp    Ys stress_intensity   STI   GMP   SSI rank_mean
#>   <chr> <chr>    <dbl> <dbl>            <dbl> <dbl> <dbl> <dbl>     <dbl>
#> 1 G3    moderate  44.8  33.1            0.272 0.742  38.5 0.959      1.67
#> 2 G1    moderate  45.4  32.7            0.272 0.745  38.6 1.03       2   
#> 3 G2    moderate  44.9  32.5            0.272 0.730  38.2 1.02       3   
#> 4 G4    moderate  43.6  31.9            0.272 0.697  37.3 0.992      3.33
#> 5 G4    severe    43.6  21.0            0.565 0.459  30.3 0.918      1   
#> 6 G1    severe    45.4  19.9            0.565 0.452  30.1 0.995      2   
#> 7 G2    severe    44.9  19.4            0.565 0.435  29.5 1.01       3   
#> 8 G3    severe    44.8  17.5            0.565 0.391  28.0 1.08       4   
#> # ℹ 1 more variable: rank_overall <int>
```

## Learn more

```r
vignette("plantstressR")
```

## Citation

```r
citation("plantstressR")
```

## License

GPL (>= 3)
