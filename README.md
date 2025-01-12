
<!-- README.md is generated from README.Rmd. Please edit that file -->

# medRCT

<!-- badges: start -->

[![MIT
license](http://img.shields.io/badge/license-MIT-brightgreen.svg)](http://opensource.org/licenses/MIT)
[![R-CMD-check](https://github.com/T0ngChen/medRCT/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/T0ngChen/medRCT/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/T0ngChen/medRCT/graph/badge.svg)](https://app.codecov.io/gh/T0ngChen/medRCT)
<!-- badges: end -->

> Causal Mediation Analysis Estimating Interventional Effects Mapped to
> A Target Trial

The R package `medRCT` for causal mediation analysis supports the
estimation of interventional effects (VanderWeele, Vansteelandt, and
Robins 2014), specifically interventional effects that are defined
explicitly in terms of a “target trial” (Hernán and Robins 2016), as
proposed recently by Moreno-Betancur et al. (2021). In the target trial,
the treatment strategies are specified to reflect hypothetical
interventions targeting and thus shifting the joint distribution of the
mediators. `medRCT` can accommodate any number of potentially correlated
mediators, including mediators that are not of primary interest but that
are intermediate (exposure-induced) mediator-outcome confounders.

## Installation

The `medRCT` package is not yet available on CRAN. You can install the
latest stable version from [GitHub](https://github.com/T0ngChen/medRCT)
using the following command:

``` r
remotes::install_github("T0ngChen/medRCT")
```

## Example

Here we illustrate how to use `medRCT` to estimate interventional
indirect effects that emulate a target trial using a simulated dataset
based on a published case study from the Longitudinal Study of
Australian Children (Goldfeld et al. 2023). Specifically, we aim to
estimate the difference in expected outcome (risk of mental health
problems) under exposure (low socioeconomic position) with versus
without a hypothetical intervention that individually shifts the
distribution of each mediator (parental mental health and preschool
attendance) to the levels in the unexposed, while accounting for
baseline confounders, an intermediate (exposure-induced)
mediator-outcome confounder and correlations amongst mediators.

``` r
# Load the medRCT package
library(medRCT)

# Set a seed for reproducibility
set.seed(2024)

# Display the first few rows of the dataset
head(LSACdata)
#> # A tibble: 6 × 10
#>   child_sex child_atsi mat_cob mat_engl mat_age   sep fam_stress parent_mh preschool_att child_mh
#>       <dbl>      <dbl>   <dbl>    <dbl>   <dbl> <dbl>      <dbl>     <dbl>         <dbl>    <dbl>
#> 1         1          0       0        0       1     1          0         0             1        0
#> 2         0          0       0        0       0     0          0         0             1        0
#> 3         1          0       1        1       0     0          0         1             0        0
#> 4         0          0       0        0       0     0          0         0             1        0
#> 5         0          0       0        0       0     0          0         0             0        0
#> 6         0          0       0        0       0     1          0         0             0        0

# Define confounders for the analysis
confounders <- c("child_sex", "child_atsi", "mat_cob", "mat_engl", "mat_age")

# Estimate interventional indirect effects for a hypothetical intervention
# that shifts the distribution of each mediator individually
med_res <- medRCT(
  dat = LSACdata,                      
  exposure = "sep",                    
  outcome = "child_mh",                
  mediators = c("parent_mh", "preschool_att"), 
  intermediate_confs = "fam_stress",  # intermediate confounders 
  confounders = confounders,           
  bootstrap = TRUE,                    
  intervention_type = "shift_k",       
  mcsim = 100                          
)

# Summarise the results
summary(med_res)
#> 
#> Estimated interventional effect: 
#> 
#>                      Estimate Std. Error  CI Lower  CI Upper p-value
#> IIE_1 (p_trt - p_1)  0.010075   0.002731  0.004813  0.015517 0.00022
#> IIE_2 (p_trt - p_2)  0.000549   0.001733 -0.002863  0.003930 0.75129
#> TCE (p_trt - p_ctr)  0.106266   0.016395  0.073678  0.137945 9.1e-11
#> 
#> Estimated expected outcome in each trial arm:
#> 
#>       Estimate Std. Error CI Lower CI Upper p-value
#> p_1    0.32338    0.01412  0.29483  0.35017  <2e-16
#> p_2    0.33290    0.01423  0.30425  0.36001  <2e-16
#> p_ctr  0.22719    0.00765  0.21187  0.24184  <2e-16
#> p_trt  0.33345    0.01406  0.30510  0.36023  <2e-16
#> 
#> Sample Size: 5107 
#> 
#> Simulations: 100
```

Based on the estimated interventional effect (IIE_1), a hypothetical
intervention improving the mental health of parents of children in low
socioeconomic position to the levels in those with high socioeconomic
position could potentially prevent 1 per 100 cases of child mental
health problems. Meanwhile, the effect of a hypothetical intervention on
preschool attendance (IIE_2) is negligible.

For detailed guidance on using the package to handle more complex
scenarios, please refer to the
[vignette](https://t0ngchen.github.io/medRCT/articles/intro.html).

## Citation

For work involving the `medRCT` R package, please cite the following
works:

    @software{Chen2024medRCT,
       author = {Tong Chen and Margarita Moreno-Betancur and S. Ghazaleh Dashti},
       title = {medRCT: Causal Mediation Analysis Estimating Interventional Effects Mapped to a Target Trial},
       year  = {2025},
       url = {https://t0ngchen.github.io/medRCT/},
       note = {R package version 0.0.0.9080}
       }
    @article{Moreno2021Mediation,
       author={Margarita Moreno-Betancur and Paul Moran and Denise Becker and George C Patton and John B Carlin},
       title={Mediation effects that emulate a target randomised trial: Simulation-based evaluation of ill-defined interventions on multiple mediators},
       journal={Statistical Methods in Medical Research},
       volume={30},
       number={6},
       pages={1395--1412},
       year={2021},
       URL={https://doi.org/10.1177/0962280221998409},
       doi={10.1177/0962280221998409},
       publisher={SAGE Publications Ltd}
       }    

## References

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-Goldfeld2023" class="csl-entry">

Goldfeld, Sharon, Margarita Moreno-Betancur, Sarah Gray, Shuaijun Guo,
Marnie Downes, Elodie O’Connor, Francisco Azpitarte, et al. 2023.
“Addressing Child Mental Health Inequities Through Parental Mental
Health and Preschool Attendance.” *Pediatrics* 151 (5): e2022057101.
<https://doi.org/10.1542/peds.2022-057101>.

</div>

<div id="ref-hernan2016using" class="csl-entry">

Hernán, Miguel A, and James M Robins. 2016. “Using Big Data to Emulate a
Target Trial When a Randomized Trial Is Not Available.” *American
Journal of Epidemiology* 183 (8): 758–64.
<https://doi.org/10.1093/aje/kwv254>.

</div>

<div id="ref-Moreno2021Mediation" class="csl-entry">

Moreno-Betancur, Margarita, Paul Moran, Denise Becker, George C Patton,
and John B Carlin. 2021. “Mediation Effects That Emulate a Target
Randomised Trial: Simulation-Based Evaluation of Ill-Defined
Interventions on Multiple Mediators.” *Statistical Methods in Medical
Research* 30 (6): 1395–1412. <https://doi.org/10.1177/0962280221998409>.

</div>

<div id="ref-vanderweele2014effect" class="csl-entry">

VanderWeele, Tyler J, Stijn Vansteelandt, and James M Robins. 2014.
“Effect Decomposition in the Presence of an Exposure-Induced
Mediator-Outcome Confounder.” *Epidemiology* 25 (2): 300–306.
<https://doi.org/10.1097/EDE.0000000000000034>.

</div>

</div>
