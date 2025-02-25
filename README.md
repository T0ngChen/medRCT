
<!-- README.md is generated from README.Rmd. Please edit that file -->

# medRCT

<!-- badges: start -->

[![MIT
license](http://img.shields.io/badge/license-MIT-brightgreen.svg)](http://opensource.org/licenses/MIT)
[![R-CMD-check](https://github.com/T0ngChen/medRCT/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/T0ngChen/medRCT/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/T0ngChen/medRCT/graph/badge.svg)](https://app.codecov.io/gh/T0ngChen/medRCT)
<!-- badges: end -->

> Causal mediation analysis estimating interventional effects mapped to
> a target trial

The R package `medRCT` for causal mediation analysis supports the
estimation of interventional effects (VanderWeele, Vansteelandt, and
Robins 2014), specifically interventional effects that are defined such
that they map explicitly to a “target trial” (Hernán and Robins 2016),
as recently proposed by Moreno-Betancur et al. (2021). In the target
trial, the treatment strategies are specified to reflect hypothetical
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

Using a simulated dataset based on a published case study from the
Longitudinal Study of Australian Children (Goldfeld et al. 2023), we
illustrate how to use `medRCT` to estimate the interventional effects
that emulate a target trial. Specifically, we aim to estimate the
difference in expected outcome (risk of child mental health problems)
under exposure (low family socioeconomic position) with versus without a
hypothetical intervention that individually shifts the distribution of
each mediator (parental mental health and preschool attendance) to the
levels in the unexposed (high family socioeconomic position), while
accounting for baseline confounders, an intermediate (exposure-induced)
mediator-outcome confounder (family stressful life events), and
correlations amongst mediators.

We begin by loading the library and dataset, and defining the confounder
vector.

``` r
# Load the medRCT package
library(medRCT)

# Set a seed for reproducibility
set.seed(2025)

# Display the first few rows of the dataset
head(LSACdata)
#>   child_sex child_atsi mat_cob mat_engl mat_age sep fam_stress parent_mh preschool_att child_mh
#> 1         0          1       0        0       1   0          0         0             1        0
#> 2        NA          0       0        0      NA   0         NA        NA             0        0
#> 3        NA          0       0        0      NA   0         NA        NA             0        1
#> 4        NA          0       0        0      NA   0         NA        NA             0        0
#> 5         1          0       0        0       1   1          0         0             0        1
#> 6         1          0       0        0       1   0          1         1             0        1
#>   child_SDQscore
#> 1       8.924660
#> 2       7.349826
#> 3      12.824643
#> 4       6.611369
#> 5      10.329341
#> 6      13.552515

# Define confounders for the analysis
confounders <- c("child_sex", "child_atsi", "mat_cob", "mat_engl", "mat_age")
```

Next we run the analyses, estimating interventional effects for a
hypothetical intervention that shifts the distribution of each mediator
individually. **Note 1:** the dataset has missing data. Incomplete
records are by default deleted before the analysis. **Note 2:** It is
recommended to perform the analysis with at least 200 Monte Carlo
simulations by setting `mcsim = 200`. For illustration purposes, we use
`mcsim = 50`, which takes approximately 90 seconds to run.

``` r
# Estimate interventional effects for a hypothetical intervention
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
  mcsim = 50                          
)
#> Conducting complete case analysis, 2499 observations were excluded due to missing data.
#> Note: It is recommended to run analysis with no fewer than 200 Monte Carlo simulations.

# Summarise the results
summary(med_res)
#> 
#> Estimated interventional indirect effect: 
#> 
#>                      Estimate Std. Error  CI Lower  CI Upper p-value
#> IIE_1 (p_trt - p_1)  0.011155   0.004181  0.002814  0.019203  0.0076
#> IIE_2 (p_trt - p_2) -0.000763   0.002501 -0.005443  0.004362  0.7604
#> TCE (p_trt - p_ctr)  0.128669   0.024554  0.082420  0.178668 1.6e-07
#> 
#> Estimated interventional direct effect: 
#> 
#>                     Estimate Std. Error CI Lower CI Upper p-value
#> IDE_1 (p_1 - p_ctr)   0.1175     0.0247   0.0712   0.1679 1.9e-06
#> IDE_2 (p_2 - p_ctr)   0.1294     0.0244   0.0833   0.1789 1.1e-07
#> 
#> Estimated expected outcome in each trial arm:
#> 
#>       Estimate Std. Error CI Lower CI Upper p-value
#> p_1     0.3302     0.0225   0.2872   0.3755  <2e-16
#> p_2     0.3421     0.0221   0.2995   0.3862  <2e-16
#> p_ctr   0.2127     0.0100   0.1922   0.2315  <2e-16
#> p_trt   0.3413     0.0223   0.2987   0.3860  <2e-16
#> 
#> Sample Size: 2608 
#> 
#> Simulations: 50
```

Based on the estimated interventional effect (IIE_1), a hypothetical
intervention improving the mental health of parents of children from
families with low socioeconomic position to the levels of those from
families with high socioeconomic position could potentially prevent 1
per 100 cases of child mental health problems. Meanwhile, the effect of
a hypothetical intervention on preschool attendance (IIE_2) is
negligible.

For detailed guidance on using the package to handle more complex
scenarios, please refer to the
[vignette](https://t0ngchen.github.io/medRCT/articles/intro.html).

## Citation

For work involving the `medRCT` R package, please cite the following:

    @software{Chen2024medRCT,
       author = {Tong Chen and Margarita Moreno-Betancur and S. Ghazaleh Dashti},
       title = {medRCT: Causal mediation analysis estimating interventional effects mapped to a target trial},
       year  = {2025},
       url = {https://t0ngchen.github.io/medRCT/},
       note = {R package version 0.1.0}
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
