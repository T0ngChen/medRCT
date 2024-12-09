
<!-- README.md is generated from README.Rmd. Please edit that file -->

# medRCT

<!-- badges: start -->

[![MIT
license](http://img.shields.io/badge/license-MIT-brightgreen.svg)](http://opensource.org/licenses/MIT)
[![R-CMD-check](https://github.com/T0ngChen/medRCT/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/T0ngChen/medRCT/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

> Causal **Med**iation Analysis for Estimating Interventional Effects
> that Emulate a **R**andomised **C**ontrolled **T**rial

The R package `medRCT` supports the estimation of interventional effects
(VanderWeele, Vansteelandt, and Robins 2014) that are defined explicitly
in terms of a target trial (Hernán and Robins 2016), where the treatment
strategies are specified to reflect the hypothetical interventions,
which are encoded by shifts in the joint distribution of the mediators
(Moreno-Betancur et al. 2021). `medRCT` can accommodates multiple
mediators in the presence of exposure-induced mediator-outcome
confounders.

## Installation

The `medRCT` package is not yet available on CRAN. You can install the
latest stable version from [GitHub](https://github.com/T0ngChen/medRCT)
using the following command:

``` r
remotes::install_github("T0ngChen/medRCT")
```

## Example

We here illustrate how to use `medRCT` to estiamte interventional
indirect effects that emulate a target trail. Specifically, we consider
a hypothetical intervention that shifts the distribution of each
mediator individually. It allows for the investigation of causal
pathways and the quantification of indirect effects to directly address
real-world research questions. Consider the following example using
simulated data based on the Longitudinal Study of Australian Children
(Sanson and Johnstone 2004):

``` r
# Load the medRCT package
library(medRCT)

# Set a seed for reproducibility
set.seed(2024)

# Display the first few rows of the dataset
head(LSACdata)
#>   child_sex child_atsi mat_cob mat_engl mat_age sep fam_stress parent_mh
#> 1         1          0       0        0       1   1          0         0
#> 2         0          0       0        0       0   0          0         0
#> 3         1          0       1        1       0   0          0         1
#> 4         0          0       0        0       0   0          0         0
#> 5         0          0       0        0       0   0          0         0
#> 6         0          0       0        0       0   1          0         0
#>   preschool_att child_mh
#> 1             1        0
#> 2             1        0
#> 3             0        0
#> 4             1        0
#> 5             0        0
#> 6             0        0

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
#> IIE_1 (p_trt - p_1)  0.010075   0.002731  0.005061  0.015409 0.00022
#> IIE_2 (p_trt - p_2)  0.000549   0.001733 -0.002301  0.004548 0.75129
#> TCE (p_trt - p_ctr)  0.106266   0.016395  0.076784  0.141251 9.1e-11
#> 
#> Estimated expected outcome in each trial arm:
#> 
#>       Estimate Std. Error CI Lower CI Upper p-value
#> p_trt  0.33345    0.01406  0.31244  0.36755  <2e-16
#> p_ctr  0.22719    0.00765  0.21490  0.24367  <2e-16
#> p_1    0.32338    0.01412  0.30067  0.35502  <2e-16
#> p_2    0.33290    0.01423  0.31124  0.36460  <2e-16
#> 
#> Sample Size: 5107 
#> 
#> Simulations: 100
```

For detailed guidance on using the package to handle more complex
scenarios, please refer to the
[vignette](https://t0ngchen.github.io/medRCT/).

## Citation

For work involving the `medRCT` R package, please cite the following:

    @software{Chen2024medRCT,
       author = {Tong Chen and Margarita Moreno-Betancur and Ghazaleh Dashti},
       title = {medRCT: Causal Mediation Analysis Estimating Interventional Effects Mapped to a Target Trial},
       year  = {2024},
       url = {https://t0ngchen.github.io/medRCT/},
       note = {R package version 0.0.0.9020}
       }
    @article{Moreno2021Mediation,
       author = {Margarita Moreno-Betancur and Paul Moran and Denise Becker and George C Patton and John B Carlin},
       title = {Mediation effects that emulate a target randomised trial: Simulation-based evaluation of ill-defined interventions on multiple mediators},
       journal = {Statistical Methods in Medical Research},
       volume = {30},
       number = {6},
       pages = {1395-1412},
       year = {2021}
       }    

## References

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-hernan2016using" class="csl-entry">

Hernán, Miguel A, and James M Robins. 2016. “Using Big Data to Emulate a
Target Trial When a Randomized Trial Is Not Available.” *American
Journal of Epidemiology* 183 (8): 758–64.

</div>

<div id="ref-Moreno2021Mediation" class="csl-entry">

Moreno-Betancur, Margarita, Paul Moran, Denise Becker, George C Patton,
and John B Carlin. 2021. “Mediation Effects That Emulate a Target
Randomised Trial: Simulation-Based Evaluation of Ill-Defined
Interventions on Multiple Mediators.” *Statistical Methods in Medical
Research* 30 (6): 1395–1412.

</div>

<div id="ref-Sanson2004GrowingUI" class="csl-entry">

Sanson, Ann V., and Robert E Johnstone. 2004. “Growing up in Australia
Takes Its First Steps.” *Family Matters* 67: 46–53.

</div>

<div id="ref-vanderweele2014effect" class="csl-entry">

VanderWeele, Tyler J, Stijn Vansteelandt, and James M Robins. 2014.
“Effect Decomposition in the Presence of an Exposure-Induced
Mediator-Outcome Confounder.” *Epidemiology* 25 (2): 300–306.

</div>

</div>
