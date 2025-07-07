
<!-- README.md is generated from README.Rmd. Please edit that file -->

# medRCT <a href="https://t0ngchen.github.io/medRCT"><img src="man/figures/logo.png" align="right" height="100" alt="medRCT website" style="max-width: 100%; height: 138;" /></a>

<!-- badges: start -->

[![MIT
license](http://img.shields.io/badge/license-MIT-brightgreen.svg)](http://opensource.org/licenses/MIT)
[![R-CMD-check](https://github.com/T0ngChen/medRCT/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/T0ngChen/medRCT/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/T0ngChen/medRCT/graph/badge.svg)](https://app.codecov.io/gh/T0ngChen/medRCT)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.08063/status.svg)](https://doi.org/10.21105/joss.08063)
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

## Statement of need:

Causal mediation analysis generally seeks to investigate the extent to
which the causal effect of an exposure on an outcome is mediated through
intermediate variables. Natural (in)direct effects (Robins and Greenland
1992; Pearl 2001) were initially proposed as the estimands of interest
in these analyses. Natural effects are defined based on cross-world
counterfactuals (Robins and Richardson 2011) and their identifiability
relies on a cross-world independence assumption. Given their reliance on
cross-world counterfactuals, these effects have been criticized for not
capturing the effects of interventions or policy measures that could be
conducted in the real world (Naimi, Kaufman, and MacLehose 2014).
Further, the independence assumption required can never be guaranteed,
even in an experiment (Robins and Richardson 2011; Didelez, Dawid, and
Geneletti 2006), and it renders the estimands unidentifiable in the
common settings of exposure-induced mediator-outcome confounding and
multiple mediators (Avin, Shpitser, and Pearl 2005; VanderWeele,
Vansteelandt, and Robins 2014; Vansteelandt and VanderWeele 2012).
However, in the context of multiple mediators, certain path-specific
natural effects, also defined in terms of cross-world counterfactuals,
can still be identified and may be of substantive interest (VanderWeele
and Vansteelandt 2014).

Interventional effects have been proposed as an alternative to address
these limitations. Firstly, these effects can be shown to map to a
hypothetical randomized trial that evaluates the impact of hypothetical
interventions shifting the distribution of the mediators
(Moreno-Betancur and Carlin 2018). Secondly, interventional effects
remain identifiable in the presence of exposure-induced mediator-outcome
confounding and multiple interrelated mediators of interest.

The `medRCT` package implements the estimation of interventional effects
that are defined explicitly as effects in a hypothetical randomized
trial (the target trial) , as proposed by Moreno-Betancur et al. (2021).
This assists with clarifying the research question and ensuring that the
study findings are meaningful and relevant to policy and practice. In
the target trial, the treatment strategies are specified to reflect
hypothetical interventions targeting and thus shifting the joint
mediator distribution. The `medRCT` package implements the estimation of
interventional effects that correspond to effects of hypothetical
interventions which:

1.  shift the joint distribution of all mediators under exposure to that
    under no exposure,

2.  shift the distribution of a specific mediator under exposure, given
    confounders, to match the corresponding distribution under no
    exposure, independent of and without considering flow-on effects on
    other mediators,

3.  shift the distribution of a specific mediator under exposure, given
    confounders, to match the corresponding distribution under no
    exposure, while considering flow-on effects on causally descendant
    mediators.

`medRCT` estimates these interventional effects using a Monte Carlo
simulation-based g-computation approach. It should be noted that this
method can be computationally intensive and is sensitive to model
misspecification, as all nuisance parameters are estimated via
restrictive parametric models.

Researchers should consider using `medRCT` when their ultimate goal for
conducting mediation analysis is to examine the effects of hypothetical
interventions targeting multiple, potentially interdependent mediators.

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
#>   child_sex child_atsi mat_cob mat_engl mat_age sep fam_stress parent_mh
#> 1         0          1       0        0       1   0          0         0
#> 2        NA          0       0        0      NA   0         NA        NA
#> 3        NA          0       0        0      NA   0         NA        NA
#> 4        NA          0       0        0      NA   0         NA        NA
#> 5         1          0       0        0       1   1          0         0
#> 6         1          0       0        0       1   0          1         1
#>   preschool_att child_mh child_SDQscore
#> 1             1        0       8.924660
#> 2             0        0       7.349826
#> 3             0        1      12.824643
#> 4             0        0       6.611369
#> 5             0        1      10.329341
#> 6             0        1      13.552515

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
#> 
#> Effect Measure: Risk Difference 
#> Results are based on all 100 bootstrap samples.
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

    @article{Chen2025medRCT,
      author = {Tong Chen and S. Ghazaleh Dashti and Margarita Moreno-Betancur},
      title = {{medRCT}: Causal mediation analysis estimating interventional effects mapped to a target trial in {R}},
      year = {2025},
      doi = {10.21105/joss.08063},
      url = {https://doi.org/10.21105/joss.08063},
      journal = {Journal of Open Source Software},
      volume = {10},
      number = {110},
      pages = {8063},
      publisher = {The Open Journal}
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

<div id="ref-Avin2005" class="csl-entry">

Avin, Chen, Ilya Shpitser, and Judea Pearl. 2005. “Identifiability of
Path-Specific Effects.” In *Proceedings of the 19th International Joint
Conference on Artificial Intelligence*, 357–63. IJCAI’05. San Francisco,
CA, USA: Morgan Kaufmann Publishers Inc.
<https://dl.acm.org/doi/10.5555/1642293.1642350>.

</div>

<div id="ref-Didelez2006" class="csl-entry">

Didelez, Vanessa, Philip Dawid, and Sara Geneletti. 2006. “Direct and
Indirect Effects of Sequential Treatments.” In *Proceedings of the
Twenty-Second Conference on Uncertainty in Artificial Intelligence*,
138–46. UAI’06. Arlington, Virginia, USA: AUAI Press.
<https://dl.acm.org/doi/10.5555/3020419.3020437>.

</div>

<div id="ref-Goldfeld2023" class="csl-entry">

Goldfeld, Sharon, Margarita Moreno-Betancur, Sarah Gray, Shuaijun Guo,
Marnie Downes, Elodie O’Connor, Francisco Azpitarte, et al. 2023.
“Addressing Child Mental Health Inequities Through Parental Mental
Health and Preschool Attendance.” *Pediatrics* 151 (5): e2022057101.
<https://doi.org/10.1542/peds.2022-057101>.

</div>

<div id="ref-hernan2016using" class="csl-entry">

Hernán, Miguel A, and James M. Robins. 2016. “Using Big Data to Emulate
a Target Trial When a Randomized Trial Is Not Available.” *American
Journal of Epidemiology* 183 (8): 758–64.
<https://doi.org/10.1093/aje/kwv254>.

</div>

<div id="ref-MorenoBetancur2018" class="csl-entry">

Moreno-Betancur, Margarita, and John B. Carlin. 2018. “Understanding
Interventional Effects: A More Natural Approach to Mediation Analysis.”
*Epidemiology* 29 (5): 614–17.
<https://doi.org/10.1097/EDE.0000000000000866>.

</div>

<div id="ref-Moreno2021Mediation" class="csl-entry">

Moreno-Betancur, Margarita, Paul Moran, Denise Becker, George C. Patton,
and John B. Carlin. 2021. “Mediation Effects That Emulate a Target
Randomised Trial: Simulation-Based Evaluation of Ill-Defined
Interventions on Multiple Mediators.” *Statistical Methods in Medical
Research* 30 (6): 1395–1412. <https://doi.org/10.1177/0962280221998409>.

</div>

<div id="ref-Naimi2014Mediation" class="csl-entry">

Naimi, Ashley I., Jay S. Kaufman, and Richard F. MacLehose. 2014.
“Mediation Misgivings: Ambiguous Clinical and Public Health
Interpretations of Natural Direct and Indirect Effects.” *International
Journal of Epidemiology* 43 (5): 1656–61.
<https://doi.org/10.1093/ije/dyu107>.

</div>

<div id="ref-Pearl2001" class="csl-entry">

Pearl, Judea. 2001. “Direct and Indirect Effects.” In *Proceedings of
the 17th Conference on Uncertainty in Artificial Intelligence*, 411–20.
UAI’01. San Francisco, CA, USA: Morgan Kaufmann Publishers Inc.
<https://dl.acm.org/doi/10.5555/2074022.2074073>.

</div>

<div id="ref-Robins1992" class="csl-entry">

Robins, James M., and Sander Greenland. 1992. “Identifiability and
Exchangeability for Direct and Indirect Effects.” *Epidemiology* 3 (2):
143–55. <https://doi.org/10.1097/00001648-199203000-00013>.

</div>

<div id="ref-Robins2011" class="csl-entry">

Robins, James M., and Thomas S. Richardson. 2011. “Alternative Graphical
Causal Models and the Identification of Direct Effects.” In *Causality
and Psychopathology: Finding the Determinants of Disorders and Their
Cures*. Oxford University Press.
<https://doi.org/10.1093/oso/9780199754649.003.0011>.

</div>

<div id="ref-VanderWeele2014" class="csl-entry">

VanderWeele, Tyler J., and Stijn Vansteelandt. 2014. “Mediation Analysis
with Multiple Mediators.” *Epidemiologic Methods* 2 (1): 95–115.
<https://doi.org/10.1515/em-2012-0010>.

</div>

<div id="ref-vanderweele2014effect" class="csl-entry">

VanderWeele, Tyler J., Stijn Vansteelandt, and James M. Robins. 2014.
“Effect Decomposition in the Presence of an Exposure-Induced
Mediator-Outcome Confounder.” *Epidemiology* 25 (2): 300–306.
<https://doi.org/10.1097/EDE.0000000000000034>.

</div>

<div id="ref-Vansteelandt2012Natural" class="csl-entry">

Vansteelandt, Stijn, and Tyler J. VanderWeele. 2012. “Natural Direct and
Indirect Effects on the Exposed: Effect Decomposition Under Weaker
Assumptions.” *Biometrics* 68 (4): 1019–27.
<https://doi.org/10.1111/j.1541-0420.2012.01777.x>.

</div>

</div>
