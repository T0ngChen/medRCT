# Causal Mediation Analysis Estimating Interventional Effects Mapped to a Target Trial

## Background

Causal mediation analysis methods are commonly used in research studies
to explore the causal pathways linking exposures to outcomes.
Traditional mediation analysis aims to decompose the total effect of an
exposure on an outcome into the natural direct and indirect effects,
with the natural indirect effect capturing the exposure effect acting
through the mediator. However, the goal of research studies, especially
in health and medical research, is often to elucidate the extent to
which hypothetical interventions targeting mediators could counter
exposure effects. This suggests that a shift in focus is warranted, from
effect decomposition to explicitly evaluating hypothetical mediator
interventions that may inform future policy, practice and intervention
development.

Interventional effects (VanderWeele, Vansteelandt, and Robins 2014) are
mediation estimands that are more aligned with this interventional focus
while addressing key identifiability challenges associated with natural
effects (Robins and Greenland 1992). Indeed, interventional effects have
been shown to implicitly emulate effects in target trials (Hernán and
Robins 2016) assessing the impact of distributional shifts in mediators
(Moreno-Betancur and Carlin 2018). Based on this concept,
Moreno-Betancur et al. (2021) proposed definitions of interventional
effects that map explicitly to a target trial assessing hypothetical
mediator interventions of specific interest, making them particularly
useful for policy-relevant evaluations in applied research as
illustrated in several published examples (Dashti et al. 2022; Goldfeld
et al. 2023; Afshar et al. 2024).

## Interventional effects: definition and estimation methods

Moreno-Betancur et al. (2021) provided definitions for interventional
effects based on three types of hypothetical mediator interventions. In
this vignette, we describe these estimands and their corresponding
mediator interventions, as implemented in the `medRCT` R package. The
type of mediator intervention can be specified using the argument
`intervention_type`. The default, `intervention_type = "all"`, estimates
all three types of interventional effects, which are:

- **`shift_all`**:  
  This estimand corresponds to an intervention that shifts the **joint
  distribution of all mediators** in the exposed, given confounders, to
  match the corresponding distribution in the unexposed.

- **`shift_k`**:  
  This estimand corresponds to an intervention that shifts the
  **distribution of a specific mediator** in the exposed, given
  confounders, to match the corresponding distribution in the unexposed,
  independent of and **without considering flow-on effects** on other
  mediators.

- **`shift_k_order`**:  
  This estimand corresponds to an intervention that shifts the
  **distribution of a specific mediator** in the exposed, given
  confounders, to match the corresponding distribution in the unexposed,
  **while considering flow-on effects** on causally descendant
  mediators.

`medRCT` estimates interventional effects based on these interventions
using a Monte Carlo simulation-based g-computation approach, as
described by Vansteelandt and Daniel (2017). The number of Monte Carlo
simulations can be specified using the argument `mcsim`, allowing users
to balance computational efficiency and estimation accuracy. The default
number of Monte Carlo simulations used is 200. Additionally, `medRCT`
automatically removes records with missing data for any of the analysis
variables, performing a complete-case analysis.

## Getting started with `medRCT`: Estimating interventional effects

The `medRCT` package is designed to handle complex mediation settings
involving **multiple interrelated mediators** and **intermediate
confounders** (exposure-induced mediator-outcome confounders). It
supports a range of variable types. Specifically, it supports
**categorical exposures** with two or more levels (Levels $\geq$ 2). For
the **outcome variable**, `medRCT` supports both **binary** and
**continuous** outcomes. Additionally, it supports any number of
mediators. Each **mediator**, including those not of primary interest
(i.e., **intermediate confounders**), can be either **binary** or
**continuous**. Any number and type of **baseline confounders** are
supported. This makes the package suitable for a broad range of
real-world applications in epidemiology, public health, and social
science research.

We begin by illustrating how to estimate the three different types of
interventional effects in the presence of multiple mediators and
intermediate confounders. Such intermediate confounders can be specified
using the argument `intermediate_confs`. We consider the following
example using simulated data inspired by a case study within the
Longitudinal Study of Australian Children (LSAC) (Goldfeld et al. 2023).
Specifically, we aim to estimate the difference in expected outcome
(risk of child mental health problems) under exposure (low socioeconomic
position) with versus without a hypothetical intervention that
individually shifts the distribution of each mediator (parental mental
health and preschool attendance) to the levels in the unexposed (high
socioeconomic position), while accounting for baseline confounders, an
intermediate (exposure-induced) mediator-outcome confounder and
correlations amongst mediators.

``` r
library(medRCT)
#> medRCT: Causal mediation analysis estimating interventional effects mapped to a target trial
#> Note: When setting intervention_type = 'shift_k_order', the order of the mediators as specified 
#> in the 'mediators' argument is important.

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

# Define intermediate confounders
intermediate_confs <- "fam_stress"

# Estimate interventional indirect effects
med_res <- medRCT(
  dat = LSACdata,                      
  exposure = "sep",                    
  outcome = "child_mh",                
  mediators = c("parent_mh", "preschool_att"), 
  intermediate_confs = intermediate_confs, 
  confounders = confounders,
  interactions_XC = "all",
  intervention_type = "all",
  bootstrap = TRUE,
  boot_args = list(R = 100, stype = "i", ci.type = "norm"),
  mcsim = 200                          
)
#> Assumed causal order for estimating effect of type 'shift_k_order': parent_mh, preschool_att
#> Conducting complete case analysis, 2499 observations were excluded due to missing data.

# Summarise the results
summary(med_res)
#> 
#> Estimated interventional indirect effect: 
#> 
#>                                  Estimate Std. Error  CI Lower  CI Upper
#> IIE_1 (p_trt - p_1)              0.011029   0.004155  0.002601  0.018890
#> IIE_2 (p_trt - p_2)             -0.000779   0.002482 -0.005452  0.004278
#> IIE_all (p_trt - p_all)          0.011012   0.004736  0.001513  0.020078
#> IIE_1_prime (p_trt - p_1_prime)  0.011824   0.003960  0.003646  0.019167
#> TCE (p_trt - p_ctr)              0.127928   0.024522  0.080967  0.177092
#>                                 p-value
#> IIE_1 (p_trt - p_1)              0.0080
#> IIE_2 (p_trt - p_2)              0.7538
#> IIE_all (p_trt - p_all)          0.0201
#> IIE_1_prime (p_trt - p_1_prime)  0.0028
#> TCE (p_trt - p_ctr)             1.8e-07
#> 
#> Estimated interventional direct effect: 
#> 
#>                                 Estimate Std. Error CI Lower CI Upper p-value
#> IDE_1 (p_1 - p_ctr)               0.1169     0.0247   0.0699   0.1666 2.2e-06
#> IDE_2 (p_2 - p_ctr)               0.1287     0.0243   0.0819   0.1773 1.2e-07
#> IDE_all (p_all - p_ctr)           0.1169     0.0245   0.0702   0.1663 1.9e-06
#> IDE_1_prime (p_1_prime - p_ctr)   0.1161     0.0247   0.0691   0.1661 2.7e-06
#> 
#> Estimated expected outcome in each trial arm:
#> 
#>           Estimate Std. Error CI Lower CI Upper p-value
#> p_1         0.3300     0.0225   0.2868   0.3751  <2e-16
#> p_2         0.3418     0.0221   0.2990   0.3855  <2e-16
#> p_all       0.3300     0.0224   0.2870   0.3748  <2e-16
#> p_ctr       0.2131     0.0101   0.1929   0.2324  <2e-16
#> p_trt       0.3410     0.0222   0.2981   0.3853  <2e-16
#> p_1_prime   0.3292     0.0226   0.2860   0.3745  <2e-16
#> 
#> Sample Size: 2608 
#> 
#> Simulations: 200
```

In the summary table:

- `IIE_all` is the estimated interventional indirect effect when
  applying the hypothetical intervention specified by `shift_all`. A
  hypothetical intervention improving both the mental health of parents
  of children and preschool attendance in low socioeconomic position to
  the levels in those with high socioeconomic position could potentially
  prevent 1 per 100 cases of child mental health problems.

- `IIE_k` is the estimated interventional indirect effect when applying
  the hypothetical intervention described by `shift_k` to a specific
  mediator `k`. A hypothetical intervention improving the mental health
  of parents of children (IIE_1) in low socioeconomic position to the
  levels in those with high socioeconomic position could potentially
  prevent 1 per 100 cases of child mental health problems. Meanwhile,
  the effect of a hypothetical intervention on preschool attendance
  (IIE_2) is negligible.

- `IIE_k_prime` is the estimated interventional indirect effect for a
  specific mediator `k`, accounting for the flow-on effects on its
  descendant mediators under the hypothetical intervention described by
  `shift_k_order`. A hypothetical intervention improving the mental
  health of parents of children (IIE_1_prime) in low socioeconomic
  position to the levels in those with high socioeconomic position,
  while accounting for flow-on effects on preschool attendance, could
  potentially prevent 1 per 100 cases of child mental health problems.

- `TCE` is the estimated total causal effect.

## Specifying interactions with `interactions_XC`

The argument `interactions_XC` in `medRCT` allows users to specify
two-way interactions amongst exposure and baseline confounders. Of note,
the algorithm automatically includes all two-way interactions between
the exposure, mediators, and intermediate confounders. These interaction
terms are fixed and cannot be modified by the user.

#### Default option:

- By default, `interactions_XC = "all"` includes **all two-way
  interactions** between the exposure and baseline confounders, while
  **baseline confounder-baseline confounder interactions are excluded**.

#### Removing interactions:

- Setting `interactions_XC = "none"` **removes all two-way
  interactions** involving the exposure and baseline confounders.

#### Custom interactions:

- **Custom specification:** Users can manually specify interaction terms
  by providing the model formula in the format
  `interactions_XC = "exposure:confounder1 + exposure:confounder2"`.
  This option provides full control over interaction terms between the
  exposure and baseline confounders in the model, enabling customisation
  for specific research questions.

## Special cases

#### Analysis without intermediate confounders:

If there are **no intermediate confounders** in the analysis, set the
argument **`intermediate_confs = NULL`** when calling the `medRCT`
function, as described in the following example:

``` r
med_res <- medRCT(
  dat = LSACdata,
  exposure = "sep",
  outcome = "child_mh",
  mediators = c("parent_mh", "preschool_att"),
  intermediate_confs = intermediate_confs,
  confounders = confounders,
  interactions_XC = "all",
  intervention_type = "all",
  intermediate_confs = NULL,
  bootstrap = TRUE,
  boot_args = list(R = 100, stype = "i", ci.type = "norm"),
  mcsim = 200
)
```

#### Single mediator case:

If there is only **one mediator**, only the `intervention_type`
**`shift_k`** can be estimated, as the other options (`shift_all` and
`shift_k_order`) require **multiple mediators** to define meaningful
interventional effects.

#### Causal ordering for the estimation of the estimand `shift_k_order`:

When estimating the estimand `shift_k_order`, the **order of mediators**
specified in the argument **`mediators`** is important, as it defines
the assumed **causal ordering** among the mediators in the analysis.
This ordering determines the **causal descendants** for the mediator of
interest, providing essential information for **accounting for flow-on
effects** after the hypothetical intervention. Specifying an incorrect
mediator ordering may result in **biased estimates** for the estimand
`shift_k_order`.

## Model assessment

We have developed a **Shiny** application for convenient model
assessment, which can be launched using the command:

``` r
medRCT_shiny(data = data)
```

Before conducting the mediation analysis, users are encouraged to assess
the models fitted by the algorithm using this Shiny application. The
Shiny app aims to provide a user-friendly interface to review model
summaries and identify potential warnings and errors, ensuring that the
models are appropriately specified before proceeding with the analysis.
If issues with model fitting are detected, users are encouraged to
adjust the exposure-confounder interaction term as needed. However,
**mediators or confounders must not be selected based on model fitting
results**.

In the Shiny application, after specifying the relevant arguments and
clicking the `Run medRCT models` button, the app will fit all regression
models required by the `medRCT` algorithm and generate detailed
summaries for each model. These summaries support model assessment,
allowing users to identify potential issues with fitting what are
usually richly-specified models. Additionally, any **warnings** or
**error messages** encountered during the model fitting process are
highlighted in the app.

## References

Afshar, Nina, S. Ghazaleh Dashti, Victoria Mar, Luc te Marvelde, Sue
Evans, Roger L. Milne, and Dallas R. English. 2024. “Do Age at
Diagnosis, Tumour Thickness and Tumour Site Explain Sex Differences in
Melanoma Survival? A Causal Mediation Analysis Using Cancer Registry
Data.” *International Journal of Cancer* 154 (5): 793–800.
<https://doi.org/10.1002/ijc.34752>.

Dashti, S. Ghazaleh, Julie A. Simpson, Vivian Viallon, Amalia
Karahalios, Margarita Moreno-Betancur, Theodore Brasky, Kathy Pan, et
al. 2022. “Adiposity and Breast, Endometrial, and Colorectal Cancer Risk
in Postmenopausal Women: Quantification of the Mediating Effects of
Leptin, c-Reactive Protein, Fasting Insulin, and Estradiol.” *Cancer
Medicine* 11 (4): 1145–59. <https://doi.org/10.1002/cam4.4434>.

Goldfeld, Sharon, Margarita Moreno-Betancur, Sarah Gray, Shuaijun Guo,
Marnie Downes, Elodie O’Connor, Francisco Azpitarte, et al. 2023.
“Addressing Child Mental Health Inequities Through Parental Mental
Health and Preschool Attendance.” *Pediatrics* 151 (5): e2022057101.
<https://doi.org/10.1542/peds.2022-057101>.

Hernán, Miguel A, and James M. Robins. 2016. “Using Big Data to Emulate
a Target Trial When a Randomized Trial Is Not Available.” *American
Journal of Epidemiology* 183 (8): 758–64.
<https://doi.org/10.1093/aje/kwv254>.

Moreno-Betancur, Margarita, and John B. Carlin. 2018. “Understanding
Interventional Effects: A More Natural Approach to Mediation Analysis.”
*Epidemiology* 29 (5): 614–17.
<https://doi.org/10.1097/EDE.0000000000000866>.

Moreno-Betancur, Margarita, Paul Moran, Denise Becker, George C. Patton,
and John B. Carlin. 2021. “Mediation Effects That Emulate a Target
Randomised Trial: Simulation-Based Evaluation of Ill-Defined
Interventions on Multiple Mediators.” *Statistical Methods in Medical
Research* 30 (6): 1395–1412. <https://doi.org/10.1177/0962280221998409>.

Robins, James M., and Sander Greenland. 1992. “Identifiability and
Exchangeability for Direct and Indirect Effects.” *Epidemiology* 3 (2):
143–55. <https://doi.org/10.1097/00001648-199203000-00013>.

VanderWeele, Tyler J., Stijn Vansteelandt, and James M. Robins. 2014.
“Effect Decomposition in the Presence of an Exposure-Induced
Mediator-Outcome Confounder.” *Epidemiology* 25 (2): 300–306.
<https://doi.org/10.1097/EDE.0000000000000034>.

Vansteelandt, Stijn, and Rhian M. Daniel. 2017. “Interventional Effects
for Mediation Analysis with Multiple Mediators.” *Epidemiology* 28 (2):
258–65. <https://doi.org/10.1097/EDE.0000000000000596>.
