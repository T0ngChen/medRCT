# Causal Mediation Analysis Estimating Interventional Effects

This function performs the actual causal mediation analysis to estimate
interventional effects mapped to a hypothetical target trial.

## Usage

``` r
medRCT.fun(
  dat,
  ind = 1:nrow(dat),
  first,
  K,
  fam_type,
  mediators,
  interactions_XC,
  use_interactions_XM,
  intervention_type,
  effect_measure,
  separation_method,
  mcsim
)
```

## Arguments

- dat:

  A `data.frame` containing the dataset for analysis.

- ind:

  A `vector` of indices specifying the subset of `dat` to use for the
  analysis. Defaults to all rows of `dat`. This parameter is
  particularly useful when using this function within the `boot()`
  function from the `boot` package, as it enables resampling by
  specifying subsets of the data.

- first:

  An `integer` specifying the index of the first mediator of interest in
  the combined list of intermediate confounders and mediators.

- K:

  An `integer` specifying the total number of mediators and intermediate
  confounders. Mediators are considered sequentially based on their
  order.

- fam_type:

  A `character` string specifying the family type for modeling. Options
  typically include `"gaussian"` for continuous variables or
  `"binomial"` for binary variables.

- mediators:

  A `character` vector including the variable names for mediators
  (including intermediate confounders).

- interactions_XC:

  A `character` string specifying the two-way interactions amongst
  exposure and baseline confounders to include in the regression models
  in the estimation procedure. The default value, `"all"`, includes all
  two-way exposure-confounder interactions but excludes
  confounder-confounder interactions. Specify `"none"` to exclude all
  two-way interactions amongst exposure and baseline confounders.

- use_interactions_XM:

  Logical. Include exposure–mediator and exposure–intermediate
  confounder interactions (default is TRUE).

- intervention_type:

  A `character` string indicating the type of interventional effect to
  be estimated. Options include:

  - `"all"` (default): Estimates all types of interventional indirect
    effects.

  - `"shift_all"`: Estimates an interventional indirect effect mapped to
    a target trial assessing the impact of shifting the joint
    distribution of all mediators in the exposed to match the
    corresponding distribution in the unexposed.

  - `"shift_k"`: Estimates an interventional indirect effect mapped to a
    target trial assessing the impact of shifting the distribution of a
    specific mediator (`k`) in the exposed to match the corresponding
    distribution in the unexposed.

  - `"shift_k_order"`: Estimates an interventional indirect effect
    mapped to a target trial assessing the impact of shifting the
    distribution of a specific mediator (`k`) in the exposed to match
    the corresponding distribution in the unexposed while accounting for
    the flow-on effects on causally descendant mediators.

- effect_measure:

  A `character` string specifying the effect measure to calculate. For a
  continuous outcome, only `"Diff"` (difference in means) is allowed.
  For a binary outcome, must be either `"RD"` (risk difference) or
  `"RR"` (risk ratio). If not specified, defaults to `"Diff"` for
  continuous outcomes and `"RD"` for binary outcomes. Only one effect
  measure can be specified at a time.

- separation_method:

  Method to handle separation, only relevant for binomial (binary
  outcome) models. Options are `"brglm"` (Logistic regression models are
  fitted using bias reduction methods for generalised linear models
  implemented in the `brglm2` package) or `"discard"` (if separation is
  detected, the function returns `NA`. If this occurs during the main
  estimation, the program stops; if it occurs during bootstrapping, the
  affected bootstrap samples are discarded).

- mcsim:

  An `integer` specifying the number of Monte Carlo simulations to
  perform.
