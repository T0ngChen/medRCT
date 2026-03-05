# Estimation of Joint Distribution of Mediators and Intermediate Confounders under Exposed and Unexposed

Estimation of Joint Distribution of Mediators and Intermediate
Confounders under Exposed and Unexposed

## Usage

``` r
joint_dist(
  k,
  K,
  data,
  dat2,
  fam_type,
  mediators,
  interactions_XC,
  use_interactions_XM,
  exposure_level,
  separation_method,
  n
)
```

## Arguments

- k:

  An integer specifying the index of the current mediator being
  processed.

- K:

  An integer specifying the total number of mediators and intermediate
  confounders.

- data:

  A `data.frame` containing the dataset for analysis. The dataset is
  used to fit the model for the mediators.

- dat2:

  A `data.frame` used to perform counterfactual predictions and random
  draws.

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
  in the estimation procedure.

- use_interactions_XM:

  Logical. Include exposure–mediator and exposure–intermediate
  confounder interactions (default is TRUE).

- exposure_level:

  A numeric vector specifying the levels of the exposure (e.g.,
  `c(0, 1)`) for which counterfactual predictions are performed.

- separation_method:

  Method to handle separation, only relevant for binomial (binary
  outcome) models. Options are `"brglm"` (Logistic regression models are
  fitted using bias reduction methods for generalised linear models
  implemented in the `brglm2` package) or `"discard"` (if separation is
  detected, the function returns `NA`. If this occurs during the main
  estimation, the program stops; if it occurs during bootstrapping, the
  affected bootstrap samples are discarded).

- n:

  An integer specifying the number of observations for `dat2`.
