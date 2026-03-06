# Estimation of Joint Distribution of All Other Mediators under Exposed

Estimation of Joint Distribution of All Other Mediators under Exposed

## Usage

``` r
joint_X_nonzero(
  MM,
  k,
  first,
  K,
  data,
  dat2,
  fam_type,
  mediators,
  interactions_XC,
  use_interactions_XM,
  lnzero,
  separation_method,
  n,
  index,
  exposure
)
```

## Arguments

- MM:

  An integer specifying the index of the mediator whose distribution
  will be shifted.

- k:

  An integer specifying the index of the current mediator being
  processed.

- first:

  An integer specifying the index of the first mediator.

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

- lnzero:

  A numeric vector specifying the non-zero levels of the exposure.

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

- index:

  An integer vector specifying the indices of all mediators, excluding
  the mediator specified by `MM`.

- exposure:

  A `character` string specifying the name of the exposure variable in
  the dataset.
