# Generate Model Formula for Mediator Models

This function generates a model formula for mediator models based on the
input parameters.

## Usage

``` r
gen_formula(
  k,
  first = NULL,
  MM = NULL,
  K = NULL,
  interactions_XC,
  use_interactions_XM,
  include_all = FALSE,
  marginal = FALSE
)
```

## Arguments

- k:

  An integer specifying the index of the mediator for which the formula
  is being generated.

- first:

  An integer (optional) specifying the index of the first mediator.
  Defaults to `NULL`, indicating that the function will generate
  formulas without explicitly considering this parameter.

- MM:

  An integer (optional) specifying the index of the mediator whose
  distribution will be shifted. Defaults to `NULL`, indicating that the
  function will generate formulas without explicitly considering this
  parameter.

- K:

  An integer (optional) specifying the total number of mediators and
  intermediate confounders. Defaults to `NULL`, indicating that the
  function will generate formulas without explicitly considering this
  parameter.

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

- include_all:

  Logical.

- marginal:

  Logical. If `TRUE`, estimating marginals under `X=0`.
