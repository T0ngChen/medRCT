# Generate Model Formulas for Shiny-Based Model Assessment

This function generates model formulas using the original variable names
as they appear in the dataset.

## Usage

``` r
gen_formula_shiny(
  mediators = mediators,
  exposure = exposure,
  k,
  first = NULL,
  MM = NULL,
  K = NULL,
  interactions_XC,
  include_all = FALSE,
  marginal = FALSE
)
```

## Arguments

- mediators:

  A `character` vector specifying the names of the variables in the
  dataset corresponding to the mediators of interest. The mediators can
  be either binary or continuous. When estimating the effect type
  `"shift_k_order"`, the order of mediators in the vector is important,
  as it interpreted as the assumed causal ordering of the mediators.

- exposure:

  A `character` string specifying the name of the exposure variable in
  the dataset. The exposure variable must be categorical, with `0`
  explicitly denoting the unexposed (or control) group, which is taken
  as the reference group. Other values represent different,
  non-reference exposure categories.

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

- include_all:

  Logical.

- marginal:

  Logical. If `TRUE`, estimating marginals under `X=0`.
