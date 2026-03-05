# Collect All Models fitted by the algorithm

`collect_models` fits and collects models required for the algorithm.

## Usage

``` r
collect_models(
  data,
  exposure,
  outcome,
  mediators,
  intermediate_confs,
  confounders,
  interactions_XC = "all",
  intervention_type = c("all", "shift_all", "shift_k", "shift_k_order")
)
```

## Arguments

- data:

  A `data.frame` containing the dataset for analysis.

- exposure:

  A `character` string specifying the name of the exposure variable in
  the dataset. The exposure variable must be categorical, with `0`
  explicitly denoting the unexposed (or control) group, which is taken
  as the reference group. Other values represent different,
  non-reference exposure categories.

- outcome:

  A `character` string specifying the name of the outcome variable in
  the dataset. The outcome variable can be either binary or continuous.

- mediators:

  A `character` vector specifying the names of the variables in the
  dataset corresponding to the mediators of interest. The mediators can
  be either binary or continuous. When estimating the effect type
  `"shift_k_order"`, the order of mediators in the vector is important,
  as it interpreted as the assumed causal ordering of the mediators.

- intermediate_confs:

  A `character` vector specifying the names of the variables in the
  dataset corresponding to intermediate confounders. The intermediate
  confounders can be either binary or continuous. If `NULL`, no
  intermediate confounders are specified.

- confounders:

  A `character` vector listing the names of the variables in the dataset
  corresponding to the baseline confounders.

- interactions_XC:

  A `character` string specifying the two-way interactions amongst
  exposure and baseline confounders to include in the regression models
  in the estimation procedure. The default value, `"all"`, includes all
  two-way exposure-confounder interactions but excludes
  confounder-confounder interactions. Specify `"none"` to exclude all
  two-way interactions amongst exposure and baseline confounders.

- intervention_type:

  A `character` string indicating the type of interventional effect to
  be estimated.
