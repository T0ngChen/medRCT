# Determine Appropriate Family Types for GLM

Identifies the appropriate family type (`"binomial()"` or
`"gaussian()"`) for a set of variables.

## Usage

``` r
family_type(data, variable_names, unique_threshold = 10)
```

## Arguments

- data:

  A `data.frame` containing the variables to be analyzed.

- variable_names:

  A `character` vector specifying the names of the variables to
  evaluate. Each variable name should correspond to a column in the
  `data.frame`.

- unique_threshold:

  An `integer` value specifying the minimum number of unique values
  required for a variable to be classified as continuous. Defaults to
  `10`.

## Value

A `list` where each element corresponds to a variable in
`variable_names`, and the value indicates the family type: either
`"binomial()"` for binary variables or `"gaussian()"` for continuous
variables.
