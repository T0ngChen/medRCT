# Counterfactual Prediction and Random Draw

Counterfactual Prediction and Random Draw

## Usage

``` r
cf_predict(fit, data, var_name, n, family)
```

## Arguments

- fit:

  A fitted glm model object

- data:

  A `data.table` containing the data to which the counterfactual
  predictions will be applied.

- var_name:

  A character string specifying the name of the variable to store the
  random draws.

- n:

  An integer specifying the number of observations in the dataset.

- family:

  A character string specifying the distribution family to use for
  generating random draws. Must be either `"binomial"` or `"gaussian"`.
