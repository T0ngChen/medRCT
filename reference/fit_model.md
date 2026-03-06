# Fit a generalized linear model with separation handling

This function fits a generalized linear model and handles data
separation in binomial models using either the `brglm2` package for
bias-reduced estimation or the `detectseparation` method to discard
separated data.

## Usage

``` r
fit_model(
  formula,
  data,
  family,
  separation_method,
  exposure_name = NULL,
  mediator_names = NULL,
  ...
)
```

## Arguments

- formula:

  A formula specifying the model.

- data:

  A data frame containing the variables in the model.

- family:

  Link function to be used in the model.

- separation_method:

  Method to handle separation, only relevant for binomial (binary
  outcome) models. Options are `"brglm"` (Logistic regression models are
  fitted using bias reduction methods for generalised linear models
  implemented in the `brglm2` package) or `"discard"` (if separation is
  detected, the function returns `NA`. If this occurs during the main
  estimation, the program stops; if it occurs during bootstrapping, the
  affected bootstrap samples are discarded).

- ...:

  Additional arguments passed to `glm`.

## Value

An object of class `glm` or `brglmFit`, depending on the method used.
