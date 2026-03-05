# Capture Warnings and Errors from Model Fitting

`catch_model_messages` fits a generalised linear model while capturing
and storing any warnings or errors generated during the fitting process.
This is useful for debugging.

## Usage

``` r
catch_model_messages(formula, data, family)
```

## Arguments

- formula:

  A formula specifying the model to be fitted.

- data:

  A data frame containing the variables referenced in the formula.

- family:

  A description of the error distribution and link function to be used
  in the model, as specified in `glm` (it can be either `binomial` or
  `gaussian`).

## Details

This function uses `tryCatch` and `withCallingHandlers` to handle
warnings and errors separately. Warnings are captured and stored in the
resulting model object under the `warnings` attribute, while errors are
stored under the `errors` attribute. If the model fits successfully
without warnings or errors, these attributes will be `NULL`.
