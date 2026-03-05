# Launch Interactive Shiny App for Model Assessment

`medRCT_shiny` launches a Shiny application to provide model summaries
of all models fitted by the algorithm. The app provides a user-friendly
interface for model assessment and facilitates adjustments of
interaction terms to improve model fit.

## Usage

``` r
medRCT_shiny(data, ...)
```

## Arguments

- data:

  A `data.frame` containing the dataset for analysis. It should include
  variables for the exposure, outcome, mediators, confounders, and
  exposure-induced mediator-outcome confounders specified in the
  analysis.

- ...:

  additional arguments for shiny

## Examples

``` r
if (interactive()) {
   medRCT_shiny(data=LSACdata)
}
```
