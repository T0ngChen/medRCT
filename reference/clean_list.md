# Clean Nested Lists

`clean_list` removes empty lists and `NULL` elements from a nested list
structure, cleaning up the list for further processing.

## Usage

``` r
clean_list(x)
```

## Arguments

- x:

  A nested list.

## Value

A cleaned list with empty lists and `NULL` elements removed.

## Details

This function performs two cleaning operations on a list:

1.  At the first level, it removes any empty lists
    ([`list()`](https://rdrr.io/r/base/list.html)).

2.  At the second level, it removes `NULL` elements from sublists while
    preserving other elements.
