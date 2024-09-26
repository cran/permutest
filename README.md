# permutest

The goal of permutest is to allow users to run permutation tests and construct associated confidence intervals for *any* test statistic and *any* set of permutations.

## Installation

You can install the development version of permutest from [GitHub](https://github.com/akglazer/permutest) with:

``` r
# install.packages("devtools")
devtools::install_github("akglazer/permutest")
```

## Example

This is a basic example which shows you how to use the permutation test function with a difference in means test statistic:

``` r
library(permutest)

data <- data.frame(group = c(rep(1, 10), rep(2, 10)),
                     outcome = c(rep(1, 10), rep(1, 10)))

output <- permutation_test(df = data, group_col = "group", outcome_col = "outcome",
                   test_stat = "diff_in_means", perm_func = permute_group,
                   alternative = "greater", shift = 0, reps = 10,
                   return_perm_dist = F, return_test_dist = T, seed = 42)
```

Alternatively, for the two-sample problem with difference in means test statistic, you can use the function `two_sample()` (and for the one-sample problem you can use `one_sample()`):

``` r
output <- two_sample(x = c(10, 9, 11), y = c(12, 11, 13), alternative = "less", seed = 42)
```

The following is an example of how to construct a confidence interval for a shift parameter in the two-sample problem:

``` r
x <- c(35.3, 35.9, 37.2, 33.0, 31.9, 33.7, 36.0, 35.0,
         33.3, 33.6, 37.9, 35.6, 29.0, 33.7, 35.7)
y <- c(32.5, 34.0, 34.4, 31.8, 35.0, 34.6, 33.5, 33.6,
         31.5, 33.8, 34.6)
df <- data.frame(outcome = c(x, y), group = c(rep(1, length(x)), rep(0, length(y))))

output <- permutation_test_ci(df, "group", "outcome", strata_col = NULL,
                      test_stat = "diff_in_means",
                      perm_func = permute_group,
                      upper_bracket = NULL,
                      lower_bracket = NULL,
                      cl = 0.95,
                      e = 0.01,
                      reps = 100000,
                      seed = 42)
```
