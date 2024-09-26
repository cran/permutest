test_that("permutation test works", {
  # outcome the same for each unit, test that p-value is 1
  data <- data.frame(group = c(rep(1, 10), rep(2, 10)),
                     outcome = c(rep(1, 10), rep(1, 10)))
  outcome <- permutation_test(df = data, group_col = "group", outcome_col = "outcome",
                   test_stat = "diff_in_means", perm_func = permute_group,
                   alternative = "greater",
                   shift = 0, reps = 10,
                   return_perm_dist = F, return_test_dist = T, seed = 42)
  expect_equal(outcome$p_value, 1)

  # with difference in medians test stat
  data <- data.frame(group = c(rep(1, 3), rep(2, 3)), out = c(rep(1, 3), rep(2, 3)))
  output <- permutation_test(df = data, group_col = "group", outcome_col = "out",
                              test_stat = "diff_in_medians", perm_func = permute_group,
                              alternative = "greater",
                              shift = 0, reps = 10^4,
                              return_perm_dist = F, return_test_dist = T, seed = 42)
  expect_equal(unique(output$test_stat_dist), c(1, -1))
  expect_equal(output$p_value, 1)

  # mice example from https://www.thoughtco.com/example-of-a-permutation-test-3997741
  data <- data.frame(experimental_group = c(1, 0, 1, 0, 1, 0),
                     race_time = c(10, 12, 9, 11, 11, 13))
  outcome <- permutation_test(df = data, group_col = "experimental_group",
                              outcome_col = "race_time", test_stat = "diff_in_means",
                              perm_func = permute_group,
                              alternative = "less",
                              shift = 0, reps = 1000, seed = 42)
  expect_equal(outcome$p_value, 0.1, tolerance = 0.1)

  # Strata and group the same, so no new group assignments
  data <- data.frame(group = c(rep(1, 10), rep(2, 10)),
                     strata = c(rep(1, 10), rep(2, 10)),
                     outcome = c(rep(4, 10), rep(1, 10)))
  outcome <- permutation_test(df = data, group_col = "group", outcome_col = "outcome",
                              strata_col = "strata",
                              test_stat = "diff_in_means", perm_func = strat_permute_group, reps = 10,
                              return_perm_dist = T, return_test_dist = T, seed = 42)
  # check that test stat is the same for all permutations
  expect_equal(outcome$test_stat_dist, rep(3, 10))
  # check that indices only permuted within strata
  expect_equal(sum(outcome$perm_indices_mat[,1:10] < 11), 100)

  # Create difference in max test statistic and run permutation test
  test_stat_func <- function(df, group_col, outcome_col){
    value <- max(df[[outcome_col]][df[group_col] == 1]) - max(df[[outcome_col]][df[group_col] == 0])

    return(value)
  }
  data <- data.frame(group = c(1, 1, 1, 0, 0, 0),
                     outcome = c(3, 7, 5, 6, 2, 1))
  # P-value should be 1/2 for anytime 7 is in group 1
  output <- permutation_test(df = data, group_col = "group", outcome_col = "outcome",
                             strata_col = NULL, test_stat = test_stat_func,
                             perm_func = permute_group, alternative = 'greater', reps = 1000, seed = 42)
  expect_equal(round(output$p_value, 1), 0.5)

  # One-sample problem
  data <- data.frame(group = rep(1, 3), outcome = c(-1, 1, 2))
  perm_set <- as.matrix(expand.grid(c(-1, 1), c(-1, 1), c(-1, 1)))
  output <- permutation_test(df = data, group_col = "group", outcome_col = "outcome",
                              test_stat = "mean", perm_func = permute_sign,
                              alternative = "greater",
                              perm_set = perm_set,
                              complete_enum = T)
  expect_equal(output$p_value, 0.375)

  # One-sample problem, 2/32 combn same
  data <- data.frame(group = rep(1, 5), out = 0:4)
  output <- permutation_test(df = data, group_col = "group", outcome_col = "out",
                             test_stat = "mean", perm_func = permute_sign,
                             alternative = "greater", seed = 42)
  expect_equal(round(output$p_value, 2), 0.06)

})

test_that("one-sample permutation test works", {
  output <- one_sample(c(-1, 1, 2), seed = 42)
  expect_equal(round(output, 2), 0.39)
})

test_that("two-sample permutation test works", {
  output <- two_sample(x = c(10, 9, 11), y = c(12, 11, 13),
                       alternative = "less", seed = 42)
  expect_equal(round(output, 1), 0.1)
})

test_that("permutation test confidence interval works", {
  x <- c(35.3, 35.9, 37.2, 33.0, 31.9, 33.7, 36.0, 35.0,
         33.3, 33.6, 37.9, 35.6, 29.0, 33.7, 35.7)
  y <- c(32.5, 34.0, 34.4, 31.8, 35.0, 34.6, 33.5, 33.6,
         31.5, 33.8, 34.6)
  df <- data.frame(outcome = c(x, y),
                   group = c(rep(1, length(x)), rep(0, length(y))))

  output <- permutation_test_ci(df, "group", "outcome", strata_col = NULL,
                      test_stat = "diff_in_means",
                      perm_func = permute_group,
                      upper_bracket = NULL,
                      lower_bracket = NULL,
                      cl = 0.95,
                      e = 0.01,
                      reps = 100000,
                      seed = 42)

  expect_equal(output$ci[1], -0.65, tolerance=.1)
  expect_equal(output$ci[2], 2.34, tolerance=.1)

  # One-sample problem CI
  x <- c(49, -67, 8, 6, 16, 23, 28, 41, 14, 29, 56, 24, 75, 60, -48)
  combinations <- as.matrix(expand.grid(rep(list(c(-1, 1)), 15)))
  df <- data.frame(outcome = x, group = rep(1, length(x)))
  output <- permutation_test_ci(df, "group", "outcome", strata_col = NULL,
                                test_stat = "mean",
                                perm_func = permute_sign,
                                upper_bracket = c(10, 75),
                                lower_bracket = c(-20, 10),
                                perm_set = combinations,
                                cl = 0.95,
                                e = 0.001,
                                seed = 42)
  expect_equal(round(output$ci[1], 2), -0.2)
  expect_equal(round(output$ci[2], 2), 41)
})
