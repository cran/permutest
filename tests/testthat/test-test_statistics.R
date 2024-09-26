test_that("difference in means function works", {
  df <- data.frame(group = c(1, 1, 2, 2),
                   outcome = c(1, 1, 4, 4))
  expect_equal(diff_in_means(df, "group", "outcome", 1), -3)
  expect_equal(diff_in_means(df, "group", "outcome", 2), 3)
})

test_that("difference in medians function works", {
  df <- data.frame(group = c(1, 1, 1, 2, 2, 2),
                   outcome = c(1, 2, 3, 4, 5, 6))
  expect_equal(diff_in_means(df, "group", "outcome", 1), -3)
  expect_equal(diff_in_means(df, "group", "outcome", 2), 3)
})

test_that("t-test statistic function works", {
  df <- data.frame(group = c(1, 1, 2, 2),
                   outcome = c(1, 2, 3, 4))
  expect_equal(ttest_stat(df, "group", "outcome"), -2.82, tolerance = 0.01)
})

test_that("one-way anova test statistic function works", {
  df <- data.frame(group = c(1, 1, 2, 2, 3, 3),
                   outcome = c(1, 2, 3, 4, 5, 6))
  expect_equal(one_way_anova_stat(df, "group", "outcome"), 16)
})
