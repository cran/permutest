test_that("permute group function works", {
  df <- data.frame(group = c(1, 1, 2, 2),
                   outcome = c(1, 1, 4, 4))
  output <- permute_group(df, "group", seed = 42)
  expect_equal(output$group, c(1, 2, 2, 1))
})

test_that("stratified permute group function works", {
  df <- data.frame(group = c(1, 1, 1, 2, 2, 2),
                   outcome = c(1, 1, 1, 4, 4, 4),
                   strata = c(1, 1, 1, 2, 2, 2))
  output <- strat_permute_group(df, "group", "strata", seed = 42)
  # test that group assignment is same since strata and group is the same
  expect_equal(output$group, c(1, 1, 1, 2, 2, 2))
})

test_that("permute sign function works", {
  df <- data.frame(group = rep(1, 6),
                   outcome = 1:6)
  output <- permute_sign(df, "group", seed=42)
  expect_equal(sum(abs(output$group)), 6)
})


