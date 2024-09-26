#' Run permutation test
#'
#' Run permutation test with user inputted data, test statistic, and permutation function
#'
#' @param df A data frame
#' @param group_col The name of the column in df that corresponds to the group label
#' @param outcome_col The name of the column in df that corresponds to the outcome variable
#' @param strata_col The name of the column in df that corresponds to the strata
#' @param test_stat Test statistic function
#' @param perm_func Function to permute group
#' @param alternative String, two-sided or one-sided (greater or less) p-value; options are 'greater', 'less', or 'two-sided'
#' @param shift Value of shift to apply in one- or two-sample problem
#' @param reps Number of iterations to use when calculating permutation p-value
#' @param perm_set Matrix of group assignments to use instead of reps iterations of perm_func
#' @param complete_enum Boolean, whether to calculate P-value under complete enumeration of permutations
#' @param return_test_dist Boolean, whether to return test statistic distribution under permutations
#' @param return_perm_dist Boolean, whether to return a matrix where each row is the group assignment under that permutation
#' @param seed An integer seed value
#' @return
#' \code{p_value}: the permutation test p-value
#'
#' \code{test_stat_dist}: array, the distribution of the test statistic under the set of permutations,
#' if return_test_dist is set to TRUE
#'
#' \code{perm_indices_mat}: matrix, each row corresponds to a permutation used
#' in the permutation test calculation
#' @examples
#' data <- data.frame(group = c(rep(1, 10), rep(2, 10)), outcome = c(rep(1, 10), rep(1, 10)))
#'
#' permutation_test(df = data, group_col = "group", outcome_col = "outcome",
#' test_stat = "diff_in_means", perm_func = permute_group, alternative = "greater",
#' shift = 0, reps = 10, return_perm_dist = TRUE, return_test_dist = TRUE, seed = 42)
#'
#' @export
permutation_test <- function(df, group_col, outcome_col, strata_col = NULL,
                             test_stat = "diff_in_means",
                             perm_func = permute_group,
                             alternative = 'two-sided',
                             shift = 0,
                             reps = 10000,
                             perm_set = NULL,
                             complete_enum = FALSE,
                             return_test_dist = FALSE, return_perm_dist = FALSE,
                             seed = NULL){
  # get unique groups
  group_names <- unique(df[[group_col]])

  # check that if nonzero shift, there's exactly 2 groups
  if((shift != 0) & (length(group_names) > 2)){
    stop("To apply shift number of unique groups must be 1 or 2.")
  }

  # set test statistic function
  if(is.character(test_stat) && test_stat == "diff_in_means"){
    warning(paste0(group_names[1], " being used as value for treatment"))
    test_stat_func <- function(df, group_col, outcome_col){
      diff_in_means(df, group_col, outcome_col, treatment_value = group_names[1])
    }
  } else if(is.character(test_stat) && test_stat == "diff_in_medians") {
    test_stat_func <- function(df, group_col, outcome_col){
      diff_in_medians(df, group_col, outcome_col, treatment_value = group_names[1])
    }
  } else if(is.character(test_stat) && test_stat == "mean") {
    test_stat_func <- function(df, group_col, outcome_col){
      one_sample_mean(df, group_col, outcome_col)
    }
  } else {
    test_stat_func <- test_stat
  }

  # if nonzero shift and 1 group, apply shift
  if((shift != 0) & length(group_names) == 1){
    df[[outcome_col]] <- df[[outcome_col]] - shift
  }

  # calculate test stat on original dataset
  og_test_stat_value <- test_stat_func(df, group_col, outcome_col)

  if(is.null(perm_set)){
    # initialize array of permutation test statistic values
    perm_stat_values <- rep(0, reps)
    # if return_perm_dist, return array or permuted indices
    if(return_perm_dist){
      perm_indices_mat <- matrix(NA, nrow = reps, ncol = nrow(df))
    }else{
      perm_indices_mat <- NULL
    }

    # set seed
    set.seed(seed)

    for(i in 1:reps){
      # shuffle data
      perm_df <- perm_func(df, group_col, strata_col, NULL)
      # if return_perm_dist, insert permutation into array
      if(return_perm_dist){
        perm_indices_mat[i, ] <- perm_df[[group_col]]
      }
      # if shift != 0 and 2 groups, modify outcome variable
      if(shift != 0 & length(group_names) == 2){
        control_to_treat <- which(perm_df[[group_col]] == group_names[1] & df[[group_col]] == group_names[2])
        treat_to_control <- which(perm_df[[group_col]] == group_names[2] & df[[group_col]] == group_names[1])
        perm_df[[outcome_col]][control_to_treat] <- perm_df[[outcome_col]][control_to_treat] + shift
        perm_df[[outcome_col]][treat_to_control] <- perm_df[[outcome_col]][treat_to_control] - shift
      }

      # calculate test statistic for that shuffling
      test_stat_value <- test_stat_func(perm_df, group_col, outcome_col)
      # insert test statistic value from permuted data into perm_stat_values array
      perm_stat_values[i] <- test_stat_value
    }
  } else {

    perm_indices_mat <- perm_set

    # initialize array of permutation test statistic values
    perm_stat_values <- rep(0, nrow(perm_indices_mat))
    reps <- length(perm_stat_values)

    for(i in 1:nrow(perm_set)){
      # apply permutation to data
      perm_df <- df
      perm_df[[group_col]] <- perm_set[i, ]

      # if shift != 0, modify outcome variable
      if(shift != 0 & length(group_names) == 2){
        control_to_treat <- which(perm_df[[group_col]] == group_names[1] & df[[group_col]] == group_names[2])
        treat_to_control <- which(perm_df[[group_col]] == group_names[2] & df[[group_col]] == group_names[1])
        perm_df[[outcome_col]][control_to_treat] <- perm_df[[outcome_col]][control_to_treat] + shift
        perm_df[[outcome_col]][treat_to_control] <- perm_df[[outcome_col]][treat_to_control] - shift
      }
      # calculate test statistic for that shuffling
      test_stat_value <- test_stat_func(perm_df, group_col, outcome_col)
      # insert test statistic value from permuted data into perm_stat_values array
      perm_stat_values[i] <- test_stat_value
    }

  }

  # calculate p-value
  if(alternative == 'two-sided'){
    if(complete_enum){
      p_value <- mean(abs(perm_stat_values) >= abs(og_test_stat_value))
    } else {
      p_value <- (sum(abs(perm_stat_values) >= abs(og_test_stat_value)) + 1) / (reps + 1)
    }
  } else if(alternative == 'greater') {
    if(complete_enum){
      p_value <- mean(perm_stat_values >= og_test_stat_value)
    } else {
      p_value <- (sum(perm_stat_values >= og_test_stat_value) + 1) / (reps + 1)
    }
  } else if(alternative == 'less') {
    if(complete_enum){
      p_value <- mean(perm_stat_values <= og_test_stat_value)
    } else {
      p_value <- (sum(perm_stat_values <= og_test_stat_value) + 1) / (reps + 1)
    }
  } else {
    stop("Enter a valid alternative: 'two-sided', 'greater' or 'less'")
  }

  # return test statistic distribution
  if(return_test_dist){
    test_stat_dist <- perm_stat_values
  } else {
    test_stat_dist <- c()
  }

  return(list(p_value = p_value,
              test_stat_dist = test_stat_dist,
              perm_indices_mat = perm_indices_mat))
}


#' One-sample permutation test
#'
#' This function runs a permutation test for the one-sample problem by calling
#' the permutation_test function using the one-sample mean test statistic.
#'
#' @param x array of data
#' @param shift Value of shift to apply in one-sample problem
#' @param alternative String, two-sided or one-sided (greater or less) p-value
#' @param reps Number of iterations to use when calculating permutation p-value
#' @param seed An integer seed value
#' @return The permutation test p-value
#' @examples
#' one_sample(x = c(-1, 1, 2), seed = 42)
#'
#' @export
one_sample <- function(x, shift = 0, alternative = "greater",
                       reps = 10^4, seed = NULL){
  # Set up dataframe for one-sample problem
  data <- data.frame(outcome = x, group = rep(1, length(x)))
  # Run permutation test
  output <- permutation_test(df = data, group_col = "group", outcome_col = "outcome",
                             strata_col = NULL,
                             test_stat = "mean",
                             perm_func = permute_sign,
                             alternative = alternative,
                             shift = shift,
                             reps = reps,
                             perm_set = NULL,
                             complete_enum = FALSE,
                             seed=seed)

  return(output$p_value)
}


#' Two-sample permutation test
#'
#' This function runs a permutation test with difference in means test statistic
#' for the two-sample problem by calling
#' the permutation_test function.
#'
#' @param x array of data for treatment group
#' @param y array of data for control group
#' @param shift Value of shift to apply in two-sample problem
#' @param alternative String, two-sided or one-sided (greater or less) p-value; options are 'greater', 'less', or 'two-sided'
#' @param reps Number of iterations to use when calculating permutation p-value
#' @param seed An integer seed value
#' @return The permutation test p-value
#' @export
#' @examples
#' two_sample(x = c(10, 9, 11), y = c(12, 11, 13), alternative = "less", seed = 42)
#'
two_sample <- function(x, y, shift = 0, alternative = "greater",
                       reps = 10^4, seed = NULL){
  # Set up dataframe for one-sample problem
  data <- data.frame(outcome = c(x, y),
                     group = c(rep(1, length(x)), rep(0, length(y)))
                    )

  # Run permutation test
  output <- permutation_test(df = data, group_col = "group",
                             outcome_col = "outcome",
                             strata_col = NULL,
                             test_stat = "diff_in_means",
                             perm_func = permute_group,
                             alternative = alternative,
                             shift = shift,
                             reps = reps,
                             perm_set = NULL,
                             complete_enum = FALSE,
                             seed=seed)

  return(output$p_value)
}


#' Construct confidence interval by inverting permutation tests
#'
#' This function constructs a confidence interval by inverting permutation tests and applying the method in Glazer and Stark, 2024.
#'
#' @param df A data frame
#' @param group_col The name of the column in df that corresponds to the group label
#' @param outcome_col The name of the column in df that corresponds to the outcome variable
#' @param strata_col The name of the column in df that corresponds to the strata
#' @param test_stat Test statistic function
#' @param perm_func Function to permute group
#' @param upper_bracket Array with 2 values that bracket upper confidence bound
#' @param lower_bracket Array with 2 values that bracket lower confidence bound
#' @param cl Confidence level, default 0.95
#' @param e Maximum distance from true confidence bound value
#' @param reps Number of iterations to use when calculating permutation p-value
#' @param perm_set Matrix of group assignments to use instead of reps iterations of perm_func
#' @param seed An integer seed value
#' @return A list containing the permutation test p-value, and the test statistic distribution if applicable
#' @export
#' @examples
#' x <- c(35.3, 35.9, 37.2, 33.0, 31.9, 33.7, 36.0, 35.0, 33.3, 33.6, 37.9, 35.6, 29.0, 33.7, 35.7)
#' y <- c(32.5, 34.0, 34.4, 31.8, 35.0, 34.6, 33.5, 33.6, 31.5, 33.8, 34.6)
#' df <- data.frame(outcome = c(x, y), group = c(rep(1, length(x)), rep(0, length(y))))
#' permutation_test_ci(df = df, group_col = "group", outcome_col = "outcome", strata_col = NULL,
#'                     test_stat = "diff_in_means", perm_func = permute_group,
#'                     upper_bracket = NULL, lower_bracket = NULL,
#'                     cl = 0.95, e = 0.01, reps = 10^3, seed = 42)
permutation_test_ci <- function(df, group_col, outcome_col, strata_col = NULL,
                                test_stat = "diff_in_means",
                                perm_func = permute_group,
                                upper_bracket = NULL,
                                lower_bracket = NULL,
                                cl = 0.95,
                                e = 0.1,
                                reps = 10000,
                                perm_set = NULL,
                                seed = 42){
  alpha <- 1-cl

  # get unique groups
  group_names <- unique(df[[group_col]])

  # initialize params for permutation_test calls
  params <- list(
    df = df,
    group_col = group_col,
    outcome_col = outcome_col,
    strata_col = strata_col,
    test_stat = test_stat,
    perm_func = perm_func,
    alternative = 'greater',
    shift = 0,
    reps = reps,
    perm_set = perm_set,
    return_test_dist = TRUE,
    return_perm_dist = TRUE,
    seed = seed
  )

  # set test statistic function
  if(is.character(test_stat) && test_stat == "diff_in_means"){
    warning(paste0(group_names[1], " being used as value for treatment"))
    test_stat_func <- function(df, group_col, outcome_col){
      diff_in_means(df, group_col, outcome_col, treatment_value = group_names[1])
    }
  } else if(is.character(test_stat) && test_stat == "diff_in_medians") {
    test_stat_func <- function(df, group_col, outcome_col){
      diff_in_medians(df, group_col, outcome_col, treatment_value = group_names[1])
    }
  } else if(is.character(test_stat) && test_stat == "mean") {
    test_stat_func <- function(df, group_col, outcome_col){
      one_sample_mean(df, group_col, outcome_col)
    }
  } else {
    test_stat_func <- test_stat
  }

  # get test statistic value from observed data
  obs_diff <- test_stat_func(df, group_col, outcome_col)
  # run permutation test with shift 0
  output <- do.call(permutation_test, params)
  # if perm_set is NULL, generate set of permutations to use
  if(is.null(perm_set)){
    params$perm_set <- output$perm_indices_mat
  }
  # make sure reps matches perm_set
  params$reps <- nrow(params$perm_set)

  # Speed up for 2 sample problem -- calculate necessary values
  if(is.character(test_stat) && test_stat == "diff_in_means"){
    nt <- sum(df[[group_col]] == group_names[1])
    nc <- nrow(df) - nt

    n_switch <- rep(0, nrow(params$perm_set))
    for(row in 1:nrow(params$perm_set)){
      n_switch[row] <- sum(df[[group_col]] == group_names[1] & params$perm_set[row, ] == group_names[2])
    }
    # Save test stat dist under shift = 0
    test_stat_dist <- output$test_stat_dist
  } else if(is.character(test_stat) && test_stat == "mean") {
    # get size of data
    n <- nrow(df)
    # get number of group values equal to -1 for each permutation
    num_neg <- rep(0, nrow(params$perm_set))
    for(row in 1:nrow(params$perm_set)){
      num_neg[row] <- sum(params$perm_set[row, ] == -1)
    }
    # Save test stat dist under shift = 0
    test_stat_dist <- output$test_stat_dist
  }

  ## Upper CI bound ##
  # Set initial values
  if(is.null(upper_bracket)){
    u_upper <- max(df[[outcome_col]])
  } else {
    u_upper <- upper_bracket[2]
  }

  if(is.null(lower_bracket)){
    u_lower <- obs_diff
  } else {
    u_lower <- upper_bracket[1]
  }

  # Check that root is bracketed
  # Check upper bound
  m <- u_upper
  if(is.character(test_stat) && test_stat == "diff_in_means"){
    new_test_dist <- test_stat_dist + m * n_switch * (1/nt + 1/nc)
    upper_p_value <- 2*min((sum(new_test_dist >= obs_diff) + 1)/ (params$reps + 1),
                           (sum(new_test_dist <= obs_diff) + 1)/ (params$reps + 1))
  } else if(is.character(test_stat) && test_stat == "mean") {
    new_test_dist <- test_stat_dist + (m * (2 * num_neg) / n)
    upper_p_value <- 2*min((sum(new_test_dist >= obs_diff) + 1)/ (params$reps + 1),
                           (sum(new_test_dist <= obs_diff) + 1)/ (params$reps + 1))
  } else {
    params$shift <- m
    # run permutation test
    output <- do.call(permutation_test, params)
    # calculate p-value
    upper_p_value <- 2*min(output$p_value,
                           (sum(output$test_stat_dist <= obs_diff) + 1) / (params$reps + 1))
  }

  if(upper_p_value > alpha){
    stop("Upper bound of upper bracket must yield p-value less than alpha.
         Choose a different value for upper_bracket.")
  }

  # Check lower bound
  m <- u_lower
  if(is.character(test_stat) && test_stat == "diff_in_means"){
    new_test_dist <- test_stat_dist + m * n_switch * (1/nt + 1/nc)
    upper_p_value <- 2*min((sum(new_test_dist >= obs_diff) + 1)/ (params$reps + 1),
                           (sum(new_test_dist <= obs_diff) + 1)/ (params$reps + 1))
  } else if(identical(test_stat_func, one_sample_mean)) {
    new_test_dist <- test_stat_dist + (m * (2 * num_neg) / n)
    upper_p_value <- 2*min((sum(new_test_dist >= obs_diff) + 1)/ (params$reps + 1),
                           (sum(new_test_dist <= obs_diff) + 1)/ (params$reps + 1))
  } else {
    params$shift <- m
    # run permutation test
    output <- do.call(permutation_test, params)
    # calculate p-value
    upper_p_value <- 2*min(output$p_value,
                           (sum(output$test_stat_dist <= obs_diff) + 1) / (params$reps + 1))
  }

  if(upper_p_value < alpha){
    stop("Lower bound of upper bracket must yield p-value greater than alpha.
         Choose a different value for upper_bracket.")
  }

  # Bisection method
  while(u_upper-u_lower > e){
    # calculate midpoint
    m <- (u_upper + u_lower)/2
    # use speed up if difference in means test statistic
    if(is.character(test_stat) && test_stat == "diff_in_means"){
      new_test_dist <- test_stat_dist + m * n_switch * (1/nt + 1/nc)
      upper_p_value <- 2*min((sum(new_test_dist >= obs_diff) + 1)/ (params$reps + 1),
                             (sum(new_test_dist <= obs_diff) + 1)/ (params$reps + 1))
    } else if(is.character(test_stat) && test_stat == "mean") {
      new_test_dist <- test_stat_dist + (m * (2 * num_neg) / n)
      upper_p_value <- 2*min((sum(new_test_dist >= obs_diff) + 1)/ (params$reps + 1),
                             (sum(new_test_dist <= obs_diff) + 1)/ (params$reps + 1))
    } else {
      params$shift <- m
      # run permutation test
      output <- do.call(permutation_test, params)
      # calculate p-value
      upper_p_value <- 2*min(output$p_value,
                             (sum(output$test_stat_dist <= obs_diff) + 1) / (params$reps + 1))
    }
    # update bounds
    if(upper_p_value > alpha){
      u_lower <- m
    } else {
      u_upper <- m
    }
  }

  ## Lower CI bound ##
  # Set initial values
  if(is.null(lower_bracket)){
    l_upper <- obs_diff
  } else {
    l_upper <- lower_bracket[2]
  }

  if(is.null(lower_bracket)){
    l_lower <- -max(df[[outcome_col]])
  } else {
    l_lower <- lower_bracket[1]
  }

  # Check that root is bracketed
  # Check upper bound
  m <- l_upper
  if(is.character(test_stat) && test_stat == "diff_in_means"){
    new_test_dist <- test_stat_dist + m * n_switch * (1/nt + 1/nc)
    lower_p_value <- 2*min((sum(new_test_dist >= obs_diff) + 1)/ (params$reps + 1),
                           (sum(new_test_dist <= obs_diff) + 1)/ (params$reps + 1))
  } else if(is.character(test_stat) && test_stat == "mean") {
    new_test_dist <- test_stat_dist + (m * (2 * num_neg) / n)
    lower_p_value <- 2*min((sum(new_test_dist >= obs_diff) + 1)/ (params$reps + 1),
                           (sum(new_test_dist <= obs_diff) + 1)/ (params$reps + 1))
  } else {
    params$shift <- m
    # run permutation test
    output <- do.call(permutation_test, params)
    # calculate p-value
    lower_p_value <- 2*min(output$p_value,
                           (sum(output$test_stat_dist <= obs_diff) + 1) / (params$reps + 1))
  }

  if(lower_p_value < alpha){
    stop("Upper bound of lower bracket must yield p-value greater than alpha.
         Choose a different value for lower_bracket.")
  }

  # Check lower bound
  m <- l_lower
  if(is.character(test_stat) && test_stat == "diff_in_means"){
    new_test_dist <- test_stat_dist + m * n_switch * (1/nt + 1/nc)
    lower_p_value <- 2*min((sum(new_test_dist >= obs_diff) + 1)/ (params$reps + 1),
                           (sum(new_test_dist <= obs_diff) + 1)/ (params$reps + 1))
  } else if(is.character(test_stat) && test_stat == "mean") {
    new_test_dist <- test_stat_dist + (m * (2 * num_neg) / n)
    lower_p_value <- 2*min((sum(new_test_dist >= obs_diff) + 1)/ (params$reps + 1),
                           (sum(new_test_dist <= obs_diff) + 1)/ (params$reps + 1))
  } else {
    params$shift <- m
    # run permutation test
    output <- do.call(permutation_test, params)
    # calculate p-value
    lower_p_value <- 2*min(output$p_value,
                           (sum(output$test_stat_dist <= obs_diff) + 1) / (params$reps + 1))
  }

  if(lower_p_value > alpha){
    stop("Lower bound of lower bracket must yield p-value less than than alpha.
         Choose a different value for lower_bracket.")
  }

  while(l_upper-l_lower > e){
    # calculate midpoint
    m <- (l_upper + l_lower)/2
    # use speed up if difference in means test statistic
    if(is.character(test_stat) && test_stat == "diff_in_means"){
      new_test_dist <- test_stat_dist + m * n_switch * (1/nt + 1/nc)
      lower_p_value <- 2*min((sum(new_test_dist >= obs_diff) + 1)/ (params$reps + 1),
                             (sum(new_test_dist <= obs_diff) + 1)/ (params$reps + 1))
    } else if(is.character(test_stat) && test_stat == "mean") {
      new_test_dist <- test_stat_dist + (m * (2 * num_neg) / n)
      lower_p_value <- 2*min((sum(new_test_dist >= obs_diff) + 1)/ (params$reps + 1),
                             (sum(new_test_dist <= obs_diff) + 1)/ (params$reps + 1))
    } else {
      params$shift <- m
      # run permutation test
      output <- do.call(permutation_test, params)
      # calculate p-value
      lower_p_value <- 2*min(output$p_value,
                             (sum(output$test_stat_dist <= obs_diff) + 1) / (params$reps + 1))
    }

    # update bounds
    if(lower_p_value >= alpha){
      l_upper <- m
    } else {
      l_lower <- m
    }
  }

  return(list(obs_diff = obs_diff, ci = c(l_lower, u_upper)))
}


