#' Fisher combining function
#'
#' This function takes an array of p-values and returns a combined p-value using fisher's combining function:
#' \eqn{-2 \sum_i \log(p_i)}
#'
#' @param pvalues Array of p-values
#' @return Combined p-value using fisher's method
#' @export
#' @examples
#' fisher(pvalues = c(.05, .1, .5))
#'
fisher <- function(pvalues){
  # compute fisher combination
  value <- -2 * log(prod(pvalues))

  return(value)
}


#' Tippett combining function
#'
#' This function takes an array of p-values and returns a combined p-value using Tippett's combining function:
#' \eqn{\max_i \{1-p_i\}}
#'
#' @param pvalues Array of p-values
#' @return Combined p-value using Tippett's method
#' @export
#' @examples
#' tippett(pvalues = c(.05, .1, .5))
#'
tippett <- function(pvalues){
  # compute tippett
  value <- max(1-pvalues)

  return(value)
}


#' Liptak combining function
#'
#' This function takes an array of p-values and returns a combined p-value using Liptak's combining function:
#' \eqn{\sum_i \Phi^{-1}(1-p_i)} where \eqn{\Phi} is the CDF of the Normal distribution
#'
#' @importFrom stats qnorm
#' @param pvalues Array of p-values
#' @return Combined p-value using Liptak's method
#' @export
#' @examples
#' liptak(pvalues = c(.05, .1, .5))
#'
liptak <- function(pvalues){

  value <- sum(qnorm(1-pvalues))

  return(value)
}


#' Run NPC
#'
#' This function takes a data frame and group and outcome column names as input
#' and returns the nonparametric combination of tests (NPC) omnibus p-value
#'
#' @param df A data frame
#' @param group_col The name of the column in df that corresponds to the group label
#' @param outcome_cols The names of the columns in df that corresponds to the outcome variable
#' @param strata_col The name of the column in df that corresponds to the strata
#' @param test_stat Test statistic function
#' @param perm_func Function to permute group, default is permute_group which randomly permutes group assignment
#' @param combn Combining function method to use, takes values 'fisher', 'tippett', or 'liptak', or a user defined function
#' @param shift Value of shift to apply in one- or two-sample problem
#' @param reps Number of iterations to use when calculating permutation p-value
#' @param perm_set Matrix of permutations to use instead of reps iterations of perm_func
#' @param complete_enum Boolean, whether to calculate P-value under complete enumeration of permutations
#' @param seed An integer seed value
#' @return The omnibus p-value
#' @export
#' @examples
#' data <- data.frame(group = c(rep(1, 4), rep(2, 4)),
#' out1 = c(0, 1, 0, 0, 1, 1, 1, 0),
#' out2 = rep(1, 8))
#' output <- npc(df = data, group_col = "group",
#'               outcome_cols = c("out1", "out2"), perm_func = permute_group,
#'               combn = "fisher", reps = 10^4, seed=42)
#'
npc <- function(df, group_col, outcome_cols, strata_col = NULL,
                test_stat = "diff_in_means",
                perm_func = permute_group,
                combn = "fisher",
                shift = 0,
                reps = 10000,
                perm_set = NULL,
                complete_enum = FALSE,
                seed = NULL){

  # Check the combination function
  if(is.character(combn) && combn == "fisher"){
    combn_func <- fisher
  } else if(is.character(combn) && combn == "tippett"){
    combn_func <- tippett
  } else if(is.character(combn) && combn == "liptak"){
    combn_func <- liptak
  } else {
    combn_func <- combn
  }

  # Get matrix of p-values
  n_out <- length(outcome_cols)
  obs_p_value <- rep(NA, n_out)
  p_value_mat <- matrix(NA, nrow = reps, ncol = n_out)

  for(i in 1:n_out){
      output <- permutation_test(df = df, group_col = group_col,
                                 outcome_col = outcome_cols[i], strata_col = strata_col,
                                 test_stat = test_stat, perm_func = perm_func,
                                 shift = shift, reps = reps, perm_set = perm_set,
                                 complete_enum = complete_enum, alternative = "greater",
                                 return_test_dist = T, return_perm_dist = T,
                                 seed = seed)

      obs_p_value[i] <- output$p_value

    if(i == 1){
      perm_set <- output$perm_indices_mat
      reps <- nrow(perm_set)
    }
      p_value_mat[,i] <- (reps - rank(output$test_stat_dist, ties.method = "min") + 1)/(reps+1)
  }
  # Get combined p-values
  combn_pvalues <- rep(NA, reps)
  obs_combn_pvalue <- combn_func(obs_p_value)

  for(j in 1:reps){
    combn_pvalues[j] <- combn_func(p_value_mat[j, ])
  }
  # Get omnibus p-values
  omnibus_p <- (sum(combn_pvalues >= obs_combn_pvalue)+1) / (reps+1)
  # Return omnibus p-value
  return(omnibus_p)
}

#' Adjust p-values for multiple testing
#'
#' This function takes an array of p-values and returns adjusted p-values using
#' user-inputted FWER or FDR correction method
#'
#' @param pvalues Array of p-values
#' @param method The FWER or FDR correction to use, either 'holm-bonferroni',
#' 'bonferroni', or 'benjamini-hochberg'
#' @return Adjusted p-values
#' @export
#' @examples
#' adjust_p_value(pvalues = c(.05, .1, .5), method='holm-bonferroni')
#'
adjust_p_value <- function(pvalues, method='holm-bonferroni'){
  # get number of p-values
  n <- length(pvalues)
  if(method == 'holm-bonferroni'){
    order <- rank(pvalues, ties.method = 'last')
    adj_pvalues <- pmin(pvalues * (n - order + 1), rep(1, n))
    prev_index <- which(order == 1)
    for (i in 1:n) {
      current_index <- which(order == i)
      adj_pvalues[current_index] <- max(adj_pvalues[prev_index], adj_pvalues[current_index])
      prev_index <- current_index
    }
  } else if (method == "bonferroni"){
    adj_pvalues <- pmin(pvalues*n, rep(1, n))
  } else if (method == "benjamini-hochberg"){
    order <- rank(pvalues, ties.method = 'last')
    adj_pvalues <- pmin(pvalues * (n / order), rep(1, n))
    prev_index <- which(order == n)
    for (i in n:1) {
      current_index <- which(order == i)
      adj_pvalues[current_index] <- min(adj_pvalues[prev_index], adj_pvalues[current_index])
      prev_index <- current_index
    }
  } else {
    stop("Method must be 'holm-bonferroni', 'bonferroni', or 'benjamini-hochberg'")
  }

  return(adj_pvalues)
}

