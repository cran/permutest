#' Unstratified group permutation
#'
#' This function takes a data frame and group column name as input
#' and returns the dataframe with the group column randomly permuted
#'
#' @param df A data frame
#' @param group_col String, the name of the column in df that corresponds to the group label
#' @param strata_col The name of the column in df that corresponds to the strata, should be NULL for unstratified permutation
#' @param seed An integer seed value
#' @return The inputted data frame with the group column randomly shuffled
#' @export
#' @examples
#' data <- data.frame(group_label = c(1, 2, 2, 1, 2, 1), outcome = 1:6)
#' permute_group(df = data, group_col = "group_label", strata_col = NULL, seed = 42)
#'
permute_group <- function(df, group_col, strata_col=NULL, seed=NULL){
  if(is.null(seed) == FALSE){
    set.seed(seed)
  }

  if(is.null(strata_col) == FALSE){
    warning("For stratified permutation, use function strat_permute_group")
  }

  perm_df <- df
  # Sample the indices of the group column
  permuted_indices <- sample(seq_along(perm_df[[group_col]]), length(perm_df[[group_col]]), replace = FALSE)

  # Reorder the group column based on the permuted indices
  perm_df[[group_col]] <- df[[group_col]][permuted_indices]

  return(perm_df)
}

#' Stratified group permutation
#'
#' This function takes a data frame and group and strata column name as input
#' and returns the dataframe with the group column randomly permuted by strata
#'
#' @param df A data frame
#' @param group_col The name of the column in df that corresponds to the group label
#' @param strata_col The name of the column in df that corresponds to the strata
#' @param seed An integer seed value
#' @return The inputted data frame with the group column randomly shuffled by strata
#' @export
#' @examples
#' data <- data.frame(group_label = c(1, 2, 2, 1, 2, 1), stratum = c(1, 1, 1, 2, 2, 2), outcome = 1:6)
#' permute_group(df = data, group_col = "group_label", strata_col = "stratum", seed = 42)
#'
strat_permute_group <- function(df, group_col, strata_col, seed=NULL){
  if(is.null(seed) == FALSE){
    set.seed(seed)
  }

  strata_values <- unique(df[[strata_col]])
  perm_df <- df

  for (s in strata_values) {
    # Get the indices of rows in the current stratum
    strata_indices <- which(df[[strata_col]] == s)

    # Sample the indices within this stratum
    permuted_indices <- sample(strata_indices, length(strata_indices), replace = FALSE)

    # Reorder the group column based on the permuted indices
    perm_df[[group_col]][strata_indices] <- df[[group_col]][permuted_indices]
  }

  return(perm_df)
}

#' Sign permutation
#'
#' This function takes a data frame and group and outcome column name as input
#' and returns the dataframe with the group column replaced with randomly assigned signs
#'
#' @param df A data frame
#' @param group_col The name of the column in df that corresponds to the group label
#' @param strata_col The name of the column in df that corresponds to the strata, should be NULL for this function
#' @param seed An integer seed value
#' @return The inputted data frame with the group column replaced with randomly assigned signs
#' @export
#' @examples
#' data <- data.frame(group_label = rep(1, 6), outcome = 1:6)
#' permute_group(df = data, group_col = "group_label", strata_col = NULL, seed = 42)
#'
permute_sign <- function(df, group_col, strata_col=NULL, seed=NULL){
  # Set seed if not NULL
  if(is.null(seed) == FALSE){
    set.seed(seed)
  }

  # Issue warning if strata inputted
  if(is.null(strata_col) == FALSE){
    warning("Stratified permutation not supported for this function")
  }

  # Permute sign
  perm_df <- df
  perm_sign <- sample(c(1, -1), size = length(df[[group_col]]), replace = T)
  perm_df[[group_col]] <- perm_sign

  return(perm_df)

}

