#' Apply superSeq model to subsampling data
#'
#' The superSeq function fits a non-linear least squares model to subsampling data to learn the relationship between statistical power and read depth.
#'
#' @param object A subSeq summary object.
#' @param control Specify convergence criteria for non-linear least squares algorithm.
#' See  \code{\link{defaultNLSControl}}.
#' @param starts A list of combinations of starting guesses to try with
#' \code{nls}. If NULL, uses \code{\link{defaultNLSStarts}}.
#' @param new_p A vector of subsampling proportions to predict using the superSeq model fits. By default, triple the subsampling
#' depth is predicted.
#' 
#' @return A superSeq object, which is a \code{data.frame}:
#'
#' \item{fits}{The fitted objected from nls function}
#' \item{subsample}{The subsampled object from subSeq}
#' \item{predictions}{A data frame with the following columns: method used, total predicted discoveries, proportion read depth, predicted number of DE genes, and estimated statistical power.}
#' 
#' @examples
#'\dontrun{
#'library(superSeq)
#'library(subSeq)
#'library(Biobase)
#'# Load bottomly data
#'data(bottomly)
#'bottomly_counts <- exprs(bottomly)
#'bottomly_design <- pData(bottomly)
#'bottomly_counts <- bottomly_counts[rowSums(bottomly_counts) >= 10, ]
#'bottomly_proportions <- 10 ^ seq(-2, 0, 0.1)
#'# Apply subsampling methodology subSeq
#'ss = subsample(counts = bottomly_counts,
#'               proportions = bottomly_proportions,
#'               treatment=bottomly_design$strain, 
#'               method=c("voomLimma"),
#'               replications = 3,
#'               seed = 12345)
#'               ss_sum <- summary(ss)
#'               
#'# apply superSeq model
#'ss_obj <- superSeq(ss_sum)
#'
#'# plot results
#'plot(ss_obj)
#'
#'# summarise results
#'summary(ss_obj)
#'}
#' 
#' @seealso \code{\link{fitnls}}
#'
#' @import dplyr subSeq purrr ggplot2
#' @export
superSeq <- function(object, control = defaultNLSControl, starts = NULL, new_p = NULL) {
  if (is.null(new_p)) new_p <- seq(0, 3, 0.05)
  # apply non-linear least squares algorithm
  fits <- object %>%
    group_by(method) %>%
    do(fit = fitnls(., control = control, starts = starts)) %>%
    ungroup()
  # use model to provide predictions
  predictions <- fits %>% 
    group_by(method) %>%
    do(data.frame(proportion = new_p,
                  predicted = predict_func(fit =.$fit[[1]][[1]], data = new_p))) %>%
    distinct()
  
  predictions <- fits %>%
    group_by(method) %>%
    do(data.frame(total_discoveries = coef(.$fit[[1]]$fit)[1])) %>%
    right_join(predictions, by = "method") %>%
    mutate(estimated_power = predicted / total_discoveries)
  out <- list(call = match.call(), subsample = object, fits = fits, predictions = predictions)
  class(out) <- "superSeq"
  out
}

# Function to get predictions from NLS model fits 
predict_func <- function(fit, data) {
  predict(fit, newdata = data.frame(proportion = data))
}
