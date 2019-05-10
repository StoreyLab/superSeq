#' Summary function of superSeq object
#' 
#' @param object superSeq object
#' @param depth Specified read depth proportion to provide predictions in the summary.
#' @param digits Number of digits to print.
#' @param ... Not currently used
#' 
#' @export
summary.superSeq <- function(object, depth = c(1, 2, 3), digits = getOption("digits"), ...) {
  # Call
  model_fit <- object$fits$fit[[1]]
  total_dis <- object$predictions$total_discoveries[1]
  cat("\nCall:\n", deparse(object$call), "\n\n", sep = "")
  # Expected number of discoveries
  cat("Asymptotic number of discoveries:", format(total_dis, digits=digits), "\n", sep="\t")
  cat("\n")
  # Number of significant values for p-value, q-value and local FDR
  new_p <- predict(model_fit[[1]], newdata = data.frame(proportion = depth))
  df <- data.frame(depth, format(new_p, digits = digits), format(new_p / total_dis, digits = digits))

  colnames(df) <- c("Read depth proportion", "Estimated number of discoveries", "Estimated experimental power")
  print(df)
  cat("\n")
}