#' Plotting for superSeq object
#'
#' @param x superSeq object
#' @param ... not used
#'
#' @keywords plot
#' @aliases plot, plot.superSeq
#' @export
plot.superSeq <- function(x, ...) {
  subseq_object <- x$subsample
  predictions <- x$predictions
  cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  p0 <- subseq_object %>% 
    ggplot(aes(proportion, significant, color = method))  +
    geom_point(size = 1) +
    geom_line(dat = predictions, aes(y = predicted, x = proportion, color = method), lty = 2, size = 1.0) +
    geom_vline(xintercept = 1, lty = 3) +
    scale_colour_manual(values=cbbPalette) +
    ylab("Number of significant genes") + 
    xlab("Read depth proportion") +
    theme_bw() 
  p0
}
