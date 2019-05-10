#' Fit significance curve using nonlinear least squares
#'
#' Fit a nonlinear least squares model explaining the number of significant
#' genes based on the proportion of reads used. Tries each of the starting
#' combinations of parameters in the given order. If all fail, returns NULL.
#'
#' @param dat a data.frame from summary.subsamples, containing columns for
#' "proportion" and "significant"
#' @param method specifies the method used to fit the model. Default is probit (superSeq). "logit" and "smoother" are possible
#' options but are mainly for use in our manuscript.
#' @param starts A list of combinations of starting guesses to try with
#' \code{nls}. If NULL, uses \code{\link{defaultNLSStarts}}
#' @param control control for \code{nls}, by default \link{defaultNLSControl}
#'
#' @seealso \code{\link{defaultNLSStarts}}
#' @export
fitnls <- function (dat, method = "probit", starts = NULL, control = defaultNLSControl) {
  nls0 <- purrr::possibly(nls, NULL, quiet = TRUE)  
  # starting parameters
  if (is.null(starts)) {
    starts <- defaultNLSStarts(max(dat$significant))
  }
  # Create upper bound for k
  pi0 <- unique(dat$pi0)
  up_lim <- pmin(unique(dat$genes) * (1 - pi0 / 1.4) / 0.95, unique(dat$genes))
  for (s in starts) {
    if (method == "probit") {
      fit = nls0(significant ~ k * plnorm(proportion + b, mu, s),
                 dat,
                 start = s,
                 control = control,
                 lower = c(.95 * max(dat$significant), -15, 0.00001, 0),
                 upper = c(up_lim, 15, 15, 10), 
                 algorithm = "port")
    } else if (method == "logit") {
      fit = nls0(significant ~ k * plogis(log10(proportion + b), location = mu, scale = s),
                 dat,
                 start = s,
                 lower = c(.95 * max(dat$significant), -15, 0.00001, 0), 
                 upper = c(unique(dat$genes), 15, 15, 10),
                 control = control,
                 algorithm = "port")
      
    }  else {
      fit = smooth.spline(x = dat$proportion, y = dat$significant, df = 4)
    }
    if (!is.null(fit)) {
      return(list(fit = fit))
    }
  }
  warning("Failed to converge. Try different starting parameters")
  return(NULL)
}


#' Default control for nonlinear least squares fitting of saturation curves
#'
#' This is a wrapper for nls.control. It is a list that can be directly appended.
#' See nls.control for more details.
#' @export
defaultNLSControl <- nls.control(maxiter=2000, minFactor=1/1e3, tol=1e-2, printEval=FALSE,
                                 warnOnly=FALSE)

#' Default combinations of starting values for nonlinear least squares fitting
#' of saturation curves
#'
#' @param maxsig the maximum number of significant genes in the curve
#'
#' @return a list of lists, each of which serves as a set of starting
#' parameter estimates for nls
#'
#' @export
defaultNLSStarts <- function(maxsig) {
  list(
    list(k=1 * maxsig, mu=0, s=2, b=0),
    list(k=1 * maxsig, mu=-2, s=1, b=5),
    list(k=1 * maxsig, mu=.5, s=.2, b=1),
    list(k=1 * maxsig, mu=0, s=20, b=10 ),
    list(k=1 * maxsig, mu=-2, s=8, b=-10),
    list(k=2 * maxsig, mu=.5, s=.2, b=1),
    list(k=2 * maxsig, mu=0, s=2, b=0),
    list(k=2 * maxsig, mu=-2, s=1, b=0),
    list(k=2 * maxsig, mu=.5, s=7, b=7),
    list(k=2 * maxsig, mu=0, s=9, b=-7),
    list(k=3 * maxsig, mu=-2, s=1, b=0),
    list(k=3 * maxsig, mu=.5, s=.2, b=1),
    list(k=3 * maxsig, mu=0, s=2, b=0),
    list(k=3 * maxsig, mu=-2, s=3.1, b=10),
    list(k=3 * maxsig, mu=.5, s=2.2, b=-20),
    list(k=4 * maxsig, mu=0, s=2, b=0),
    list(k=4 * maxsig, mu=-2, s=1, b=0),
    list(k=4 * maxsig, mu=.5, s=.2, b=1),
    list(k=4 * maxsig, mu=0, s=8, b=2),
    list(k=4 * maxsig, mu=-2, s=4, b=10),
    list(k=1.5 * maxsig, mu=1.5, s=.2, b=0),
    list(k=1.5 * maxsig, mu=10, s=2.6, b=0),
    list(k=1.5 * maxsig, mu=-20, s=1.2, b=0),
    list(k=1.5 * maxsig, mu=5, s=5.2, b=0),
    list(k=1.5 * maxsig, mu=2, s=22, b=0),
    list(k=2.5 * maxsig, mu=-2, s=1, b=0),
    list(k=2.5 * maxsig, mu=.5, s=.2, b=1),
    list(k=2.5 * maxsig, mu=0.1, s=2, b=0),
    list(k=2.5 * maxsig, mu=-.2, s=0.01, b=0),
    list(k=3.5 * maxsig, mu=1.5, s=.2, b=1),
    list(k=3.5 * maxsig, mu=20, s=2, b=0),
    list(k=3.5 * maxsig, mu=-2.3, s=1, b=0),
    list(k=3.5 * maxsig, mu=2.5, s=.2, b=1),
    list(k = 1 * maxsig, mu = 0, s = 3, b = 0),
    list(k = 2 * maxsig, mu = -1, s = 2, b = 0),
    list(k = 3 * maxsig, mu = -1, s = 0.001, b = 0),
    list(k = 5 * maxsig, mu = -2, s = 6, b = 0),
    list(k = 3 * maxsig, mu = -6, s = 0.2, b = 0),
    list(k = 6 * maxsig, mu = 0, s = 0, b = 1),
    list(k = 3 * maxsig, mu = -40, s = -10, b = 0),
    list(k = 2 * maxsig, mu = -2, s = -5, b = 0),
    list(k = 1 * maxsig, mu = 0, s = 2, b = 0),
    list(k = 1 *maxsig,  mu = -2, s = 1, b = 5),
    list(k = 1 * maxsig, mu = 0.5, s = 0.2, b = 1),
    list(k = 1 * maxsig, mu = 0, s = 20, b = 10),
    list(k = 1 * maxsig, mu = -2, s = 8, b = -10),
    list(k = 2 * maxsig, mu = 0.5, s = 0.2, b = 1),
    list(k = 2 * maxsig, mu = 0, s = 2, b = 0),
    list(k = 2 * maxsig, mu = -2, s = 1, b = 0),
    list(k = 2 * maxsig, mu = 0.5, s = 7, b = 7),
    list(k = 2 * maxsig, mu = 0, s = 9, b = -7),
    list(k = 3 * maxsig, mu = -2, s = 1, b = 0),
    list(k = 3 * maxsig, mu = 0.5, s = 0.2, b = 1),
    list(k = 3 * maxsig, mu = 0, s = 2, b = 0),
    list(k = 3 * maxsig, mu = -2, s = 3.1, b = 10),
    list(k = 3 * maxsig, mu = 0.5, s = 2.2, b = -20),
    list(k = 4 * maxsig, mu = 0, s = 2, b = 0),
    list(k = 4 * maxsig, mu = -2, s = 1, b = 0),
    list(k = 4 * maxsig, mu = 0.5, s = 0.2, b = 1),
    list(k = 4 * maxsig, mu = 0, s = 8, b = 2),
    list(k = 4 * maxsig, mu = -2, s = 4, b = 10),
    list(k = 1.5 * maxsig, mu = 1.5, s = 0.2, b = 0),
    list(k = 1.5 * maxsig, mu = 10, s = 2.6, b = 0),
    list(k = 1.5 * maxsig, mu = -20, s = 1.2, b = 0),
    list(k = 1.5 * maxsig, mu = 5, s = 5.2, b = 0),
    list(k = 1.5 * maxsig, mu = 2, s = 22, b = 0),
    list(k = 2.5 * maxsig, mu = -2, s = 1, b = 0),
    list(k = 2.5 * maxsig, mu = 0.5, s = 0.2, b = 1),
    list(k = 2.5 * maxsig, mu = 0.1, s = 2, b = 0),
    list(k = 2.5 * maxsig, mu = -0.2, s = 0.01, b = 0),
    list(k = 3.5 * maxsig, mu = 1.5, s = 0.2, b = 1),
    list(k = 3.5 * maxsig, mu = 20, s = 2, b = 0),
    list(k = 3.5 * maxsig, mu = -2.3, s = 1, b = 0),
    list(k = 3.5 * maxsig, mu = 2.5, s = 0.2, b = 1),
    list(k = 1 * maxsig, mu = 0, s = 1/3, b = 0),
    list(k = 2 * maxsig, mu = -1, s = 1/2, b = 0),
    list(k = 3 * maxsig, mu = -1, s = 1/0.001, b = 0),
    list(k = 5 * maxsig, mu = -2, s = 6, b = 0),
    list(k = 3 * maxsig, mu = -6, s = 1/0.2, b = 0),
    list(k = 6 * maxsig, mu = 0.00010, s = 1, b = 1),
    list(k = 3 * maxsig, mu = -40, s = -10, b = 0),
    list(k = 2 * maxsig, mu = -2, s = -5, b = 0),
    list(k = 1 * maxsig, mu = 0, s = 12, b = 0),
    list(k = 1 * maxsig, mu = -2, s = 11, b = 5),
    list(k = 1 * maxsig, mu = 0.05, s = 0.2, b = 1),
    list(k = 1 * maxsig, mu = 0, s = 1/20, b = 10),
    list(k = 1 * maxsig, mu = -2, s = 1/8, b = -10),
    list(k = 2 * maxsig, mu = 0.5, s = 1/0.2, b = 1),
    list(k = 2 * maxsig, mu = 0, s = 2, b = 0),
    list(k = 2 * maxsig, mu = -0.002, s = 1, b = 0),
    list(k = 2 * maxsig, mu = 0.5, s = 3/7, b = 7),
    list(k = 2 * maxsig, mu = 0, s = 1/9, b = -7),
    list(k = 3 * maxsig, mu = -0.002, s = .00011, b = 0),
    list(k = 3 * maxsig, mu = 0.5, s = 0.2, b = 1),
    list(k = 3 * maxsig, mu = 0, s = 2, b = 0),
    list(k = 3 * maxsig, mu = -2, s = 3.1, b = 10),
    list(k = 3 * maxsig, mu = 0.5, s = 2.2, b = -20),
    list(k = 3 * maxsig, mu = 0, s = 2, b = 0),
    list(k = 4 * maxsig, mu = -2, s = 31, b = 0),
    list(k = 4 * maxsig, mu = 0.5, s = 10.2, b = 1),
    list(k = 4 * maxsig, mu = 0, s = 8, b = 2),
    list(k = 4 * maxsig, mu = -2, s = 42, b = 0),
    list(k = 1.5 * maxsig, mu = 1.5, s = 0.2, b = 0),
    list(k = 1.5 *  maxsig, mu = 2, s = 2.6, b = 0),
    list(k = 1.5 * maxsig, mu = -2, s = 1.2, b = 0),
    list(k = 1.5 * maxsig, mu = -1, s =3.2, b = 0),
    list(k = 1.5 * maxsig, mu = 1, s = 2, b = 0),
    list(k = 2.5 * maxsig, mu = -2, s = 0.4, b = 0),
    list(k = 2.5 * maxsig, mu = 0.5, s = 0.2, b = 0),
    list(k = 2.5 * maxsig, mu = 0.1, s = 2, b = 0),
    list(k = 2.5 * maxsig, mu = -0.2, s = 0.01, b = 0),
    list(k = 3.5 * maxsig, mu = 1.5, s = 0.2, b = 1),
    list(k = 2 *  maxsig, mu = 2, s = 2, b = .0005),
    list(k = 2 * maxsig, mu = -4, s = 1, b = 0),
    list(k = 2 * maxsig, mu = 4, s = 2, b = -.5)
  )
}
