% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/superSeq.R
\name{superSeq}
\alias{superSeq}
\title{Apply superSeq model to subsampling data}
\usage{
superSeq(object, control = defaultNLSControl, starts = NULL,
  new_p = NULL)
}
\arguments{
\item{object}{A subSeq summary object.}

\item{control}{Specify convergence criteria for non-linear least squares algorithm.
See  \code{\link{defaultNLSControl}}.}

\item{starts}{A list of combinations of starting guesses to try with
\code{nls}. If NULL, uses \code{\link{defaultNLSStarts}}.}

\item{new_p}{A vector of subsampling proportions to predict using the superSeq model fits. By default, triple the subsampling
depth is predicted.}
}
\value{
A superSeq object, which is a \code{data.frame}:

\item{fits}{The fitted objected from nls function}
\item{subsample}{The subsampled object from subSeq}
\item{predictions}{A data frame with the following columns: method used, total predicted discoveries, proportion read depth, predicted number of DE genes, and estimated statistical power.}
}
\description{
The superSeq function fits a non-linear least squares model to subsampling data to learn the relationship between statistical power and read depth.
}
\examples{
\dontrun{
library(superSeq)
library(subSeq)
library(Biobase)
# Load bottomly data
data(bottomly)
bottomly_counts <- exprs(bottomly)
bottomly_design <- pData(bottomly)
bottomly_counts <- bottomly_counts[rowSums(bottomly_counts) >= 10, ]
bottomly_proportions <- 10 ^ seq(-2, 0, 0.1)
# Apply subsampling methodology subSeq
ss = subsample(counts = bottomly_counts,
              proportions = bottomly_proportions,
              treatment=bottomly_design$strain, 
              method=c("voomLimma"),
              replications = 3,
              seed = 12345)
              ss_sum <- summary(ss)
              
# apply superSeq model
ss_obj <- superSeq(ss_sum)

# plot results
plot(ss_obj)

# summarise results
summary(ss_obj)
}

}
\seealso{
\code{\link{fitnls}}
}
