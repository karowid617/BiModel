#' Function to perform M step in EM algorithm.
#'
#' Function description
#'
#'
#' @param X Matrix of binary data to cluster where in rows there are features and in columns observations.
#' @param p_nk Posterior probabilities of the hidden variables.
#'
#' @returns Function returns a \code{dlist} which contains: \describe{
#'  \item{p}{Matrix of probabilites for each observation and cluster}
#'  \item{a}{Vector of each cluster weight in Bernoulli mixture model.}
#' }
#'
#' @importFrom matrixStats rowMaxs
#' 
#'
#' @export
Bernoulli_M_step <- function(X, p_nk) {
  p <- crossprod(p_nk, X) + 1e-10
  p <- p/matrixStats::rowMaxs(p + 2*1e-10)
  a <- colSums(p_nk)/nrow(X)
  return(list(p = p, a = a))
}