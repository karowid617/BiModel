#' Function to perform initialization for EM step
#'
#' Function description
#'
#'
#' @param xdata Matrix of binary data to cluster where in rows there are features and in columns observations.
#' @param k Number of clusters to divide your data.
#' @param ini Method for parameters initialization. There are three options "random","kmeans" or "kmeanspp".
#'
#' @returns Function returns a \code{dlist} which contains: \describe{
#'  \item{p}{Matrix of probabilites for each observation and cluster.}
#'  \item{a}{Vecotr of each cluster weight in Bernoulli mixture model.}
#' }
#'
#' 
#' @seealso \code{\link{Bernoulli_M_step}}
#'
#' @import stats ClusterR
#'
#' @export
bernoulliEM_ini <- function(xdata, k, ini){
  switch(ini,
         "random" = {
           p <- matrix(ncol = ncol(xdata), nrow = k)
           for (i in 1:k) {
             p_random <- runif(ncol(xdata), min = 0.05)
             p[i,] <- p_random
           }
           a <- runif(k, min = 0.05)
           a <- a/sum(a)
         },
         "kmeans" = {
           kmn <- kmeans(xdata, centers = k, nstart = 20)
           p <- abs(kmn$centers + 1e-4)
           a <- kmn$size/sum(kmn$size)
         },
         "kmeanspp" = {
           kmnpp <- ClusterR::KMeans_rcpp(xdata, clusters = k, num_init = 20)
           p <- abs(kmnpp$centroids + 1e-4)
           a <- kmnpp$obs_per_cluster/sum(kmnpp$obs_per_cluster)
         }
  )
  out <- list("p" = p,
              "a" = a)
  return(out)
}