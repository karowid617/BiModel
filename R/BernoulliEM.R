#' Function to perform binomial mixture model (BMM)
#'
#' Function description
#'
#'
#' @param X Matrix of binary data to cluster where in rows there are features and in columns observations.
#' @param K Number of clusters to divide your data. Default K=2.
#' @param start_ini Number of starting initialization.
#' @param ini Method for parameters initialization. There are three options 'random' (default), 'kmeans' or "kmeanspp".
#' @param m_iter Maximum number of iteration for EM algorithm. Default value 3000.
#' @param eps Minimum delta of model LL to stop EM algorithm.Default value 1e-40.
#'
#' @returns Function returns a \code{list} which contains: \describe{
#'  \item{probs}{Matrix of finall probabilities for each observation and cluster.}
#'  \item{alphas}{Finall vector of each cluster weight in Bernoulli mixture model.}
#'  \item{clusters}{Cluster assigment of each observation.}
#'  \item{iter}{Number of EM interations.}
#'  \item{delta}{Final EM step delta fo algorithm stop.}
#'  \item{bic}{Bayesian information criterion (BIC) value for fitted model.}
#'  \item{ll}{Log-Likelihood of fitted model in clustering.}
#' }
#' 
#' @importFrom matrixStats rowLogSumExps
#' @importFrom methods hasArg
#' 
#' @examples
#' \dontrun{
#' data(example)
#' res<-BernoulliEM(example$nouli_data, 4, start_ini = 20, ini = "random", m_iter=3000, eps=1e-40)
#' }
#' 
#' @seealso \code{\link{bernoulliEM_ini}}
#'
#' @export
BernoulliEM <- function(X, K=2, start_ini = 20, ini = "random", m_iter=3000, eps=1e-40){
  
    # Out definition ----
    res <- list('probs' = NULL, 'alphas' = NULL, 'clusters' = NULL, 'itr' = NULL, 
                'delta' = NULL, 'bic' = NULL, 'll' = NULL)
    
    # Arguments check ----
    if (!hasArg("X")){
      stop("No data included.")}
    
    
    
    
    # Param initialization ----
    if(ini == 'random'){
      logLik_ini <- matrix(nrow = start_ini, ncol = 1)
      params_ini <- list()
      for(i in 1:start_ini){
        params <- bernoulliEM_ini(xdata = X, k = K, ini = ini)
        logL <- tcrossprod(X,log(params$p))+tcrossprod((1-X),(log(1-params$p)))
        logL <- sweep(logL, MARGIN = 2, log(params$a), FUN = '+')
        logL <- sum(matrixStats::rowLogSumExps(logL))
        
        logLik_ini[i,] <- logL 
        params_ini[[length(params_ini)+1]] <- params
      }
      best_logL <- which.min(logLik_ini)
      p <- params_ini[[best_logL]]$p
      a <- params_ini[[best_logL]]$a
      
      rm(logLik_ini, params_ini, params)
    }else if(ini == 'kmeans'){
      params <- bernoulliEM_ini(xdata = X, k = K, ini = ini)
      p <- params$p
      a <- params$a
      rm(params)
    }else if(ini == 'kmeanspp'){
      params <- bernoulliEM_ini(xdata = X, k = K, ini = ini)
      p <- params$p
      a <- params$a
      rm(params)
    }
  
    res$alphas <- a
    res$probs <- p
    
    # EM  ----
    max_iter <- m_iter; n_iter <- 0; delta = 1
    while(delta>eps && n_iter != max_iter) {
      n_iter = n_iter + 1
      
      p.old <- p
      a.old <- a
      
      ##--- E step ---##
      
      p <- (p+p.old)/2
      
      p <- replace(p, p > 1, 1-1e-5)
      p <- replace(p, p == 0, 0+1e-5)
      
      p.not.stand <- tcrossprod(X, log(p)) + tcrossprod((1-X), (log(1-p)))
      p.not.stand <- sweep(p.not.stand, MARGIN = 2, log(a), FUN = '+')
      
      p_nk <- exp(sweep(p.not.stand, 1, rowLogSumExps(p.not.stand), '-'))
      
      ##--- M step ---##
      parameters <- Bernoulli_M_step(X, p_nk)
      p <- parameters$p
      a <- parameters$a
      
      delta <- abs(sum(p.old - p)) + abs(sum(a.old - a))
      
    }
    
    
    
    # Output prep ----
    profil <- apply(p_nk, 1, which.max)
    
    
    logL <- tcrossprod(X,log(p))+tcrossprod((1-X),(log(1-p)))
    logL <- sweep(logL, MARGIN = 2, log(a), FUN = '+')
    logL <- sum(rowLogSumExps(logL))
    
    W = dim(p)[2]+1 # 1 stends for alpha for each k
    n = nrow(X)
    
    BIC = -2*logL+((W*K)-1)*log(n)
    
    res$clusters <- profil
    res$bic <- BIC
    res$ll <- logL
    res$delta <- delta
    res$itr <- n_iter

   return(results = res)
} 
