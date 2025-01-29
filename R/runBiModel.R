#' Function to perform binomial mixture model (BMM)
#'
#' Function description
#'
#'
#' @param X Matrix of binary data to cluster where in rows there are features and in columns observations.
#' @param K Number of clusters to divide your data. Default K=2.
#' @param start_ini Number of starting initialization.
#' @param ini Method for parameters initialization. There are three options "random" (default), "kmeans" or "kmeanspp".
#' @param m_iter Maximum number of iteration for EM algorithm. Default value 1000.
#' @param eps Minimum delta of model LL to stop EM algorithm.Default value 1e-6.
#' @param IC Information criterion used to select the number of model components. Possible methods are "AIC","AICc", "BIC" (default).
#' @param quick_stop Logical value. Determines if stop searching of the number of components earlier based on the Likelihood Ratio Test. Used to speed up the function (TRUE, by default).
#' @param signi Significance level set for Likelihood Ratio Test (0.05, by default).
#' @param fixed Logical value. Fit BMM for selected number of components given by K (FALSE, by default).
#' @param plot Logical value. If TRUE, the IC plot will be displayed (FALSE, by default).
#' 
#' @returns Function returns a \code{list} which contains: \describe{
#'  \item{probs}{Matrix of finall probabilities for each observation and cluster.}
#'  \item{alphas}{Finall vector of each cluster weight in Bernoulli mixture model.}
#'  \item{clusters}{Cluster assigment of each observation.}
#'  \item{iter}{Number of EM interations.}
#'  \item{delta}{Final EM step delta fo algorithm stop.}
#'  \item{bic}{Bayesian information criterion (BIC) value for fitted model.}
#'  \item{ll}{Log-Likelihood of fitted model in clustering.}
#'  \item{n_params}{Number of parameters of fitted model.}
#' }
#' 
#' @import ggplot2
#' @importFrom methods hasArg
#' 
#' @examples
#' \dontrun{
#' data(example)
#' res<-BernoulliEM(example$nouli_data, 4, start_ini = 20, ini = "random", m_iter=3000, eps=1e-40)
#' }
#' 
#' @seealso \code{\link{BernoulliEM}}
#'
#' @export



# library(dplyr)
# library(purrr)
# library(ggplot2)

runBiModel <- function(X, K=2, start_ini = 20, ini = "random", 
                       m_iter=1000, eps=1e-6, IC = "BIC", 
                       signi = 0.05, quick_stop = TRUE,
                       fixed = FALSE, plot = FALSE){
  print("test")
  if (!hasArg("X")){
    stop("No data.")}
  
  if (length(X) < 2){
    stop("Not enough data.")}
  
  
  IC_list <- c("AIC","AICc", "BIC")
  if (!IC %in% IC_list) {
    stop("Criterion not implemented. Please use AIC, AICc, BIC.")
  }
  
  if(fixed == FALSE){
    
    if (K < 2){
      stop("K (No. of components) must be larger than 1.")}
    
    res <- list()
    k <- 2
    stop <- FALSE
    
    while(k<=K && stop == FALSE){
      cat("Component:", k, '\n')
      res[[k-1]] <- BernoulliEM(X, k, start_ini, ini, m_iter, eps, IC) 
      
      if(quick_stop == TRUE && length(res)>1){
        Loglik.test <- -2*(res[[k-2]]$ll - res[[k-1]]$ll)
        p.val <- pchisq(Loglik.test, df = res[[k-1]]$n_params - res[[k-2]]$n_params, lower.tail = FALSE)
        if(p.val<signi){
          stop = TRUE
        }
      }
      k = k+1
    }

    names(res) <- paste0("K.", 2:(k-1))
    
    ic_vec <- c()
    # ll_vec <- c()
    for(i in names(res)){
      ic_vec <- c(ic_vec, res[[i]]$ic)
      # ll_vec <- c(ll_vec, res[[i]]$ll)
    }
    
    # df <- data.frame("K" = 2:(k-1), "IC" = ic_vec, "LogLik" = ll_vec)
    df <- data.frame("K" = 2:(k-1), "IC" = ic_vec)
    
    
    plt <- ggplot(df, aes(x = K, y = IC))+
      geom_line(color = 'steelblue2', linewidth = 1.5)+
      geom_point(color = 'steelblue4', size = 3)+
      geom_vline(xintercept = df$K[which.min(df$IC)], color = 'red3', linetype = 'dashed', linewidth = 1)+
      scale_x_continuous(breaks=seq(2, k, 1))+
      xlab("Number of components (K)")+
      ylab(IC)+
      theme_bw()+
      ggtitle(paste0(IC," values for different number of clusters"))
    
    if(plot == TRUE){
      print(plt)
    }
    
    res <- c(res, plot = plt)
    
  }else{
    res <- BernoulliEM(X, K, start_ini, ini, m_iter, eps, IC)
  }
  
  return(results = res)
}
