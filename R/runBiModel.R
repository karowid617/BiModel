
runBiModel <- function(X, fixed = TRUE, K=2, start_ini = 20, ini = "random", 
                       m_iter=1000, eps=1e-6, IC = "BIC", 
                       quick_stop = FALSE, signi = 0.05){
  if (!hasArg("X")){
    stop("No data.")}
  
  if (length(X) < 2){
    stop("Not enough data.")}
  
  if (opts$KS < 2){
    stop("KS (no of components) must be larger than 1.")}
  
  IC_list <- c("AIC","AICc", "BIC")
  if (!opts$IC %in% IC_list) {
    stop("Criterion not implemented. Please use AIC, AICc, BIC")
  }
  
  if(fixed == FALSE){
    res <- list()
    k <- 2
    stop <- FALSE
    # for(k in 2:K){
    while(k<=K && stop == FALSE){
      res[[k-1]] <- BernoulliEM(X, k, start_ini, ini, m_iter, eps, IC) 
      
      if(quick_stop == TRUE && length(res)>1){
        Loglik.test <- -2*(res[[k-1]]$ll - res[[k-2]]$ll)
        p.val <- pchisq(Loglik.test, df = res[[k-2]]$n_params - res[[k-1]]$n_params, lower.tail = FALSE)
        if(p.val<signi){
          stop = TRUE
        }
        k = k+1
      }
    }
    names(res) <- paste0("K.", 2:k)
  }else{
    res <- BernoulliEM(X, K, start_ini, ini, m_iter, eps, IC)
  }
  
  return(results = res)
}