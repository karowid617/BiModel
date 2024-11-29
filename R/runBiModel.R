library(dplyr)
library(purrr)
library(ggplot2)

runBiModel <- function(X, fixed = TRUE, K=2, start_ini = 20, ini = "random", 
                       m_iter=1000, eps=1e-6, IC = "BIC", 
                       quick_stop = FALSE, signi = 0.05,
                       plot = FALSE){
  if (!hasArg("X")){
    stop("No data.")}
  
  if (length(X) < 2){
    stop("Not enough data.")}

  
  IC_list <- c("AIC","AICc", "BIC")
  if (!opts$IC %in% IC_list) {
    stop("Criterion not implemented. Please use AIC, AICc, BIC.")
  }
  
  if(fixed == FALSE){
    
    if (K < 2){
      stop("K (no of components) must be larger than 1.")}
    
    res <- list()
    k <- 2
    stop <- FALSE
    
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
    
    if(plot == TRUE){
      tmp <- res %>% bind_rows() %>% split.default(names(.)) %>% map(na.omit)
      tmp$ic$ic
      
      df <- data.frame("K" = 1:k, "IC" = tmp$ic$ic)
      
      
      plt <- ggplot(df, aes(x = K, y = IC))+
        geom_line(color = 'grey', linewidth = 1)+
        geom_point(color = 'red', size = 2)+
        geom_vline(xintercept = min(df$IC), color = 'blue', linetype = 'dashed')+
        scale_x_continuous(breaks=seq(2, k, 1))+
        xlab("Number of components (K)")+
        ylab(IC)+
        theme_bw()+
        ggtitle(paste0(IC," values for different number of clusters"))
      
      print(plt)
    }
  }else{
    res <- BernoulliEM(X, K, start_ini, ini, m_iter, eps, IC)
  }
  
  return(results = res)
}