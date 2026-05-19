setwd("~/projects/wangxinfeng/CML/BIC")





source("~/projects/TCGA/Death_LUAD/count/base.R")


library(pbapply)
library(parallel)
library(ggplot2)
library(mvtnorm)
library(reshape2)



biQ_function <- function(par, prob_log, omega_log, X, k, n1, n2,n3,
                         times1, times2, times3){
  n = dim(X)[1]; d = dim(X)[2];
  X1 = X[,1:n1]
  X2 = X[,(n1+1):(n1+n2)]
  X3 = X[,(n1+n2+1):(n1+n2+n3)]
  
  par.cov <- par[1:6]
  par.mu <- array(par[-c(1:6)], dim = c(k,2,3))
  
  cov1 = get_SAD1_covmatrix(par.cov[1:2], n1)
  cov2 = get_SAD1_covmatrix(par.cov[3:4], n2)
  cov3 = get_SAD1_covmatrix(par.cov[5:6], n3)
  
  mu1 <- power_equation(times1, par.mu[,,1][1:k,])
  mu2 <- power_equation(times2, par.mu[,,2][1:k,])
  mu3 <- power_equation(times3, par.mu[,,3][1:k,])
  
  mvn_log1 <- sapply(1:k, function(c) dmvnorm(X1, mu1[c,], cov1, log = T))
  mvn_log2 <- sapply(1:k, function(c) dmvnorm(X2, mu2[c,], cov2, log = T))
  mvn_log3 <- sapply(1:k, function(c) dmvnorm(X3, mu3[c,], cov3, log = T))
  
  
  mvn.log = mvn_log1 + mvn_log2 + mvn_log3
  tmp = sweep(mvn.log, 2, FUN = "+", STATS = prob_log) - omega_log
  Q = -sum(tmp*exp(omega_log))
  return(Q)
}


biget_par_int <- function(X, k, times1, times2,times3,
                          n1, n2, n3){
  n = dim(X)[1]; d = dim(X)[2];
  X1 = X[,1:n1]
  X2 = X[,(n1+1):(n1+n2)]
  X3 = X[,(n1+n2+1):(n1+n2+n3)]
  
  # cov.int = c(0.3,IQR(diag(cov(X1))), 0.3,IQR(diag(cov(X2))))
  
  cov.int = c(1e-5,mean(diag(cov(X1))),
              1e-5,mean(diag(cov(X2))),
              1e-5,mean(diag(cov(X3))))
  
  
  init.cluster <- kmeans(X,centers = k)
  prob <- table(init.cluster$cluster)/n
  
  fit1 <- lapply(1:k,function(c) power_equation_all(times1,init.cluster$centers[c,1:n1]))
  fit2 <- lapply(1:k,function(c) power_equation_all(times2,init.cluster$centers[c,(n1+1):(n1+n2)]))
  fit3 <- lapply(1:k,function(c) power_equation_all(times3,init.cluster$centers[c,(n1+n2+1):(n1+n2+n3)]))
  
  mu.par.int1 <- t(sapply(fit1, coef))
  mu.par.int2 <- t(sapply(fit2, coef))
  mu.par.int3 <- t(sapply(fit3, coef))
  
  return_obj <- list(initial_cov_params = cov.int,
                     initial_mu_params = c(mu.par.int1,mu.par.int2,mu.par.int3),
                     initial_probibality = prob)
  return(return_obj)
}


bifun_clu <- function(data1, data2,data3,
                      k, Time1 = NULL, Time2 = NULL,Time3 = NULL,trans=trans,
                      initial.pars = NULL, iter.max = 1e2, parscale = 1e-3){
  #sort from low to high
  data1 = data1[,order(colSums(data1))]
  data2 = data2[,order(colSums(data2))]
  data3 = data3[,order(colSums(data3))]
  
  
  #attribute
  data = as.matrix(cbind(data1,data2,data3)); n = dim(data)[1]; d = dim(data)[2]
  n1 = dim(data1)[2] 
  n2 = dim(data2)[2]
  n3 = dim(data3)[2]
  
  eplison = 1; iter = 0;
  
  
  if (is.null(Time1)|is.null(Time2)|is.null(Time3) ) {
    
    if (is.null(trans)) {
      times1 = as.numeric(colSums(data1))
      times2 = as.numeric(colSums(data2))
      times3 = as.numeric(colSums(data3))
      X = data
    } else{
      times1 = as.numeric(trans(colSums(data1)+1));
      times2 = as.numeric(trans(colSums(data2)+1));
      times3 = as.numeric(trans(colSums(data3)+1));
      X = trans(data+1)
    }
    
  } else{
    times1 = as.numeric(Time1)
    times2 = as.numeric(Time2)
    times3 = as.numeric(Time3)
  }
  
  # initial pars
  if (is.null(initial.pars)) {
    initial.pars = biget_par_int(X, k, times1, times2,times3,
                                 n1, n2, n3)
  }
  
  X1 = X[,1:n1]
  X2 = X[,(n1+1):(n1+n2)]
  X3 = X[,(n1+n2+1):(n1+n2+n3)]
  
  par.int <- c(initial.pars$initial_cov_params, initial.pars$initial_mu_params)
  prob_log <- log(initial.pars$initial_probibality)
  
  #parscale = par.int[2]*par.int[4]*parscale
  
  while( abs(eplison) > 1e-3 && iter <= iter.max ){
    #E step
    par.mu <- array(par.int[-c(1:6)], dim = c(k,2,3))
    par.cov = par.int[1:6]
    cov1 = get_SAD1_covmatrix(par.cov[1:2], n1)
    cov2 = get_SAD1_covmatrix(par.cov[3:4], n2)
    cov3 = get_SAD1_covmatrix(par.cov[5:6], n3)
    
    mu1 <- power_equation(times1, par.mu[,,1][1:k,])
    mu2 <- power_equation(times2, par.mu[,,2][1:k,])
    mu3 <- power_equation(times3, par.mu[,,3][1:k,])
    
    
    mvn_log1 <- sapply(1:k, function(c) dmvnorm(X1, mu1[c,], cov1, log = T))
    mvn_log2 <- sapply(1:k, function(c) dmvnorm(X2, mu2[c,], cov2, log = T))
    mvn_log3 <- sapply(1:k, function(c) dmvnorm(X3, mu3[c,], cov3, log = T))
    
    
    mvn_log = mvn_log1 + mvn_log2 + mvn_log3
    
    mvn = sweep(mvn_log, 2, FUN = '+', STATS =  prob_log )
    
    omega_log = t(sapply(1:n, function(c) mvn[c,] - logsumexp(mvn[c,]) ))
    omega = exp(omega_log)
    
    LL.mem <- biQ_function(par = par.int, prob_log = prob_log, omega_log, X, k, n1, n2,n3,
                           times1, times2, times3)
    
    #M step
    prob_exp = apply(omega_log, 2, logsumexp)
    prob_log = prob_exp - log(n)
    
    Q.maximization <- try(optim(par = par.int, biQ_function,
                                prob_log = prob_log,
                                omega_log = omega_log,
                                X = X,
                                k = k,
                                n1 = n1,
                                n2 = n2,
                                n3 = n3,
                                
                                times1 = times1,
                                times2 = times2,
                                times3 = times3,
                                method = "BFGS",
                                lower = c(-10,-10,-10,-10,-10,-10,rep(-Inf,6*k)),
                                upper = c(10,10,10,10,10,10,rep(Inf,6*k)),
                                control = list(trace = TRUE,
                                               parscale = c(rep(parscale,6),rep(1,6*k)),
                                               maxit = 1e2
                                )))
    if ('try-error' %in% class(Q.maximization))
      break
    par.hat <- Q.maximization$par
    par.int = par.hat
    LL.next <- biQ_function(par = par.int, prob_log = prob_log, omega_log, X, k, n1, n2,n3,
                            times1, times2, times3)
    eplison <-  LL.next - LL.mem
    LL.mem <- LL.next
    iter = iter + 1
    
    cat("\n", "iter =", iter, "\n", "Log-Likelihood = ", LL.next, "\n")
  }
  AIC = 2*(LL.next) + 2*(length(par.hat)+k-1)
  BIC = 2*(LL.next) + log(n)*(length(par.hat)+k-1)
  
  omega = exp(omega_log)
  X.clustered <- data.frame(X, apply(omega,1,which.max),check.names = F)
  
  return_obj <- list(cluster_number = k,
                     Log_likelihodd = LL.mem,
                     AIC = AIC,
                     BIC = BIC,
                     cov_par = par.cov,
                     mu_par = par.mu,
                     probibality = exp(prob_log),
                     omega = omega,
                     cluster = X.clustered,
                     cluster2 = data.frame(data, apply(omega,1,which.max), check.names = F),
                     Time1 = times1,
                     Time2 = times2,
                     Time3 = times3,
                     original_data = data)
  return(return_obj)
}



bifun_clu_parallel <- function(data1, data2, data3, 
                               Time1 = NULL, Time2 = NULL,Time3 = NULL,trans = trans, start, end, iter.max = 100, thread = 2){
  core.number <- thread
  cl <- makeCluster(getOption("cl.cores", core.number))
  clusterEvalQ(cl, {library(mvtnorm)})
  clusterExport(cl,c(c("bifun_clu","biget_par_int","biQ_function","logsumexp","power_equation_base",
                       "power_equation_all","power_equation","get_biSAD1","get_SAD1_covmatrix"),ls()),
                envir=environment())
  result <- parLapply(cl=cl, start:end, function(c)bifun_clu(data1 = data1,
                                                             data2 = data2,
                                                             data3 = data3,
                                                             
                                                             k = c,
                                                             Time1 = Time1,
                                                             Time2 = Time2,
                                                             Time3 = Time3,
                                                             trans = trans,
                                                             iter.max = iter.max))
  stopCluster(cl)
  
  return(result)
}


result_fit <- readRDS("~/projects/wangxinfeng/CML/result_fit_Normal.rds")
dat_N <- result_fit


result_fit <- readRDS("~/projects/wangxinfeng/CML/result_fit_CML.rds")
dat_CML <- result_fit


result_fit <- readRDS("~/projects/wangxinfeng/CML/result_fit_remission.rds")
dat_R <- result_fit



data1 = dat_N$original_data
data2 = dat_CML$original_data
data3 = dat_R$original_data



set.seed(1234)

res_5 = bifun_clu_parallel(data1 = data1, 
                           data2 = data2,
                           data3 = data3,
                           Time1 = NULL,
                           Time2 = NULL,
                           Time3 = NULL,
                           trans = log10,
                           start = 5, 
                           end = 5, 
                           thread = 9,
                           iter.max = 2)


save(res_5,file="res_5.RData")





