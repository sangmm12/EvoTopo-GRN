setwd("~/projects/wangxinfeng/CML/BIC/M19/sM4/gene")




library(pbapply)
library(parallel)
library(ggplot2)





power_equation <- function(x, power_par){ t(sapply(1:nrow(power_par),
                                                   function(c) power_par[c,1]*x^power_par[c,2] ) )}


power_equation_base <- function(x, y){
  x1 <- x
  y1 <- y
  nn <- which(is.na(y1)==T)
  
  if(length(nn)==0){
    x <- as.numeric(x1)
    y <- as.numeric(y1)
  }else{
    x <- as.numeric(x1[-nn])
    y <- as.numeric(y1[-nn])
  }
  
  min_value = min(y[y!=0])
  
  lmFit <- lm( log( y + runif(1, min = 0, max = min_value))  ~ log(x))
  coefs <- coef(lmFit)
  a <- exp(coefs[1])
  b <- coefs[2]
  
  model <- try(nls(y~a*x^b,start = list(a = a, b = b),
                   control = nls.control(maxiter = 10000000, minFactor = 1e-200)))
  if( 'try-error' %in% class(model)) {
    result = NULL
  }
  else{
    result = model
  }
  return(result)
}




power_equation_all <- function(x,y, maxit=1e4){
  result <- power_equation_base(x,y)
  iter <- 1
  while( is.null(result) && iter <= maxit) {
    iter <- iter + 1
    try(result <- power_equation_base(x,y))
  }
  return(result)
}

#x=X
#y=trans_data[51,]
power_equation_fit <- function(data, n=30, trans = log10, thread = 2) {
  data = data[,order(colSums(data))]
  if ( is.null(trans)) {
    X = rowSums(data)
    trans_data = data
  } else{
    X = trans(colSums(data)+1)
    trans_data = trans(data+1)
  }
  
  colnames(trans_data) = X
  
  core.number <- thread
  cl <- makeCluster(getOption("cl.cores", core.number))
  #cl <- parallel::makeCluster(2, setup_strategy = "sequential")
  clusterExport(cl, c("power_equation_all", "power_equation_base", "trans_data", "X"), envir = environment())
  all_model = parLapply(cl = cl, 1:nrow(data), function(c) power_equation_all(X, trans_data[c,]))
  stopCluster(cl)
  
  
  
  names(all_model) = rownames(data)
  no = which(sapply(all_model, length)>=1)
  all_model2 = all_model[no]
  data2 = data[no,]
  trans_data2 = trans_data[no,]
  
  new_x = seq(min(X), max(X), length = n)
  power_par = t(vapply(all_model2, coef, FUN.VALUE = numeric(2), USE.NAMES = TRUE))
  power_fit = t(vapply(all_model2, predict, newdata = data.frame(x=new_x),
                       FUN.VALUE = numeric(n), USE.NAMES = TRUE))
  
  colnames(power_fit) = new_x
  result = list(original_data = data2, trans_data = trans_data2,
                power_par = power_par, power_fit = power_fit,
                Time = X)
  return(result)
}



dat1_N <- readRDS("~/projects/wangxinfeng/CML/BIC/M19/sM4/data_N_sM4.rds")
dat1_CML <- readRDS("~/projects/wangxinfeng/CML/BIC/M19/sM4/data_CML_sM4.rds")
dat1_R <- readRDS("~/projects/wangxinfeng/CML/BIC/M19/sM4/data_R_sM4.rds")

dat_N <- as.matrix(dat1_N)
dat_CML <- as.matrix(dat1_CML)
dat_R <- as.matrix(dat1_R)


res_N <- power_equation_fit(dat_N, n=20, trans = log10,  thread = 2)
res_CML <- power_equation_fit(dat_CML, n=20, trans = log10,  thread = 2)
res_R <- power_equation_fit(dat_R, n=20, trans = log10,  thread = 2)


saveRDS(res_N,file="res_N.rds")
saveRDS(res_CML,file="res_CML.rds")
saveRDS(res_R,file="res_R.rds")


