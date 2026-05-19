setwd("~/projects/wangxinfeng/CML")




dat_CML <- readRDS("data_CML.rds")
dat_CML <- as.matrix(dat_CML)
dat_CML <- t(dat_CML)

dat_remission <- readRDS("data_remission.rds")
dat_remission <- as.matrix(dat_remission)
dat_remission <- t(dat_remission)


dat_Normal <- readRDS("data_Normal.rds")
dat_Normal <- as.matrix(dat_Normal)
dat_Normal <- t(dat_Normal)



del_num_Normal <- c()
for(i in 1:dim(dat_Normal)[2]){
  n1 <- length(which(dat_Normal[,i]==0))
  if(n1 > dim(dat_Normal)[1]*0.7){
    del_num_Normal <- c( del_num_Normal,i)
  }
}


del_num_CML <- c()
for(i in 1:dim(dat_CML)[2]){
  n1 <- length(which(dat_CML[,i]==0))
  if(n1 > dim(dat_CML)[1]*0.7){
    del_num_CML <- c( del_num_CML,i)
  }
}


del_num_remission <- c()
for(i in 1:dim(dat_remission)[2]){
  n1 <- length(which(dat_remission[,i]==0))
  if(n1 > dim(dat_remission)[1]*0.7){
    del_num_remission <- c( del_num_remission,i)
  }
}




del_num <- unique(c(del_num_Normal,del_num_CML,del_num_remission))

dat1_Normal <- dat_Normal[,-del_num]
dat1_CML <- dat_CML[,-del_num]
dat1_remission <- dat_remission[,-del_num]




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




data_Normal <- t(dat1_Normal)

cl <- parallel::makeCluster(2, setup_strategy = "sequential")

result_fit <- power_equation_fit(data_Normal, n=30, trans = log10, thread = cl)

saveRDS(result_fit,file="result_fit_Normal.rds")


data_CML <- t(dat1_CML)

cl <- parallel::makeCluster(2, setup_strategy = "sequential")

result_fit <- power_equation_fit(data_CML, n=30, trans = log10, thread = cl)

saveRDS(result_fit,file="result_fit_CML.rds")




data_remission <- t(dat1_remission)

cl <- parallel::makeCluster(2, setup_strategy = "sequential")

result_fit <- power_equation_fit(data_remission, n=30, trans = log10, thread = cl)

saveRDS(result_fit,file="result_fit_remission.rds")

