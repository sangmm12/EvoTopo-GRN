
setwd("~/projects/wangxinfeng/CML/BIC")



sum_fun <- function(x){
  dd <- c()
  for(i in 1:dim(x)[2]){
    d1 <- sum(na.omit(x[,i]))
    dd <- c(dd,d1)
  }
  dd
}


power_equation <- function(par,x){
  y=par[1]*x^par[2]
  y
}


dat_N <- readRDS("~/projects/wangxinfeng/CML/result_fit_Normal.rds")
dat_CML <- readRDS("~/projects/wangxinfeng/CML/result_fit_CML.rds")
dat_R <- readRDS("~/projects/wangxinfeng/CML/result_fit_remission.rds")



data_N <- dat_N$original_data
data_CML <- dat_CML$original_data
data_R <- dat_R$original_data


Time_N <- sum_fun(data_N)
Time_CML <- sum_fun(data_CML)
Time_R <- sum_fun(data_R)


TT_N <- seq(min(log10(Time_N)),max(log10(Time_N)),length = 20)
TT_CML <- seq(min(log10(Time_CML)),max(log10(Time_CML)),length = 20) 
TT_R <- seq(min(log10(Time_R)),max(log10(Time_R)),length = 20) 



load("~/projects/wangxinfeng/CML/BIC/res_25.RData")



cluster <- res_25[[1]]$cluster2$`apply(omega, 1, which.max)`
uni_cluster <- sort(unique(res_25[[1]]$cluster2$`apply(omega, 1, which.max)` ))




res_line_N <- c()
res_line_CML <- c()
res_line_R <- c()

for(i in 1:length(uni_cluster)){
  
  dat_N <- data_N[cluster==uni_cluster[i],]
  dat_CML <- data_CML[cluster==uni_cluster[i],]
  dat_R <- data_R[cluster==uni_cluster[i],]
  
  par1 <- res_25[[1]]$mu_par[,,1]
  par2 <- res_25[[1]]$mu_par[,,2]
  par3 <- res_25[[1]]$mu_par[,,3]
  
  
  line_N <- power_equation(par1[uni_cluster[i],],TT_N)
  line_CML <- power_equation(par2[uni_cluster[i],], TT_CML)
  line_R <- power_equation(par3[uni_cluster[i],], TT_R)
  
  
  res_line_N <- rbind(res_line_N,line_N)
  res_line_CML <- rbind(res_line_CML,line_CML)
  res_line_R <- rbind(res_line_R,line_R)
}



colnames(res_line_N) <- TT_N
rownames(res_line_N) <- paste0("M",1:25)

colnames(res_line_CML) <- TT_CML
rownames(res_line_CML) <- paste0("M",1:25)

colnames(res_line_R) <- TT_R
rownames(res_line_R) <- paste0("M",1:25)


save(res_line_N,file="res_line_N.RData")
save(res_line_CML,file="res_line_CML.RData")
save(res_line_R,file="res_line_R.RData")

res_par_N <- res_25[[1]]$mu_par[,,1]
res_par_CML <- res_25[[1]]$mu_par[,,2]
res_par_R <- res_25[[1]]$mu_par[,,3]

res_par_N <- res_par_N[uni_cluster,]
res_par_CML <- res_par_CML[uni_cluster,]
res_par_R <- res_par_R[uni_cluster,]

rownames(res_par_N) <- paste0("M",1:25)
rownames(res_par_CML) <- paste0("M",1:25)
rownames(res_par_R) <- paste0("M",1:25)


save(res_par_N,file="res_par_N.RData")
save(res_par_CML,file="res_par_CML.RData")
save(res_par_R,file="res_par_R.RData")


