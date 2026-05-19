setwd("~/projects/wangxinfeng/CML/BIC")


library(patchwork)
library(glmnet)
library(parallel)
library(orthopolynom)
library(deSolve)
library(igraph)
library(qgraph)


source("~/projects/Brain/PCNSL/base.R")
source("~/projects/Brain/MPNST/BIC2/result/ode_1.R")






load("res_line_N.RData")
load("res_line_CML.RData")
load("res_line_R.RData")

load("res_par_N.RData")
load("res_par_CML.RData")
load("res_par_R.RData")


set.seed(1234)
result1 = list(original_data=res_line_N,  
               power_par = res_par_N,
               power_fit = res_line_N, 
               Time = as.numeric(colnames(res_line_N)))
cl <- parallel::makeCluster(2, setup_strategy = "sequential")
res_ode_N = try(qdODE_parallel(result = result1, thread = cl, reduction = T))

save(res_ode_N,file="res_ode_N.RData")



result1 = list(original_data=res_line_CML,  
               power_par = res_par_CML,
               power_fit = res_line_CML, 
               Time = as.numeric(colnames(res_line_CML)))
cl <- parallel::makeCluster(2, setup_strategy = "sequential")
res_ode_CML = try(qdODE_parallel(result = result1, thread = cl, reduction = T))


save(res_ode_CML,file="res_ode_CML.RData")




result1 = list(original_data=res_line_R,  
               power_par = res_par_R,
               power_fit = res_line_R, 
               Time = as.numeric(colnames(res_line_R)))
cl <- parallel::makeCluster(2, setup_strategy = "sequential")
res_ode_R = try(qdODE_parallel(result = result1, thread = cl, reduction = T))


save(res_ode_R,file="res_ode_R.RData")





