setwd("~/projects/wangxinfeng/CML/BIC/M19/sM4/gene")


library(patchwork)
library(glmnet)
library(parallel)
library(orthopolynom)
library(deSolve)
library(igraph)
library(qgraph)


source("~/projects/Brain/PCNSL/base.R")
source("~/projects/Brain/MPNST/BIC2/result/ode_1.R")





set.seed(1234)

res <- readRDS("res_N.rds")
result1 = res
cl <- parallel::makeCluster(2, setup_strategy = "sequential")
res_ode_N = try(qdODE_parallel(result = result1, thread = cl, reduction = T))

saveRDS(res_ode_N,file="res_ode_N.rds")


res <- readRDS("res_CML.rds")
result1 = res
cl <- parallel::makeCluster(2, setup_strategy = "sequential")
res_ode_CML = try(qdODE_parallel(result = result1, thread = cl, reduction = T))


saveRDS(res_ode_CML,file="res_ode_CML.rds")



res <- readRDS("res_R.rds")
result1 = res
cl <- parallel::makeCluster(2, setup_strategy = "sequential")
res_ode_R = try(qdODE_parallel(result = result1, thread = cl, reduction = T))


saveRDS(res_ode_R,file="res_ode_R.rds")



