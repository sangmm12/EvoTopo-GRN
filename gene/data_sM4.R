
setwd("~/projects/wangxinfeng/CML/BIC/M19/sM4")


load("~/projects/wangxinfeng/CML/BIC/M19/res_3.RData")



cluster <- res_3[[2]]$cluster2$`apply(omega, 1, which.max)`
uni_cluster <- sort(unique(res_3[[2]]$cluster2$`apply(omega, 1, which.max)` ))



data_N <- readRDS("~/projects/wangxinfeng/CML/BIC/M19/data_N_M19.rds")
data_CML <- readRDS("~/projects/wangxinfeng/CML/BIC/M19/data_CML_M19.rds")
data_R <- readRDS("~/projects/wangxinfeng/CML/BIC/M19/data_R_M19.rds")

data_N_sM4 <- data_N[cluster==uni_cluster[4],]
data_CML_sM4 <- data_CML[cluster==uni_cluster[4],]
data_R_sM4 <- data_R[cluster==uni_cluster[4],]


saveRDS(data_N_sM4,file="data_N_sM4.rds")
saveRDS(data_CML_sM4,file="data_CML_sM4.rds")
saveRDS(data_R_sM4,file="data_R_sM4.rds")



