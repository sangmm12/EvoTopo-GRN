
gene_cluster <- data.frame(gene=rownames(res_25[[1]]$cluster2),
                           cluster=res_25[[1]]$cluster2$`apply(omega, 1, which.max)`)


write.csv(gene_cluster,file="gene_cluster.csv",q=T)
saveRDS(gene_cluster,file="gene_cluster.rds")
