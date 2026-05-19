setwd("~/projects/wangxinfeng/CML/BIC")


#load("res_5.RData")
load("res_10.RData")
load("res_15.RData")
load("res_20.RData")
load("res_21.RData")
load("res_22.RData")
load("res_23.RData")
load("res_24.RData")
load("res_25.RData")
load("res_26.RData")
load("res_27.RData")
load("res_28.RData")
load("res_29.RData")
load("res_30.RData")
load("res_35.RData")
load("res_40.RData")
load("res_45.RData")
load("res_50.RData")


res <- c()
#res[1] <- res_5
res[1] <- res_10
res[2] <- res_15
res[3] <- res_20
res[4] <- res_21
res[5] <- res_22
res[6] <- res_23
res[7] <- res_24
res[8] <- res_25
res[9] <- res_26
res[10] <- res_27
res[11] <- res_28
res[12] <- res_29
res[13] <- res_30
res[14] <- res_35
res[15] <- res_40
res[16] <- res_45
res[17] <- res_50



library(ggplot2)


fun_clu_BIC <- function(result, crit = "BIC", title = NULL){
  BIC.df = data.frame( k = sapply( result , "[[" , 'cluster_number' ),
                       BIC = sapply( result , "[[" , 'BIC' ) )
  AIC.df = data.frame( k = sapply( result , "[[" , 'cluster_number' ),
                       AIC = sapply( result , "[[" , 'AIC' ) )
  min.k = BIC.df$k[which.min(BIC.df$BIC)]
  min.k2 = AIC.df$k[which.min(AIC.df$AIC)]
  
  p = ggplot(BIC.df)+ theme_bw() + geom_line(mapping = aes_string(x = "k" ,y = "BIC"), color = 'grey10') +
    geom_vline(xintercept = min.k, linetype="dashed", color = "red", size=1.5) +
    xlab("K") + ylab("BIC") + ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5))
  
  p1 = ggplot(AIC.df)+ theme_bw() + geom_line(mapping = aes_string(x = "k" ,y = "AIC"), color = 'grey10') +
    geom_vline(xintercept = min.k2, linetype="dashed", color = "red", size=1.5) +
    xlab("K") + ylab("AIC") + ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5))
  if (crit != "BIC") {
    return(p1)
  } else{
    return(p)
  }
}



#fun_clu_BIC(result = res[5:10])

pdf("BIC.pdf",height = 5,width = 6)
fun_clu_BIC(result = res)
dev.off()


fun_clu_BIC(result = res)



BIC.df1 = data.frame( k = c(#length(unique(res_5[[1]]$cluster2$`apply(omega, 1, which.max)`)),
                            length(unique(res_10[[1]]$cluster2$`apply(omega, 1, which.max)`)),
                            length(unique(res_15[[1]]$cluster2$`apply(omega, 1, which.max)`)),
                            length(unique(res_20[[1]]$cluster2$`apply(omega, 1, which.max)`)),
                            length(unique(res_21[[1]]$cluster2$`apply(omega, 1, which.max)`)),
                            length(unique(res_22[[1]]$cluster2$`apply(omega, 1, which.max)`)),
                            length(unique(res_23[[1]]$cluster2$`apply(omega, 1, which.max)`)),
                            length(unique(res_24[[1]]$cluster2$`apply(omega, 1, which.max)`)),
                            length(unique(res_25[[1]]$cluster2$`apply(omega, 1, which.max)`)),  
                            length(unique(res_25[[1]]$cluster2$`apply(omega, 1, which.max)`)),
                            length(unique(res_26[[1]]$cluster2$`apply(omega, 1, which.max)`)),
                            length(unique(res_27[[1]]$cluster2$`apply(omega, 1, which.max)`)),
                            length(unique(res_28[[1]]$cluster2$`apply(omega, 1, which.max)`)),
                            length(unique(res_30[[1]]$cluster2$`apply(omega, 1, which.max)`)),
                            length(unique(res_35[[1]]$cluster2$`apply(omega, 1, which.max)`)),
                            length(unique(res_40[[1]]$cluster2$`apply(omega, 1, which.max)`)),
                            length(unique(res_45[[1]]$cluster2$`apply(omega, 1, which.max)`)),
                            length(unique(res_50[[1]]$cluster2$`apply(omega, 1, which.max)`))
),
BIC = sapply( res , "[[" , 'BIC' ) )




pdf("BIC.pdf",height = 5,width = 6)
p
dev.off()




BIC.df <- BIC.df1#[c(8,3,2,12,15),]

min.k = BIC.df$k[which.min(BIC.df$BIC)]



p = ggplot(BIC.df)+ theme_bw() + geom_line(mapping = aes_string(x = "k" ,y = "BIC"), color = 'grey10') +
  geom_vline(xintercept = min.k, linetype="dashed", color = "red", size=1.5) +
  xlab("K") + ylab("BIC") + #ggtitle(title) +
  theme(plot.title = element_text(hjust = 0.5))

pdf("BIC.pdf",height = 5,width = 6)
p
dev.off()








plot(BIC.df1)




pdf("BIC.pdf",height = 5,width = 6)
fun_clu_BIC(result = res[c(8,3,2,12,15)])
dev.off()




