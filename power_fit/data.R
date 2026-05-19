setwd("~/projects/wangxinfeng/CML")

data <- read.table("GSE144119_raw_counts_GRCh38.p13_NCBI.tsv",header = T,row.names = 1)
sample <- read.csv("sample.csv")

sample_CML <- sample$X[which(sample$Stage=="chronic phase")]
sample_remission <- sample$X[which(sample$Stage=="remission")]
sample_Normal <- sample$X[which(sample$Stage=="n/a")]



data_CML <- data[,sapply(sample_CML,function(x){
  grep(x,colnames(data))
})]


data_remission <- data[,sapply(sample_remission,function(x){
  grep(x,colnames(data))
})]


data_Normal <- data[,sapply(sample_Normal,function(x){
  grep(x,colnames(data))
})]


library(data.table)
anno <- fread("Human.GRCh38.p13.annot.tsv",header=T)

anno <- anno[-c(2032,2152),]

which(rownames(data_CML)=="107985614")
which(rownames(data_CML)=="107985615")
which(rownames(data_CML)=="107985753")

data_CML <- data_CML[-c(2032,2152),]
data_remission <- data_remission[-c(2032,2152),]
data_Normal <- data_Normal[-c(2032,2152),]


rownames(data_CML) <-  anno$Symbol[match(rownames(data_CML),anno$GeneID)]
rownames(data_remission) <-  anno$Symbol[match(rownames(data_remission),anno$GeneID)]
rownames(data_Normal) <-  anno$Symbol[match(rownames(data_Normal),anno$GeneID)]


saveRDS(data_CML,file="data_CML.rds")
saveRDS(data_remission,file="data_remission.rds")
saveRDS(data_Normal,file="data_Normal.rds")





