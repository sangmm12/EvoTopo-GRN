
rm(list=ls())

# install.packages("readxl")
setwd("~/projects/wangxinfeng/CML/BIC")    #设置工作目录

# 读取CSV文件
data <- read.csv("fig2_Normal.csv")  

# 读取txt文件中的基因名
genes_txt <- paste0("M",1:25)#readLines("20PCD.txt")  
# tail(genes_txt)

# 提取两列基因名数据
From1 <- data$from  
To2 <- data$to  

# 计算txt文件中基因在两列基因名中出现的次数
gene_From <- table(factor(From1, levels = genes_txt))
gene_To <- table(factor(To2, levels = genes_txt))

# 将结果转换为数据框
gene_From_df1 <- as.data.frame(gene_From)
gene_To_df2 <- as.data.frame(gene_To)
print(gene_From_df1)
print(gene_To_df2)

mixed <- merge(gene_From_df1, gene_To_df2, by = "Var1")
print(mixed)

# 修改列名 
# 将第二列行名“Var1”改为“gene”;
# 将第三列行名“Freq.x”改为“from”。;
# 将第四列行名“Freq.y”改为“to”。
colnames(mixed) <- c("gene", "from", "to")
# 添加第五列内容,计算from to 的总和
mixed$total <- mixed$from + mixed$to
# 删除total为0的行
mixed_filtered <- mixed[mixed$total != 0, ]
# 输出为CSV表格
write.csv(mixed_filtered, file = "fig2_Normal_FromTo.csv", row.names = FALSE)




