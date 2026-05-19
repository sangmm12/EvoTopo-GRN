
setwd("~/projects/wangxinfeng/CML/BIC/M1")

dat_N <- readRDS("data_N_M1.rds")
data_N <- dat_N

data1_N <- data_N




load("res.RData")



cluster <- res[[19]]$cluster2$`apply(omega, 1, which.max)`
uni_cluster <- sort(unique(res[[19]]$cluster2$`apply(omega, 1, which.max)` ))






num_N <- c()

for(i in 1:21){
  dat_N_1 <- data_N[which(cluster==uni_cluster[i]),]
  
  if(is.null(dim(dat_N_1)[1])==F){
    num_N <- c(num_N,dim(dat_N_1)[1])
  }else{
    num_N <- c(num_N,1)
  }
  
}


load("res_line_N.RData")
load("res_line_CML.RData")
load("res_line_R.RData")

load("res_par_N.RData")
load("res_par_CML.RData")
load("res_par_R.RData")


text_name <- data.frame(value = paste0(paste0("sM",1:21)," (",num_N,")"))



data_N <-  readRDS("data_N_M1.rds")
data_N <- as.matrix(data_N)


data_CML <-  readRDS("data_CML_M1.rds")
data_CML <- as.matrix(data_CML)

data_R <-  readRDS("data_R_M1.rds")
data_R <- as.matrix(data_R)




Time_N <- sum_fun(data_N)
Time_CML <- sum_fun(data_CML)
Time_R <- sum_fun(data_R)



colnames(data_N) <- Time_N
colnames(data_CML) <- Time_CML
colnames(data_R) <- Time_R



dat_raw_plot_N <- c()
dat_raw_plot_CML <- c()
dat_raw_plot_R <- c()


for(i in 1:21){
  
  dat_N_1 <- data_N[cluster==uni_cluster[i],]
  dat_CML_1 <- data_CML[cluster==uni_cluster[i],]
  dat_R_1 <- data_R[cluster==uni_cluster[i],]
  
  if(is.null(dim(dat_N_1))==T){
    
    dat_N_2 <- matrix(dat_N_1,nrow=1)
    colnames(dat_N_2) <- as.numeric(colnames(data_N))
    rownames(dat_N_2) <- rownames(data_N)[cluster==uni_cluster[i]]
    
    
    dat_CML_2 <- matrix(dat_CML_1,nrow=1)
    colnames(dat_CML_2) <- as.numeric(colnames(data_CML))
    rownames(dat_CML_2) <- rownames(data_CML)[cluster==uni_cluster[i]]
    
    
    dat_R_2 <- matrix(dat_R_1,nrow=1)
    colnames(dat_R_2) <- as.numeric(colnames(data_R))
    rownames(dat_R_2) <- rownames(data_R)[cluster==uni_cluster[i]]
    
    dat_1_plot.N <- reshape2::melt(dat_N_2)
    dat_1_plot.CML <- reshape2::melt(dat_CML_2)
    dat_1_plot.R <- reshape2::melt(dat_R_2)
  }else{
    dat_1_plot.N <- reshape2::melt(dat_N_1)
    dat_1_plot.CML <- reshape2::melt(dat_CML_1)
    dat_1_plot.R <- reshape2::melt(dat_R_1)
  }
  
  
  M_name <- i
  
  dat_2_plot.N <- data.frame("Time" = log10(dat_1_plot.N$Var2),
                             "Value" = log10(dat_1_plot.N$value+1),
                             "Cluster" = rep(M_name,dim(dat_1_plot.N)[1]),
                             "Age" = rep("1Normal",dim(dat_1_plot.N)[1]),
                             "Genus"=dat_1_plot.N$Var1)
  
  
  dat_2_plot.CML <- data.frame("Time" = log10(dat_1_plot.CML$Var2),
                               "Value" = log10(dat_1_plot.CML$value+1),
                               "Cluster" = rep(M_name,dim(dat_1_plot.CML)[1]),
                               "Age" = rep("2CML",dim(dat_1_plot.CML)[1]),
                               "Genus"=dat_1_plot.CML$Var1)
  
  
  dat_2_plot.R <- data.frame("Time" = log10(dat_1_plot.R$Var2),
                             "Value" = log10(dat_1_plot.R$value+1),
                             "Cluster" = rep(M_name,dim(dat_1_plot.R)[1]),
                             "Age" = rep("3Remission",dim(dat_1_plot.R)[1]),
                             "Genus"=dat_1_plot.R$Var1)
  
  dat_raw_plot_N <- rbind(dat_raw_plot_N,dat_2_plot.N)
  dat_raw_plot_CML <- rbind(dat_raw_plot_CML,dat_2_plot.CML)
  dat_raw_plot_R <- rbind(dat_raw_plot_R,dat_2_plot.R)
}


dat_raw_plot <- rbind(dat_raw_plot_N,dat_raw_plot_CML,dat_raw_plot_R)




dat_line_plot_N <- c()
dat_line_plot_CML <- c()
dat_line_plot_R <- c()

for(i in 1:21){
  dat_N_1 <- res_line_N[i,]
  dat_CML_1 <- res_line_CML[i,]
  dat_R_1 <- res_line_R[i,]
  dat_1_plot.N <- reshape2::melt(dat_N_1)
  dat_1_plot.CML <- reshape2::melt(dat_CML_1)
  dat_1_plot.R <- reshape2::melt(dat_R_1)
  
  M_name <- i
  
  dat_2_plot.N <- data.frame("Time" = as.numeric(rownames(dat_1_plot.N)),
                             "Value" = dat_1_plot.N$value,
                             "Cluster" = rep(M_name,dim(dat_1_plot.N)[1]),
                             "Age" = rep("1Normal",dim(dat_1_plot.N)[1]))
  
  
  dat_2_plot.CML <- data.frame("Time" = as.numeric(rownames(dat_1_plot.CML)),
                               "Value" = dat_1_plot.CML$value,
                               "Cluster" = rep(M_name,dim(dat_1_plot.CML)[1]),
                               "Age" = rep("2CML",dim(dat_1_plot.CML)[1]))
  
  dat_2_plot.R <- data.frame("Time" = as.numeric(rownames(dat_1_plot.R)),
                             "Value" = dat_1_plot.R$value,
                             "Cluster" = rep(M_name,dim(dat_1_plot.R)[1]),
                             "Age" = rep("3Remission",dim(dat_1_plot.R)[1]))
  
  dat_line_plot_N <- rbind(dat_line_plot_N,dat_2_plot.N)
  dat_line_plot_CML <- rbind(dat_line_plot_CML,dat_2_plot.CML)
  dat_line_plot_R <- rbind(dat_line_plot_R,dat_2_plot.R)
  
}


dat_line_plot <- rbind(dat_line_plot_N,dat_line_plot_CML,dat_line_plot_R)


library(ggplot2)
library(reshape2)

load("res.RData")


result <- res#[[9]]$cluster_number
best.k <- 21

cluster.result = result[[which(sapply( result , "[[" , 'cluster_number' )==best.k)]]
kk = length(table(cluster.result$cluster$`apply(omega, 1, which.max)`))
if ( kk!= best.k) stop("Please use a smaller k or rerun functional clustering")

times1 = cluster.result$Time1
times2 = cluster.result$Time2
times3 = cluster.result$Time3
times1_new = seq(min(times1),max(times1),length = 20)
times2_new = seq(min(times2),max(times2),length = 20)
times3_new = seq(min(times3),max(times3),length = 20)

n1 = length(times1);n2 = length(times2);n3 = length(times3)
timesall = c(min(c(times1,times2,times3)),max(c(times1,times2,times3)))



par.mu = cluster.result$mu_par
k = cluster.result$cluster_number
alpha = as.numeric(table(cluster.result$cluster$`apply(omega, 1, which.max)`))

power_equation <- function(x, power_par){ t(sapply(1:nrow(power_par),
                                                   function(c) power_par[c,1]*x^power_par[c,2] ) )}


mu.fit1 = power_equation(times1_new, par.mu[,,1][1:k,])
mu.fit2 = power_equation(times2_new, par.mu[,,2][1:k,])
mu.fit3 = power_equation(times3_new, par.mu[,,2][1:k,])

colnames(mu.fit1) = times1_new
mu.fit1 = melt(as.matrix(mu.fit1))
colnames(mu.fit1) = c("cluster","x","y")
mu.fit1$x = as.numeric(as.character(mu.fit1$x))

colnames(mu.fit2) = times2_new
mu.fit2 = melt(as.matrix(mu.fit2))
colnames(mu.fit2) = c("cluster","x","y")
mu.fit2$x = as.numeric(as.character(mu.fit2$x))

colnames(mu.fit3) = times3_new
mu.fit3 = melt(as.matrix(mu.fit3))
colnames(mu.fit3) = c("cluster","x","y")
mu.fit3$x = as.numeric(as.character(mu.fit3$x))

#mu.fit1$type = "L"
#mu.fit2$type = "R"
#mu.fit = rbind(mu.fit1,mu.fit2)

plot.df1 = cluster.result$cluster[,c(1:n1,ncol(cluster.result$cluster))]
X1 = plot.df1[,-ncol(plot.df1)]
plot.df1$name = rownames(plot.df1)
colnames(plot.df1) = c(times1,"cluster","name")

plot.df1 = melt(plot.df1, id.vars = c('cluster',"name"))
colnames(plot.df1) = c("cluster","name", "x","y")
plot.df1$x = as.numeric(as.character(plot.df1$x))
plot.df1$name = paste0("a_",plot.df1$name)

plot.df2 = cluster.result$cluster[,((n1+1):ncol(cluster.result$cluster))]
X2 = plot.df2[,-ncol(plot.df2)]
plot.df2$name = rownames(plot.df2)
colnames(plot.df2) = c(times2,"cluster","name")

plot.df2 = melt(plot.df2, id.vars = c('cluster',"name"))
colnames(plot.df2) = c("cluster","name", "x","y")
plot.df2$x = as.numeric(as.character(plot.df2$x))
plot.df2$name = paste0("b_",plot.df2$name)


plot.df3 = cluster.result$cluster[,((n1+n2+1):ncol(cluster.result$cluster))]
X3 = plot.df3[,-ncol(plot.df3)]
plot.df3$name = rownames(plot.df3)
colnames(plot.df3) = c(times3,"cluster","name")

plot.df3 = melt(plot.df3, id.vars = c('cluster',"name"))
colnames(plot.df3) = c("cluster","name", "x","y")
plot.df3$x = as.numeric(as.character(plot.df3$x))
plot.df3$name = paste0("c_",plot.df3$name)



alpha = as.numeric(table(cluster.result$cluster$`apply(omega, 1, which.max)`))


name.df1 = data.frame(label = paste0("M",1:best.k,"(",alpha ,")"),
                      x = mean(timesall), y = max(X1,X2)*0.9, Cluster = 1:best.k)




name.df <- name.df1#[1:6,]



label <- 10


color1 = "#38E54D"
color2 = "#FF8787"
color3 = "#80B1D3"


p<- ggplot() +  geom_point(dat_raw_plot,  mapping = aes_string(x = "Time", y = "Value",
                                                               colour = factor(dat_raw_plot$Age),
                                                               alpha = factor(rep(0.1,length(dat_raw_plot$Cluster)))),size=1.25) +
  geom_line(dat_line_plot, mapping = aes_string(x = "Time", y = "Value", 
                                                colour =  factor(dat_line_plot$Age)),size=1.25) +
  facet_wrap(~Cluster)  +
  xlab("Expression Index") + ylab("Individual Expression of Each Gene") + theme(axis.title=element_text(size=18)) +
  geom_text(name.df, mapping = aes_string(x = "x", y = "y",label = "label"),
            check_overlap = TRUE, size = 4) +
  theme_bw() + scale_color_manual(name = "Cluster",values = c("1Normal" = color1,
                                                              "2CML" = color2,
                                                              "3Remission" = color3))+
  guides(alpha = "none") +
  theme(axis.title=element_text(size=15),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10,hjust = 0),
        panel.spacing = unit(0.0, "lines"),
        plot.margin = unit(c(1,1,1,1), "lines"),
        strip.background = element_blank(),
        plot.background = element_blank(),
        strip.text = element_blank())

if (is.null(label)) {
  p = p
} else {
  xlabel = ggplot_build(p)$layout$panel_params[[1]]$x.sec$breaks
  ylabel = ggplot_build(p)$layout$panel_params[[1]]$y.sec$breaks
  
  xlabel2 = parse(text= paste(label,"^", xlabel, sep="") )
  
  #if (ylabel[1] == 0) {
  ylabel2 = parse(text=c(0,paste(label,"^", ylabel[2:length(ylabel)], sep="")))
  # } else{
  # ylabel2 = parse(text= paste(label,"^", ylabel, sep="") )
  # }
  p = p + scale_x_continuous(labels = xlabel2) + scale_y_continuous(labels = ylabel2)
}
#return(p)


pdf("fig1.pdf",height = 7.2,width = 10)

p

dev.off()


