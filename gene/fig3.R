setwd("~/projects/wangxinfeng/CML/BIC/M19/sM4/gene")

library(ggrepel)

dat <- readRDS("~/projects/wangxinfeng/CML/BIC/M19/sM4/data_N_sM4.rds")


anno <- rownames(dat)

text_name2 <- c()
for(i in 1:length(anno)){
  anno1 <- anno[i]
  an1 <- paste0(anno1," (",i,")")
  text_name2 <- c(text_name2,an1)
}

text_name1 <- data.frame(value=text_name2,
                         raw=anno,
                         num=1:length(anno))






res_ode_N <- readRDS("res_ode_N.rds")


res_ode_CML <- readRDS("res_ode_CML.rds")


res_ode_R <- readRDS("res_ode_R.rds")



region <- c(1:21)


result <-  res_ode_N

raw_plot <-  c()
ind_plot <- c()
dep_plot <- c() 
text_name <- c()
for(i in 1:21){
  d1 <- result$ode_result[[i]]$fit
  
  d_raw <- data.frame(
    value=d1$y,
    times=d1$x
  )
  
  name1 <- colnames(d1)[4]
  
  n1 <- which(text_name1$raw==name1)
  
  
  
  d_raw$group <- region[i]
  d_raw$cluster <- region[i]
  
  raw_plot <-  rbind(raw_plot,d_raw)
  
  n2 <- dim(d1)[2]
  
  if((n2-4) >= 2){
    dd <- as.matrix(d1[,5:dim(d1)[2]])
  }else{
    dd <- as.matrix(d1[,5])
    colnames(dd) <- colnames(d1)[5]
  }
  
  
  dd <- t(dd)
  colnames(dd) <- d1$x
  d_dep <- reshape2::melt(dd)
  colnames(d_dep) <- c("group","times","value")
  d_dep$cluster <- region[i]
  dep_plot <-  rbind(dep_plot,d_dep)
  
  
  y1 <- dd[,20]
  x1 <- rep(d1$x[16],length(y1))
  cell_na <- rownames(dd)
  
  num <- sapply(cell_na,function(x){text_name1$num[text_name1$raw==x]})
  clu_na <- rep(text_name1$num[n1],dim(dd)[1])
  
  dd_name <- data.frame(x=x1,y=y1,Value=cell_na,num=num,cluster=clu_na)
  
  text_name <- rbind(text_name,dd_name)
  
  d_ind <- data.frame(
    value=d1[,4],
    times=d1$x
  )
  d_ind$group <- region[i]
  d_ind$cluster <- region[i]
  
  ind_plot <-  rbind(ind_plot,d_ind)
  
}

raw_plot_Normal <- raw_plot
ind_plot_Normal <- ind_plot
dep_plot_Normal <- dep_plot
text_name_Normal <- text_name

raw_plot_Normal$gender <- "1Normal"
ind_plot_Normal$gender <- "1Normal"
dep_plot_Normal$gender <- "1Normal"
text_name_Normal$gender <- "1Normal"


result <-  res_ode_CML

raw_plot <-  c()
ind_plot <- c()
dep_plot <- c() 
text_name <- c()
for(i in 1:21){
  d1 <- result$ode_result[[i]]$fit
  
  d_raw <- data.frame(
    value=d1$y,
    times=d1$x
  )
  
  name1 <- colnames(d1)[4]
  
  n1 <- which(text_name1$raw==name1)
  
  
  
  d_raw$group <- region[i]
  d_raw$cluster <- region[i]
  
  raw_plot <-  rbind(raw_plot,d_raw)
  
  n2 <- dim(d1)[2]
  
  if((n2-4) >= 2){
    dd <- as.matrix(d1[,5:dim(d1)[2]])
  }else{
    dd <- as.matrix(d1[,5])
    colnames(dd) <- colnames(d1)[5]
  }
  
  
  dd <- t(dd)
  colnames(dd) <- d1$x
  d_dep <- reshape2::melt(dd)
  colnames(d_dep) <- c("group","times","value")
  d_dep$cluster <- region[i]
  dep_plot <-  rbind(dep_plot,d_dep)
  
  
  y1 <- dd[,20]
  x1 <- rep(d1$x[16],length(y1))
  cell_na <- rownames(dd)
  
  num <- sapply(cell_na,function(x){text_name1$num[text_name1$raw==x]})
  clu_na <- rep(text_name1$num[n1],dim(dd)[1])
  
  dd_name <- data.frame(x=x1,y=y1,Value=cell_na,num=num,cluster=clu_na)
  
  text_name <- rbind(text_name,dd_name)
  
  d_ind <- data.frame(
    value=d1[,4],
    times=d1$x
  )
  d_ind$group <- region[i]
  d_ind$cluster <- region[i]
  
  ind_plot <-  rbind(ind_plot,d_ind)
  
}



raw_plot_Tumor <- raw_plot
ind_plot_Tumor <- ind_plot
dep_plot_Tumor <- dep_plot
text_name_Tumor <- text_name

raw_plot_Tumor$gender <- "2CML"
ind_plot_Tumor$gender <- "2CML"
dep_plot_Tumor$gender <- "2CML"
text_name_Tumor$gender <- "2CML"




result <-  res_ode_R

raw_plot <-  c()
ind_plot <- c()
dep_plot <- c() 
text_name <- c()
for(i in 1:21){
  d1 <- result$ode_result[[i]]$fit
  
  d_raw <- data.frame(
    value=d1$y,
    times=d1$x
  )
  
  name1 <- colnames(d1)[4]
  
  n1 <- which(text_name1$raw==name1)
  
  
  
  d_raw$group <- region[i]
  d_raw$cluster <- region[i]
  
  raw_plot <-  rbind(raw_plot,d_raw)
  
  n2 <- dim(d1)[2]
  
  if((n2-4) >= 2){
    dd <- as.matrix(d1[,5:dim(d1)[2]])
  }else{
    dd <- as.matrix(d1[,5])
    colnames(dd) <- colnames(d1)[5]
  }
  
  
  dd <- t(dd)
  colnames(dd) <- d1$x
  d_dep <- reshape2::melt(dd)
  colnames(d_dep) <- c("group","times","value")
  d_dep$cluster <- region[i]
  dep_plot <-  rbind(dep_plot,d_dep)
  
  
  y1 <- dd[,20]
  x1 <- rep(d1$x[16],length(y1))
  cell_na <- rownames(dd)
  
  num <- sapply(cell_na,function(x){text_name1$num[text_name1$raw==x]})
  clu_na <- rep(text_name1$num[n1],dim(dd)[1])
  
  dd_name <- data.frame(x=x1,y=y1,Value=cell_na,num=num,cluster=clu_na)
  
  text_name <- rbind(text_name,dd_name)
  
  d_ind <- data.frame(
    value=d1[,4],
    times=d1$x
  )
  d_ind$group <- region[i]
  d_ind$cluster <- region[i]
  
  ind_plot <-  rbind(ind_plot,d_ind)
  
}



raw_plot_R <- raw_plot
ind_plot_R <- ind_plot
dep_plot_R <- dep_plot
text_name_R <- text_name

raw_plot_R$gender <- "3Remission"
ind_plot_R$gender <- "3Remission"
dep_plot_R$gender <- "3Remission"
text_name_R$gender <- "3Remission"






raw_plot <- rbind(raw_plot_Tumor,raw_plot_Normal,raw_plot_R)
ind_plot <- rbind(ind_plot_Tumor,ind_plot_Normal,ind_plot_R)
dep_plot <- rbind(dep_plot_Tumor,dep_plot_Normal,dep_plot_R)
text_name <- rbind(text_name_Tumor,text_name_Normal,text_name_R)




library(ggplot2)





theme_zg1 <- function(..., bg='white'){
  require(grid)
  theme_classic(...,base_family="serif") +
    theme(rect=element_rect(fill=bg),
          plot.margin=unit(c(0.5,0.5,0.5,1.0), 'lines'),
          panel.background=element_rect(fill='#FFF5F0', color='black'),
          panel.border=element_rect(fill='transparent', color='transparent'),
          panel.grid=element_blank(),
          plot.title=element_text(size=25,vjust=2),
          axis.title = element_text(color='black',size=20),
          axis.title.x = element_text(size=20,vjust=0),
          axis.title.y = element_text(size=20,vjust=1.5,margin = margin(r = 8.8)),
          axis.ticks.length = unit(-0.4,"lines"),
          axis.text.x=element_text(hjust=0.5,size=20,margin = margin(t = 8.8)),
          axis.text.y=element_text(angle = 0,hjust=0.5,size=20,margin = margin(r = 8.8)),
          axis.ticks = element_line(color='black'),
          legend.position='none',
          legend.title=element_blank(),
          legend.key=element_rect(fill='transparent', color='transparent'),
          strip.background=element_rect(fill='transparent', color='transparent'),
          strip.text=element_text(color='white',vjust=-1.8,hjust=0.06,size=20),
          strip.text.x=element_text(color='white',vjust=-1.8,hjust=0.06,size=20)
    )
  
}


dep_plot$cluster <- factor(dep_plot$cluster,levels=region)
ind_plot$cluster <- factor(ind_plot$cluster,levels=region)
raw_plot$cluster <- factor(raw_plot$cluster,levels=region)
text_name$cluster <- factor(text_name$cluster,levels=region)




data_label <- data.frame(x=6.5,y=3,label="Gene")



label=10


p1 <- list()
for(i in 1:length(region)){
  region1 <- region[i]
  
  text_name_2 <- text_name1[i,1]
  
  scaleFUN <- function(x) sprintf("%.0f", x)
  
  dep_plot_1 <- subset(dep_plot,subset = cluster ==region1)
  ind_plot_1 <- subset(ind_plot,subset = cluster ==region1)
  raw_plot_1 <- subset(raw_plot,subset = cluster ==region1)
  text_name_1 <- subset(text_name,subset = cluster ==region1)
  
  p1[[i]] <- ggplot() + geom_line(dep_plot_1, mapping = aes(x = times, y = value, group= group),col="green3", size=1.25)+#alpha = 0.25,
    geom_line(ind_plot_1, mapping = aes(x = times, y = value, group= group),col="blue3", size=1.25)+
    geom_line(raw_plot_1, mapping = aes(x = times, y = value, group= group),col="red3", size=1.25)+
    geom_text(text_name_1,mapping=aes(x,y, label= num),col="black",alpha = 0.75,size=8)+
    facet_wrap(~gender, nrow = 1, ncol =3
               ,scales = "free_x"
    )+ coord_cartesian(ylim=c(-5, 5))+
    xlab("Expression Index") + ylab("Expression of Individual Gene")+ 
    #geom_text_repel(data=data_label,aes(x= 3.6,y= 5,label="Gene"),
      #             force=1,point.padding=unit(1,"lines"),
        #            vjust=1,direction="y",nudge_x = 0.1)+
  #  geom_text_repel(data=data_label,aes(x= 3.6,y= -5,label="Gene"),
          #          force=1,point.padding=unit(1,"lines"),
         #           vjust=1,direction="y",nudge_x = 0.1)+
    
    ggtitle(text_name_2)+theme_zg1()
  
  
  # if (is.null(label)) {
  # p1[[i]] = p1[[i]]
  # } else {
  xlabel = c(3.2,3.6,3.7,4,4.4,4.8)#ggplot_build(p1[[i]])$layout$panel_params[[1]]$x.sec$breaks
  ylabel = ggplot_build(p1[[i]])$layout$panel_params[[1]]$y.sec$breaks
  
  xlabel2 = parse(text= paste(label,"^", xlabel, sep="") )
  
  #if (ylabel[1] == 0) {
  #ylabel2 = parse(text=c(0,paste(label,"^", ylabel[2:length(ylabel)], sep="")))
  # } else{
  ylabel2 = parse(text= paste(label,"^", ylabel, sep="") )
  # }
  p1[[i]] = p1[[i]] + scale_x_continuous(breaks=xlabel,labels = xlabel2) + scale_y_continuous(breaks=ylabel,labels = ylabel2)
  #  }
  
}










library(cowplot)
library(dplyr)
p <- plot_grid(p1[[1]], p1[[2]],p1[[3]],
               p1[[4]],p1[[5]],p1[[6]],
               p1[[7]],p1[[8]],p1[[9]],
               p1[[10]],p1[[11]],p1[[12]],
               p1[[13]],p1[[14]],p1[[15]],
               p1[[16]],p1[[17]],p1[[18]],
               p1[[19]],p1[[20]],p1[[21]],
               ncol=4)

pdf("fig3-1.pdf",height = 1.5*2.4*4,width=8*6)

p
dev.off()





