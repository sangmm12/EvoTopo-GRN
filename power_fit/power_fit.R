setwd("~/projects/wangxinfeng/CML")



result_fit <- readRDS("result_fit_Normal.rds")

dat_fit <- result_fit$power_fit

dat_fit_N = melt(dat_fit,id=c("cluster"))

colnames(dat_fit_N)[1:2] <- c("gene","time")
dat_fit_N$group <- "1Normal"


dat_N <- result_fit$trans_data
dat_N <- as.matrix(dat_N)

dat_raw_N = melt(dat_N,id=c("cluster"))

colnames(dat_raw_N)[1:2] <- c("gene","time")
dat_raw_N$group <- "1Normal"





result_fit <- readRDS("result_fit_CML.rds")

dat_fit <- result_fit$power_fit

dat_fit_T = melt(dat_fit,id=c("cluster"))

colnames(dat_fit_T)[1:2] <- c("gene","time")

dat_fit_T$group <- "2CML"


dat_T <- result_fit$trans_data
dat_T <- as.matrix(dat_T)

dat_raw_T = melt(dat_T,id=c("cluster"))

colnames(dat_raw_T)[1:2] <- c("gene","time")

dat_raw_T$group <- "2CML"





result_fit <- readRDS("result_fit_remission.rds")

dat_fit <- result_fit$power_fit

dat_fit_R = melt(dat_fit,id=c("cluster"))

colnames(dat_fit_R)[1:2] <- c("gene","time")

dat_fit_R$group <- "3Remission"


dat_R <- result_fit$trans_data
dat_R <- as.matrix(dat_R)

dat_raw_R = melt(dat_R,id=c("cluster"))

colnames(dat_raw_R)[1:2] <- c("gene","time")

dat_raw_R$group <- "3Remission"





dat_fit <- rbind(dat_fit_N,dat_fit_T,dat_fit_R)
dat_raw <- rbind(dat_raw_N,dat_raw_T,dat_raw_R)





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





label=10


p1 <- list()
for(i in c(20,100,200,
           300,410,520,
           610,800)){
  
  region1 <- unique(dat_raw$gene)[i]
  
  
  dat_raw_plot <- subset(dat_raw,subset = gene ==region1)
  dat_line_plot <- subset(dat_fit,subset = gene ==region1)
  
  
  p1[[i]] <-   ggplot() +  
    geom_point(dat_raw_plot, mapping = aes(x = time, y = value, colour=group,fill=group,alpha = 1/10),size=3,shape=1) +
    
    geom_line(dat_line_plot, mapping = aes(x = time, y = value,  colour=group),size=2) +
    
    facet_wrap(~group,scales = "free_x",ncol=3
    ) + scale_fill_manual(values= c( "#38E54D",   "#FF8787",   "#80B1D3"),aesthetics = c("colour", "fill"))+
    coord_cartesian(ylim=c(0, 3)) +#scale_y_continuous(limits = c(0,15)) + # ,labels = ylabels
    xlab("Expression Index") + ylab("Individual Expression of Each Gene")+ 
    ggtitle(region1)+theme_zg1()
  
  
  # if (is.null(label)) {
  # p1[[i]] = p1[[i]]
  # } else {c(10,12,14)#
  xlabel = c(6.8,7,7.2,7.4)
  #ggplot_build(p1[[i]])$layout$panel_params[[1]]$x.sec$breaks
  ylabel = ggplot_build(p1[[i]])$layout$panel_params[[1]]$y.sec$breaks#c(2,2.5,3,3.5,4,4.5,5,5.5,6,7,8)
  #ggplot_build(p1[[i]])$layout$panel_params[[1]]$y.sec$breaks
  
  xlabel2 = parse(text= paste(label,"^", xlabel, sep="") )
  
  # if (ylabel[1] == 0) {
  # ylabel2 = parse(text=c(0,paste(label,"^", ylabel[2:length(ylabel)], sep="")))
  # } else{
  ylabel2 = parse(text= paste(label,"^", ylabel, sep="") )
  # }
  p1[[i]] = p1[[i]] + scale_x_continuous(breaks=xlabel,labels = xlabel2) + scale_y_continuous(breaks=ylabel,labels = ylabel2)
  #  }
  
}






library(cowplot)
library(dplyr)
p <- plot_grid(p1[[20]], p1[[100]],
               p1[[200]],p1[[300]],
               p1[[410]], p1[[520]],
               p1[[610]],p1[[800]],
               ncol=2)

pdf("power_fit.pdf",height = 1.5*3*3,width=5*3)

p
dev.off()


