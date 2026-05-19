
setwd("~/projects/wangxinfeng/CML/BIC")



library(igraph)

load("res_ode_N.RData")

result <-  res_ode_N

ode_result <-  c()
for(i in 1:length(result$ode_result)){
  d1 <- result$ode_result[[i]]$fit
  d2 <- d1[,4:c(dim(d1)[2])]
  ode_result[[i]] <- list(fit.df=d2)
}




network_conversion <- function(result, j = NULL){
  n = ncol(result$fit.df)
  lop_order = length(result$parameters)/n
  
  if ( is.null(j) == TRUE ) {
    effect.mean = apply(result$fit.df,2,mean)
  } else {
    effect.mean = unlist(result$fit.df[j,])
  }
  output <- list(ind.name = colnames(result$fit.df)[1],
                 dep.name = colnames(result$fit.df)[-1],
                 #ind.par = result$parameters[1:lop_order],
                 #dep.par = matrix(result$parameters[-(1:lop_order)],nrow = (n-1)),
                 effect.mean = effect.mean,
                 effect.all = result$fit.df2)
  
  return(output)
}


net_all = list(mean.net = lapply(ode_result,network_conversion))



replace_character <- function(i){
  require(stringr)
  tmp <- str_remove_all(i, "]")
  tmp1 <- str_remove(tmp, "\\(")
  tmp2 <- as.numeric(unlist(strsplit(tmp1,split='\\,')))
  return(tmp2)
}

get_after <- function(i){
  temp <- matrix(NA,nrow = length(i[[2]]),ncol=3)
  temp[,1] <- i[[2]]
  temp[,2] <- i[[1]]
  temp[,3] <- i[[3]][2:(length(i[[2]])+1)]
  
  colnames(temp) <- c('from','to','dep_effect')
  temp <- data.frame(temp)
  temp[,3] <- as.numeric(as.character(temp[,3]))
  return(temp)
}

get_max_effect <- function(k){
  after <- do.call(rbind,lapply(k, get_after))
  max_dep_effect <- max(abs(after$dep_effect))
  
  temp <- aggregate(dep_effect ~ to, data = after, sum)
  all_dep_effect <- max(abs(temp$dep_effect))
  return(c(max_dep_effect,all_dep_effect))
}

max_effect = t(sapply(net_all,get_max_effect))

k=net_all$mean.net
color_size=2



extra <- sapply(k,function(c) c[[3]][1])

after <- do.call(rbind,lapply(k, get_after))

colfunc <- colorRampPalette(c("#0095EF",
                              "#ffdead",
                              "#FE433C"))

links_col = data.frame(color = colfunc(color_size),
                       y = cut(seq(-max(max_effect[,1])-0.1,max(max_effect[,1])+0.1,length=color_size),color_size))

links_interval = t(sapply(1:color_size,function(c) replace_character(links_col$y[c])))

links <- after
colnames(links) <- c("from","to","weight")
for (i in 1:nrow(links)) {
  links$edge.colour[i] <- links_col$color[which(sapply(1:color_size,function(c)
    findInterval(links$weight[i],c(links_interval[c,])))==1)] #add colour for links
}


#links$edge.colour[grep("H",links$from)] <- "blue"
#links$edge.colour[grep("N",links$from)] <- "green3"


node_col <- data.frame(color = colfunc(color_size),
                       y = cut(seq(-max(max_effect[,2])-0.1,max(max_effect[,2])+0.1,length=color_size),color_size))
node_interval = t(sapply(1:color_size,function(c) replace_character(node_col$y[c])))
nodes <- data.frame(unique(links[,2]),unique(links[,2]),extra)
colnames(nodes) <- c("id","name","ind_effect")
node2 <- aggregate(weight ~ to, data = links, sum)
nodes$influence  <- node2[match(nodes$id,node2[,1]),2]



#nodes$colour[grep("H",nodes$name)] <- "#3DB2FF"
#nodes$colour[grep("N",nodes$name)] <- "green3"


for (i in 1:nrow(nodes)) {
  nodes$colour[i] <- node_col$color[which(sapply(1:color_size,function(c)
    findInterval(nodes$influence[i],c(node_interval[c,])))==1)] #add colour for links
}

#nodes$name = c("H1.a","H1.b","N1.a","N1.b","N1.c","N1.d")


#normalization
normalization <- function(x){(x-min(x))/(max(x)-min(x))*1.2+0.8}

#final plot
links[,3] <- normalization(abs(links[,3]))
nodes[,3:4] <- normalization(abs(nodes[,3:4]))



net <- graph_from_data_frame( d=links,vertices = nodes,directed = T )




pdf("fig2_Normal.pdf")
plot.igraph(net,
            rescale = T,
            #edge.arrow.size=E(net)$weight,
            vertex.label=V(net)$name,
            vertex.label.color="black",
            vertex.shape="circle",
            vertex.label.cex=1.5,#V(net)$ind_effect*1,
            vertex.size=V(net)$ind_effect*10,
            edge.curved=0.2,
            edge.color=E(net)$edge.colour,
            edge.frame.color=E(net)$edge.colour,
            edge.width=E(net)$weight*2,
            vertex.color="grey",#V(net)$colour,
            vertex.frame.color="grey",#V(net)$colour,
            layout=layout_in_circle,
            #main=title,
            margin=c(-0.05,-0.05,-0.05,-0.05)
)
dev.off()

write.csv(after,file="fig2_Normal.csv",quote=F,row.names = F)

