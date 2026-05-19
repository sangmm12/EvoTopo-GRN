setwd("~/projects/wangxinfeng/CML/BIC")



#library(showtext)
#showtext_auto(enable=TRUE)
#font_add("Times New Roman","/System/Library/Fonts/Supplemental/Times New Roman.ttf")
#font_add("Times New Roman1",regular = "/System/Library/Fonts/Supplemental/Times New Roman Italic.ttf")



links_N <- read.csv("fig2_Normal_FromTo.csv",row.names = 1)
colnames(links_N) <- c("incoming","outgoing","all")


name <- paste0("M",1:25)#sort(unique(c(rownames(links_N),rownames(links_T))))


n3 <- 1:25

links_N_1 <- links_N[match(name,rownames(links_N)),]



pdf("links_layer1_Normal.pdf",width = 4.8,height = 5)

ab <-5

num <- 25

par(mar=c(0,0,0,0),oma=c(0,0,0,0))
plot(c(0,0), c(0,0), type="n",xaxt="n",yaxt="n",frame=FALSE,xlab="",ylab="",xlim=c(5,25),ylim=c(-1, num+5))


segments(5,num+1,20,num+1,lwd=2)
segments(20,0,20,num+1,lwd=2)


M_name <-  name

for(i in 1:num){
  text(20.2,num+0.5-i,M_name[n3][i],cex=1,font=1,adj=0)
}




for(i in 1:num){
  sub_rc <-c(20,num-i+1,20-links_N_1[n3[i],3], num-i)
  rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1,col="#F0FFF0")
  
  sub_rc <-c(20,num-i+1,20-links_N_1[n3[i],2], num-i)
  rect(sub_rc[1],sub_rc[2],sub_rc[3],sub_rc[4],border="black",lwd=1,col="#FFC0CB")
  
}




for(i in 0:4){
  text(20-3*i,num+2,i*3,cex=1.2,font=1)
  segments(20-3*i,num+1,20-3*i,num+1-0.5)
}



text(20,-1,"Normal",cex=1.2,adj=1,font=1)

dev.off()




