setwd("~/Documents/scientific/CAR/analysis")


library("cluster")
library("factoextra")



A549_merged_GeneCount<-read.delim(file = "A549_merged_GeneCount.txt",sep = "\t",header = T)
Hela_merged_GeneCount<-read.delim(file = "Hela_merged_GeneCount.txt",sep = "\t",header = T)
merged_GeneCount<-read.delim(file = "merged_GeneCount.txt",sep = "\t",header = T)

data<-merged_GeneCount[2:17]
row.names(data)<-merged_GeneCount$gene
#filt genes that RPM>10 in at least 5 samples
log=rep(x =0,times=nrow(data))
for(r in 1:nrow(data)){
  count=0
  for(c in 1:ncol(data)){
    if((data[r,c])/sum(data[[c]])>=10/1000000){
      count=count+1
    }
  }
  if(count>=4){
    log[r]=1
    next
  }
}

log<-as.logical(log)
my_data<-data[log,]
t_data<-t(my_data)



res.dist <- get_dist(t_data, stand = F, method = "spearman")
fviz_dist(res.dist, gradient = list(low = "#0000FF", high = "white"))

fviz_nbclust(t_data, kmeans, method = "gap_stat")

km.res <- kmeans(t_data, 2, nstart = 25)
#pam.res <- pam(t_data, 6)
fviz_cluster(km.res, data = t_data, ellipse.type = "convex")+  theme_minimal()

