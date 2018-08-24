gene_table<-read.delim("gene_info_20170509.txt",header = T,sep = "\t")


keep_mNAS<-gene_table$mNas<1e-2
summary(HeLa_mNas_Low_keep)

summary(gene_table$mNas)
summary(gene_table$iNas)

library(ggplot2)
g2<-ggplot(gene_table[IM_keep,],mapping =aes(log(iNas+1e-5,base = 10),log(mNas+1e-5,base = 10),color=IM_HeLa))
g2<-g2+stat_density2d(aes(fill=..density..),geom="raster", contour=FALSE)
g2


g2<-ggplot(gene_table[IM_keep,],mapping =aes(log(iNas+1e-5,base = 10),log(mNas+1e-5,base = 10),color=IM_HeLa))
g2<-g2+geom_point()+stat_density2d(aes(alpha=..density..),geom="raster", contour=FALSE)
g2



g2<-ggplot(gene_table,mapping =aes(log(iNas+1e-5,base = 10),log(mNas+1e-5,base = 10),color=IM_HeLa))
g2<-g2+geom_point()+stat_density2d(aes(alpha=..density..),geom="raster", contour=F)
g2


g2<-ggplot(gene_table,mapping =aes(log(mean_MP+1e-5,base = 10),log(mNas+1e-5,base = 10),color=IM_HeLa))
g2<-g2+geom_point()+stat_density2d(aes(alpha=..density..),geom="raster", contour=F)
g2<-g2+facet_grid(.~IM_HeLa)
g2



gene_table$mean_MP


HeLa_mNas_Low_keep<-gene_table$mNas<1e-3
IM_keep<-gene_table$IM_HeLa=="IM_CAR"

HeLa_mNas_Low<-gene_table[HeLa_mNas_Low_keep&IM_keep,]



g2<-ggplot(HeLa_mNas_Low,mapping =aes(IM_CAR,fill=gene_type))
g2<-g2+geom_bar()
g2




