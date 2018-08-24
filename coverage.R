

# exon coverage analysis 


options(width=40)
library(CoverageView)

#read IM——CAR information from bed file
A549_IM<-system.file("extdata","A549_IM.bed",package="CoverageView")




BC1Bam<-"~/IonXpressRNA_001/BC1.STARBowtie2.bam"
BC2Bam<-"~/IonXpressRNA_002/BC2.STARBowtie2.bam"
BC3Bam<-"~/IonXpressRNA_003/BC3.STARBowtie2.bam"
BC4Bam<-"~/IonXpressRNA_004/BC4.STARBowtie2.bam"

BC5Bam<-"~/IonXpressRNA_005/IonXpressRNA_005.STARBowtie2.bam"
BC6Bam<-"~/IonXpressRNA_006/IonXpressRNA_006.STARBowtie2.bam"
BC7Bam<-"~/IonXpressRNA_007/IonXpressRNA_007.STARBowtie2.bam"
BC8Bam<-"~/IonXpressRNA_008/IonXpressRNA_008.STARBowtie2.bam"

BC1Bam<-CoverageBamFile(BC1Bam,reads_mapped=27753502)
BC2Bam<-CoverageBamFile(BC2Bam,reads_mapped = 44945452)
BC3Bam<-CoverageBamFile(BC3Bam,reads_mapped = 36169044)
BC4Bam<-CoverageBamFile(BC4Bam,reads_mapped=48672937)

BC5Bam<-CoverageBamFile(BC5Bam,reads_mapped=44636078)
BC6Bam<-CoverageBamFile(BC6Bam,reads_mapped=44818427)
BC7Bam<-CoverageBamFile(BC7Bam,reads_mapped=33485718)
BC8Bam<-CoverageBamFile(BC8Bam,reads_mapped=36169044)



BC8_HeLaIMCOV<-cov.matrix(BC8Bam,coordfile=A549_IM,no_windows=100, offset=10,num_cores=6,normalization="rpm")



BC4_A549MOCOV<-as.data.frame(BC4_A549MOCOV)
BC8_HeLaMOCOV<-as.data.frame(BC8_HeLaMOCOV)

data<-BC4_A549MOCOV
data<-BC8_HeLaMOCOV

for(i in 1:nrow(data)){
  meani<-mean(as.numeric(data[i,(1:ncol(data))]))
  if(i==1){
    mean<-meani
  }
  else{
    mean<-c(mean,meani)
  }
}


BC4_A549MOCOV_mean<-cbind(bin=(1:nrow(BC4_A549MOCOV)),mean,sample="A549_exoncunt_multi")
BC8_HeLaMOCOV_mean<-cbind(bin=(1:nrow(BC8_HeLaMOCOV)),mean,sample="HeLa_exoncunt_multi")

write.table(BC4_A549MOCOV_mean, file="BC4_coverage on A549 MO_CAR multiexon.txt",sep = "\t",quote = F)
write.table(BC8_HeLaMOCOV_mean, file="BC8_coverage on HeLa MO_CAR multiexon.txt",sep = "\t",quote = F)

BC4_A549MOCOV_mean<-read.delim(file="BC4_coverage on A549 MO_CAR multiexon.txt",header = T,sep = "\t")

BC8_HeLaMOCOV_mean<-read.delim(file="BC8_coverage on HeLa MO_CAR multiexon.txt",header = T,sep = "\t")


data<-rbind(BC4_A549MOCOV_mean,BC8_HeLaMOCOV_mean)

library(ggplot2)
g2 <- ggplot(data,mapping = aes(bin,mean, group=sample))
#g2 <- ggplot(data,mapping = aes(bin,mean))
g2 <- g2 + geom_line(aes(colour = sample))
#g2 <- g2 + geom_line()
g2<-g2+geom_vline(xintercept = c(11,21,31,41,51,61,71,81,91,101,111))
g2 <- g2 + theme(legend.direction = 'horizontal', 
                 legend.position = 'top',
                 legend.text = element_text(size = 8))
g2<-g2+labs(y = "total number of genes")
g2<-g2+ggtitle("gene type of common and specific  CARs")
print(g2)


