
# this script is for defining CARs from A549 and hela cells. first, all gene level count data were merged and filter out low abundant genes.



A549_merged_GeneCount<-read.delim(file = "A549_merged_GeneCount.txt",sep = "\t",header = T)
Hela_merged_GeneCount<-read.delim(file = "Hela_merged_GeneCount.txt",sep = "\t",header = T)

data<-cbind(A549_merged_GeneCount,Hela_merged_GeneCount[2:9])

row.names(data)<-Hela_merged_GeneCount$gene

#filt genes that RPM>10 in at least 4 samples
log=rep(x =0,times=nrow(data))
for(r in 1:nrow(data)){
  count=0
  for(c in 1:ncol(data)){
    if((data[r,c])/sum(data[[c]])>=10/1000000){
      count=count+1
    }
  }
  if(count>=2){
    log[r]=1
    next
  }
}


log<-as.logical(log)
my_data<-data[log,]


########## plot sample matrix
t_data<-t(my_data)
library("cluster")
library("factoextra")
res.dist <- get_dist(t_data, stand = T, method = "spearman")
fviz_dist(res.dist, gradient = list(low = "#0000FF", high = "white"))

######### EdgeR
library("edgeR")

A549_data<-my_data[1:8]
Hela_data<-my_data[9:16]
cell_data<-Hela_data

#group <- factor(c("IS","IP","MS","MP","IS","IP","MS","MP","IS","IP","MS","MP","IS","IP","MS","MP"))
group <- factor(c("IS","IP","MS","MP","IS","IP","MS","MP"))


y<- DGEList(counts=cell_data,group=group)
design<-model.matrix(~0+group,data=y$samples)
colnames(design) <- levels(y$samples$group)
y <- estimateGLMCommonDisp(y,design)
y <- estimateGLMTrendedDisp(y,design)
y <- estimateGLMTagwiseDisp(y,design)
fit<- glmFit(y, design)


########plot MDS and BCV
plotMDS(y,pch=10)
plotMDS(y,cex = 0.5)
plotBCV(y)
maPlot(fit$fitted.values[,1],fit$fitted.values[,2])


#Plot sample matrix based on fitted 
res.dist <- get_dist(t(fit$fitted.values), stand = T, method = "spearman")
fviz_dist(res.dist, gradient = list(low = "#0000FF", high = "white"))



### make contrast
 IPvsIS <- makeContrasts(IP-IS, levels=design)
 MPvsMS <- makeContrasts(MP-MS, levels=design)
 M_I_P<-MPvsMS-IPvsIS
 PvsS<-MPvsMS+IPvsIS
 
  lrt.IpvsIs <- glmLRT(fit, contrast=IPvsIS )
  lrt.MpvsMs <- glmLRT(fit, contrast=MPvsMS)
  lrt.M_I_P<-glmLRT(fit, contrast=M_I_P)
  lrt.PvsS<-glmLRT(fit, contrast=PvsS)
  

# decide 

Is.CAR_P1<-((lrt.PvsS$table$PValue<0.05)&(lrt.PvsS$table$logFC>0))
Is.CAR_P2<-((lrt.IpvsIs$table$PValue<0.05)&(lrt.IpvsIs$table$logFC>0))
Is.CAR_P3<-((lrt.MpvsMs$table$PValue<0.05)&(lrt.MpvsMs$table$logFC>0))
Is.CAR<-Is.CAR_P1|Is.CAR_P2|Is.CAR_P3
summary(Is.CAR)


detags <- rownames(y)[as.logical(Is.CAR)]
plotSmear(lrt.PvsS, de.tags=detags)
abline(h=c(-1, 1), col="blue")


colnames(lrt.IpvsIs$table)<-c("logFC_IpvsIs_Hela","logCPM_IpvsIs" ,"LR_IpvsIs_Hela","PValue_IpvsIs_Hela")
colnames(lrt.MpvsMs$table)<-c("logFC_MpvsMs_Hela","logCPM_MpvsMs" ,"LR_MpvsMs_Hela","PValue_MpvsMs_Hela")
colnames(lrt.M_I_P$table)<-c("logFC_MvsI_Hela","logCPM_MvsI" ,"LR_MvsI_Hela","PValue_MvsI_Hela")
colnames(lrt.PvsS$table)<-c("logFC_PvsS_Hela","logCPM_PvsS" ,"LR_PvsS_Hela","PValue_PvsS_Hela")

# annotate Class
CLass_hela<-rep ("Non",times=nrow(cell_data))
length(CLass_hela)

CLass_hela[!Is.CAR]<-"Non"



Hela_table<-cbind(my_data,lrt.IpvsIs$table,lrt.MpvsMs$table,lrt.M_I_P$table,lrt.PvsS$table,CLass_hela)
A549_table<-cbind(my_data,lrt.IpvsIs$table,lrt.MpvsMs$table,lrt.M_I_P$table,lrt.PvsS$table,CLass_A549)

levels(A549_table$CLass_A549)
levels(Hela_table$CLass_hela)
merged_table2<-cbind(A549_table,Hela_table)


# write.table(Hela_table,file="Hela_table2.txt",sep="\t")

#Venn plot for cell type specificity
Is.CAR_A549<-merged_table2$CLass_A549!="Non"
Is.CAR_Hela<-merged_table2$CLass_hela!="Non"

vennDiagram(cbind(Is.IMCAR_A549,Is.IMCAR_Hela))




merged_table3<-cbind(merged_table2,CellType)

write.table(merged_table3,file="merged_table_defined_CAR.txt",sep="\t")




#ggplot2 ploting

library(ggplot2)

anno$logFC_IpvsIs
g2 <- ggplot(anno,mapping = aes(logFC_IpvsIs,logFC_MpvsMs),colour =CLass_A549)

g2 <- g2 + geom_point(aes(colour =CLass_A549),show.legend = T)
g2 <- g2 + theme(legend.direction = 'horizontal', 
                 legend.position = 'top')
#g2<-g2+ geom_vline(xintercept = 1,col="blue")+ geom_hline(yintercept = 1,col="blue")
print(g2)






x<-'perl anno_from_HGNC.pl'

setwd("~/redefinedIM")

anno=read.delim("/Users/yan/Documents/scientific/CAR/analysis/Redefined_CAR/anno.txt",header=T,sep="\t")
row.names(anno)<-anno$X


# finding overlap of i- and m-CAR, identifying IO,IM and IM CARs 
levels(anno$CLass_hela)

Is.Hela_CAR<-anno$CLass_hela!="Non"
Is.A549_CAR<-anno$CLass_A549!="Non"

summary(Is.Hela_CAR)
Is.Hela_iCAR<-Is.Hela_CAR & anno$logFC_IpvsIs_Hela>1
summary(Is.Hela_iCAR)
Is.Hela_mCAR<-Is.Hela_CAR & anno$logFC_MpvsMs_Hela>1
summary(Is.Hela_mCAR)

vennDiagram(cbind(Is.A549_iCAR,Is.A549_mCAR))

# define and anno IO,IM,and MO-CARs
IM_Hela<-rep ("Non",times=nrow(anno))
length(IM_Hela)

for(i in 1:nrow(anno)){
  if(Is.Hela_iCAR[i] & Is.Hela_mCAR[i]){
    IM_Hela[i]="IM_CAR"
    next
  }
  
  if(Is.Hela_iCAR[i] & !Is.Hela_mCAR[i]){
    IM_Hela[i]="IO_CAR"
    next
  }
  if(!Is.Hela_iCAR[i] & Is.Hela_mCAR[i]){
    IM_Hela[i]="MO_CAR"
  }
  else{
    IM_Hela[i]="nonCAR"
  }
}

anno2<-cbind(anno,IM_A549, IM_Hela)

levels(anno2$IM_A549)
# venn plot of IM_CARs in 2 cells
Is.Hela_IOCAR<-anno2$IM_Hela=="IO_CAR"
Is.Hela_IMCAR<-anno2$IM_Hela=="IM_CAR"
Is.Hela_MOCAR<-anno2$IM_Hela=="MO_CAR"

vennDiagram(cbind(Is.A549_MOCAR,Is.Hela_MOCAR))

anno2$logFC_MpvsMs_Hela

######### plot dynamic CARs on 
library(ggplot2)
g2 <- ggplot(anno2,mapping = aes(logFC_IpvsIs_Hela,logFC_MpvsMs_Hela),color = IM_Hela)
g2 <- g2 + geom_point(aes(color = IM_Hela))
g2 <- g2 + theme(legend.direction = 'horizontal', 
                 legend.position = 'top' )
g2<-g2
g2<-g2+ggtitle("chromatin/soluble Foldchange of Hela-S3 dynamic CARs")
print(g2)





### annotate with HGNC and genecode V19 information,to get back chromatin coordinate and gene type information

x<-"perl anno_from_HGNC.pl"
system(x)





#########ploting RNA types
# plot RNA type of different class of CARs
g2 <- ggplot(anno[anno$gene_type!="Nd",],mapping = aes(x=factor(CLass_A549), y=1,fill = gene_type))
g2 <- g2 + geom_bar(stat = 'identity',position="fill")
g2 <- g2 + theme(legend.direction = 'horizontal', 
                 legend.position = 'top' )
g2<-g2
g2<-g2+ggtitle("RNA types of A549 dynamic CARs")
print(g2)



#########ploting RNA types
# plot RNA type of different class of CARs in two cell lines
library(ggplot2)
g2 <- ggplot(anno2[anno2$gene_type!="Nd",],mapping = aes(x=factor(IM_A549), y=1,fill = gene_type))
g2 <- g2 + geom_bar(stat = 'identity',position="fill")
g2 <- g2 + theme(legend.direction = 'horizontal', 
                 legend.position = 'top' )
g2<-g2+facet_grid(.~IM_Hela)
g2<-g2+ggtitle("RNA types of A549 dynamic CARs")
print(g2)

# based on count
g2 <- ggplot(anno2[anno2$gene_type!="Nd",],mapping = aes(factor(IM_A549),fill = gene_type))
g2 <- g2 + geom_bar(position = "dodge")
g2 <- g2 + theme(legend.direction = 'horizontal', 
                 legend.position = 'top' )
g2<-g2+scale_y_sqrt()
g2<-g2+ggtitle("RNA types of A549 dynamic CARs")
print(g2)



g2 <- ggplot(anno2[(anno2$gene_type!="Nd") &(anno3$Class_cell!="nonCAR | nonCAR"),],mapping = aes(factor(IM_A549),fill = gene_type))
g2 <- g2 + geom_bar(position = "stack")
g2 <- g2 + theme(legend.direction = 'horizontal', 
                 legend.position = 'top' )
g2<-g2+facet_grid(.~IM_Hela)
g2<-g2+ggtitle("RNA types of HeLa-S3 dynamic CARs")
print(g2)

levels(anno2$IM_A549)
#compare CARs from two cell lines
g2 <- ggplot(anno2[(anno2$gene_type!="Nd") &(anno3$Class_cell!="nonCAR | nonCAR"),],mapping = aes(x=0,fill = gene_type))
g2 <- g2 + geom_bar(position = "stack")
g2 <- g2 + theme(legend.direction = 'horizontal', 
                 legend.position = 'top' )
g2<-g2+facet_grid(ordered(IM_Hela,levels = c("nonCAR","MO_CAR","IO_CAR","IM_CAR"))~IM_A549)
g2<-g2+ggtitle("Barplot of dynamic CARs in both cells")
print(g2)


g2 <- ggplot(anno2[(anno2$gene_type!="Nd") &(anno3$Class_cell!="nonCAR | nonCAR"),],mapping = aes(IM_A549,IM_Hela))
g2 <- g2 + geom_count(aes(size=..n.., color=..n..), max_size = 400)
g2 <- g2 + theme(legend.direction = 'horizontal', 
                 legend.position = 'top' )
g2<-g2+scale_color_gradient(low = "#56B1F7",high ="#132B43",guide = "colourbar")
g2<-g2+ggtitle("Barplot of dynamic CARs in both cells")
print(g2)

intersect(anno2$IM_A549,anno2$IM_Hela)


# compare expression level of different class of CARs
library(ggplot2)
g2 <- ggplot(anno2,mapping = aes(IM_A549,logCPM_MpvsMs))
g2 <- g2 + geom_boxplot(aes(color=IM_A549))
g2 <- g2 + theme(legend.direction = 'horizontal', 
                 legend.position = 'top' )
g2<-g2+ggtitle("over all abundance of A549 dynamic CARs")
print(g2)


## performing wilcox.test for expression
cpm_IMCAR<-as.numeric(anno2[Is.Hela_IMCAR,"logCPM_MpvsMs"])
x=list(as.numeric(anno2[Is.Hela_IMCAR,"logCPM_MpvsMs"]),
       as.numeric(anno2[Is.Hela_IOCAR,"logCPM_MpvsMs"]),
       as.numeric(anno2[Is.Hela_MOCAR,"logCPM_MpvsMs"]),
       as.numeric(anno2[!Is.Hela_CAR,"logCPM_MpvsMs"])
       )
x2=x=list(as.numeric(anno2[Is.A549_IOCAR,"logCPM_MpvsMs"]),

          as.numeric(anno2[!Is.A549_CAR,"logCPM_MpvsMs"])

)

wilcox.test(as.numeric(anno2[Is.A549_IOCAR,"logCPM_MpvsMs"]),as.numeric(anno2[!Is.A549_CAR,"logCPM_MpvsMs"]),alternative = "less")



############## compare public peaks with IM_CAR TSS region
#load libraries
library(ChIPseeker)
library(clusterProfiler)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene



#prepare all of the TSS 
TSS_all<-cbind(anno2[53],as.numeric(anno2[[55]])-3000,as.numeric(anno2[[55]])+3000,anno2$X)
colnames(TSS_all) = c("chrom","start","end","gene")

TSS_Hela_nonCAR<-as(TSS_all[!Is.Hela_CAR,],"GRanges")
TSS_Hela_iCAR<-as(TSS_all[Is.Hela_iCAR,],"GRanges")
TSS_Hela_mCAR<-as(TSS_all[Is.Hela_mCAR,],"GRanges")

TSS_Hela_IMCAR<-as(TSS_all[Is.Hela_IMCAR,],"GRanges")
TSS_Hela_IOCAR<-as(TSS_all[Is.Hela_IOCAR,],"GRanges")
TSS_Hela_MOCAR<-as(TSS_all[Is.Hela_MOCAR,],"GRanges")



TSS_commonIMCAR<-as(TSS_all[Is.commonIMCAR,],"GRanges")
TSS_commonNonIMCAR<-as(TSS_all[!Is.A549_CAR &!Is.Hela_CAR,],"GRanges")
TSS_AspecIMCAR<-as(TSS_all[Is.AspecIMCAR,],"GRanges")
TSS_HspecIMCAR<-as(TSS_all[Is.HspecIMCAR,],"GRanges")
TSS_AspecIMCAR_l<-as(TSS_all[Is.AspecIMCAR_l,],"GRanges")
TSS_HspecIMCAR_l<-as(TSS_all[Is.HspecIMCAR_l,],"GRanges")


#load ChIP-seq peaks

A549H3k27ac  <- readPeakFile("/Volumes/Elements/Useful Public Data/ENCODE data/histone marks/A549/wgEncodeBroadHistoneA549H3k27acEtoh02Pk.broadPeak.gz")
A549H3k04me1  <- readPeakFile("/Volumes/Elements/Useful Public Data/ENCODE data/histone marks/A549/wgEncodeBroadHistoneA549H3k04me1Etoh02Pk.broadPeak.gz")
A549H3k04me3  <-readPeakFile( "/Volumes/Elements/Useful Public Data/ENCODE data/histone marks/A549/wgEncodeBroadHistoneA549H3k04me3Etoh02Pk.broadPeak.gz")
A549H3k27me3  <-readPeakFile( "/Volumes/Elements/Useful Public Data/ENCODE data/histone marks/A549/wgEncodeBroadHistoneA549H3k27me3Etoh02Pk.broadPeak.gz")
A549Pol2<-readPeakFile( "/Volumes/Elements/Useful Public Data/ENCODE data/histone marks/A549/wgEncodeAwgTfbsHaibA549Pol2Pcr2xEtoh02UniPk.narrowPeak.gz")





tagMatrix_A549H3k27me3_nonCAR <- getTagMatrix(A549H3k27me3, windows=TSS_A549_nonCAR)
tagMatrix_A549H3k27me3_IMCAR <- getTagMatrix(A549H3k27me3, windows=TSS_A549_IMCAR)
tagMatrix_A549H3k27me3_IOCAR <- getTagMatrix(A549H3k27me3, windows=TSS_A549_IOCAR)
tagMatrix_A549H3k27me3_MOCAR <- getTagMatrix(A549H3k27me3, windows=TSS_A549_MOCAR)

tagMatrixList=list(A549H3k27me3_nonCAR=tagMatrix_A549H3k27me3_nonCAR,
                   A549H3k27me3_IMCAR=tagMatrix_A549H3k27me3_IMCAR,
                   A549H3k27me3_IOCAR=tagMatrix_A549H3k27me3_IOCAR,
                   A549H3k27me3_MOCAR=tagMatrix_A549H3k27me3_MOCAR)



########analysis overlap with HEK293 cheRNA

HEK293_cheRNA  <- readPeakFile("/Users/yan/Documents/scientific/CAR/downloaded_data/HEKcheRNA.bed")




tagMatrix_HEK293_cheRNA_nonCAR <- getTagMatrix(HEK293_cheRNA, windows=TSS_A549_nonCAR)
tagMatrix_HEK293_cheRNA_IMCAR <- getTagMatrix(HEK293_cheRNA, windows=TSS_A549_IMCAR)
tagMatrix_HEK293_cheRNA_IOCAR <- getTagMatrix(HEK293_cheRNA, windows=TSS_A549_IOCAR)
tagMatrix_HEK293_cheRNA_MOCAR <- getTagMatrix(HEK293_cheRNA, windows=TSS_A549_MOCAR)


tagMatrixList=list(A549IM_FIRE_nonCAR=tagMatrix_HEK293_cheRNA_nonCAR,
                   A549IM_FIRE_IMCAR=tagMatrix_HEK293_cheRNA_IMCAR,
                   A549IM_FIRE_IOCAR=tagMatrix_HEK293_cheRNA_IOCAR)


plotAvgProf(tagMatrixList, xlim=c(-3000, 3000),facet = "column")
tagHeatmap(tagMatrix_A549H3k04me1_IOCAR,xlim=c(-3000, 3000),xlab = "TSS")

  peaks<-list(HEK293_cheRNA=HEK293_cheRNA,Hela_nonCAR=TSS_Hela_nonCAR,A549_Non=TSS_A549_nonCAR)
vennplot(peaks)


####


############Functinal annotation
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

#dynamic
peakAnnoList <- lapply(list(IMCAR=TSS_Hela_IMCAR,IOCAR=TSS_Hela_IOCAR,MOCAR=TSS_Hela_MOCAR,nonCAR=TSS_Hela_nonCAR
), seq2gene, TxDb=txdb, 
tssRegion=c(-3000, 3000),  flankDistance = 5000)

compKEGG <- compareCluster(geneCluster   = peakAnnoList, 
                           fun           = "enrichKEGG",
                           pvalueCutoff  = 0.05, 
                           pAdjustMethod = "none")
plot(compKEGG, showCategory = 15, title = "Hela KEGG Pathway Enrichment Analysis")


compGO <- compareCluster(geneCluster   = peakAnnoList, 
                         fun           = "enrichGO",
                         pvalueCutoff  = 0.05, 
                         pAdjustMethod = "none",
                         OrgDb="org.Hs.eg.db")
plot(compGO, showCategory = 15, title = "Hela Gene Ontology Enrichment Analysis")


# cell type specific
peakAnnoList <- lapply(list(spec_A549_IMCAR=TSS_AspecIMCAR,spec_Hela_IMCAR=TSS_HspecIMCAR,
                            common_IMCAR=TSS_commonIMCAR), seq2gene, TxDb=txdb,
                       tssRegion=c(-3000, 3000),  flankDistance = 5000)

compKEGG <- compareCluster(geneCluster   = peakAnnoList, 
                           fun           = "enrichKEGG",
                           pvalueCutoff  = 0.05, 
                           pAdjustMethod = "none")
plot(compKEGG, showCategory = 15, title = "KEGG Pathway Enrichment Analysis")


compGO <- compareCluster(geneCluster   = peakAnnoList, 
                         fun           = "enrichGO",
                         pvalueCutoff  = 0.05, 
                         pAdjustMethod = "none",
                         OrgDb="org.Hs.eg.db")
plot(compGO, showCategory = 15, title = "Gene Ontology Enrichment Analysis")



#analysis superenhancers
#load super and non-super enhancer peaks
hela_nonsuper  <- readPeakFile("/Users/yan/Documents/scientific/CAR/downloaded_data/superenhancers/hela_nonsuper.bed.txt")
hela_super  <- readPeakFile("/Users/yan/Documents/scientific/CAR/downloaded_data/superenhancers/hela_super.bed.txt")
IMR90_nonsuper  <- readPeakFile("/Users/yan/Documents/scientific/CAR/downloaded_data/superenhancers/IMR90_nonsuper.txt")
IMR90_super<-readPeakFile("/Users/yan/Documents/scientific/CAR/downloaded_data/superenhancers/IMR90_super.txt")

tagMatrix_IMR90_super_nonCAR <- getTagMatrix(IMR90_super, windows=TSS_A549_nonCAR)
tagMatrix_IMR90_super_IMCAR <- getTagMatrix(IMR90_super, windows=TSS_A549_IMCAR)
tagMatrix_IMR90_super_IOCAR <- getTagMatrix(IMR90_super, windows=TSS_A549_IOCAR)
tagMatrix_IMR90_super_MOCAR <- getTagMatrix(IMR90_super, windows=TSS_A549_MOCAR)

tagMatrixList=list(IMR90_super_nonCAR=tagMatrix_IMR90_super_nonCAR,IMR90_super_IMCAR=tagMatrix_IMR90_super_IMCAR,IMR90_super_IOCAR=tagMatrix_IMR90_super_IOCAR,IMR90_super_MOCAR=tagMatrix_IMR90_super_MOCAR)
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000),facet = "column")
tagHeatmap(tagMatrix_A549H3k04me1_IOCAR,xlim=c(-3000, 3000),xlab = "TSS")

peaks<-list(hela_nonsuper=hela_nonsuper,hela_super=hela_super,Hela_IMCAR=TSS_Hela_IMCAR,Hela_IOCAR=TSS_Hela_IOCAR,Hela_MOCAR=TSS_Hela_MOCAR)
peaks<-list(IMR_nonsuper=IMR90_nonsuper,A549_IMCAR=TSS_A549_nonCAR)

vennplot(peaks)

mat<-matrix( c(6,19,1,1,0,1,119,612),nrow=2,ncol=4)
fisher.test(mat)


#analysis TFBS
#load super and non-super enhancer peaks
TfbsCluster  <- readPeakFile("~/UCSC_data/wgEncodeRegTfbsClusteredWithCellsV3.bed")

tagMatrix_TfbsCluster_nonCAR <- getTagMatrix(TfbsCluster, windows=TSS_A549_nonCAR)
tagMatrix_TfbsCluster_IMCAR <- getTagMatrix(TfbsCluster, windows=TSS_A549_IMCAR)
tagMatrix_TfbsCluster_IOCAR <- getTagMatrix(TfbsCluster, windows=TSS_A549_IOCAR)
tagMatrix_TfbsCluster_MOCAR <- getTagMatrix(TfbsCluster, windows=TSS_A549_MOCAR)

tagMatrixList=list(TfbsCluster_nonCAR=tagMatrix_TfbsCluster_nonCAR,TfbsCluster_IMCAR=tagMatrix_TfbsCluster_IMCAR,TfbsCluster_IOCAR=tagMatrix_TfbsCluster_IOCAR,TfbsCluster_MOCAR=tagMatrix_TfbsCluster_MOCAR)
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000),facet = "column")
tagHeatmap(tagMatrixList,xlim=c(-3000, 3000),xlab = "TSS")

peaks<-list(TfbsCluster=TfbsCluster,Hela_IMCAR=TSS_Hela_IMCAR,Hela_IOCAR=TSS_Hela_IOCAR,Hela_MOCAR=TSS_Hela_MOCAR)
peaks<-list(TfbsCluster=TfbsCluster,A549_nonCAR=TSS_A549_nonCAR)

vennplot(peaks)

mat<-matrix( c(6,19,1,1,0,1,119,612),nrow=2,ncol=4)
fisher.test(mat)

#export tss to bed files
write.table(TSS_all[Is.HspecIMCAR,],file="TSS_HspecIMCAR.bed",quote = F,sep = "\t")





# meta-exon coverage analysis 



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



BC4_A549IMCOV<-cov.matrix(BC4Bam,coordfile=A549_IM,no_windows=100, offset=10,num_cores=6,normalization="rpm")



BC4_A549IMCOV<-as.data.frame(BC4_A549IMCOV)

data<-BC4_A549IMCOV

for(i in 1:nrow(data)){
  meani<-mean(as.numeric(data[i,(1:ncol(data))]))
  if(i==1){
    mean<-meani
  }
  else{
    mean<-c(mean,meani)
  }
}


BC4_A549IMCOV_mean<-cbind(bin=(1:nrow(BC4_A549IMCOV)),mean,sample="A549_exoncunt_multi")

write.table(BC4_A549MOCOV_mean, file="BC4_coverage on A549 IM_CAR multiexon.txt",sep = "\t",quote = F)

BC4_A549IMCOV_mean<-read.delim(file="BC4_coverage on A549 IM_CAR multiexon.txt",header = T,sep = "\t")



data<-rbind(BC4_A549IMCOV_mean,BC8_HeLaIMCOV_mean)

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





#### chip-seq peak intergration


files<-c("TSS_A549Ctcf.bed", "TSS_A549H2az.bed", "TSS_A549H3k04me1.bed", "TSS_A549H3k04me2.bed", "TSS_A549H3k04me3.bed", "TSS_A549H3k09ac.bed", "TSS_A549H3k09me3.bed", "TSS_A549H3k27me3.bed", "TSS_A549H3k36me3.bed", "TSS_A549H3k79me2.bed", "TSS_Helas3Ctcf.bed", "TSS_Helas3Ctcf.bed", "TSS_Ezh2.bed", "TSS_Helas3H2az.bed", "TSS_Helas3H3k04me1.bed", "TSS_Helas3H3k04me2.bed", "TSS_Helas3H3k04me3.bed", "TSS_Helas3H3k9ac.bed", "TSS_Helas3H3k27me3.bed", "TSS_Helas3H3k36me3.bed", "TSS_Helas3H3k79me2.bed", "TSS_Helas3H4k20me1.bed", "TSS_Helas3Pol2.bed")


TSS_A549Ctcf<-read.delim("TSS_A549Ctcf.bed",header = F,sep = "\t")
colnames(TSS_A549Ctcf)


import<-function (file){
  temp<-read.delim(file,header = F,sep = "\t")
  colnames(temp)<-colnames(A549_H3k27Ac_TSS)
return (temp)
}

TSS_A549Ctcf<-import("TSS_A549H2az.bed")
TSS_Helas3Pol2<-import("TSS_Helas3Pol2.bed")

TSS_A549A549H3k27me3_poli<-intersected_polish_noscore(data =TSS_A549A549H3k27me3,gene="gene",col_obs=15 )


# intergrate data and count mean score of overlaped intervals
intersected_polish <- function(data = NULL, gene, col_obs,score) {
  scalar<-length(levels(factor(data[[gene]])))
  genes<-levels(factor(data[[gene]]))
  for (i in 1:scalar) {
    temp <- data[(data[[gene]] == genes[i]),]
    col_obs<-as.integer(col_obs)
    if (as.character(temp[1, col_obs]) == ".") {
      N_record <- 0
      score_mean <- 0
    }
    else{
      N_record <- nrow(temp)
      score_mean <- sum(as.numeric(as.character(temp[[score]]))) / N_record
    }
    temp <- cbind(temp[1, ], N_record, score_mean)
    if (i == 1) {
      polished <- temp[1,]
    }
    else{
      polished <- rbind(polished, temp)
    }
  }
  return (polished)
}

# intergrate data and count mean score of overlaped intervals
intersected_polish_noscore <- function(data = NULL, gene, col_obs) {
  scalar<-length(levels(factor(data[[gene]])))
  genes<-levels(factor(data[[gene]]))
  for (i in 1:scalar) {
    temp <- data[(data[[gene]] == genes[i]),]
    col_obs<-as.integer(col_obs)
    if (as.character(temp[1, col_obs]) == ".") {
      N_record <- 0
      # score_mean <- 0
    }
    else{
      N_record <- nrow(temp)
      # score_mean <- sum(as.numeric(temp[[score]])) / N_record
    }
    temp <- cbind(temp[1, ], N_record)
    if (i == 1) {
      polished <- temp[1,]
    }
    else{
      polished <- rbind(polished, temp)
    }
  }
  return (polished)
}


poli<-function (data=NULL){
  temp<-intersected_polish_noscore(data =data,gene="gene",col_obs=15 )
  return (temp)
}


integrat<-function(files,data=NULL){
  files<-as.character(files)
  for (i in 1:length(files)){
    temp<-read.delim(files[i],header = F,sep = "\t")
    colnames(temp)<-colnames(A549_H3k27Ac_TSS)
    
    temp2<-poli(data =temp )
    
    data[[files[i]]][temp2$N_record==0]="negative"
    data[[files[i]]][temp2$N_record>0]="positive"
  }
  return (data)
}



#### coilin iCLIP data integration


superCAR<-read.delim("super_commomCAR_info.bed",header = F)
colin_ChIP<-read.delim("colin_Chip.txt",header=T)
colin_iCLIP<-read.delim("colin_iCLIP.txt",header = T)

colin_ChIP$chr_start<-colin_ChIP$start
colin_ChIP$chr_start[colin_ChIP<-colin_ChIP$start

colin_iCLIP$start<-colin_iCLIP$gene.start
colin_iCLIP$start[colin_iCLIP$strand=="-"]<-colin_iCLIP$gene.end[colin_iCLIP$strand=="-"]

colin_iCLIP$end<-colin_iCLIP$gene.end
colin_iCLIP$end[colin_iCLIP$strand=="-"]<-colin_iCLIP$gene.start[colin_iCLIP$strand=="-"]



                                      
 write.table(cbind(as.character(colin_ChIP$chr),colin_ChIP$start,colin_ChIP$end),quote = F,row.names = F,col.names = F,sep = "\t", file = "colin_Chip.bed")                    

 write.table(cbind(as.character(colin_iCLIP$X.chromosome),colin_iCLIP$start,colin_iCLIP$end),quote = F,row.names = F,col.names = F,sep = "\t", file = "colin_iCLIP.bed")                    
 
 
 
 
x<-"cd ~/Documents/scientific/CAR/analysis/Redefined_CAR/redefinedIM/chr_fixed"
system(x)

x<-"sort -k1,1 -k2,2n colin_Chip.bed >colin_Chip_sorted.bed
bedtools merge -i colin_Chip_sorted.bed -c 1 -o count >colin_Chip_merged.bed"
system(x)
x<-"bedtools intersect -wa -loj -a super_commomCAR_info.bed -b colin_Chip_merged.bed >super_commomCAR_CHIP.bed"
system(x)
super_ChIP<-read.delim("super_commomCAR_CHIP.bed",header = F,sep = "\t")



x<-"sort -k1,1 -k2,2n colin_iCLIP.bed >colin_iCLIP.bed_sorted.bed
bedtools merge -i colin_iCLIP.bed_sorted.bed -c 1 -o count >colin_iCLIP_merged.bed"
system(x)
x<-"bedtools intersect -wa -loj -a super_commomCAR_info.bed -b colin_iCLIP_merged.bed >super_commomCAR_iCLIP.bed"
system(x)
super_iCLIP<-read.delim("super_commomCAR_iCLIP.bed",header = F,sep = "\t")


super_iCLIP_polish<-intersected_polish_noscore(data = super_iCLIP,gene = "V4",col_obs = 9)
super_ChIP_polish<-intersected_polish_noscore(data = super_ChIP,gene = "V4",col_obs = 9)


super_iCLIP_polish$V4==super_ChIP_polish$V4


superCAR_colin<-super_iCLIP_polish[1:8]
superCAR_colin$Colin_CLIP<-"positive"
superCAR_colin$Colin_CLIP[super_iCLIP_polish$N_record==0]<-"negative"
superCAR_colin$Colin_ChIP<-"positive"
superCAR_colin$Colin_ChIP[super_ChIP_polish$N_record==0]<-"negative"


library(ggplot2)
g2<-ggplot(data = superCAR_colin, aes(Colin_ChIP,Colin_CLIP))
g2 <- g2 + geom_count(aes(size=(..n..)^0.5, color=..prop..))+scale_size_continuous()
g2


superCAR_colin$V7

g2<-ggplot(data = superCAR_colin, aes(V7,fill=Colin_CLIP))
g2 <- g2 + geom_bar(position = "fill")
g2

gene_info<-read.delim("gene_info_20170509.txt",header = T,sep = "\t")
levels(gene_info$IM_CAR)
write.table(cbind(gene_info[1:10],gene_info$IM_CAR), quote = F,
            col.names = F,row.names = F, sep = "\t", file="gene_info.20171219.bed")


x<-"sort -k1,1 -k2,2n gene_info.20171219.bed >gene_info.20171219_sort.bed
bedtools merge -i gene_info.20171219_sort.bed -c 1 -o count >gene_info.20171219_merged.bed"
system(x)
x<-"bedtools intersect -wa -loj -a gene_info.20171219.bed -b colin_iCLIP_merged.bed >gene_iCLIP.bed"
system(x)
gene_iCLIP<-read.delim("gene_iCLIP.bed",header = F,sep = "\t")

gene_iCLIP_POLISH<-intersected_polish_noscore(data = gene_iCLIP,gene = "V4",col_obs = 12)

gene_iCLIP_POLISH$iCLIP<-"positive"
gene_iCLIP_POLISH$iCLIP[gene_iCLIP_POLISH$N_record==0]<-"negative"
g2<-ggplot(gene_iCLIP_POLISH,mapping = aes(V9,fill=iCLIP))
g2+geom_bar(position = "fill")



g2<-ggplot(gene_iCLIP_POLISH[gene_iCLIP_POLISH$V7=="RNA, small nucleolar",],mapping = aes(V11,fill=iCLIP))
g2+geom_bar(position = "fill")+xlim(c("A_com4","B_comsomatic","C_comNF","D_comDi","E_Specific","Aspecific","F_non"))


g2<-ggplot(gene_iCLIP_POLISH,mapping = aes(V11,fill=iCLIP))
g2+geom_bar(position = "fill")+xlim(c("A_com4","B_comsomatic","C_comNF","D_comDi","E_Specific","Aspecific","F_non"))



sno<-gene_iCLIP_POLISH[gene_iCLIP_POLISH$V7=="RNA, small nucleolar",]

sno$CAR<-"negative"
sno$CAR[(grepl("A_",sno$V11)|grepl("B_",sno$V11)|grepl("C_",sno$V11))]<-"comCAR"
sno$CAR[grepl("E_",sno$V11)|grepl("F_",sno$V11)]<-"speficic"       
sno$CAR[grepl("D_",sno$V11)]<-"MARGI"

g2<-ggplot(sno,mapping = aes(CAR,fill=iCLIP))
g2+geom_bar(position = "fill")+xlim(c("comCAR","speficic","MARGI","negative"))+
  scale_fill_manual(values = c("#3399cc", "#ff6666"))

superCAR_colin[superCAR_colin$V4=="RMRP","V7"]<-"RNA, small nucleolar"
superCAR_colin$V9<-superCAR_colin$V7
superCAR_colin$V9[grepl("SNHG",superCAR_colin$V4)]<-"snoRNA_host"
superCAR_colin$V8[l]<-"RNU_pseudogene"


g2<-ggplot(superCAR_colin,mapping = aes(V7,fill=Colin_CLIP))
g2+geom_bar()+
  scale_fill_manual(values = c("#3399cc", "#ff6666"))

sno_info<-read.delim("snoRNA_info.txt",header = T)
back<-superCAR_colin
superCAR_colin<-back

write.table(superCAR_colin, file="superCAR_inCajal.txt",quote = F,sep = "\t",row.names = F)
superCAR_inCajal<-read.delim("superCAR_inCajal.txt",header = T)



back<-gene_iCLIP_POLISH


gene_iCLIP_POLISH$CAR<-"negative"
gene_iCLIP_POLISH$CAR[(grepl("A_",gene_iCLIP_POLISH$V11)|grepl("B_",gene_iCLIP_POLISH$V11)|grepl("C_",gene_iCLIP_POLISH$V11))]<-"comCAR"
gene_iCLIP_POLISH$CAR[grepl("E_",gene_iCLIP_POLISH$V11)|grepl("Aspecific",gene_iCLIP_POLISH$V11)]<-"speficic"       
gene_iCLIP_POLISH$CAR[grepl("D_",gene_iCLIP_POLISH$V11)]<-"MARGI"


g2<-ggplot(gene_iCLIP_POLISH,mapping = aes(V9,fill=iCLIP))
g2+geom_bar(position = "fill")+
  scale_fill_manual(values = c("#3399cc", "#ff6666"))+xlim(c("comCAR","MARGI","speficic","negative"))
  facet_grid(V9~.)





##### compare phastconsv

x<-"perl symbol2exon.pl"

system(x)

x<-"bedtools intersect -loj -wa -wb -a gene_exon.txt -b /Volumes/ZANY64G/UCSC_data/phaseConsElements.bed >exon_phaseCons.bed" 
system(x)

exon_phastconsEl<-read.delim("exon_phaseCons.bed",header=F,sep="\t")
colnames(gene_phastconsEl)<-c("chrom","start","end","gene","strand","gene_type","gene_fam","IM_A549","IM_HeLa","Class_cell", "EL_chrom","El_start","El_end","El_lod","El_score")


exon_phastconsEl_polish<-intersected_polish(data = exon_phastconsEl[exon_phastconsEl$V9!=".",],gene = "V7",col_obs =9,score = 13 )



gene_info<-read.delim("gene_info_20170509.txt",header = T,sep = "\t")

A549_IM<-as.character(gene_info$gene[gene_info$IM_A549=="IM_CAR"])
A549_IO<-as.character(gene_info$gene[gene_info$IM_A549=="IO_CAR"])
A549_MO<-as.character(gene_info$gene[gene_info$IM_A549=="MO_CAR"])


ACD<-read.delim("ACD_comIM_CAR.bed",header = F,sep = "\t")


A<-is.element(exon_phastconsEl_polish$V7,as.character(ACD$V4[ACD$V9=="A_com4"]))
C<-is.element(exon_phastconsEl_polish$V7,as.character(ACD$V4[ACD$V9=="C_comNF"]))
D<-is.element(exon_phastconsEl_polish$V7,as.character(ACD$V4[ACD$V9=="D_comDi"]))


exon_phastconsEl_polish$CAR<-"non"
exon_phastconsEl_polish$CAR[A]<-"A"
exon_phastconsEl_polish$CAR[C]<-"C"
exon_phastconsEl_polish$CAR[D]<-"D"

colnames(phastconsEl_polished)<-c("chrom","start","end","gene","strand","gene_type","gene_fam","IM_A549","IM_HeLa","Class_cell", "N_record","El_score")

library(ggplot2)

exon_phastconsEl_polish$score_mean
#plot peaks number over promoter CARs

g2 <- ggplot(exon_phastconsEl_polish,mapping = aes(factor(CAR),score_mean, fill = CAR))
g2 <- g2 + geom_boxplot()+scale_y_log10()
g2 <- g2 + theme(legend.direction = 'horizontal', 
                 legend.position = 'top' )
g2<-g2+facet_grid(.~IM_A549)
g2<-g2+ggtitle("No. phastcons element peaks per promoter of in dynamic CARs")
print(g2)

rt<-Ranksumtest(exon_phastconsEl_polish,group = "CAR",observ = "score_mean")

#plot peaks number over promoter CARs

g2 <- ggplot(phastconsEl_polished,mapping = aes(factor(IM_A549),1000*N_record/(end-start), fill = IM_A549))
g2 <- g2 + geom_boxplot()+ylim(0,100)
g2 <- g2 + theme(legend.direction = 'horizontal', 
                 legend.position = 'top' )
g2<-g2+facet_grid(.~IM_HeLa)
g2<-g2+ggtitle("No. phastcons element peaks per promoter of in dynamic CARs in both cells")
print(g2)


#plot total peaks number over promoter CARs
g2 <- ggplot(polished,mapping = aes(factor(IM_HeLa),score_sum/N_record, fill = IM_HeLa))
g2 <- g2 + geom_boxplot()
g2 <- g2 + theme(legend.direction = 'horizontal', 
                 legend.position = 'top' )
g2<-g2+facet_grid(.~IM_A549)
g2<-g2+ggtitle("mean H3K27ac peaks score per promoter of in dynamic CARs")
print(g2)









