# from polished data 

enrichment.Fishertest <- function (data = NULL,
                                   group = NULL,
                                   observ = NULL,
                                   cut=NULL, alternative="two.sided") {
  # first make enrichment matrix
  groups<-levels(factor(data[[group]]))
  for (i in 1:length(groups)) {
    groupi <- groups[i]
    temp <- data[data[[group]] == groupi, ]
    Negative <- nrow(temp[temp[[observ]] <= cut, ])
    positive <- nrow(temp[temp[[observ]] > cut, ])
    r <- cbind(positive,Negative)
    if (i == 1) {
      matrix_enrich <- r
    }
    else
      matrix_enrich <- rbind(matrix_enrich, r)
  }
row.names(matrix_enrich) <-groups

# then perform fisher test
P <-
  matrix(
    rep(0, nrow(matrix_enrich) * nrow(matrix_enrich)),
    nrow  = nrow(matrix_enrich),
    dimnames = list(groups , groups )
  )
OR <-
  matrix(
    rep(0, nrow(matrix_enrich) * nrow(matrix_enrich)),
    nrow  = nrow(matrix_enrich),
    dimnames = list(groups ,groups )
  )


for (i in 1:nrow(matrix_enrich) ) {
  for (j in i:nrow(matrix_enrich)) {
    test <- fisher.test(matrix_enrich[c(i, j), 1:2],alternative = alternative)
    P[i, j] <- test$p.value
    OR[i, j] <- test$estimate
  }
}

return (list(
  count = matrix_enrich,
  method=list(test$method,test$alternative),
  Pvalue = P,
  OddRatio = OR
))

}

enrich_test <- function(data=NULL,
                        test_group=NULL,
                        target_group=NULL,target_level=NULL) {
  test_groups <- levels(factor(data[[test_group]]))
  target_groups <- levels(factor(data[[target_group]]))
  for (i in 1:length(test_groups)) {
    pos_pos <-
      nrow(data[(data[[test_group]] == test_groups[i]) &
                  (data[[target_group]] == target_level), ])
    neg_pos <-
      nrow(data[(data[[test_group]] != test_groups[i]) &
                  (data[[target_group]] == target_level), ])
    pos_neg <-
      nrow(data[(data[[test_group]] == test_groups[i]) &
                  (data[[target_group]] != target_level), ])
    neg_neg <-
      nrow(data[(data[[test_group]] != test_groups[i]) &
                  (data[[target_group]] != target_level), ])
    
    t <-
      fisher.test(matrix(c(pos_pos, neg_pos, pos_neg, neg_neg), ncol = 2))
    temp <-
      cbind(
        test_group = test_groups[i],
        method = t$method,
        OR = t$estimate,
        p_value = t$p.value,
        alternative = t$alternative
      )
    
    if (i == 1) {
      out <- temp
    }
    else{
      out <- rbind(out, temp)
    }
    
  }  
  rownames(out)<-1:length(test_groups)
  return(out)
}

Ranksumtest <- function (data = NULL,
                         group = NULL,
                         observ = NULL) {
  # first make enrichment matrix
  groups<-levels(factor(data[[group]]))
  for (i in 1:length(groups)) {
    groupi <- groups[i]
    mean <- mean(data[data[[group]] == groupi, observ])
    median <- median(data[data[[group]] == groupi, observ])
    
    r <- cbind(mean, median)
    if (i == 1) {
      matrix_enrich <- r
    }
    else
      matrix_enrich <- rbind(matrix_enrich, r)
  }
  row.names(matrix_enrich) <- groups
  
  
  
  
  # then perform fisher test
  P <-
    matrix(
      rep(NA, nrow(matrix_enrich) * nrow(matrix_enrich)),
      nrow  = nrow(matrix_enrich),
      dimnames = list(groups, groups)
    )
  
  P2 <-
    matrix(
      rep(NA, nrow(matrix_enrich) * nrow(matrix_enrich)),
      nrow  = nrow(matrix_enrich),
      dimnames = list(groups, groups)
    )
  
  for (i in 1:nrow(matrix_enrich)) {
    for (j in 1:nrow(matrix_enrich)) {
      tempi <- data[data[[group]] == groups[i], observ]
      tempj <- data[data[[group]] == groups[j], observ]
      test <- wilcox.test(tempi, tempj)
      P[i, j] <- test$p.value
      test2 <- wilcox.test(tempi, tempj, alternative = "greater")
      P2[i, j] <- test2$p.value
    }
  }
  
  return (list(
    count = matrix_enrich,
    Pvalue = P,
    Pgreater = P2
    
  ))
  
}


enrich_test <- function(data=NULL,
                        test_group=NULL,
                        target_group=NULL,target_level=NULL) {
  test_groups <- levels(factor(data[[test_group]]))
  target_groups <- levels(factor(data[[target_group]]))
  for (i in 1:length(test_groups)) {
    pos_pos <-
      nrow(data[(data[[test_group]] == test_groups[i]) &
                  (data[[target_group]] == target_level), ])
    neg_pos <-
      nrow(data[(data[[test_group]] != test_groups[i]) &
                  (data[[target_group]] == target_level), ])
    pos_neg <-
      nrow(data[(data[[test_group]] == test_groups[i]) &
                  (data[[target_group]] != target_level), ])
    neg_neg <-
      nrow(data[(data[[test_group]] != test_groups[i]) &
                  (data[[target_group]] != target_level), ])
    
    t <-
      fisher.test(matrix(c(pos_pos, neg_pos, pos_neg, neg_neg), ncol = 2))
    temp <-
      cbind(
        test_group = test_groups[i],
        method = t$method,
        OR = t$estimate,
        p_value = t$p.value,
        alternative = t$alternative
      )
    
    if (i == 1) {
      out <- temp
    }
    else{
      out <- rbind(out, temp)
    }
    
  }  
  rownames(out)<-1:length(test_groups)
  return(out)
}


enrich_test1<-enrich_test(data = diMARGI_H9_degrees,test_group = "gene_type",target_group = "IM_A549",target_level = "IM_CAR")

enrich_test2<-enrich_test(data = diMARGI_H9_degrees,test_group = "gene_type",target_group = "IM_HeLa",target_level = "IM_CAR")

enrich_testb<-cbind(enrich_test1,enrich_test2)

enrichment.Fishertest1 <-enrichment.Fishertest(data = data,group = "IM_HeLa",N_record ="N_record",cut = 0,alternative = "two.sided" )

Ranksumtest1<-Ranksumtest(data = data,group = "IM_A549",observ = "m6A_density")

write.table(enrichment.Fishertest1 ,file="whether m6A CITS sites per RNA Fishertest, HeLa.txt",sep = "\t",quote = F)

write.table(Ranksumtest1,file="m6A CITS sites per 1kb Ranksumtest, A549.txt",sep = "\t",quote = F)

m6A_density<-m6A_CITS_polished$N_record/(as.integer(m6A_CITS_polished$end)-as.integer(m6A_CITS_polished$start))
m6A_CITS_polished2<-cbind(m6A_CITS_polished,m6A_density)

library(ggplot2)
data = m6A_CITS_polished2[m6A_CITS_polished2$N_record>0,]
data$score_mean
g2 <- ggplot(data = data,mapping = aes(factor(IM_HeLa),m6A_density, fill =factor(IM_HeLa)))
g2 <- g2 + geom_boxplot()+scale_y_log10()
g2 <- g2 + theme(legend.direction = 'horizontal', 
                 legend.position = 'top' )
g2<-g2+facet_grid(.~IM_A549)
g2<-g2+ggtitle("m6A CITS sites  per kilobase  in both(Hela,A549) cells")
print(g2)


g2 <- ggplot(data[data$N_record>0,],mapping = aes(IM_A549, fill =factor(IM_A549)))
g2 <- g2 + geom_bar()
g2 <- g2 + theme(legend.direction = 'horizontal', 
                 legend.position = 'top' )
#g2<-g2+facet_grid(.~IM_A549)
g2<-g2+ggtitle("overall abounance of  RNA in HeLa-S3 ")
print(g2)

data$

g2 <- ggplot(geneType_testHeLa,mapping = aes(X.lnP, OR,size=X.lnP))
g2 <- g2 + geom_point(aes(color=test_group))
g2 <- g2 + theme(legend.direction = 'horizontal', 
                 legend.position = 'top' )
#g2<-g2+facet_grid(.~IM_A549)
g2<-g2+scale_y_log10()
g2<-g2+ggtitle("sigificance of gene_type in  IM_CARs in HeLa-S3 ")
print(g2)

geneType_testA549$test_group

g2 <- ggplot(data = geneType_testA549,mapping = aes(factor(X.lnP),OR, fill =factor(test_group)))
g2 <- g2 + geom_dotplot()+scale_y_log10()
g2 <- g2 + theme(legend.direction = 'horizontal', 
                 legend.position = 'top' )
#g2<-g2+facet_grid(.~IM_A549)
g2<-g2+ggtitle("m6A CITS sites  per kilobase  in both(Hela,A549) cells")
print(g2)
