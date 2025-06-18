rm(list=ls()) 

#没有请先安装这个包
# BiocManager::install("sva")

#导入数据
count1 <- read.table("data/mRNA.counts.matrix",
                     header = T, row.names = 1)

count2 <- read.table("data/mRNA.TMM.EXPR.matrix",
                     header = T, row.names = 1)



#过滤掉低表达的基因
Expr0 <- count1[rowSums(count1)>1,]

Expr <- count2[rowSums(count2)>0.055,]


##载入样品信息，含批次
data <- read.table("data/data.txt",header = T)
data[,2] <- as.factor(data$type)
data[,3] <- as.factor(data$batch)
row.names(data) <- data$sample


############不做校正的PCA分析
library(ggplot2)

#install.packages("FactoMineR")
#install.packages("factoextra")
library("FactoMineR")
library("factoextra")

pre.pca <- PCA(t(Expr),graph = FALSE)

fviz_pca_ind(pre.pca,
             geom = "point",
             col.ind = data$batch,
             addEllipses = TRUE,
             legend.title ="Group")


############DESeq2包加入批次信息进行差异分析
library(DESeq2)

dds <- DESeqDataSetFromMatrix(countData = Expr0,
                              colData = data1,
#                              design = ~ type + batch
                              )

dds <- DESeq(dds)
resultsNames(dds)

##type_purple_vs_green
res <- results(dds,
               name="type_purple_vs_green") 

summary(res)
head(res)

resOdered <- res[order(res$padj),] 
deg <- as.data.frame(resOdered)
deg <- na.omit(deg) 

allDiff4 <- deg
diffLab4 <- subset(allDiff4,
                   abs(log2FoldChange)>1 & padj<0.05)

###查看差异基因个数
dim(diffLab4)


############使用limma包自带的removeBatchEffect函数消除批次效应的PCA分析
#if(!require("limma"))
#  BiocManager::install("limma",
#                       update = F,
#                       ask = F)

library(limma)


##消除批次效应
design <- model.matrix(~type,
                     data = data)
rB_Expr0 <- removeBatchEffect(Expr0,
                           batch = data$batch,
                           design = design)


rB_Expr <- removeBatchEffect(Expr,
                              batch = data$batch,
                              design = design)


##差异分析
#fit <- lmFit(rB_Expr,design)
#fit2 <- eBayes(fit)
#allDiff3 = topTable(fit2,
#                    adjust='fdr',
#                    coef=2,
#                    number=Inf)
#diffLab3 <- subset(allDiff3,
#                  abs(logFC)>1 & adj.P.Val<0.05)

#dim(diffLab3)###查看差异基因的个数


##PCA分析
af2.pca <- PCA(t(rB_Expr0),
               graph = FALSE
               )

af2.pca <- PCA(t(rB_Expr),
               graph = FALSE)


af2.pca <- PCA(t(rrr1),
               graph = FALSE)


fviz_pca_ind(af2.pca,
             geom = "point",
             col.ind = data$batch,
             addEllipses = TRUE,
             legend.title ="Group"
             )


##removeBatchEffect后数据取整并将小于0的数改为0
rrr0 <- round(rB_Expr0,0)

rrr0[rrr0<0] <- 0

rrr1 <- round(rB_Expr,3)

rrr1[rrr1<0] <- 0

write.table(rrr0,
            file = 'mRNA.counts.removeBatchEffect.matrix', 
            sep = '\t', 
            quote = F, 
            row.names = T)

write.table(rrr1,
            file = 'mRNA.TMM.EXPR.removeBatchEffect.matrix', 
            sep = '\t', 
            quote = F, 
            row.names = T)
