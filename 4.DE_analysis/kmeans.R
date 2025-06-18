#====================kmeans==========================

if(!requireNamespace('pheatmap',quietly = TRUE))  install.packages("pheatmap",update=F)
if(!requireNamespace('RColorBrewer',quietly = TRUE))  install.packages('RColorBrewer',update=F)
if(!requireNamespace('fpc',quietly = TRUE))  install.packages('fpc',update=F)
if(!requireNamespace('NbClust',quietly = TRUE))  install.packages('NbClust',update=F)

rm(list=ls())

#读入数据
expr = read.delim('mRNA.TMM.EDSall.removeBE_spike_means.matrix',
                header=T,row.names=1,sep="\t"
                )



expr <- round(expr,3)

expr = expr[rowSums(expr)!=0,] 


expr = scale(expr,center=T,scale=T)

#===选K值，多个指标估计分类数==========
library(NbClust)

system.time(nc<-NbClust(expr,min.nc = 2,max.nc = 20,method = "kmeans"))
barplot(table(nc$Best.nc[1,]),xlab = "聚类数",ylab = "支持指标数")
Selected_K=unname(which.max(table(nc$Best.nc[1,])))
