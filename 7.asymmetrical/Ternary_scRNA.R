rm(list=ls())


library(tidyverse)
library(ggtern)
library(hrbrthemes)
library(ggtext)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(scales)
library(tidyr)
library(matrixStats)
library(ggrepel)
require(pals)
library(rdist)
library(rpart)
library(ggthemes)
library(ComplexHeatmap)
library(ggdendro)
library(ggsci)
library(pals)
library(reshape2)

#####################################

raw.expr <- read.delim("data/Meristem_scRNA.ABD.txt",
                       header = T)


raw.expr.1 <- mutate_all(raw.expr, ~replace(., is.na(.), 0))

# my_data <- mutate_at(my_data, c(`C1`, `C4`), ~replace(., is.na(.), 0))
# my_data[is.na(my_data)] <- 0
# raw.expr.2  <- mutate_all(raw.expr.1 ,  ~replace(., .<0, 0))

C0<-raw.expr.1

###1.draw the ternary diagram for the expression of triads in bulk cells of wheat root, and group the triads into 7 groups, i.e,
#Central (Balance), A.dominant, B.dominant, D.dominant, A.suppressed, B.suppressed, D.suppressed.
#load files
#C0<-read.csv("C0_abd_triads_R1.csv",header=T)
#count the ratio for each triad
cluster0<-cbind(C0[,2],C0[,3],C0[,4])
cluster0  <- round(cluster0,2)
# cluster0<- dplyr::mutate_all(cluster0,as.numeric)
head(cluster0)
sc0<-rowSums(cluster0)
# max0 <- rowMaxs(cluster0)
max0 <- rowSums(cluster0)
head(max0)
aai<-cluster0[,1]
bbi<-cluster0[,2]
ddi<-cluster0[,3]
clu0<-cbind(aai/sc0,bbi/sc0,ddi/sc0, max0)
head(clu0)
row.names(clu0)<-C0[,1]
colnames(clu0)<-c("AA","BB","DD", "size")
# clu0 <- mutate_all(clu0, ~replace(., is.na(.), 0))
# clu0 <- clu0[rowSums(clu0) > 0.5,]
clu0 <- clu0[clu0[,4] > 0.5,]
# clu0[is.na(clu0)] <- 0
# clu0 <- clu0[rowSums(clu0) > 0.5,]
head(clu0)
write.csv(clu0,"celltype_ABD/Meristem.sanyuan.expr.out.csv")
#Calculate the Euclidean distance
centers<-t(matrix(c(0.33,0.33,0.33,1,0,0,0,1,0,0,0,1,0,0.5,0.5,0.5,0,0.5,0.5,0.5,0), nrow=3))
colnames(centers)<-c("A","B","D")
rownames(centers)<-c("Central","A.dominant","B.dominant","D.dominant","A.suppressed","B.suppressed","D.suppressed")
head(centers)
tridat<-clu0[,1:3]
#rownames(tridat) <- rownames(tridat)
tridat <- data.frame(tridat)
tridat$group <- apply(tridat, 1, function(x){
  rownames(centers)[which.min(cdist(t(x), centers))]
})
tridat$size <- clu0[,4]
head(tridat)

tridat01<-matrix(unlist(tridat[rownames(tridat),c('AA','BB','DD')]), ncol=3)  %>%
  as.data.frame() %>% 
  magrittr::set_colnames(c('AA','BB','DD')) %>% 
  mutate(group=tridat$group, size=tridat$size) %>% 
  magrittr::set_rownames(rownames(tridat))
#tridat01<-na.omit(tridat01)
rowname<-rownames(tridat01)
head(rowname)
#tridat01 <- apply(tridat01,2,as.character)
write.csv(tridat01,"celltype_ABD/Meristem.tridat_cluster.csv")
# tridat01<-read.csv("C0_tridat_clusterR1.csv",header=T,row.names=1)
row.names(tridat01)<-rowname
dim(tridat01)
head(tridat01)
mypal1 = pal_lancet("lanonc", alpha = 0.7)(8)
mypal2 = pal_startrek("uniform", alpha=0.8)(7)
mypal3 = pal_npg("nrc", alpha=0.8)(10)
mypal4 = pal_nejm("default", alpha=0.8)(7)
mypal5 = pal_uchicago("default", alpha=0.8)(9)
mypal=c(mypal1, mypal2,mypal3,mypal4,mypal5)
show_col(mypal)
#require(pals)
library('pals')
cubebasis <- pal.compress(brewer.accent,n=8)
#cubebasis <- brewer.accent(8)
show_col(cubebasis)
coll<-c("#F781BF","#A65628","#4DAF4A","#984EA3","#A6CEE3","#FF7F00","#E41A1C","#377EB8")

#tridat01<-tridat01[which(tridat01[,4]!="character(0)"),]
#plot_data<-tridat01[which(tridat01[,4]=="Central"),]
#plot_data1<-tridat01[which(tridat01[,4]!="Central"),]

tridat01<-as.data.frame(tridat01)

tridat01$AA <- round(as.numeric(tridat01$AA),2)
tridat01$BB <- round(as.numeric(tridat01$BB),2)
tridat01$DD <- round(as.numeric(tridat01$DD),2)
# tridat01$size <- round((as.numeric(tridat01$size)),2)
tridat01$size <- round(log2(as.numeric(tridat01$size)),2)

p <- ggtern(as.data.frame(tridat01), mapping = aes(AA, BB, DD, size, group)) +
  geom_mask() +
  #geom_point(data=as.data.frame(tridat01), aes(size=size), color="grey", show.legend=F)  +
  #geom_point(data=plot_data, aes(size=size), color="grey", alpha=0.8, show.legend=F) +
  geom_point(data=as.data.frame(tridat01), aes(size=size, color=group), show.legend=T)  +
  scale_colour_manual(values=coll) +
  scale_size(range=c(0, 4)) +
  guides(size="none") +
  theme_bw() +
  theme(axis.text=element_blank(),
        axis.ticks=element_blank()) +
  ggtitle(paste0("Bulk_RNA_seq")) +
  theme(plot.title = element_text(hjust = 0.5, vjust = -2))

pdf("celltype_ABD/Meristem.sanyuan.pdf",width = 10, height = 10)
print(p)
dev.off()
# length(tridat01$group)
write.csv(as.data.frame(table(tridat01$group)),"celltype_ABD/Meristem.sanyuan.abd.csv")

as.data.frame(table(tridat01$group)/length(tridat01$group))
write.csv(as.data.frame(table(tridat01$group)/length(tridat01$group)),"celltype_ABD/Meristem.sanyuan.abd.bili.csv")

#####################################

raw.expr <- read.delim("data/Pith_scRNA.ABD.txt",
                       header = T)


raw.expr.1 <- mutate_all(raw.expr, ~replace(., is.na(.), 0))

# my_data <- mutate_at(my_data, c(`C1`, `C4`), ~replace(., is.na(.), 0))
# my_data[is.na(my_data)] <- 0
# raw.expr.2  <- mutate_all(raw.expr.1 ,  ~replace(., .<0, 0))

C0<-raw.expr.1

###1.draw the ternary diagram for the expression of triads in bulk cells of wheat root, and group the triads into 7 groups, i.e,
#Central (Balance), A.dominant, B.dominant, D.dominant, A.suppressed, B.suppressed, D.suppressed.
#load files
#C0<-read.csv("C0_abd_triads_R1.csv",header=T)
#count the ratio for each triad
cluster0<-cbind(C0[,2],C0[,3],C0[,4])
cluster0  <- round(cluster0,2)
# cluster0<- dplyr::mutate_all(cluster0,as.numeric)
head(cluster0)
sc0<-rowSums(cluster0)
# max0 <- rowMaxs(cluster0)
max0 <- rowSums(cluster0)
head(max0)
aai<-cluster0[,1]
bbi<-cluster0[,2]
ddi<-cluster0[,3]
clu0<-cbind(aai/sc0,bbi/sc0,ddi/sc0, max0)
head(clu0)
row.names(clu0)<-C0[,1]
colnames(clu0)<-c("AA","BB","DD", "size")
# clu0 <- mutate_all(clu0, ~replace(., is.na(.), 0))
# clu0 <- clu0[rowSums(clu0) > 0.5,]
clu0 <- clu0[clu0[,4] > 0.5,]
# clu0[is.na(clu0)] <- 0
# clu0 <- clu0[rowSums(clu0) > 0.5,]
head(clu0)
write.csv(clu0,"celltype_ABD/Pith.sanyuan.expr.out.csv")
#Calculate the Euclidean distance
centers<-t(matrix(c(0.33,0.33,0.33,1,0,0,0,1,0,0,0,1,0,0.5,0.5,0.5,0,0.5,0.5,0.5,0), nrow=3))
colnames(centers)<-c("A","B","D")
rownames(centers)<-c("Central","A.dominant","B.dominant","D.dominant","A.suppressed","B.suppressed","D.suppressed")
head(centers)
tridat<-clu0[,1:3]
#rownames(tridat) <- rownames(tridat)
tridat <- data.frame(tridat)
tridat$group <- apply(tridat, 1, function(x){
  rownames(centers)[which.min(cdist(t(x), centers))]
})
tridat$size <- clu0[,4]
head(tridat)

tridat01<-matrix(unlist(tridat[rownames(tridat),c('AA','BB','DD')]), ncol=3)  %>%
  as.data.frame() %>% 
  magrittr::set_colnames(c('AA','BB','DD')) %>% 
  mutate(group=tridat$group, size=tridat$size) %>% 
  magrittr::set_rownames(rownames(tridat))
#tridat01<-na.omit(tridat01)
rowname<-rownames(tridat01)
head(rowname)
#tridat01 <- apply(tridat01,2,as.character)
write.csv(tridat01,"celltype_ABD/Pith.tridat_cluster.csv")
# tridat01<-read.csv("C0_tridat_clusterR1.csv",header=T,row.names=1)
row.names(tridat01)<-rowname
dim(tridat01)
head(tridat01)
mypal1 = pal_lancet("lanonc", alpha = 0.7)(8)
mypal2 = pal_startrek("uniform", alpha=0.8)(7)
mypal3 = pal_npg("nrc", alpha=0.8)(10)
mypal4 = pal_nejm("default", alpha=0.8)(7)
mypal5 = pal_uchicago("default", alpha=0.8)(9)
mypal=c(mypal1, mypal2,mypal3,mypal4,mypal5)
show_col(mypal)
#require(pals)
library('pals')
cubebasis <- pal.compress(brewer.accent,n=8)
#cubebasis <- brewer.accent(8)
show_col(cubebasis)
coll<-c("#F781BF","#A65628","#4DAF4A","#984EA3","#A6CEE3","#FF7F00","#E41A1C","#377EB8")

#tridat01<-tridat01[which(tridat01[,4]!="character(0)"),]
#plot_data<-tridat01[which(tridat01[,4]=="Central"),]
#plot_data1<-tridat01[which(tridat01[,4]!="Central"),]

tridat01<-as.data.frame(tridat01)

tridat01$AA <- round(as.numeric(tridat01$AA),2)
tridat01$BB <- round(as.numeric(tridat01$BB),2)
tridat01$DD <- round(as.numeric(tridat01$DD),2)
# tridat01$size <- round((as.numeric(tridat01$size)),2)
tridat01$size <- round(log2(as.numeric(tridat01$size)),2)

p <- ggtern(as.data.frame(tridat01), mapping = aes(AA, BB, DD, size, group)) +
  geom_mask() +
  #geom_point(data=as.data.frame(tridat01), aes(size=size), color="grey", show.legend=F)  +
  #geom_point(data=plot_data, aes(size=size), color="grey", alpha=0.8, show.legend=F) +
  geom_point(data=as.data.frame(tridat01), aes(size=size, color=group), show.legend=T)  +
  scale_colour_manual(values=coll) +
  scale_size(range=c(0, 4)) +
  guides(size="none") +
  theme_bw() +
  theme(axis.text=element_blank(),
        axis.ticks=element_blank()) +
  ggtitle(paste0("Bulk_RNA_seq")) +
  theme(plot.title = element_text(hjust = 0.5, vjust = -2))

pdf("celltype_ABD/Pith.sanyuan.pdf",width = 10, height = 10)
print(p)
dev.off()
# length(tridat01$group)
write.csv(as.data.frame(table(tridat01$group)),"celltype_ABD/Pith.sanyuan.abd.csv")

as.data.frame(table(tridat01$group)/length(tridat01$group))
write.csv(as.data.frame(table(tridat01$group)/length(tridat01$group)),"celltype_ABD/Pith.sanyuan.abd.bili.csv")


#####################################

raw.expr <- read.delim("data/Protoxylem_scRNA.ABD.txt",
                       header = T)


raw.expr.1 <- mutate_all(raw.expr, ~replace(., is.na(.), 0))

# my_data <- mutate_at(my_data, c(`C1`, `C4`), ~replace(., is.na(.), 0))
# my_data[is.na(my_data)] <- 0
# raw.expr.2  <- mutate_all(raw.expr.1 ,  ~replace(., .<0, 0))

C0<-raw.expr.1

###1.draw the ternary diagram for the expression of triads in bulk cells of wheat root, and group the triads into 7 groups, i.e,
#Central (Balance), A.dominant, B.dominant, D.dominant, A.suppressed, B.suppressed, D.suppressed.
#load files
#C0<-read.csv("C0_abd_triads_R1.csv",header=T)
#count the ratio for each triad
cluster0<-cbind(C0[,2],C0[,3],C0[,4])
cluster0  <- round(cluster0,2)
# cluster0<- dplyr::mutate_all(cluster0,as.numeric)
head(cluster0)
sc0<-rowSums(cluster0)
# max0 <- rowMaxs(cluster0)
max0 <- rowSums(cluster0)
head(max0)
aai<-cluster0[,1]
bbi<-cluster0[,2]
ddi<-cluster0[,3]
clu0<-cbind(aai/sc0,bbi/sc0,ddi/sc0, max0)
head(clu0)
row.names(clu0)<-C0[,1]
colnames(clu0)<-c("AA","BB","DD", "size")
# clu0 <- mutate_all(clu0, ~replace(., is.na(.), 0))
# clu0 <- clu0[rowSums(clu0) > 0.5,]
clu0 <- clu0[clu0[,4] > 0.5,]
# clu0[is.na(clu0)] <- 0
# clu0 <- clu0[rowSums(clu0) > 0.5,]
head(clu0)
write.csv(clu0,"celltype_ABD/Protoxylem.sanyuan.expr.out.csv")
#Calculate the Euclidean distance
centers<-t(matrix(c(0.33,0.33,0.33,1,0,0,0,1,0,0,0,1,0,0.5,0.5,0.5,0,0.5,0.5,0.5,0), nrow=3))
colnames(centers)<-c("A","B","D")
rownames(centers)<-c("Central","A.dominant","B.dominant","D.dominant","A.suppressed","B.suppressed","D.suppressed")
head(centers)
tridat<-clu0[,1:3]
#rownames(tridat) <- rownames(tridat)
tridat <- data.frame(tridat)
tridat$group <- apply(tridat, 1, function(x){
  rownames(centers)[which.min(cdist(t(x), centers))]
})
tridat$size <- clu0[,4]
head(tridat)

tridat01<-matrix(unlist(tridat[rownames(tridat),c('AA','BB','DD')]), ncol=3)  %>%
  as.data.frame() %>% 
  magrittr::set_colnames(c('AA','BB','DD')) %>% 
  mutate(group=tridat$group, size=tridat$size) %>% 
  magrittr::set_rownames(rownames(tridat))
#tridat01<-na.omit(tridat01)
rowname<-rownames(tridat01)
head(rowname)
#tridat01 <- apply(tridat01,2,as.character)
write.csv(tridat01,"celltype_ABD/Protoxylem.tridat_cluster.csv")
# tridat01<-read.csv("C0_tridat_clusterR1.csv",header=T,row.names=1)
row.names(tridat01)<-rowname
dim(tridat01)
head(tridat01)
mypal1 = pal_lancet("lanonc", alpha = 0.7)(8)
mypal2 = pal_startrek("uniform", alpha=0.8)(7)
mypal3 = pal_npg("nrc", alpha=0.8)(10)
mypal4 = pal_nejm("default", alpha=0.8)(7)
mypal5 = pal_uchicago("default", alpha=0.8)(9)
mypal=c(mypal1, mypal2,mypal3,mypal4,mypal5)
show_col(mypal)
#require(pals)
library('pals')
cubebasis <- pal.compress(brewer.accent,n=8)
#cubebasis <- brewer.accent(8)
show_col(cubebasis)
coll<-c("#F781BF","#A65628","#4DAF4A","#984EA3","#A6CEE3","#FF7F00","#E41A1C","#377EB8")

#tridat01<-tridat01[which(tridat01[,4]!="character(0)"),]
#plot_data<-tridat01[which(tridat01[,4]=="Central"),]
#plot_data1<-tridat01[which(tridat01[,4]!="Central"),]

tridat01<-as.data.frame(tridat01)

tridat01$AA <- round(as.numeric(tridat01$AA),2)
tridat01$BB <- round(as.numeric(tridat01$BB),2)
tridat01$DD <- round(as.numeric(tridat01$DD),2)
# tridat01$size <- round((as.numeric(tridat01$size)),2)
tridat01$size <- round(log2(as.numeric(tridat01$size)),2)

p <- ggtern(as.data.frame(tridat01), mapping = aes(AA, BB, DD, size, group)) +
  geom_mask() +
  #geom_point(data=as.data.frame(tridat01), aes(size=size), color="grey", show.legend=F)  +
  #geom_point(data=plot_data, aes(size=size), color="grey", alpha=0.8, show.legend=F) +
  geom_point(data=as.data.frame(tridat01), aes(size=size, color=group), show.legend=T)  +
  scale_colour_manual(values=coll) +
  scale_size(range=c(0, 4)) +
  guides(size="none") +
  theme_bw() +
  theme(axis.text=element_blank(),
        axis.ticks=element_blank()) +
  ggtitle(paste0("Bulk_RNA_seq")) +
  theme(plot.title = element_text(hjust = 0.5, vjust = -2))

pdf("celltype_ABD/Protoxylem.sanyuan.pdf",width = 10, height = 10)
print(p)
dev.off()
# length(tridat01$group)
write.csv(as.data.frame(table(tridat01$group)),"celltype_ABD/Protoxylem.sanyuan.abd.csv")

as.data.frame(table(tridat01$group)/length(tridat01$group))
write.csv(as.data.frame(table(tridat01$group)/length(tridat01$group)),"celltype_ABD/Protoxylem.sanyuan.abd.bili.csv")

#####################################

raw.expr <- read.delim("data/Epidermis_scRNA.ABD.txt",
                       header = T)


raw.expr.1 <- mutate_all(raw.expr, ~replace(., is.na(.), 0))

# my_data <- mutate_at(my_data, c(`C1`, `C4`), ~replace(., is.na(.), 0))
# my_data[is.na(my_data)] <- 0
# raw.expr.2  <- mutate_all(raw.expr.1 ,  ~replace(., .<0, 0))

C0<-raw.expr.1

###1.draw the ternary diagram for the expression of triads in bulk cells of wheat root, and group the triads into 7 groups, i.e,
#Central (Balance), A.dominant, B.dominant, D.dominant, A.suppressed, B.suppressed, D.suppressed.
#load files
#C0<-read.csv("C0_abd_triads_R1.csv",header=T)
#count the ratio for each triad
cluster0<-cbind(C0[,2],C0[,3],C0[,4])
cluster0  <- round(cluster0,2)
# cluster0<- dplyr::mutate_all(cluster0,as.numeric)
head(cluster0)
sc0<-rowSums(cluster0)
# max0 <- rowMaxs(cluster0)
max0 <- rowSums(cluster0)
head(max0)
aai<-cluster0[,1]
bbi<-cluster0[,2]
ddi<-cluster0[,3]
clu0<-cbind(aai/sc0,bbi/sc0,ddi/sc0, max0)
head(clu0)
row.names(clu0)<-C0[,1]
colnames(clu0)<-c("AA","BB","DD", "size")
# clu0 <- mutate_all(clu0, ~replace(., is.na(.), 0))
# clu0 <- clu0[rowSums(clu0) > 0.5,]
clu0 <- clu0[clu0[,4] > 0.5,]
# clu0[is.na(clu0)] <- 0
# clu0 <- clu0[rowSums(clu0) > 0.5,]
head(clu0)
write.csv(clu0,"celltype_ABD/Epidermis.sanyuan.expr.out.csv")
#Calculate the Euclidean distance
centers<-t(matrix(c(0.33,0.33,0.33,1,0,0,0,1,0,0,0,1,0,0.5,0.5,0.5,0,0.5,0.5,0.5,0), nrow=3))
colnames(centers)<-c("A","B","D")
rownames(centers)<-c("Central","A.dominant","B.dominant","D.dominant","A.suppressed","B.suppressed","D.suppressed")
head(centers)
tridat<-clu0[,1:3]
#rownames(tridat) <- rownames(tridat)
tridat <- data.frame(tridat)
tridat$group <- apply(tridat, 1, function(x){
  rownames(centers)[which.min(cdist(t(x), centers))]
})
tridat$size <- clu0[,4]
head(tridat)

tridat01<-matrix(unlist(tridat[rownames(tridat),c('AA','BB','DD')]), ncol=3)  %>%
  as.data.frame() %>% 
  magrittr::set_colnames(c('AA','BB','DD')) %>% 
  mutate(group=tridat$group, size=tridat$size) %>% 
  magrittr::set_rownames(rownames(tridat))
#tridat01<-na.omit(tridat01)
rowname<-rownames(tridat01)
head(rowname)
#tridat01 <- apply(tridat01,2,as.character)
write.csv(tridat01,"celltype_ABD/Epidermis.tridat_cluster.csv")
# tridat01<-read.csv("C0_tridat_clusterR1.csv",header=T,row.names=1)
row.names(tridat01)<-rowname
dim(tridat01)
head(tridat01)
mypal1 = pal_lancet("lanonc", alpha = 0.7)(8)
mypal2 = pal_startrek("uniform", alpha=0.8)(7)
mypal3 = pal_npg("nrc", alpha=0.8)(10)
mypal4 = pal_nejm("default", alpha=0.8)(7)
mypal5 = pal_uchicago("default", alpha=0.8)(9)
mypal=c(mypal1, mypal2,mypal3,mypal4,mypal5)
show_col(mypal)
#require(pals)
library('pals')
cubebasis <- pal.compress(brewer.accent,n=8)
#cubebasis <- brewer.accent(8)
show_col(cubebasis)
coll<-c("#F781BF","#A65628","#4DAF4A","#984EA3","#A6CEE3","#FF7F00","#E41A1C","#377EB8")

#tridat01<-tridat01[which(tridat01[,4]!="character(0)"),]
#plot_data<-tridat01[which(tridat01[,4]=="Central"),]
#plot_data1<-tridat01[which(tridat01[,4]!="Central"),]

tridat01<-as.data.frame(tridat01)

tridat01$AA <- round(as.numeric(tridat01$AA),2)
tridat01$BB <- round(as.numeric(tridat01$BB),2)
tridat01$DD <- round(as.numeric(tridat01$DD),2)
# tridat01$size <- round((as.numeric(tridat01$size)),2)
tridat01$size <- round(log2(as.numeric(tridat01$size)),2)

p <- ggtern(as.data.frame(tridat01), mapping = aes(AA, BB, DD, size, group)) +
  geom_mask() +
  #geom_point(data=as.data.frame(tridat01), aes(size=size), color="grey", show.legend=F)  +
  #geom_point(data=plot_data, aes(size=size), color="grey", alpha=0.8, show.legend=F) +
  geom_point(data=as.data.frame(tridat01), aes(size=size, color=group), show.legend=T)  +
  scale_colour_manual(values=coll) +
  scale_size(range=c(0, 4)) +
  guides(size="none") +
  theme_bw() +
  theme(axis.text=element_blank(),
        axis.ticks=element_blank()) +
  ggtitle(paste0("Bulk_RNA_seq")) +
  theme(plot.title = element_text(hjust = 0.5, vjust = -2))

pdf("celltype_ABD/Epidermis.sanyuan.pdf",width = 10, height = 10)
print(p)
dev.off()
# length(tridat01$group)
write.csv(as.data.frame(table(tridat01$group)),"celltype_ABD/Epidermis.sanyuan.abd.csv")

as.data.frame(table(tridat01$group)/length(tridat01$group))
write.csv(as.data.frame(table(tridat01$group)/length(tridat01$group)),"celltype_ABD/Epidermis.sanyuan.abd.bili.csv")

#####################################

raw.expr <- read.delim("data/Protophloem_scRNA.ABD.txt",
                       header = T)


raw.expr.1 <- mutate_all(raw.expr, ~replace(., is.na(.), 0))

# my_data <- mutate_at(my_data, c(`C1`, `C4`), ~replace(., is.na(.), 0))
# my_data[is.na(my_data)] <- 0
# raw.expr.2  <- mutate_all(raw.expr.1 ,  ~replace(., .<0, 0))

C0<-raw.expr.1

###1.draw the ternary diagram for the expression of triads in bulk cells of wheat root, and group the triads into 7 groups, i.e,
#Central (Balance), A.dominant, B.dominant, D.dominant, A.suppressed, B.suppressed, D.suppressed.
#load files
#C0<-read.csv("C0_abd_triads_R1.csv",header=T)
#count the ratio for each triad
cluster0<-cbind(C0[,2],C0[,3],C0[,4])
cluster0  <- round(cluster0,2)
# cluster0<- dplyr::mutate_all(cluster0,as.numeric)
head(cluster0)
sc0<-rowSums(cluster0)
# max0 <- rowMaxs(cluster0)
max0 <- rowSums(cluster0)
head(max0)
aai<-cluster0[,1]
bbi<-cluster0[,2]
ddi<-cluster0[,3]
clu0<-cbind(aai/sc0,bbi/sc0,ddi/sc0, max0)
head(clu0)
row.names(clu0)<-C0[,1]
colnames(clu0)<-c("AA","BB","DD", "size")
# clu0 <- mutate_all(clu0, ~replace(., is.na(.), 0))
# clu0 <- clu0[rowSums(clu0) > 0.5,]
clu0 <- clu0[clu0[,4] > 0.5,]
# clu0[is.na(clu0)] <- 0
# clu0 <- clu0[rowSums(clu0) > 0.5,]
head(clu0)
write.csv(clu0,"celltype_ABD/Protophloem.sanyuan.expr.out.csv")
#Calculate the Euclidean distance
centers<-t(matrix(c(0.33,0.33,0.33,1,0,0,0,1,0,0,0,1,0,0.5,0.5,0.5,0,0.5,0.5,0.5,0), nrow=3))
colnames(centers)<-c("A","B","D")
rownames(centers)<-c("Central","A.dominant","B.dominant","D.dominant","A.suppressed","B.suppressed","D.suppressed")
head(centers)
tridat<-clu0[,1:3]
#rownames(tridat) <- rownames(tridat)
tridat <- data.frame(tridat)
tridat$group <- apply(tridat, 1, function(x){
  rownames(centers)[which.min(cdist(t(x), centers))]
})
tridat$size <- clu0[,4]
head(tridat)

tridat01<-matrix(unlist(tridat[rownames(tridat),c('AA','BB','DD')]), ncol=3)  %>%
  as.data.frame() %>% 
  magrittr::set_colnames(c('AA','BB','DD')) %>% 
  mutate(group=tridat$group, size=tridat$size) %>% 
  magrittr::set_rownames(rownames(tridat))
#tridat01<-na.omit(tridat01)
rowname<-rownames(tridat01)
head(rowname)
#tridat01 <- apply(tridat01,2,as.character)
write.csv(tridat01,"celltype_ABD/Protophloem.tridat_cluster.csv")
# tridat01<-read.csv("C0_tridat_clusterR1.csv",header=T,row.names=1)
row.names(tridat01)<-rowname
dim(tridat01)
head(tridat01)
mypal1 = pal_lancet("lanonc", alpha = 0.7)(8)
mypal2 = pal_startrek("uniform", alpha=0.8)(7)
mypal3 = pal_npg("nrc", alpha=0.8)(10)
mypal4 = pal_nejm("default", alpha=0.8)(7)
mypal5 = pal_uchicago("default", alpha=0.8)(9)
mypal=c(mypal1, mypal2,mypal3,mypal4,mypal5)
show_col(mypal)
#require(pals)
library('pals')
cubebasis <- pal.compress(brewer.accent,n=8)
#cubebasis <- brewer.accent(8)
show_col(cubebasis)
coll<-c("#F781BF","#A65628","#4DAF4A","#984EA3","#A6CEE3","#FF7F00","#E41A1C","#377EB8")

#tridat01<-tridat01[which(tridat01[,4]!="character(0)"),]
#plot_data<-tridat01[which(tridat01[,4]=="Central"),]
#plot_data1<-tridat01[which(tridat01[,4]!="Central"),]

tridat01<-as.data.frame(tridat01)

tridat01$AA <- round(as.numeric(tridat01$AA),2)
tridat01$BB <- round(as.numeric(tridat01$BB),2)
tridat01$DD <- round(as.numeric(tridat01$DD),2)
# tridat01$size <- round((as.numeric(tridat01$size)),2)
tridat01$size <- round(log2(as.numeric(tridat01$size)),2)

p <- ggtern(as.data.frame(tridat01), mapping = aes(AA, BB, DD, size, group)) +
  geom_mask() +
  #geom_point(data=as.data.frame(tridat01), aes(size=size), color="grey", show.legend=F)  +
  #geom_point(data=plot_data, aes(size=size), color="grey", alpha=0.8, show.legend=F) +
  geom_point(data=as.data.frame(tridat01), aes(size=size, color=group), show.legend=T)  +
  scale_colour_manual(values=coll) +
  scale_size(range=c(0, 4)) +
  guides(size="none") +
  theme_bw() +
  theme(axis.text=element_blank(),
        axis.ticks=element_blank()) +
  ggtitle(paste0("Bulk_RNA_seq")) +
  theme(plot.title = element_text(hjust = 0.5, vjust = -2))

pdf("celltype_ABD/Protophloem.sanyuan.pdf",width = 10, height = 10)
print(p)
dev.off()
# length(tridat01$group)
write.csv(as.data.frame(table(tridat01$group)),"celltype_ABD/Protophloem.sanyuan.abd.csv")

as.data.frame(table(tridat01$group)/length(tridat01$group))
write.csv(as.data.frame(table(tridat01$group)/length(tridat01$group)),"celltype_ABD/Protophloem.sanyuan.abd.bili.csv")

#####################################

raw.expr <- read.delim("data/Cortex_scRNA.ABD.txt",
                       header = T)


raw.expr.1 <- mutate_all(raw.expr, ~replace(., is.na(.), 0))

# my_data <- mutate_at(my_data, c(`C1`, `C4`), ~replace(., is.na(.), 0))
# my_data[is.na(my_data)] <- 0
# raw.expr.2  <- mutate_all(raw.expr.1 ,  ~replace(., .<0, 0))

C0<-raw.expr.1

###1.draw the ternary diagram for the expression of triads in bulk cells of wheat root, and group the triads into 7 groups, i.e,
#Central (Balance), A.dominant, B.dominant, D.dominant, A.suppressed, B.suppressed, D.suppressed.
#load files
#C0<-read.csv("C0_abd_triads_R1.csv",header=T)
#count the ratio for each triad
cluster0<-cbind(C0[,2],C0[,3],C0[,4])
cluster0  <- round(cluster0,2)
# cluster0<- dplyr::mutate_all(cluster0,as.numeric)
head(cluster0)
sc0<-rowSums(cluster0)
# max0 <- rowMaxs(cluster0)
max0 <- rowSums(cluster0)
head(max0)
aai<-cluster0[,1]
bbi<-cluster0[,2]
ddi<-cluster0[,3]
clu0<-cbind(aai/sc0,bbi/sc0,ddi/sc0, max0)
head(clu0)
row.names(clu0)<-C0[,1]
colnames(clu0)<-c("AA","BB","DD", "size")
# clu0 <- mutate_all(clu0, ~replace(., is.na(.), 0))
# clu0 <- clu0[rowSums(clu0) > 0.5,]
clu0 <- clu0[clu0[,4] > 0.5,]
# clu0[is.na(clu0)] <- 0
# clu0 <- clu0[rowSums(clu0) > 0.5,]
head(clu0)
write.csv(clu0,"celltype_ABD/Cortex.sanyuan.expr.out.csv")
#Calculate the Euclidean distance
centers<-t(matrix(c(0.33,0.33,0.33,1,0,0,0,1,0,0,0,1,0,0.5,0.5,0.5,0,0.5,0.5,0.5,0), nrow=3))
colnames(centers)<-c("A","B","D")
rownames(centers)<-c("Central","A.dominant","B.dominant","D.dominant","A.suppressed","B.suppressed","D.suppressed")
head(centers)
tridat<-clu0[,1:3]
#rownames(tridat) <- rownames(tridat)
tridat <- data.frame(tridat)
tridat$group <- apply(tridat, 1, function(x){
  rownames(centers)[which.min(cdist(t(x), centers))]
})
tridat$size <- clu0[,4]
head(tridat)

tridat01<-matrix(unlist(tridat[rownames(tridat),c('AA','BB','DD')]), ncol=3)  %>%
  as.data.frame() %>% 
  magrittr::set_colnames(c('AA','BB','DD')) %>% 
  mutate(group=tridat$group, size=tridat$size) %>% 
  magrittr::set_rownames(rownames(tridat))
#tridat01<-na.omit(tridat01)
rowname<-rownames(tridat01)
head(rowname)
#tridat01 <- apply(tridat01,2,as.character)
write.csv(tridat01,"celltype_ABD/Cortex.tridat_cluster.csv")
# tridat01<-read.csv("C0_tridat_clusterR1.csv",header=T,row.names=1)
row.names(tridat01)<-rowname
dim(tridat01)
head(tridat01)
mypal1 = pal_lancet("lanonc", alpha = 0.7)(8)
mypal2 = pal_startrek("uniform", alpha=0.8)(7)
mypal3 = pal_npg("nrc", alpha=0.8)(10)
mypal4 = pal_nejm("default", alpha=0.8)(7)
mypal5 = pal_uchicago("default", alpha=0.8)(9)
mypal=c(mypal1, mypal2,mypal3,mypal4,mypal5)
show_col(mypal)
#require(pals)
library('pals')
cubebasis <- pal.compress(brewer.accent,n=8)
#cubebasis <- brewer.accent(8)
show_col(cubebasis)
coll<-c("#F781BF","#A65628","#4DAF4A","#984EA3","#A6CEE3","#FF7F00","#E41A1C","#377EB8")

#tridat01<-tridat01[which(tridat01[,4]!="character(0)"),]
#plot_data<-tridat01[which(tridat01[,4]=="Central"),]
#plot_data1<-tridat01[which(tridat01[,4]!="Central"),]

tridat01<-as.data.frame(tridat01)

tridat01$AA <- round(as.numeric(tridat01$AA),2)
tridat01$BB <- round(as.numeric(tridat01$BB),2)
tridat01$DD <- round(as.numeric(tridat01$DD),2)
# tridat01$size <- round((as.numeric(tridat01$size)),2)
tridat01$size <- round(log2(as.numeric(tridat01$size)),2)

p <- ggtern(as.data.frame(tridat01), mapping = aes(AA, BB, DD, size, group)) +
  geom_mask() +
  #geom_point(data=as.data.frame(tridat01), aes(size=size), color="grey", show.legend=F)  +
  #geom_point(data=plot_data, aes(size=size), color="grey", alpha=0.8, show.legend=F) +
  geom_point(data=as.data.frame(tridat01), aes(size=size, color=group), show.legend=T)  +
  scale_colour_manual(values=coll) +
  scale_size(range=c(0, 4)) +
  guides(size="none") +
  theme_bw() +
  theme(axis.text=element_blank(),
        axis.ticks=element_blank()) +
  ggtitle(paste0("Bulk_RNA_seq")) +
  theme(plot.title = element_text(hjust = 0.5, vjust = -2))

pdf("celltype_ABD/Cortex.sanyuan.pdf",width = 10, height = 10)
print(p)
dev.off()
# length(tridat01$group)
write.csv(as.data.frame(table(tridat01$group)),"celltype_ABD/Cortex.sanyuan.abd.csv")

as.data.frame(table(tridat01$group)/length(tridat01$group))
write.csv(as.data.frame(table(tridat01$group)/length(tridat01$group)),"celltype_ABD/Cortex.sanyuan.abd.bili.csv")

#####################################

raw.expr <- read.delim("data/Bundle_sheath_scRNA.ABD.txt",
                       header = T)


raw.expr.1 <- mutate_all(raw.expr, ~replace(., is.na(.), 0))

# my_data <- mutate_at(my_data, c(`C1`, `C4`), ~replace(., is.na(.), 0))
# my_data[is.na(my_data)] <- 0
# raw.expr.2  <- mutate_all(raw.expr.1 ,  ~replace(., .<0, 0))

C0<-raw.expr.1

###1.draw the ternary diagram for the expression of triads in bulk cells of wheat root, and group the triads into 7 groups, i.e,
#Central (Balance), A.dominant, B.dominant, D.dominant, A.suppressed, B.suppressed, D.suppressed.
#load files
#C0<-read.csv("C0_abd_triads_R1.csv",header=T)
#count the ratio for each triad
cluster0<-cbind(C0[,2],C0[,3],C0[,4])
cluster0  <- round(cluster0,2)
# cluster0<- dplyr::mutate_all(cluster0,as.numeric)
head(cluster0)
sc0<-rowSums(cluster0)
# max0 <- rowMaxs(cluster0)
max0 <- rowSums(cluster0)
head(max0)
aai<-cluster0[,1]
bbi<-cluster0[,2]
ddi<-cluster0[,3]
clu0<-cbind(aai/sc0,bbi/sc0,ddi/sc0, max0)
head(clu0)
row.names(clu0)<-C0[,1]
colnames(clu0)<-c("AA","BB","DD", "size")
# clu0 <- mutate_all(clu0, ~replace(., is.na(.), 0))
# clu0 <- clu0[rowSums(clu0) > 0.5,]
clu0 <- clu0[clu0[,4] > 0.5,]
# clu0[is.na(clu0)] <- 0
# clu0 <- clu0[rowSums(clu0) > 0.5,]
head(clu0)
write.csv(clu0,"celltype_ABD/Bundle_sheath.sanyuan.expr.out.csv")
#Calculate the Euclidean distance
centers<-t(matrix(c(0.33,0.33,0.33,1,0,0,0,1,0,0,0,1,0,0.5,0.5,0.5,0,0.5,0.5,0.5,0), nrow=3))
colnames(centers)<-c("A","B","D")
rownames(centers)<-c("Central","A.dominant","B.dominant","D.dominant","A.suppressed","B.suppressed","D.suppressed")
head(centers)
tridat<-clu0[,1:3]
#rownames(tridat) <- rownames(tridat)
tridat <- data.frame(tridat)
tridat$group <- apply(tridat, 1, function(x){
  rownames(centers)[which.min(cdist(t(x), centers))]
})
tridat$size <- clu0[,4]
head(tridat)

tridat01<-matrix(unlist(tridat[rownames(tridat),c('AA','BB','DD')]), ncol=3)  %>%
  as.data.frame() %>% 
  magrittr::set_colnames(c('AA','BB','DD')) %>% 
  mutate(group=tridat$group, size=tridat$size) %>% 
  magrittr::set_rownames(rownames(tridat))
#tridat01<-na.omit(tridat01)
rowname<-rownames(tridat01)
head(rowname)
#tridat01 <- apply(tridat01,2,as.character)
write.csv(tridat01,"celltype_ABD/Bundle_sheath.tridat_cluster.csv")
# tridat01<-read.csv("C0_tridat_clusterR1.csv",header=T,row.names=1)
row.names(tridat01)<-rowname
dim(tridat01)
head(tridat01)
mypal1 = pal_lancet("lanonc", alpha = 0.7)(8)
mypal2 = pal_startrek("uniform", alpha=0.8)(7)
mypal3 = pal_npg("nrc", alpha=0.8)(10)
mypal4 = pal_nejm("default", alpha=0.8)(7)
mypal5 = pal_uchicago("default", alpha=0.8)(9)
mypal=c(mypal1, mypal2,mypal3,mypal4,mypal5)
show_col(mypal)
#require(pals)
library('pals')
cubebasis <- pal.compress(brewer.accent,n=8)
#cubebasis <- brewer.accent(8)
show_col(cubebasis)
coll<-c("#F781BF","#A65628","#4DAF4A","#984EA3","#A6CEE3","#FF7F00","#E41A1C","#377EB8")

#tridat01<-tridat01[which(tridat01[,4]!="character(0)"),]
#plot_data<-tridat01[which(tridat01[,4]=="Central"),]
#plot_data1<-tridat01[which(tridat01[,4]!="Central"),]

tridat01<-as.data.frame(tridat01)

tridat01$AA <- round(as.numeric(tridat01$AA),2)
tridat01$BB <- round(as.numeric(tridat01$BB),2)
tridat01$DD <- round(as.numeric(tridat01$DD),2)
# tridat01$size <- round((as.numeric(tridat01$size)),2)
tridat01$size <- round(log2(as.numeric(tridat01$size)),2)

p <- ggtern(as.data.frame(tridat01), mapping = aes(AA, BB, DD, size, group)) +
  geom_mask() +
  #geom_point(data=as.data.frame(tridat01), aes(size=size), color="grey", show.legend=F)  +
  #geom_point(data=plot_data, aes(size=size), color="grey", alpha=0.8, show.legend=F) +
  geom_point(data=as.data.frame(tridat01), aes(size=size, color=group), show.legend=T)  +
  scale_colour_manual(values=coll) +
  scale_size(range=c(0, 4)) +
  guides(size="none") +
  theme_bw() +
  theme(axis.text=element_blank(),
        axis.ticks=element_blank()) +
  ggtitle(paste0("Bulk_RNA_seq")) +
  theme(plot.title = element_text(hjust = 0.5, vjust = -2))

pdf("celltype_ABD/Bundle_sheath.sanyuan.pdf",width = 10, height = 10)
print(p)
dev.off()
# length(tridat01$group)
write.csv(as.data.frame(table(tridat01$group)),"celltype_ABD/Bundle_sheath.sanyuan.abd.csv")

as.data.frame(table(tridat01$group)/length(tridat01$group))
write.csv(as.data.frame(table(tridat01$group)/length(tridat01$group)),"celltype_ABD/Bundle_sheath.sanyuan.abd.bili.csv")

#####################################

raw.expr <- read.delim("data/Vascular_initial_cells_scRNA.ABD.txt",
                       header = T)


raw.expr.1 <- mutate_all(raw.expr, ~replace(., is.na(.), 0))

# my_data <- mutate_at(my_data, c(`C1`, `C4`), ~replace(., is.na(.), 0))
# my_data[is.na(my_data)] <- 0
# raw.expr.2  <- mutate_all(raw.expr.1 ,  ~replace(., .<0, 0))

C0<-raw.expr.1

###1.draw the ternary diagram for the expression of triads in bulk cells of wheat root, and group the triads into 7 groups, i.e,
#Central (Balance), A.dominant, B.dominant, D.dominant, A.suppressed, B.suppressed, D.suppressed.
#load files
#C0<-read.csv("C0_abd_triads_R1.csv",header=T)
#count the ratio for each triad
cluster0<-cbind(C0[,2],C0[,3],C0[,4])
cluster0  <- round(cluster0,2)
# cluster0<- dplyr::mutate_all(cluster0,as.numeric)
head(cluster0)
sc0<-rowSums(cluster0)
# max0 <- rowMaxs(cluster0)
max0 <- rowSums(cluster0)
head(max0)
aai<-cluster0[,1]
bbi<-cluster0[,2]
ddi<-cluster0[,3]
clu0<-cbind(aai/sc0,bbi/sc0,ddi/sc0, max0)
head(clu0)
row.names(clu0)<-C0[,1]
colnames(clu0)<-c("AA","BB","DD", "size")
# clu0 <- mutate_all(clu0, ~replace(., is.na(.), 0))
# clu0 <- clu0[rowSums(clu0) > 0.5,]
clu0 <- clu0[clu0[,4] > 0.5,]
# clu0[is.na(clu0)] <- 0
# clu0 <- clu0[rowSums(clu0) > 0.5,]
head(clu0)
write.csv(clu0,"celltype_ABD/Vascular_initial_cells.sanyuan.expr.out.csv")
#Calculate the Euclidean distance
centers<-t(matrix(c(0.33,0.33,0.33,1,0,0,0,1,0,0,0,1,0,0.5,0.5,0.5,0,0.5,0.5,0.5,0), nrow=3))
colnames(centers)<-c("A","B","D")
rownames(centers)<-c("Central","A.dominant","B.dominant","D.dominant","A.suppressed","B.suppressed","D.suppressed")
head(centers)
tridat<-clu0[,1:3]
#rownames(tridat) <- rownames(tridat)
tridat <- data.frame(tridat)
tridat$group <- apply(tridat, 1, function(x){
  rownames(centers)[which.min(cdist(t(x), centers))]
})
tridat$size <- clu0[,4]
head(tridat)

tridat01<-matrix(unlist(tridat[rownames(tridat),c('AA','BB','DD')]), ncol=3)  %>%
  as.data.frame() %>% 
  magrittr::set_colnames(c('AA','BB','DD')) %>% 
  mutate(group=tridat$group, size=tridat$size) %>% 
  magrittr::set_rownames(rownames(tridat))
#tridat01<-na.omit(tridat01)
rowname<-rownames(tridat01)
head(rowname)
#tridat01 <- apply(tridat01,2,as.character)
write.csv(tridat01,"celltype_ABD/Vascular_initial_cells.tridat_cluster.csv")
# tridat01<-read.csv("C0_tridat_clusterR1.csv",header=T,row.names=1)
row.names(tridat01)<-rowname
dim(tridat01)
head(tridat01)
mypal1 = pal_lancet("lanonc", alpha = 0.7)(8)
mypal2 = pal_startrek("uniform", alpha=0.8)(7)
mypal3 = pal_npg("nrc", alpha=0.8)(10)
mypal4 = pal_nejm("default", alpha=0.8)(7)
mypal5 = pal_uchicago("default", alpha=0.8)(9)
mypal=c(mypal1, mypal2,mypal3,mypal4,mypal5)
show_col(mypal)
#require(pals)
library('pals')
cubebasis <- pal.compress(brewer.accent,n=8)
#cubebasis <- brewer.accent(8)
show_col(cubebasis)
coll<-c("#F781BF","#A65628","#4DAF4A","#984EA3","#A6CEE3","#FF7F00","#E41A1C","#377EB8")

#tridat01<-tridat01[which(tridat01[,4]!="character(0)"),]
#plot_data<-tridat01[which(tridat01[,4]=="Central"),]
#plot_data1<-tridat01[which(tridat01[,4]!="Central"),]

tridat01<-as.data.frame(tridat01)

tridat01$AA <- round(as.numeric(tridat01$AA),2)
tridat01$BB <- round(as.numeric(tridat01$BB),2)
tridat01$DD <- round(as.numeric(tridat01$DD),2)
# tridat01$size <- round((as.numeric(tridat01$size)),2)
tridat01$size <- round(log2(as.numeric(tridat01$size)),2)

p <- ggtern(as.data.frame(tridat01), mapping = aes(AA, BB, DD, size, group)) +
  geom_mask() +
  #geom_point(data=as.data.frame(tridat01), aes(size=size), color="grey", show.legend=F)  +
  #geom_point(data=plot_data, aes(size=size), color="grey", alpha=0.8, show.legend=F) +
  geom_point(data=as.data.frame(tridat01), aes(size=size, color=group), show.legend=T)  +
  scale_colour_manual(values=coll) +
  scale_size(range=c(0, 4)) +
  guides(size="none") +
  theme_bw() +
  theme(axis.text=element_blank(),
        axis.ticks=element_blank()) +
  ggtitle(paste0("Bulk_RNA_seq")) +
  theme(plot.title = element_text(hjust = 0.5, vjust = -2))

pdf("celltype_ABD/Vascular_initial_cells.sanyuan.pdf",width = 10, height = 10)
print(p)
dev.off()
# length(tridat01$group)
write.csv(as.data.frame(table(tridat01$group)),"celltype_ABD/Vascular_initial_cells.sanyuan.abd.csv")

as.data.frame(table(tridat01$group)/length(tridat01$group))
write.csv(as.data.frame(table(tridat01$group)/length(tridat01$group)),"celltype_ABD/Vascular_initial_cells.sanyuan.abd.bili.csv")

#####################################

raw.expr <- read.delim("data/Promeristem_scRNA.ABD.txt",
                       header = T)


raw.expr.1 <- mutate_all(raw.expr, ~replace(., is.na(.), 0))

# my_data <- mutate_at(my_data, c(`C1`, `C4`), ~replace(., is.na(.), 0))
# my_data[is.na(my_data)] <- 0
# raw.expr.2  <- mutate_all(raw.expr.1 ,  ~replace(., .<0, 0))

C0<-raw.expr.1

###1.draw the ternary diagram for the expression of triads in bulk cells of wheat root, and group the triads into 7 groups, i.e,
#Central (Balance), A.dominant, B.dominant, D.dominant, A.suppressed, B.suppressed, D.suppressed.
#load files
#C0<-read.csv("C0_abd_triads_R1.csv",header=T)
#count the ratio for each triad
cluster0<-cbind(C0[,2],C0[,3],C0[,4])
cluster0  <- round(cluster0,2)
# cluster0<- dplyr::mutate_all(cluster0,as.numeric)
head(cluster0)
sc0<-rowSums(cluster0)
# max0 <- rowMaxs(cluster0)
max0 <- rowSums(cluster0)
head(max0)
aai<-cluster0[,1]
bbi<-cluster0[,2]
ddi<-cluster0[,3]
clu0<-cbind(aai/sc0,bbi/sc0,ddi/sc0, max0)
head(clu0)
row.names(clu0)<-C0[,1]
colnames(clu0)<-c("AA","BB","DD", "size")
# clu0 <- mutate_all(clu0, ~replace(., is.na(.), 0))
# clu0 <- clu0[rowSums(clu0) > 0.5,]
clu0 <- clu0[clu0[,4] > 0.5,]
# clu0[is.na(clu0)] <- 0
# clu0 <- clu0[rowSums(clu0) > 0.5,]
head(clu0)
write.csv(clu0,"celltype_ABD/Promeristem.sanyuan.expr.out.csv")
#Calculate the Euclidean distance
centers<-t(matrix(c(0.33,0.33,0.33,1,0,0,0,1,0,0,0,1,0,0.5,0.5,0.5,0,0.5,0.5,0.5,0), nrow=3))
colnames(centers)<-c("A","B","D")
rownames(centers)<-c("Central","A.dominant","B.dominant","D.dominant","A.suppressed","B.suppressed","D.suppressed")
head(centers)
tridat<-clu0[,1:3]
#rownames(tridat) <- rownames(tridat)
tridat <- data.frame(tridat)
tridat$group <- apply(tridat, 1, function(x){
  rownames(centers)[which.min(cdist(t(x), centers))]
})
tridat$size <- clu0[,4]
head(tridat)

tridat01<-matrix(unlist(tridat[rownames(tridat),c('AA','BB','DD')]), ncol=3)  %>%
  as.data.frame() %>% 
  magrittr::set_colnames(c('AA','BB','DD')) %>% 
  mutate(group=tridat$group, size=tridat$size) %>% 
  magrittr::set_rownames(rownames(tridat))
#tridat01<-na.omit(tridat01)
rowname<-rownames(tridat01)
head(rowname)
#tridat01 <- apply(tridat01,2,as.character)
write.csv(tridat01,"celltype_ABD/Promeristem.tridat_cluster.csv")
# tridat01<-read.csv("C0_tridat_clusterR1.csv",header=T,row.names=1)
row.names(tridat01)<-rowname
dim(tridat01)
head(tridat01)
mypal1 = pal_lancet("lanonc", alpha = 0.7)(8)
mypal2 = pal_startrek("uniform", alpha=0.8)(7)
mypal3 = pal_npg("nrc", alpha=0.8)(10)
mypal4 = pal_nejm("default", alpha=0.8)(7)
mypal5 = pal_uchicago("default", alpha=0.8)(9)
mypal=c(mypal1, mypal2,mypal3,mypal4,mypal5)
show_col(mypal)
#require(pals)
library('pals')
cubebasis <- pal.compress(brewer.accent,n=8)
#cubebasis <- brewer.accent(8)
show_col(cubebasis)
coll<-c("#F781BF","#A65628","#4DAF4A","#984EA3","#A6CEE3","#FF7F00","#E41A1C","#377EB8")

#tridat01<-tridat01[which(tridat01[,4]!="character(0)"),]
#plot_data<-tridat01[which(tridat01[,4]=="Central"),]
#plot_data1<-tridat01[which(tridat01[,4]!="Central"),]

tridat01<-as.data.frame(tridat01)

tridat01$AA <- round(as.numeric(tridat01$AA),2)
tridat01$BB <- round(as.numeric(tridat01$BB),2)
tridat01$DD <- round(as.numeric(tridat01$DD),2)
# tridat01$size <- round((as.numeric(tridat01$size)),2)
tridat01$size <- round(log2(as.numeric(tridat01$size)),2)

p <- ggtern(as.data.frame(tridat01), mapping = aes(AA, BB, DD, size, group)) +
  geom_mask() +
  #geom_point(data=as.data.frame(tridat01), aes(size=size), color="grey", show.legend=F)  +
  #geom_point(data=plot_data, aes(size=size), color="grey", alpha=0.8, show.legend=F) +
  geom_point(data=as.data.frame(tridat01), aes(size=size, color=group), show.legend=T)  +
  scale_colour_manual(values=coll) +
  scale_size(range=c(0, 4)) +
  guides(size="none") +
  theme_bw() +
  theme(axis.text=element_blank(),
        axis.ticks=element_blank()) +
  ggtitle(paste0("Bulk_RNA_seq")) +
  theme(plot.title = element_text(hjust = 0.5, vjust = -2))

pdf("celltype_ABD/Promeristem.sanyuan.pdf",width = 10, height = 10)
print(p)
dev.off()
# length(tridat01$group)
write.csv(as.data.frame(table(tridat01$group)),"celltype_ABD/Promeristem.sanyuan.abd.csv")

as.data.frame(table(tridat01$group)/length(tridat01$group))
write.csv(as.data.frame(table(tridat01$group)/length(tridat01$group)),"celltype_ABD/Promeristem.sanyuan.abd.bili.csv")

#####################################

raw.expr <- read.delim("data/Protoderm_cells_scRNA.ABD.txt",
                       header = T)


raw.expr.1 <- mutate_all(raw.expr, ~replace(., is.na(.), 0))

# my_data <- mutate_at(my_data, c(`C1`, `C4`), ~replace(., is.na(.), 0))
# my_data[is.na(my_data)] <- 0
# raw.expr.2  <- mutate_all(raw.expr.1 ,  ~replace(., .<0, 0))

C0<-raw.expr.1

###1.draw the ternary diagram for the expression of triads in bulk cells of wheat root, and group the triads into 7 groups, i.e,
#Central (Balance), A.dominant, B.dominant, D.dominant, A.suppressed, B.suppressed, D.suppressed.
#load files
#C0<-read.csv("C0_abd_triads_R1.csv",header=T)
#count the ratio for each triad
cluster0<-cbind(C0[,2],C0[,3],C0[,4])
cluster0  <- round(cluster0,2)
# cluster0<- dplyr::mutate_all(cluster0,as.numeric)
head(cluster0)
sc0<-rowSums(cluster0)
# max0 <- rowMaxs(cluster0)
max0 <- rowSums(cluster0)
head(max0)
aai<-cluster0[,1]
bbi<-cluster0[,2]
ddi<-cluster0[,3]
clu0<-cbind(aai/sc0,bbi/sc0,ddi/sc0, max0)
head(clu0)
row.names(clu0)<-C0[,1]
colnames(clu0)<-c("AA","BB","DD", "size")
# clu0 <- mutate_all(clu0, ~replace(., is.na(.), 0))
# clu0 <- clu0[rowSums(clu0) > 0.5,]
clu0 <- clu0[clu0[,4] > 0.5,]
# clu0[is.na(clu0)] <- 0
# clu0 <- clu0[rowSums(clu0) > 0.5,]
head(clu0)
write.csv(clu0,"celltype_ABD/Protoderm_cells.sanyuan.expr.out.csv")
#Calculate the Euclidean distance
centers<-t(matrix(c(0.33,0.33,0.33,1,0,0,0,1,0,0,0,1,0,0.5,0.5,0.5,0,0.5,0.5,0.5,0), nrow=3))
colnames(centers)<-c("A","B","D")
rownames(centers)<-c("Central","A.dominant","B.dominant","D.dominant","A.suppressed","B.suppressed","D.suppressed")
head(centers)
tridat<-clu0[,1:3]
#rownames(tridat) <- rownames(tridat)
tridat <- data.frame(tridat)
tridat$group <- apply(tridat, 1, function(x){
  rownames(centers)[which.min(cdist(t(x), centers))]
})
tridat$size <- clu0[,4]
head(tridat)

tridat01<-matrix(unlist(tridat[rownames(tridat),c('AA','BB','DD')]), ncol=3)  %>%
  as.data.frame() %>% 
  magrittr::set_colnames(c('AA','BB','DD')) %>% 
  mutate(group=tridat$group, size=tridat$size) %>% 
  magrittr::set_rownames(rownames(tridat))
#tridat01<-na.omit(tridat01)
rowname<-rownames(tridat01)
head(rowname)
#tridat01 <- apply(tridat01,2,as.character)
write.csv(tridat01,"celltype_ABD/Protoderm_cells.tridat_cluster.csv")
# tridat01<-read.csv("C0_tridat_clusterR1.csv",header=T,row.names=1)
row.names(tridat01)<-rowname
dim(tridat01)
head(tridat01)
mypal1 = pal_lancet("lanonc", alpha = 0.7)(8)
mypal2 = pal_startrek("uniform", alpha=0.8)(7)
mypal3 = pal_npg("nrc", alpha=0.8)(10)
mypal4 = pal_nejm("default", alpha=0.8)(7)
mypal5 = pal_uchicago("default", alpha=0.8)(9)
mypal=c(mypal1, mypal2,mypal3,mypal4,mypal5)
show_col(mypal)
#require(pals)
library('pals')
cubebasis <- pal.compress(brewer.accent,n=8)
#cubebasis <- brewer.accent(8)
show_col(cubebasis)
coll<-c("#F781BF","#A65628","#4DAF4A","#984EA3","#A6CEE3","#FF7F00","#E41A1C","#377EB8")

#tridat01<-tridat01[which(tridat01[,4]!="character(0)"),]
#plot_data<-tridat01[which(tridat01[,4]=="Central"),]
#plot_data1<-tridat01[which(tridat01[,4]!="Central"),]

tridat01<-as.data.frame(tridat01)

tridat01$AA <- round(as.numeric(tridat01$AA),2)
tridat01$BB <- round(as.numeric(tridat01$BB),2)
tridat01$DD <- round(as.numeric(tridat01$DD),2)
# tridat01$size <- round((as.numeric(tridat01$size)),2)
tridat01$size <- round(log2(as.numeric(tridat01$size)),2)

p <- ggtern(as.data.frame(tridat01), mapping = aes(AA, BB, DD, size, group)) +
  geom_mask() +
  #geom_point(data=as.data.frame(tridat01), aes(size=size), color="grey", show.legend=F)  +
  #geom_point(data=plot_data, aes(size=size), color="grey", alpha=0.8, show.legend=F) +
  geom_point(data=as.data.frame(tridat01), aes(size=size, color=group), show.legend=T)  +
  scale_colour_manual(values=coll) +
  scale_size(range=c(0, 4)) +
  guides(size="none") +
  theme_bw() +
  theme(axis.text=element_blank(),
        axis.ticks=element_blank()) +
  ggtitle(paste0("Bulk_RNA_seq")) +
  theme(plot.title = element_text(hjust = 0.5, vjust = -2))

pdf("celltype_ABD/Protoderm_cells.sanyuan.pdf",width = 10, height = 10)
print(p)
dev.off()
# length(tridat01$group)
write.csv(as.data.frame(table(tridat01$group)),"celltype_ABD/Protoderm_cells.sanyuan.abd.csv")

as.data.frame(table(tridat01$group)/length(tridat01$group))
write.csv(as.data.frame(table(tridat01$group)/length(tridat01$group)),"celltype_ABD/Protoderm_cells.sanyuan.abd.bili.csv")

#####################################

raw.expr <- read.delim("data/Stem_cell_niche_scRNA.ABD.txt",
                       header = T)


raw.expr.1 <- mutate_all(raw.expr, ~replace(., is.na(.), 0))

# my_data <- mutate_at(my_data, c(`C1`, `C4`), ~replace(., is.na(.), 0))
# my_data[is.na(my_data)] <- 0
# raw.expr.2  <- mutate_all(raw.expr.1 ,  ~replace(., .<0, 0))

C0<-raw.expr.1

###1.draw the ternary diagram for the expression of triads in bulk cells of wheat root, and group the triads into 7 groups, i.e,
#Central (Balance), A.dominant, B.dominant, D.dominant, A.suppressed, B.suppressed, D.suppressed.
#load files
#C0<-read.csv("C0_abd_triads_R1.csv",header=T)
#count the ratio for each triad
cluster0<-cbind(C0[,2],C0[,3],C0[,4])
cluster0  <- round(cluster0,2)
# cluster0<- dplyr::mutate_all(cluster0,as.numeric)
head(cluster0)
sc0<-rowSums(cluster0)
# max0 <- rowMaxs(cluster0)
max0 <- rowSums(cluster0)
head(max0)
aai<-cluster0[,1]
bbi<-cluster0[,2]
ddi<-cluster0[,3]
clu0<-cbind(aai/sc0,bbi/sc0,ddi/sc0, max0)
head(clu0)
row.names(clu0)<-C0[,1]
colnames(clu0)<-c("AA","BB","DD", "size")
# clu0 <- mutate_all(clu0, ~replace(., is.na(.), 0))
# clu0 <- clu0[rowSums(clu0) > 0.5,]
clu0 <- clu0[clu0[,4] > 0.5,]
# clu0[is.na(clu0)] <- 0
# clu0 <- clu0[rowSums(clu0) > 0.5,]
head(clu0)
write.csv(clu0,"celltype_ABD/Stem_cell_niche.sanyuan.expr.out.csv")
#Calculate the Euclidean distance
centers<-t(matrix(c(0.33,0.33,0.33,1,0,0,0,1,0,0,0,1,0,0.5,0.5,0.5,0,0.5,0.5,0.5,0), nrow=3))
colnames(centers)<-c("A","B","D")
rownames(centers)<-c("Central","A.dominant","B.dominant","D.dominant","A.suppressed","B.suppressed","D.suppressed")
head(centers)
tridat<-clu0[,1:3]
#rownames(tridat) <- rownames(tridat)
tridat <- data.frame(tridat)
tridat$group <- apply(tridat, 1, function(x){
  rownames(centers)[which.min(cdist(t(x), centers))]
})
tridat$size <- clu0[,4]
head(tridat)

tridat01<-matrix(unlist(tridat[rownames(tridat),c('AA','BB','DD')]), ncol=3)  %>%
  as.data.frame() %>% 
  magrittr::set_colnames(c('AA','BB','DD')) %>% 
  mutate(group=tridat$group, size=tridat$size) %>% 
  magrittr::set_rownames(rownames(tridat))
#tridat01<-na.omit(tridat01)
rowname<-rownames(tridat01)
head(rowname)
#tridat01 <- apply(tridat01,2,as.character)
write.csv(tridat01,"celltype_ABD/Stem_cell_niche.tridat_cluster.csv")
# tridat01<-read.csv("C0_tridat_clusterR1.csv",header=T,row.names=1)
row.names(tridat01)<-rowname
dim(tridat01)
head(tridat01)
mypal1 = pal_lancet("lanonc", alpha = 0.7)(8)
mypal2 = pal_startrek("uniform", alpha=0.8)(7)
mypal3 = pal_npg("nrc", alpha=0.8)(10)
mypal4 = pal_nejm("default", alpha=0.8)(7)
mypal5 = pal_uchicago("default", alpha=0.8)(9)
mypal=c(mypal1, mypal2,mypal3,mypal4,mypal5)
show_col(mypal)
#require(pals)
library('pals')
cubebasis <- pal.compress(brewer.accent,n=8)
#cubebasis <- brewer.accent(8)
show_col(cubebasis)
coll<-c("#F781BF","#A65628","#4DAF4A","#984EA3","#A6CEE3","#FF7F00","#E41A1C","#377EB8")

#tridat01<-tridat01[which(tridat01[,4]!="character(0)"),]
#plot_data<-tridat01[which(tridat01[,4]=="Central"),]
#plot_data1<-tridat01[which(tridat01[,4]!="Central"),]

tridat01<-as.data.frame(tridat01)

tridat01$AA <- round(as.numeric(tridat01$AA),2)
tridat01$BB <- round(as.numeric(tridat01$BB),2)
tridat01$DD <- round(as.numeric(tridat01$DD),2)
# tridat01$size <- round((as.numeric(tridat01$size)),2)
tridat01$size <- round(log2(as.numeric(tridat01$size)),2)

p <- ggtern(as.data.frame(tridat01), mapping = aes(AA, BB, DD, size, group)) +
  geom_mask() +
  #geom_point(data=as.data.frame(tridat01), aes(size=size), color="grey", show.legend=F)  +
  #geom_point(data=plot_data, aes(size=size), color="grey", alpha=0.8, show.legend=F) +
  geom_point(data=as.data.frame(tridat01), aes(size=size, color=group), show.legend=T)  +
  scale_colour_manual(values=coll) +
  scale_size(range=c(0, 4)) +
  guides(size="none") +
  theme_bw() +
  theme(axis.text=element_blank(),
        axis.ticks=element_blank()) +
  ggtitle(paste0("Bulk_RNA_seq")) +
  theme(plot.title = element_text(hjust = 0.5, vjust = -2))

pdf("celltype_ABD/Stem_cell_niche.sanyuan.pdf",width = 10, height = 10)
print(p)
dev.off()
# length(tridat01$group)
write.csv(as.data.frame(table(tridat01$group)),"celltype_ABD/Stem_cell_niche.sanyuan.abd.csv")

as.data.frame(table(tridat01$group)/length(tridat01$group))
write.csv(as.data.frame(table(tridat01$group)/length(tridat01$group)),"celltype_ABD/Stem_cell_niche.sanyuan.abd.bili.csv")


##################################结果数据统计#############################

###################
rm(list=ls())

cs.3 <- read.csv("celltype_ABD/Promeristem.tridat_cluster.csv")
f1.3 <- read.csv("celltype_ABD/Protoxylem.tridat_cluster.csv" )

cs.samyuan <- cs.3[,c(1,5)]
colnames(cs.samyuan) <- c("group","cs")

f1.samyuan <-  f1.3[,c(1,5)]
colnames(f1.samyuan) <- c("group","f1")


cs.f1 <- dplyr::left_join(cs.samyuan,f1.samyuan  )
write.csv(table(cs.f1[-1]),"dif_scRNA_ABD/Promeristem.Protoxylem.bianhua.统计表格.csv")
biaoge <- read.csv("dif_scRNA_ABD/Promeristem.Protoxylem.bianhua.统计表格.csv",row.names = 1,header = T)
file.remove("dif_scRNA_ABD/Promeristem.Protoxylem.bianhua.统计表格.csv")
biaoge$bubian <- apply(biaoge,2,max)
biaoge$Sum <- apply(biaoge[,c(1:7)],2,sum)
biaoge$binahua <- biaoge$Sum  - biaoge$bubian
biaoge$binahua.bili <- biaoge$binahua/biaoge$Sum 
biaoge$jizubianhua  <- rep(0,length(biaoge$binahua.bili))
biaoge$jizubianhua[2] <- sum(biaoge$binahua.bili[1]  +  biaoge$binahua.bili[2]) *100
biaoge$jizubianhua[4] <- sum(biaoge$binahua.bili[3]  +  biaoge$binahua.bili[4]) *100
biaoge$jizubianhua[7] <- sum(biaoge$binahua.bili[6]  +  biaoge$binahua.bili[7]) *100

write.csv(biaoge,"dif_scRNA_ABD/Promeristem.Protoxylem.bianhua.统计表格.比例.csv")


###################
rm(list=ls())

cs.3 <- read.csv("celltype_ABD/Promeristem.tridat_cluster.csv")
f1.3 <- read.csv("celltype_ABD/Protophloem.tridat_cluster.csv" )

cs.samyuan <- cs.3[,c(1,5)]
colnames(cs.samyuan) <- c("group","cs")

f1.samyuan <-  f1.3[,c(1,5)]
colnames(f1.samyuan) <- c("group","f1")


cs.f1 <- dplyr::left_join(cs.samyuan,f1.samyuan  )
write.csv(table(cs.f1[-1]),"dif_scRNA_ABD/Promeristem.Protophloem.bianhua.统计表格.csv")
biaoge <- read.csv("dif_scRNA_ABD/Promeristem.Protophloem.bianhua.统计表格.csv",row.names = 1,header = T)
file.remove("dif_scRNA_ABD/Promeristem.Protophloem.bianhua.统计表格.csv")
biaoge$bubian <- apply(biaoge,2,max)
biaoge$Sum <- apply(biaoge[,c(1:7)],2,sum)
biaoge$binahua <- biaoge$Sum  - biaoge$bubian
biaoge$binahua.bili <- biaoge$binahua/biaoge$Sum 
biaoge$jizubianhua  <- rep(0,length(biaoge$binahua.bili))
biaoge$jizubianhua[2] <- sum(biaoge$binahua.bili[1]  +  biaoge$binahua.bili[2]) *100
biaoge$jizubianhua[4] <- sum(biaoge$binahua.bili[3]  +  biaoge$binahua.bili[4]) *100
biaoge$jizubianhua[7] <- sum(biaoge$binahua.bili[6]  +  biaoge$binahua.bili[7]) *100

write.csv(biaoge,"dif_scRNA_ABD/Promeristem.Protophloem.bianhua.统计表格.比例.csv")



