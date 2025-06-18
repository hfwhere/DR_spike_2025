rm(list=ls())

## 从 AnnotionHub下载 OrgDb

#BiocManager::install("AnnotationHub")

#library(AnnotationHub)
#ah <- AnnotationHub()

#查询支持的物种
#ah_orgdb = ah[ah$rdataclass=="OrgDb"]
#query(ah_orgdb, "Triticum aestivum")


## 使用AnnotationForge构建

#Rscript /pub/software/emcp/emapperx.R out.emapper.annotations proteins.fa

#install.packages('org.My.eg.db_1.0.tar.gz', 
#                 repos = NULL, #从本地安装
#                 lib = 'R_Library') # 安装文件夹
#


library(readr)
library(org.My.eg.db, lib = "R_Library/")
library(GOplot)
library(ggplot2)
library(enrichplot)
library(tidyverse)

## 准备基因列表

deg_result <- read_delim("data/de_result", 
                        delim = "\t", escape_double = FALSE, 
                        trim_ws = TRUE)


my_deg_result <- filter(deg_result, 
                        sampleA == "DR",
                        sampleB == "CK") %>%
  arrange(desc(abs(log2FoldChange)))



gene <- filter(my_deg_result, 
               direction != "NS") %>%
  pull(gene_id)



geneList <- my_deg_result$log2FoldChange
names(geneList) <- my_deg_result$gene_id
geneList <- sort(geneList, decreasing = T)


#=========================GO 富集及可视化

##=========GO 分类

library(clusterProfiler)

#ggo <- groupGO(gene     = gene,
#               OrgDb    = org.My.eg.db,
#               keyType  = "GID",
#               ont      = "BP",
#               level    = 3,
#               readable = FALSE)
#head(ggo)

## =========GO over-representation 富集

ego <- enrichGO(gene          = gene,
                #universe      = names(geneList),
                OrgDb         = org.My.eg.db,
                keyType       = "GID",
                ont           = "ALL", #MF | CC |BP
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = FALSE #true
                )

###============ 富集结果优化

# 去除部分 term

# 去除 level 1,2,3
#go <- dropGO(ego, level = 1:3) %>%
  # 去除包含 drug 的GO
# filter(!str_detect(Description, "drug"))

#去除冗余 term

ego <- enrichplot::pairwise_termsim(ego)
ego <- clusterProfiler::simplify(ego, cutoff=0.7, 
                                 by="p.adjust", 
                                 select_fun=min)

ego_df <- as.data.frame(ego)

write.table(ego_df,
            file = 'out/enrich_GO_DR_vs_CK.txt', 
            sep = '\t', 
            quote = F, 
            row.names = F)


###=========== 可视化

p2 = barplot(ego, showCategory = 15)
pdf(file = 'fig0/barplot_GO_DR_vs_CK.pdf')
p2
dev.off()


p3 = dotplot(ego, showCategory = 15)
pdf(file = 'fig0/dotplot_GO_DR_vs_CK.pdf')
p3
dev.off()


p4 = cnetplot(ego, foldChange=geneList)
pdf(file = 'fig0/cnetplot_GO_DR_vs_CK.pdf')
p4
dev.off()


## 
ekpWY <- enrichGO(gene          = gene[1:2000],
                  #universe      = names(geneList),
                  OrgDb         = org.My.eg.db,
                  keyType       = "GID",
                  ont           = "BP", #MF | CC |BP
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.05,
                  readable      = FALSE #true
)

p1 = goplot(ekpWY)
pdf(file = 'fig0/goplot_GO_DR_vs_CK.pdf')
p1
dev.off()

