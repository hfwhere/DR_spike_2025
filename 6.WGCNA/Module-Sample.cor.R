
rm(list=ls())


## ----模块与表型关联---------------------------------------
library(corrplot) 
library(vegan) 
library(ggcor)
library(readr)
library(tidyverse)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(WGCNA)

# 导入数据
load(file = 'data/WGCNA.RData')
datTraits <- read.table('data/pheno.txt',
                        header = T, row.names = 1)

# 计算相关性
moduleTraitCor <- cor(
  net$MEs,
  datTraits,
  use = "p",
  method = 'spearman'#'pearson' 
)

# 计算 Pvalue
moduleTraitPvalue <- corPvalueStudent(
  moduleTraitCor, 
  nrow(datExpr))



###模块和表型相关性COR图==========
# 转置

moduleTraitCor0 <- t(moduleTraitCor)

moduleTraitPvalue0 <- t(moduleTraitPvalue)


pdf("Module-Sample.cor.pdf", width = 18, height = 10)


link_cor <- correlate(net$MEs, datTraits, cor.test = T) %>%
  as_cor_tbl() %>%
  select(.col.names, .row.names, r, p.value) %>%
  mutate(
    r = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
            labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
    p.value = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))


quickcor(moduleTraitCor0, type = "upper") +
  geom_square() +
  scale_fill_gradient2(low = '#00A087FF', high = '#DC0000FF',
                       mid = '#F39B7FFF', midpoint = 0) +
  ggcor::anno_link(data = link_cor, aes(color = p.value, size = r)) +
  scale_size_manual(values = c(0.5, 1, 2)) +
  scale_colour_manual(values = c("#D95F02", "#1B9E77", "#A2A2A288")) 


dev.off()











