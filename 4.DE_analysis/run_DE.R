###################################################考虑已知的批次因素进行差异基因分析

rm(list=ls()) 

suppressMessages(library(DESeq2))
suppressMessages(library("RColorBrewer"))
suppressMessages(library("gplots"))
suppressMessages(library("amap"))
suppressMessages(library("ggplot2"))
suppressMessages(library("BiocParallel"))
suppressMessages(library("YSX"))
suppressMessages(library(sva))
suppressMessages(library(ggfortify))
suppressMessages(library(patchwork))
suppressMessages(library(ggbeeswarm))
suppressMessages(library(limma))


# Prefix for all output file
output_prefix = "spike.simpler"

# pipelineStar.sh生成的reads count 文件，行为基因，列为样品
file = "~/work_r/spike/data/mRNA.spike.counts.selected.matrix"

# 分组信息表
sampleFile = "~/work_r/spike/data/sampleFile.selected.txt"

# 分组信息所在列名字
covariate = NULL
# covariate = "batch"
design="conditions"

# 输入数据类型，salmon结果或reads count 矩阵
type="readscount"

# 差异基因参数
padj=0.05
log2FC=1


##数据读入和标准化

dds <- readscount2deseq(file, sampleFile, design=design, covariate = covariate)


normexpr <- deseq2normalizedExpr(dds, output_prefix=output_prefix)



## 检查数据标准化效果: 标准化后基因在不同样品的表达分布越均一越好。从下图看不出存在批次效应的影响。

# normalizedExpr2DistribBoxplot(normexpr,
#   saveplot=paste(output_prefix, "DESeq2.normalizedExprDistrib.pdf", sep="."))
# normalizedExpr2DistribBoxplot(normexpr)


### 样本聚类查看样品相似性，trt组和untrt组区分明显 (聚类采用的不同基因数目、聚类参数都可能影响聚类结果)

# clusterSampleHeatmap2(normexpr$rlog,
#                       cor_file=paste(output_prefix, "DESeq2.sampleCorrelation.txt", sep="."),
#                       saveplot=paste(output_prefix, "DESeq2.sampleCorrelation.pdf", sep="."))
# 根据前5000个表达变化幅度最大的基因进行聚类分析
# clusterSampleHeatmap2(normexpr$rlog[1:5000,], cor_file=paste(output_prefix, "DESeq2.sampleCorrelation.txt", sep="."))

## [1] "Performing sample clustering"

# clusterSampleUpperTriPlot(normexpr$rlog[1:5000,], cor_file=paste(output_prefix, "DESeq2.sampleCorrelation.txt", sep="."))

## [1] "Performing sample clustering"


### 主成分分析PCA查看样品相似性，发现在PC1轴上，样品按处理条件区分开；在PC2轴上，样品按个体区分开，这个结果也与我们前面不考虑批次因素的结果是一样的。
#是不是批次变量加错了呢，还是添加的批次变量未生效？可以说都不是，操作没问题，只是DESeq2处理时只在差异分析模型中考虑批次效应信息，而不会直接校正表达矩阵。

#  metadata = as.data.frame(colData(dds))
#  sp_pca(normexpr$rlog[1:5000,], metadata, color_variable="conditions", shape_variable = "individual") + aes(size=1) + guides(size = F)


###鉴定出差异基因，获得差异基因文件spike.simpler.batch.DESeq2.all.DE和其它可视化图表

multipleGroupDEgenes(dds, design=design, output_prefix=output_prefix, padj=padj, log2FC=log2FC)

