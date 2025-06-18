#############################################################
#################普通小麦幼穗发育转录组分析###################
#############################################################



########################准备参考基因组################
##以中国春为参考基因组
#pwd:~/spike/data/genome

#基因组序列文件：序列内容改成短行，60bp，已是短的，并改名为genome.fa
#构建index，已有genome.fa.fai


#基因注释文件：gff3改为gtf格式
gffread -T -o genes_gffread.gtf Triticum_aestivum.IWGSC.55.gff3



#mRNA分析用数据处理
##提取最长转录本
###是否注释了可变剪接
gtftk count -i genes.gtf

###筛选最长转录本
gtftk short_long -l -i genes.gtf | gtftk select_by_key -k feature -v transcript | gtftk tabulate -k transcript_id,gene_id > longest_mapid.txt

###提取最长转录本对应的蛋白序列
sed '1d' longest_mapid.txt | cut -f 1 | seqtk subseq Triticum_aestivum.IWGSC.pep.all.fa -> longest_transcript.pep.fa
sed '1d' longest_mapid.txt | cut -f 1 | seqtk subseq Triticum_aestivum.IWGSC.cds.all.fa -> longest_transcript.cds.fa


#LncRNA分析用数据处理

## 已有基因组注释, mRNA
awk '{print $0" transcript_biotype \"protein_coding\";"}' IWGSC_v1.1_HC_20170706.gtf > known_mRNA.gtf


########################准备测序数据################
#pwd:~/spike/data/mRNA or LncRNA
#下载数据
nohup sh ../run_kf_mrna.sh 1>kf.log 2>&1 &
nohup sh ../run_kf_lncrna.sh 1>kf.log 2>&1 &

#.sra转为fq格式
nohup sh ../run_fq_mrna.sh 1>fq.log 2>&1 &
nohup sh ../run_fq_lncrna.sh 1>fq.log 2>&1 &

#测序数据的过滤质控，以下面的为例，相应修改文件名
fastp -i ../KN_DR_R1_1.fastq -I ../KN_DR_R1_2.fastq -o KN_DR_R1_1.fastq -O KN_DR_R1_2.fastq -h KN_DR_R1.html -j KN_DR_R1.json

fastp -i ../CS_AM_R1.fastq -o CS_AM_R1.fastq -h CS_AM_R1.html -j CS_AM_R1.json



######################比对###############################

#为参考序列构建index
hisat2-build -p 20 \ # 线程数，记得修改
../ref/genome.fa \ # 输入：基因组序列
../ref/genome \ # 输出
1>hisat2-build.log 2>&1


#Hisat2 比对
hisat2 --new-summary --dta \
-p 20 \
-x ../ref/genome \ # 上部构建好的参考序列
-1 /home/zhxd2/workspace/lncRNA/14.data/clean/SRR10556684_1.fastq \ # read1
-2 /home/zhxd2/workspace/lncRNA/14.data/clean/SRR10556684_2.fastq \ # read2
-S SRR10556684.sam \
--rna-strandness RF \
1>SRR10556684.log 2>&1 


#比对结果压缩和排序
samtools sort \
-o SRR10556688.bam \ # 输出的bam 文件
SRR10556688.sam # 输入的sam 文件


#bam构建index
samtools index SRR10556684.bam


#比对率统计,注意修改文件名
Rscript step5.statistics.R


######################LncRNA转录本重构###############################
#单样本转录本重构，pwd:~/spike/data/LncRNA
singularity exec ~/spike/software/lnctools.sif stringtie -p 36 --rf -G ../genome/known_mRNA.gtf -o gtf/ZM_AM_R1.gtf bam/ZM_AM_R1.bam 1>ZM_AM_R1.log 2>&1 &


#重构转录本合并
singularity exec ~/spike/software/lnctools.sif stringtie --merge -o gtf/merged.gtf -G ../genome/known_mRNA.gtf gtf/*.gtf


######################LncRNA鉴定###############################
#在genek服务器, pwd:~/spike/data/LncRNA

#过滤
singularity exec ~/spike/software/lnctools.sif \
FEELnc_filter.pl -i gtf/merged.gtf -a ../genome/known_mRNA.gtf --monoex=-1 -s 200 -f 0 -p 40 > candidate_lncRNA.gtf 2>FEELnc_filter.log

# 提取CDS
singularity exec ~/spike/software/lnctools.sif \
gffread -w candidate_lncRNA.fa -g ../genome/genome.fa candidate_lncRNA.gtf

# 提取候选lncRNA ID
singularity exec ~/spike/software/lnctools.sif \
gtftk get_attr_value_list -i candidate_lncRNA.gtf -k transcript_id -o candidate_lncRNA.txt


# 编码能力预测

## 基于算法
### CPC2
singularity exec ~/spike/software/lnctools.sif \
CPC2.py -i candidate_lncRNA.fa -o cpc2output
### PLEK
singularity exec ~/spike/software/lnctools.sif \
python /opt/PLEK.1.2/PLEK.py -fasta candidate_lncRNA.fa -out plekoutput.txt -thread 36
### CNCI
singularity exec ~/spike/software/lnctools.sif \
CNCI.py -f candidate_lncRNA.fa -o cnci -m pl -p 6

##基于数据库比对
### Blast，uniprot下载数据，5.28：00
singularity exec ~/spike/software/lnctools.sif \
diamond makedb --in uniprot_sprot.fasta -d uniprot_sprot

singularity exec ~/spike/software/lnctools.sif \
diamond blastx -d uniprot_sprot -q candidate_lncRNA.fa -o diamondoutput.txt


#求交集
Rscript Lnc_intersect.R


#提取序列和GTF
singularity exec ~/spike/software/lnctools.sif \
gtftk select_by_key -i candidate_lncRNA.gtf -k transcript_id -f lncRNA.txt -o lncRNA.gtf

singularity exec ~/spike/software/lnctools.sif \
gffread -w lncRNA.fa -g ../genome/genome.fa lncRNA.gtf

awk '{print $1}' ../lncRNA.fa >lncRNA.fa


#分类
singularity exec ~/spike/software/lnctools.sif \
FEELnc_classifier.pl -i lncRNA.gtf -a ../genome/known_mRNA.gtf > lncRNA_classes.txt 2>FEELnc_classifier.log



######################定量###############################

#######################LncRNA表达定量, mRNA和LncRNA一起定量

##合并结果
awk '{print $0" transcript_biotype \"lncRNA\";"}' ../lncRNA.gtf >lncRNA.gtf
cat ../genome/known_mRNA.gtf lncRNA.gtf > genes.gtf

## LncRNA表达定量
##singularity exec ~/spike/software/lnctools.sif run-featurecounts.R -a transcript_id --strandSpecific 2 -b ../bam/ZM_AM_R1.bam -g genes.gtf -o ZM_AM_R1 1>ZM_AM_R1.log 2>&1 &
bash run_Quantification.sh

## 合并成表达矩阵并进行TMM 标准化
## pwd:~/spike/data/LncRNA/Quantification
ls *.count > quant_files.txt

singularity exec ~/spike/software/lnctools.sif \
abundance_estimates_to_matrix.pl --est_method featureCounts --quant_files quant_files.txt --out_prefix

## 结果转移到 ~/spike/data/LncRNA/matrix/
cd ../matrix/
mv ../Quantification/lncRNA.* .


####################### mRNA表达定量

##singularity exec ~/spike/software/lnctools.sif run-featurecounts.R -a transcript_id -b spike_mRNA/KN_DR_R1.bam -g genes.gtf -o Quantification/KN_DR_R1 1>KN_DR_R1.log 2>&1 &
##singularity exec ~/spike/software/lnctools.sif run-featurecounts.R -a transcript_id -i false -b spike_mRNA/CS_AM_R1.bam -g genes.gtf -o Quantification/CS_AM_R1 1>CS_AM_R1.log 2>&1 &
bash run_Quantification.sh

##合并叶片转录组数据
##pwd: ~/spike/data/mRNA/leaf_mRNA

# samtools index -c LE_DR_R1.bam
#samtools merge LE_DR_PP.bam LE_DR_R1.bam LE_DR_R2.bam LE_DR_R3.bam LE_PP_R1.bam LE_PP_R2.bam LE_PP_R3.bam

samtools merge LE_DR_PP_R1.bam LE_DR_R1.bam LE_PP_R1.bam
samtools merge LE_DR_PP_R2.bam LE_DR_R2.bam LE_PP_R2.bam
samtools merge LE_DR_PP_R3.bam LE_DR_R3.bam LE_PP_R3.bam

## 叶片mRNA表达定量
## pwd:~/spike/data/mRNA
#singularity exec ~/spike/software/lnctools.sif \
#run-featurecounts.R -a transcript_id -b leaf_mRNA/LE_DR_PP.bam -g genes.gtf -o Quantification/LE_DR_PP 1>Quantification/LE_DR_PP.log 2>&1 &

singularity exec ~/spike/software/lnctools.sif \
run-featurecounts.R -a transcript_id -b leaf_mRNA/LE_DR_PP_R1.bam -g genes.gtf -o Quantification/LE_DR_PP_R1 1>Quantification/LE_DR_PP_R1.log 2>&1 &

singularity exec ~/spike/software/lnctools.sif \
run-featurecounts.R -a transcript_id -b leaf_mRNA/LE_DR_PP_R2.bam -g genes.gtf -o Quantification/LE_DR_PP_R2 1>Quantification/LE_DR_PP_R2.log 2>&1 &

singularity exec ~/spike/software/lnctools.sif \
run-featurecounts.R -a transcript_id -b leaf_mRNA/LE_DR_PP_R3.bam -g genes.gtf -o Quantification/LE_DR_PP_R3 1>Quantification/LE_DR_PP_R3.log 2>&1 &



## 合并成表达矩阵并进行TMM 标准化
## pwd:~/spike/data/mRNA/Quantification
### 将LncRNA count导入到此文件夹

ls *.count > quant_files.txt

singularity exec ~/spike/software/lnctools.sif \
abundance_estimates_to_matrix.pl --est_method featureCounts --quant_files quant_files.txt --out_prefix

### 结果转移到 ~/spike/data/LncRNA/matrix/
cd ../matrix/
mv ../Quantification/mRNA.* .


###################### 差异表达 ###############################

###################### 获取可用数据
## LncRNA数据拆分
Rscript split_lncRNA.R


# mRNA 数据拆分
Rscript split_mRNA.R


###################### 去除批次效应 BatchEffect
####H:\Doing\spike\3_sample_dif
#运行batch.R

## 利用ComBat和limma矫正批次效应后，保存为spike.mRNA.counts.Batch.ComBat.matrix和spike.mRNA.counts.Batch.limma.matrix


######################样品差异分析###############################

# pwd:~/spike/data/mRNA/4.Quantification_Merge/sample_dif
# 运行PCA_sample.R


######################差异表达###############################
#手动创建contrasts.txt，samples.txt
# nohup sh run_DE.sh 1>DE.log 2>&1 &

nohup Rscript run_batched_DE.R 1>DE.log 2>&1 & 



###################### 功能注释及富集 ###############################
# pwd:~/spike/data/mRNA/6.enrich

#注释，在线注释eggNOGmapper

Rscript /pub/software/emcp/emapperx.R out.emapper.annotations proteins.fa

######################WGCNA分析###############################
## pwd: spike\7.WGCNA
nohup Rscript spike_wgcna.R 1>WGCNA.log 2>&1 &

































###后面暂未做
#筛选差异基因，|log2FoldChange| > 1 and padj < 0.05

sed '1d' gene.counts_alt.matrix.Y1805T1_vs_CSCK1.DESeq2.DE_results | awk 'sqrt($5*$5)>1 && $9<0.05 {print $1"\t"$5"\t"$9}' | sort -k 2n > DES_Y1805T1_vs_CSCK1.txt

sed '1d' gene.counts_alt.matrix.Y1805T1_vs_CSCK1.DESeq2.DE_results | awk 'sqrt($5*$5)>1 && $9<0.05 {print $1"\t"$5"\t"$9}' | sort -k 2n > DES_Y1805T1_vs_CSCK1.txt

sed '1d' gene.counts_alt.matrix.Y1805T1_vs_CST1.DESeq2.DE_results | awk 'sqrt($5*$5)>1 && $9<0.05 {print $1"\t"$5"\t"$9}' | sort -k 2n > DES_Y1805T1_vs_CST1.txt

sed '1d' gene.counts_alt.matrix.Y1805T1_vs_Y1805CK1.DESeq2.DE_results | awk 'sqrt($5*$5)>1 && $9<0.05 {print $1"\t"$5"\t"$9}' | sort -k 2n > DES_Y1805T1_vs_Y1805CK1.txt

sed '1d' gene.counts_alt.matrix.Y1805T1_vs_Y1805T2.DESeq2.DE_results | awk 'sqrt($5*$5)>1 && $9<0.05 {print $1"\t"$5"\t"$9}' | sort -k 2n > DES_Y1805T1_vs_Y1805T2.txt

sed '1d' gene.counts_alt.matrix.Y1805T2_vs_CSCK2.DESeq2.DE_results | awk 'sqrt($5*$5)>1 && $9<0.05 {print $1"\t"$5"\t"$9}' | sort -k 2n > DES_Y1805T2_vs_CSCK2.txt

sed '1d' gene.counts_alt.matrix.Y1805T2_vs_CST2.DESeq2.DE_results | awk 'sqrt($5*$5)>1 && $9<0.05 {print $1"\t"$5"\t"$9}' | sort -k 2n > DES_Y1805T2_vs_CST2.txt

sed '1d' gene.counts_alt.matrix.Y1805T2_vs_Y1805CK2.DESeq2.DE_results | awk 'sqrt($5*$5)>1 && $9<0.05 {print $1"\t"$5"\t"$9}' | sort -k 2n > DES_Y1805T2_vs_Y1805CK2.txt



















######################功能注释###############################
#cs_ee pep序列基因名改到和DNA序列一致，按每个文件99999个序列切割成分文件，在网站http://eggnog5.embl.de/#/app/home 进行功能注释

#手动将除第一个文件外的其余文件抬头和尾部#行去除，合并所有文件为一个文件

cd /mnt/f/data_results/mRNA_8/5_annotation

cat query_seqs.fa.emapper1.annotations query_seqs.fa.emapper2_renamed.annotations query_seqs.fa.emapper3_renamed.annotations query_seqs.fa.emapper4_renamed.annotations > query_seqs.fa.emapper_all.annotations


######################引物开发###############################
#qPCR用CDS开发引物，用cDNA或mRNA全序列验证
nohup sh qPCR.sh 1>qPCR.log 2>&1 &

#qPCR.shS
eprimer32 -sequence selected_CDS.fa -outfile qPCR.primer \
-optsize 20 -numreturn 3 \
-minsize 15 -maxsize 25 \
-opttm 50 -mintm 45 -maxtm 55 \
-psizeopt 100 -prange 80-130

awk '{if($0~/EPRIMER32/) {seq_name=$5;count=1;} else \
if($0~/FORWARD PRIMER/) forward=$7; else if ($0~/REVERSE PRIMER/) \
{reverse=$7; printf("%s@%d\t%s\t%s\n", seq_name,count,forward, reverse); \
count+=1;} }' qPCR.primer > qPCR_primer_file

primersearch -seqall cs_and_ee_cdna.fa -infile qPCR_primer_file -mismatchpercent 5 -outfile qPCR.database.primerSearch

#PCR用gene或mRNA开发引物，用基因组序列验证
nohup sh PCR.sh 1>PCR.log 2>&1 &


#PCR.sh
eprimer32 -sequence selected_gene.fa -outfile PCR.primer \
-optsize 20 -numreturn 3 \
-minsize 15 -maxsize 25 \
-opttm 50 -mintm 45 -maxtm 55 \
-psizeopt 0 -prange 100-500

awk '{if($0~/EPRIMER32/) {seq_name=$5;count=1;} else \
if($0~/FORWARD PRIMER/) forward=$7; else if ($0~/REVERSE PRIMER/) \
{reverse=$7; printf("%s@%d\t%s\t%s\n", seq_name,count,forward, reverse); \
count+=1;} }' PCR.primer >PCR_primer_file


primersearch -seqall cs_and_ee_genome.fa -infile PCR_primer_file -mismatchpercent 5 -outfile PCR.database.primerSearch




