###get target gene id for express gene###
setwd("E:/liftoff")
gtf = rtracklayer::import("Gambusia_affinis.ASM309773v1.109.gtf")
gtf_df=as.data.frame(gtf)

setwd("E:/miRNA")
target=read.csv('miRNA_target_genes_hitsOnly.txt',header=F,sep="\t")
# 加载 tidyr 包
library(tidyr)

# 将 V1 列拆分成三列
target <- target %>%
  separate(V1, into = c("seqnames", "start", "end"), sep = "[:-]", remove = FALSE)

# 使用 merge() 函数按 seqnames、start 和 end 列进行合并
matched_gtf_df <- merge(target, gtf_df, by = c("seqnames", "start", "end"))
# 筛选出 type 列值为 "three_prime_utr" 的行
matched_gtf_df <- matched_gtf_df[matched_gtf_df$type == "three_prime_utr", ]
# 提取特定列并创建新的数据框
final_df <- matched_gtf_df[, c("seqnames", "start", "end", "V2", "gene_id")]
#去除重复
final_df_1 <- final_df[!duplicated(final_df[c("V2", "gene_id")]), ]

library(dplyr)
target_uni <- final_df_1 %>% select(-1:-3)
#target_uni <- final_df_1 %>% dplyr::select(-(1:3))
names(target_uni)[1] <- "miRNA"       # 修改第一列的列名
names(target_uni)[2] <- "Gene_id"     # 修改第二列的列名


###pearson analysis for express gene###

#miRNA和mRNA

#rlog miRNA count
library(ggplot2)
library(pheatmap)
library(reshape2)
library(DESeq2)
setwd("E:/miRNA")
miRNA_count=read.table("filter_miRNA_table.txt",header = T,check.names = F,row.names =1)
sampleTable=read.csv("Sample.csv",header=T,sep = ",",check.names = F)
miRNA_count=DESeqDataSetFromMatrix(countData=miRNA_count, colData=sampleTable, design = ~Sex+Salinity+Time)
miRNA_count$Sex=relevel(miRNA_count$Sex, ref ="Female")
miRNA_count$Sex
miRNA_count=rlog(miRNA_count,blind=F)
miRNA_count=miRNA_count@assays@data@listData[[1]]
miRNA_count=as.data.frame(miRNA_count)

Female_target_count=miRNA_count[,grep("F", colnames(miRNA_count))]
Male_target_count=miRNA_count[,grep("M", colnames(miRNA_count))]


#rlog mRNA count
setwd("E:/filter")
allcount=read.csv('all_sample.txt',header=F,sep="\t")
colnames(allcount)=c('sample','gene.ID','count')
allcount=dcast(allcount,formula=gene.ID~sample,value.var="count")#gene.ID as row-,sample as col-,count as the variate
rownames(allcount)=allcount$gene.ID
noint = rownames(allcount) %in% c("__alignment_not_unique","__ambiguous","__no_feature","__not_aligned","__too_low_aQual")
head(noint,6)
count_table=allcount[!noint,]
count_table=count_table[,-1]
rm=grep("G",colnames(count_table))
count_table=count_table[,-rm]
#count_filter: gene >1 at all samples
count_table=subset(count_table,apply(count_table,1, function(x) sum(x >= 1))>= 32)
sampleTable=read.csv("Sample.csv",header=T,sep = ",",check.names = F)
rownames(sampleTable)=sampleTable$Sample
colnames(count_table)==rownames(sampleTable) #check if the count table has the same order as sample table
str(sampleTable)
id=colnames(count_table)
id=data.frame(id)
sampleTable=sampleTable[id$id,]
colnames(count_table)==rownames(sampleTable)

rlog_count_table=DESeqDataSetFromMatrix(countData=count_table, colData=sampleTable, design = ~Sex+Salinity+Time)
rlog_count_table$Sex=relevel(rlog_count_table$Sex, ref ="Female")
rlog_count_table$Sex
rlog_count_table=rlog(rlog_count_table,blind=F)
rlog_count_table=rlog_count_table@assays@data@listData[[1]]
rlog_count_table=as.data.frame(rlog_count_table)

rlog_count_table_Female=rlog_count_table[,grep("F", colnames(rlog_count_table))]
rlog_count_table_Male=rlog_count_table[,grep("M", colnames(rlog_count_table))]


###1 Female counts
# 获取 rlog_count_table_Female 的列名
new_colnames <- colnames(rlog_count_table_Female)

# 将这些列名赋值给 Female_target_count
colnames(Female_target_count) <- new_colnames
list_female=unique(target_uni$miRNA)
female<-list()
for( i in 1:length(list_female)){
  t3<-list_female[i]
  tmp<-subset(target_uni,miRNA==t3)
  female[[t3]]<-tmp
}

list_female=names(female)
female_p<-list()

for( i in 1:length(list_female)){
  t1<-list_female[i]
  #miRNA
  deLNC <- t1
  #mRNA
  dePC <- female[[t1]][["Gene_id"]]
  #rbind miRNA and mRNA count table
  miRNA=Female_target_count[deLNC,]
  mRNA=rlog_count_table_Female[dePC,]
  mRNA=mRNA %>% drop_na()
  rnaExpr <- rbind(miRNA,mRNA)
  dePC=rownames(mRNA)
  
  combination <- expand.grid(deLNC, dePC)
  names(combination)=c("miR","mR")

  
  #the relationship of miRNA-mRNA pairs
  cor_result=apply(combination,1,function(x){
    miR=as.character(x[1])
    mR=as.character(x[2])
    result=cor.test(as.numeric(rnaExpr[miR,]), as.numeric(rnaExpr[mR,]))
    score=c(pval=result$p.value,result$estimate)
    return(score)
  })
  
  result=cbind(combination,t(cor_result))
  tmp=result %>% drop_na()
  female_p[[t1]]<-tmp
}

library(data.table)
female_pearson=rbindlist(female_p, use.names=TRUE)
# 筛选条件：pval < 0.05 且 cor < 0
filtered_female_pearson <- female_pearson[female_pearson$pval < 0.05 & female_pearson$cor < 0, ]

library(dplyr)
#1.1 Female parallel
setwd("E:/miRNA/express")
all_pa_gene=read.csv("pa_express_gene.csv",header=T,sep = ",",check.names = F)
#相关（受到miRNA调控的parallel gene数量）
cor_pa_gene <- filtered_female_pearson %>%
  filter(mR %in% all_pa_gene$x)#56
write.csv(cor_pa_gene,"./cor_pa_gene.csv",row.names = F)
#相关（靶基因是parallel gene的miRNA数量）
cor_pa_miR <- data.frame(miR = unique(cor_pa_gene$miR), stringsAsFactors = FALSE)#40

#1.2 Female non-parallel
setwd("E:/miRNA/express")
all_nopa_gene=read.csv("non_pa_express_gene.csv",header=T,sep = ",",check.names = F)
#相关（受到miRNA调控的non parallel gene数量）
cor_nopa_gene <- filtered_female_pearson %>%
  filter(mR %in% all_nopa_gene$value)#749
#相关（靶基因是non parallel gene的miRNA数量）
cor_nopa_miR <- data.frame(miR = unique(cor_nopa_gene$miR), stringsAsFactors = FALSE)#130


###2 Male counts
# 获取 rlog_count_table_Male 的列名
new_colnames <- colnames(rlog_count_table_Male)

# 将这些列名赋值给 Male_target_count
colnames(Male_target_count) <- new_colnames
list_male=unique(target_uni$miRNA)
male<-list()
for( i in 1:length(list_male)){
  t3<-list_male[i]
  tmp<-subset(target_uni,miRNA==t3)
  male[[t3]]<-tmp
}

list_male=names(male)
male_p<-list()

for( i in 1:length(list_male)){
  t1<-list_male[i]
  #miRNA
  deLNC <- t1
  #mRNA
  dePC <- male[[t1]][["Gene_id"]]
  #rbind miRNA and mRNA count table
  miRNA=Male_target_count[deLNC,]
  mRNA=rlog_count_table_Male[dePC,]
  mRNA=mRNA %>% drop_na()
  rnaExpr <- rbind(miRNA,mRNA)
  dePC=rownames(mRNA)
  
  combination <- expand.grid(deLNC, dePC)
  names(combination)=c("miR","mR")
  
  
  #the relationship of miRNA-mRNA pairs
  cor_result=apply(combination,1,function(x){
    miR=as.character(x[1])
    mR=as.character(x[2])
    result=cor.test(as.numeric(rnaExpr[miR,]), as.numeric(rnaExpr[mR,]))
    score=c(pval=result$p.value,result$estimate)
    return(score)
  })
  
  result=cbind(combination,t(cor_result))
  tmp=result %>% drop_na()
  male_p[[t1]]<-tmp
}

library(data.table)
male_pearson=rbindlist(male_p, use.names=TRUE)
# 筛选条件：pval < 0.05 且 cor < 0
filtered_male_pearson <- male_pearson[male_pearson$pval < 0.05 & male_pearson$cor < 0, ]

library(dplyr)
#2.1 Male parallel
setwd("E:/miRNA/express")
all_pa_gene=read.csv("pa_express_gene.csv",header=T,sep = ",",check.names = F)
#相关（受到miRNA调控的parallel gene数量）
cor_pa_gene <- filtered_male_pearson %>%
  filter(mR %in% all_pa_gene$x)#98
write.csv(cor_pa_gene,"./cor_pa_gene.csv",row.names = F)
#相关（靶基因是parallel gene的miRNA数量）
cor_pa_miR <- data.frame(miR = unique(cor_pa_gene$miR), stringsAsFactors = FALSE)#58

#2.2 Male non-parallel
setwd("E:/miRNA/express")
all_nopa_gene=read.csv("non_pa_express_gene.csv",header=T,sep = ",",check.names = F)
#相关（受到miRNA调控的non parallel gene数量）
cor_nopa_gene <- filtered_male_pearson %>%
  filter(mR %in% all_nopa_gene$value)#1364
#相关（靶基因是non parallel gene的miRNA数量）
cor_nopa_miR <- data.frame(miR = unique(cor_nopa_gene$miR), stringsAsFactors = FALSE)#121

#bar plot没有用这里的code，在miRNA_barplot.R#
###bar plot of parallel express gene
library(ggplot2)
# 创建数据框
data <- data.frame(
  Group = rep(c("Female", "Male"), each = 2),
  Category = rep(c("genes", "miRNAs"), 2),
  Value = c(56/458*100, 40/166*100, 98/458*100, 58/166*100),
  Total = 100
)

# 自定义颜色
fill_colors <- c("genes" = "#E6A064", "miRNAs" = "#990033")

# 创建图形
p <- ggplot(data, aes(x = Group)) +
  # 绘制100%背景柱状图（浅灰色）
  geom_col(aes(y = Total, group = Category),
           position = position_dodge(width = 0.7),
           fill = "gray90", width = 0.5) +
  # 绘制实际比例柱状图
  geom_col(aes(y = Value, fill = Category),
           position = position_dodge(width = 0.7),
           width = 0.4) +
  # 设置颜色和标签
  scale_fill_manual(values = fill_colors) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 100),
                     labels = scales::percent_format(scale = 1)) +
  # 添加标签和主题
  labs(y = "Percentage", x = "") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    legend.position = "none",
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    axis.text = element_text(color = "black", size = 10),
    axis.title.y = element_text(margin = margin(r = 15))
  )
setwd("E:/miRNA/express")
ggsave("pa_plot.png", plot = p, width = 8, height = 6, dpi = 300)

###bar plot of non parallel express gene
# 创建数据框
data <- data.frame(
  Group = rep(c("Female", "Male"), each = 2),
  Category = rep(c("genes", "miRNAs"), 2),
  Value = c(749/4211*100, 130/166*100, 1364/4211*100, 121/166*100),
  Total = 100
)

# 自定义颜色
fill_colors <- c("genes" = "#0047AB", "miRNAs" = "#990033")

# 创建图形
p <- ggplot(data, aes(x = Group)) +
  # 绘制100%背景柱状图（浅灰色）
  geom_col(aes(y = Total, group = Category),
           position = position_dodge(width = 0.7),
           fill = "gray90", width = 0.5) +
  # 绘制实际比例柱状图
  geom_col(aes(y = Value, fill = Category),
           position = position_dodge(width = 0.7),
           width = 0.4) +
  # 设置颜色和标签
  scale_fill_manual(values = fill_colors) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 100),
                     labels = scales::percent_format(scale = 1)) +
  # 添加标签和主题
  labs(y = "Percentage", x = "") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    legend.position = "none",
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    axis.text = element_text(color = "black", size = 10),
    axis.title.y = element_text(margin = margin(r = 15))
  )
ggsave("non_pa_plot.png", plot = p, width = 8, height = 6, dpi = 300)

###get target gene id for splice gene###
setwd("E:/liftoff")
gtf = rtracklayer::import("Gambusia_affinis.ASM309773v1.109.gtf")
gtf_df=as.data.frame(gtf)

setwd("E:/miRNA/splice")
target=read.csv('miRNA_target_genes_hitsOnly.txt',header=F,sep="\t")
# 加载 tidyr 包
library(tidyr)

# 将 V1 列拆分成三列
target <- target %>%
  separate(V1, into = c("seqnames", "start", "end"), sep = "[:-]", remove = FALSE)

# 使用 merge() 函数按 seqnames、start 和 end 列进行合并
matched_gtf_df <- merge(target, gtf_df, by = c("seqnames", "start", "end"))
# 筛选出 type 列值为 "three_prime_utr" 的行
matched_gtf_df <- matched_gtf_df[matched_gtf_df$type == "three_prime_utr", ]
# 提取特定列并创建新的数据框
final_df <- matched_gtf_df[, c("seqnames", "start", "end", "V2", "gene_id")]
#去除重复
final_df_1 <- final_df[!duplicated(final_df[c("V2", "gene_id")]), ]

library(dplyr)
target_uni <- final_df_1 %>% select(-1:-3)
#target_uni <- final_df_1 %>% dplyr::select(-(1:3))
names(target_uni)[1] <- "miRNA"       # 修改第一列的列名
names(target_uni)[2] <- "Gene_id"     # 修改第二列的列名


###pearson analysis for splice gene###

#miRNA和mRNA

#rlog miRNA count
library(ggplot2)
library(pheatmap)
library(reshape2)
library(DESeq2)
setwd("E:/miRNA")
miRNA_count=read.table("filter_miRNA_table.txt",header = T,check.names = F,row.names =1)
sampleTable=read.csv("Sample.csv",header=T,sep = ",",check.names = F)
miRNA_count=DESeqDataSetFromMatrix(countData=miRNA_count, colData=sampleTable, design = ~Sex+Salinity+Time)
miRNA_count$Sex=relevel(miRNA_count$Sex, ref ="Female")
miRNA_count$Sex
miRNA_count=rlog(miRNA_count,blind=F)
miRNA_count=miRNA_count@assays@data@listData[[1]]
miRNA_count=as.data.frame(miRNA_count)

Female_target_count=miRNA_count[,grep("F", colnames(miRNA_count))]
Male_target_count=miRNA_count[,grep("M", colnames(miRNA_count))]


#rlog mRNA count
setwd("E:/filter")
allcount=read.csv('all_sample.txt',header=F,sep="\t")
colnames(allcount)=c('sample','gene.ID','count')
allcount=dcast(allcount,formula=gene.ID~sample,value.var="count")#gene.ID as row-,sample as col-,count as the variate
rownames(allcount)=allcount$gene.ID
noint = rownames(allcount) %in% c("__alignment_not_unique","__ambiguous","__no_feature","__not_aligned","__too_low_aQual")
head(noint,6)
count_table=allcount[!noint,]
count_table=count_table[,-1]
rm=grep("G",colnames(count_table))
count_table=count_table[,-rm]
#count_filter: gene >1 at all samples
count_table=subset(count_table,apply(count_table,1, function(x) sum(x >= 1))>= 32)
sampleTable=read.csv("Sample.csv",header=T,sep = ",",check.names = F)
rownames(sampleTable)=sampleTable$Sample
colnames(count_table)==rownames(sampleTable) #check if the count table has the same order as sample table
str(sampleTable)
id=colnames(count_table)
id=data.frame(id)
sampleTable=sampleTable[id$id,]
colnames(count_table)==rownames(sampleTable)

rlog_count_table=DESeqDataSetFromMatrix(countData=count_table, colData=sampleTable, design = ~Sex+Salinity+Time)
rlog_count_table$Sex=relevel(rlog_count_table$Sex, ref ="Female")
rlog_count_table$Sex
rlog_count_table=rlog(rlog_count_table,blind=F)
rlog_count_table=rlog_count_table@assays@data@listData[[1]]
rlog_count_table=as.data.frame(rlog_count_table)

rlog_count_table_Female=rlog_count_table[,grep("F", colnames(rlog_count_table))]
rlog_count_table_Male=rlog_count_table[,grep("M", colnames(rlog_count_table))]


###1 Female counts
# 获取 rlog_count_table_Female 的列名
new_colnames <- colnames(rlog_count_table_Female)

# 将这些列名赋值给 Female_target_count
colnames(Female_target_count) <- new_colnames
list_female=unique(target_uni$miRNA)
female<-list()
for( i in 1:length(list_female)){
  t3<-list_female[i]
  tmp<-subset(target_uni,miRNA==t3)
  female[[t3]]<-tmp
}

list_female=names(female)
female_p<-list()

for( i in 1:length(list_female)){
  t1<-list_female[i]
  #miRNA
  deLNC <- t1
  #mRNA
  dePC <- female[[t1]][["Gene_id"]]
  #rbind miRNA and mRNA count table
  miRNA=Female_target_count[deLNC,]
  mRNA=rlog_count_table_Female[dePC,]
  mRNA=mRNA %>% drop_na()
  rnaExpr <- rbind(miRNA,mRNA)
  dePC=rownames(mRNA)
  
  combination <- expand.grid(deLNC, dePC)
  names(combination)=c("miR","mR")
  
  
  #the relationship of miRNA-mRNA pairs
  cor_result=apply(combination,1,function(x){
    miR=as.character(x[1])
    mR=as.character(x[2])
    result=cor.test(as.numeric(rnaExpr[miR,]), as.numeric(rnaExpr[mR,]))
    score=c(pval=result$p.value,result$estimate)
    return(score)
  })
  
  result=cbind(combination,t(cor_result))
  tmp=result %>% drop_na()
  female_p[[t1]]<-tmp
}

library(data.table)
female_pearson=rbindlist(female_p, use.names=TRUE)
# 筛选条件：pval < 0.05 且 cor < 0
filtered_female_pearson <- female_pearson[female_pearson$pval < 0.05 & female_pearson$cor < 0, ]

library(dplyr)
#1.1 Female parallel
setwd("E:/miRNA/splice")
all_pa_gene=read.csv("pa_splice_gene.csv",header=T,sep = ",",check.names = F)
#相关（受到miRNA调控的parallel gene数量）
cor_pa_gene <- filtered_female_pearson %>%
  filter(mR %in% all_pa_gene$x)#105
#相关（靶基因是parallel gene的miRNA数量）
cor_pa_miR <- data.frame(miR = unique(cor_pa_gene$miR), stringsAsFactors = FALSE)#61

#1.2 Female non-parallel
setwd("E:/miRNA/splice")
all_nopa_gene=read.csv("non_pa_splice_gene.csv",header=T,sep = ",",check.names = F)
#相关（受到miRNA调控的non parallel gene数量）
cor_nopa_gene <- filtered_female_pearson %>%
  filter(mR %in% all_nopa_gene$value)#707
#相关（靶基因是non parallel gene的miRNA数量）
cor_nopa_miR <- data.frame(miR = unique(cor_nopa_gene$miR), stringsAsFactors = FALSE)#136

###2 Male counts
# 获取 rlog_count_table_Male 的列名
new_colnames <- colnames(rlog_count_table_Male)

# 将这些列名赋值给 Male_target_count
colnames(Male_target_count) <- new_colnames
list_male=unique(target_uni$miRNA)
male<-list()
for( i in 1:length(list_male)){
  t3<-list_male[i]
  tmp<-subset(target_uni,miRNA==t3)
  male[[t3]]<-tmp
}

list_male=names(male)
male_p<-list()

for( i in 1:length(list_male)){
  t1<-list_male[i]
  #miRNA
  deLNC <- t1
  #mRNA
  dePC <- male[[t1]][["Gene_id"]]
  #rbind miRNA and mRNA count table
  miRNA=Male_target_count[deLNC,]
  mRNA=rlog_count_table_Male[dePC,]
  mRNA=mRNA %>% drop_na()
  rnaExpr <- rbind(miRNA,mRNA)
  dePC=rownames(mRNA)
  
  combination <- expand.grid(deLNC, dePC)
  names(combination)=c("miR","mR")
  
  
  #the relationship of miRNA-mRNA pairs
  cor_result=apply(combination,1,function(x){
    miR=as.character(x[1])
    mR=as.character(x[2])
    result=cor.test(as.numeric(rnaExpr[miR,]), as.numeric(rnaExpr[mR,]))
    score=c(pval=result$p.value,result$estimate)
    return(score)
  })
  
  result=cbind(combination,t(cor_result))
  tmp=result %>% drop_na()
  male_p[[t1]]<-tmp
}

library(data.table)
male_pearson=rbindlist(male_p, use.names=TRUE)
# 筛选条件：pval < 0.05 且 cor < 0
filtered_male_pearson <- male_pearson[male_pearson$pval < 0.05 & male_pearson$cor < 0, ]

library(dplyr)
#2.1 Male parallel
setwd("E:/miRNA/splice")
all_pa_gene=read.csv("pa_splice_gene.csv",header=T,sep = ",",check.names = F)
#相关（受到miRNA调控的parallel gene数量）
cor_pa_gene <- filtered_male_pearson %>%
  filter(mR %in% all_pa_gene$x)#198
#相关（靶基因是parallel gene的miRNA数量）
cor_pa_miR <- data.frame(miR = unique(cor_pa_gene$miR), stringsAsFactors = FALSE)#75

#2.2 Male non-parallel
setwd("E:/miRNA/splice")
all_nopa_gene=read.csv("non_pa_splice_gene.csv",header=T,sep = ",",check.names = F)
#相关（受到miRNA调控的non parallel gene数量）
cor_nopa_gene <- filtered_male_pearson %>%
  filter(mR %in% all_nopa_gene$value)#1213
#相关（靶基因是non parallel gene的miRNA数量）
cor_nopa_miR <- data.frame(miR = unique(cor_nopa_gene$miR), stringsAsFactors = FALSE)#119

###bar plot of parallel splice gene
# 创建数据框
data <- data.frame(
  Group = rep(c("Female", "Male"), each = 2),
  Category = rep(c("genes", "miRNAs"), 2),
  Value = c(105/782*100, 61/166*100, 198/782*100, 75/166*100),
  Total = 100
)

# 自定义颜色
fill_colors <- c("genes" = "#FFC0CB", "miRNAs" = "#990033")

# 创建图形
p <- ggplot(data, aes(x = Group)) +
  # 绘制100%背景柱状图（浅灰色）
  geom_col(aes(y = Total, group = Category),
           position = position_dodge(width = 0.7),
           fill = "gray90", width = 0.5) +
  # 绘制实际比例柱状图
  geom_col(aes(y = Value, fill = Category),
           position = position_dodge(width = 0.7),
           width = 0.4) +
  # 设置颜色和标签
  scale_fill_manual(values = fill_colors) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 100),
                     labels = scales::percent_format(scale = 1)) +
  # 添加标签和主题
  labs(y = "Percentage", x = "") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    legend.position = "none",
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    axis.text = element_text(color = "black", size = 10),
    axis.title.y = element_text(margin = margin(r = 15))
  )
setwd("E:/miRNA/splice")
ggsave("pa_plot.png", plot = p, width = 8, height = 6, dpi = 300)

###bar plot of non parallel splice gene
# 创建数据框
data <- data.frame(
  Group = rep(c("Female", "Male"), each = 2),
  Category = rep(c("genes", "miRNAs"), 2),
  Value = c(707/4853*100, 136/166*100, 1213/4853*100, 119/166*100),
  Total = 100
)

# 自定义颜色
fill_colors <- c("genes" = "#438D74", "miRNAs" = "#990033")

# 创建图形
p <- ggplot(data, aes(x = Group)) +
  # 绘制100%背景柱状图（浅灰色）
  geom_col(aes(y = Total, group = Category),
           position = position_dodge(width = 0.7),
           fill = "gray90", width = 0.5) +
  # 绘制实际比例柱状图
  geom_col(aes(y = Value, fill = Category),
           position = position_dodge(width = 0.7),
           width = 0.4) +
  # 设置颜色和标签
  scale_fill_manual(values = fill_colors) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 100),
                     labels = scales::percent_format(scale = 1)) +
  # 添加标签和主题
  labs(y = "Percentage", x = "") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    legend.position = "none",
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    axis.text = element_text(color = "black", size = 10),
    axis.title.y = element_text(margin = margin(r = 15))
  )
ggsave("non_pa_plot.png", plot = p, width = 8, height = 6, dpi = 300)


### fisher's exact test: compare female vs. male (delete) ###
# express: parallel gene
table <- matrix(c(56, 402, 98, 360), nrow = 2, byrow = TRUE,
                dimnames = list(c("Female", "Male"), c("cor", "non-cor")))
table
result <- fisher.test(table)
result#p-value = 0.0002733

# express: parallel miRNA
table <- matrix(c(40, 126, 58, 108), nrow = 2, byrow = TRUE,
                dimnames = list(c("Female", "Male"), c("cor", "non-cor")))
table
result <- fisher.test(table)
result#p-value = 0.04052

# express: non-parallel gene
table <- matrix(c(749, 3462, 1364, 2847), nrow = 2, byrow = TRUE,
                dimnames = list(c("Female", "Male"), c("cor", "non-cor")))
table
result <- fisher.test(table)
result#p-value < 2.2e-16

# express: non-parallel miRNA
table <- matrix(c(130, 36, 121, 45), nrow = 2, byrow = TRUE,
                dimnames = list(c("Female", "Male"), c("cor", "non-cor")))
table
result <- fisher.test(table)
result#p-value = 0.3066

# splice: parallel gene
table <- matrix(c(105, 677, 198, 584), nrow = 2, byrow = TRUE,
                dimnames = list(c("Female", "Male"), c("cor", "non-cor")))
table
result <- fisher.test(table)
result#p-value = 3.156e-09

# splice: parallel miRNA
table <- matrix(c(61, 105, 75, 91), nrow = 2, byrow = TRUE,
                dimnames = list(c("Female", "Male"), c("cor", "non-cor")))
table
result <- fisher.test(table)
result#p-value = 0.1467

# splice: non-parallel gene
table <- matrix(c(707, 4146, 1213, 3640), nrow = 2, byrow = TRUE,
                dimnames = list(c("Female", "Male"), c("cor", "non-cor")))
table
result <- fisher.test(table)
result#p-value < 2.2e-16

# splice: non-parallel miRNA
table <- matrix(c(136, 30, 119, 47), nrow = 2, byrow = TRUE,
                dimnames = list(c("Female", "Male"), c("cor", "non-cor")))
table
result <- fisher.test(table)
result#p-value = 0.03706


### fisher's exact test: compare parallel vs. nonparallel ###
# express: female miRNA
table <- matrix(c(40, 126, 130, 36), nrow = 2, byrow = TRUE,
                dimnames = list(c("parallel", "nonparallel"), c("cor", "non-cor")))
table
result <- fisher.test(table)
result#p-value < 2.2e-16

# express: male miRNA
table <- matrix(c(58, 108, 121, 45), nrow = 2, byrow = TRUE,
                dimnames = list(c("parallel", "nonparallel"), c("cor", "non-cor")))
table
result <- fisher.test(table)
result#p-value = 4.813e-12

# splice: female miRNA
table <- matrix(c(61, 105, 136, 30), nrow = 2, byrow = TRUE,
                dimnames = list(c("parallel", "nonparallel"), c("cor", "non-cor")))
table
result <- fisher.test(table)
result#p-value < 2.2e-16

# splice: male miRNA
table <- matrix(c(75, 91, 119, 47), nrow = 2, byrow = TRUE,
                dimnames = list(c("parallel", "nonparallel"), c("cor", "non-cor")))
table
result <- fisher.test(table)
result#p-value = 1.444e-06