#
library(ggplot2)
library(pheatmap)
library(reshape2)
library(DESeq2)
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
write.table(count_table, file="filter_table.txt",sep="\t",quote=FALSE,row.names=T)


#Female across time DEseq2
rm=grep("M", colnames(count_table))
count_table_female=count_table[,-rm]


sampleTable=read.csv("Sample_2.csv",header=T,sep = ",",check.names = F)
rownames(sampleTable)=sampleTable$Sample
colnames(count_table_female)==rownames(sampleTable) #check if the count table has the same order as sample table
str(sampleTable)
id=colnames(count_table_female)
id=data.frame(id)
sampleTable=sampleTable[id$id,]
colnames(count_table_female)==rownames(sampleTable)
#-----------------#

dds_female_count=DESeqDataSetFromMatrix(countData=count_table_female, colData=sampleTable, design = ~Time+Salinity)
dds_female_count=DESeq(dds_female_count)
female_result=results(dds_female_count, contrast = c("Salinity", "12ppt", "0ppt"))
DEgenes_female=as.data.frame(female_result)


write.csv(DEgenes_female,file="DEgenes_female_all.csv")

resSig_female=female_result[which(female_result$padj<0.05),]
head(resSig_female[order(resSig_female$log2FoldChange, decreasing = TRUE), ])
head(resSig_female[order(resSig_female$log2FoldChange, decreasing = FALSE), ])
table(resSig_female$padj < 0.05)


write.csv(resSig_female,file="DEgenes_female_0.05_all.csv")
#---------------#

nt <- normTransform(dds_female_count, f = log2, pc = 1) # defaults to log2(x+1)
log2.norm.counts <- assay(nt)[rownames(as.data.frame(resSig_female)),]
coldata_female=sampleTable[,c("Sample","Salinity")]
coldata_female=as.data.frame(coldata_female)
#---------------#

ann_colors=list(
  Salinity=c('0ppt'='#00FFFF','12ppt'='#FF0000')
)
rownames(coldata_female)=coldata_female$Sample
coldata_female$Sample=NULL

png("./female_heatmap.png", res=600, height =5, width = 5,units="in")
pheatmap(log2.norm.counts, 
         cluster_rows=TRUE, 
         show_rownames=FALSE,
         cluster_cols=T,
         scale = "row",
         clustering_distance_cols="euclidean",
         clustering_distance_rows="euclidean",
         clustering_method = "complete",
         annotation_col=coldata_female,
         annotation_colors = ann_colors,
         border_color=NA)
dev.off()


#female 1h DEseq2
rm=grep("24h",colnames(count_table_female))
count_table_female_1h=count_table_female[,-rm]
rm=grep("48h",colnames(count_table_female_1h))
count_table_female_1h=count_table_female_1h[,-rm]
rm=grep("72h",colnames(count_table_female_1h))
count_table_female_1h=count_table_female_1h[,-rm]


sampleTable=read.csv("Sample_female_1h.csv",header=T,sep = ",",check.names = F)
rownames(sampleTable)=sampleTable$Sample
colnames(count_table_female_1h)==rownames(sampleTable) #check if the count table has the same order as sample table
str(sampleTable)
id=colnames(count_table_female_1h)
id=data.frame(id)
sampleTable=sampleTable[id$id,]
colnames(count_table_female_1h)==rownames(sampleTable)
#-----------------#

dds_female_1h_count=DESeqDataSetFromMatrix(countData=count_table_female_1h, colData=sampleTable, design = ~Salinity)
dds_female_1h_count=DESeq(dds_female_1h_count)
female_1h_result=results(dds_female_1h_count, contrast = c("Salinity", "12ppt", "0ppt"))
DEgenes_female_1h=as.data.frame(female_1h_result)


write.csv(DEgenes_female_1h,file="DEgenes_female_1h.csv")


resSig_female_1h=female_1h_result[which(female_1h_result$padj<0.05),]
head(resSig_female_1h[order(resSig_female_1h$log2FoldChange, decreasing = TRUE), ])
head(resSig_female_1h[order(resSig_female_1h$log2FoldChange, decreasing = FALSE), ])
table(resSig_female_1h$padj < 0.05)


write.csv(resSig_female_1h,file="DEgenes_female_0.05_1h.csv")
#---------------#

nt <- normTransform(dds_female_1h_count, f = log2, pc = 1) # defaults to log2(x+1)
log2.norm.counts <- assay(nt)[rownames(as.data.frame(resSig_female_1h)),]
coldata_female_1h=sampleTable[,c("Sample","Salinity")]
coldata_female_1h=as.data.frame(coldata_female_1h)
#---------------#

ann_colors=list(
  Salinity=c('0ppt'='#00FFFF','12ppt'='#FF0000')
)
rownames(coldata_female_1h)=coldata_female_1h$Sample
coldata_female_1h$Sample=NULL

png("./female_1h_heatmap.png", res=600, height =5, width = 5,units="in")
pheatmap(log2.norm.counts, 
         cluster_rows=TRUE, 
         show_rownames=FALSE,
         cluster_cols=T,
         scale = "row",
         clustering_distance_cols="euclidean",
         clustering_distance_rows="euclidean",
         clustering_method = "complete",
         annotation_col=coldata_female_1h,
         annotation_colors = ann_colors,
         border_color=NA)
dev.off()


#female 24h DEseq2
rm=grep("1h",colnames(count_table_female))
count_table_female_24h=count_table_female[,-rm]
rm=grep("48h",colnames(count_table_female_24h))
count_table_female_24h=count_table_female_24h[,-rm]
rm=grep("72h",colnames(count_table_female_24h))
count_table_female_24h=count_table_female_24h[,-rm]


sampleTable=read.csv("Sample_female_24h.csv",header=T,sep = ",",check.names = F)
rownames(sampleTable)=sampleTable$Sample
colnames(count_table_female_24h)==rownames(sampleTable) #check if the count table has the same order as sample table
str(sampleTable)
id=colnames(count_table_female_24h)
id=data.frame(id)
sampleTable=sampleTable[id$id,]
colnames(count_table_female_24h)==rownames(sampleTable)
#-----------------#

dds_female_24h_count=DESeqDataSetFromMatrix(countData=count_table_female_24h, colData=sampleTable, design = ~Salinity)
dds_female_24h_count=DESeq(dds_female_24h_count)
female_24h_result=results(dds_female_24h_count, contrast = c("Salinity", "12ppt", "0ppt"))
DEgenes_female_24h=as.data.frame(female_24h_result)


write.csv(DEgenes_female_24h,file="DEgenes_female_24h.csv")


resSig_female_24h=female_24h_result[which(female_24h_result$padj<0.05),]
head(resSig_female_24h[order(resSig_female_24h$log2FoldChange, decreasing = TRUE), ])
head(resSig_female_24h[order(resSig_female_24h$log2FoldChange, decreasing = FALSE), ])
table(resSig_female_24h$padj < 0.05)


write.csv(resSig_female_24h,file="DEgenes_female_0.05_24h.csv")
#---------------#

nt <- normTransform(dds_female_24h_count, f = log2, pc = 1) # defaults to log2(x+1)
log2.norm.counts <- assay(nt)[rownames(as.data.frame(resSig_female_24h)),]
coldata_female_24h=sampleTable[,c("Sample","Salinity")]
coldata_female_24h=as.data.frame(coldata_female_24h)
#---------------#

ann_colors=list(
  Salinity=c('0ppt'='#00FFFF','12ppt'='#FF0000')
)
rownames(coldata_female_24h)=coldata_female_24h$Sample
coldata_female_24h$Sample=NULL

png("./female_24h_heatmap.png", res=600, height =5, width = 5,units="in")
pheatmap(log2.norm.counts, 
         cluster_rows=TRUE, 
         show_rownames=FALSE,
         cluster_cols=T,
         scale = "row",
         clustering_distance_cols="euclidean",
         clustering_distance_rows="euclidean",
         clustering_method = "complete",
         annotation_col=coldata_female_24h,
         annotation_colors = ann_colors,
         border_color=NA)
dev.off()


#female 48h DEseq2
rm=grep("1h",colnames(count_table_female))
count_table_female_48h=count_table_female[,-rm]
rm=grep("24h",colnames(count_table_female_48h))
count_table_female_48h=count_table_female_48h[,-rm]
rm=grep("72h",colnames(count_table_female_48h))
count_table_female_48h=count_table_female_48h[,-rm]


sampleTable=read.csv("Sample_female_48h.csv",header=T,sep = ",",check.names = F)
rownames(sampleTable)=sampleTable$Sample
colnames(count_table_female_48h)==rownames(sampleTable) #check if the count table has the same order as sample table
str(sampleTable)
id=colnames(count_table_female_48h)
id=data.frame(id)
sampleTable=sampleTable[id$id,]
colnames(count_table_female_48h)==rownames(sampleTable)
#-----------------#

dds_female_48h_count=DESeqDataSetFromMatrix(countData=count_table_female_48h, colData=sampleTable, design = ~Salinity)
dds_female_48h_count=DESeq(dds_female_48h_count)
female_48h_result=results(dds_female_48h_count, contrast = c("Salinity", "12ppt", "0ppt"))
DEgenes_female_48h=as.data.frame(female_48h_result)


write.csv(DEgenes_female_48h,file="DEgenes_female_48h.csv")


resSig_female_48h=female_48h_result[which(female_48h_result$padj<0.05),]
head(resSig_female_48h[order(resSig_female_48h$log2FoldChange, decreasing = TRUE), ])
head(resSig_female_48h[order(resSig_female_48h$log2FoldChange, decreasing = FALSE), ])
table(resSig_female_48h$padj < 0.05)


write.csv(resSig_female_48h,file="DEgenes_female_0.05_48h.csv")
#---------------#

nt <- normTransform(dds_female_48h_count, f = log2, pc = 1) # defaults to log2(x+1)
log2.norm.counts <- assay(nt)[rownames(as.data.frame(resSig_female_48h)),]
coldata_female_48h=sampleTable[,c("Sample","Salinity")]
coldata_female_48h=as.data.frame(coldata_female_48h)
#---------------#

ann_colors=list(
  Salinity=c('0ppt'='#00FFFF','12ppt'='#FF0000')
)
rownames(coldata_female_48h)=coldata_female_48h$Sample
coldata_female_48h$Sample=NULL

png("./female_48h_heatmap.png", res=600, height =5, width = 5,units="in")
pheatmap(log2.norm.counts, 
         cluster_rows=TRUE, 
         show_rownames=FALSE,
         cluster_cols=T,
         scale = "row",
         clustering_distance_cols="euclidean",
         clustering_distance_rows="euclidean",
         clustering_method = "complete",
         annotation_col=coldata_female_48h,
         annotation_colors = ann_colors,
         border_color=NA)
dev.off()


#female 72h DEseq2
rm=grep("1h",colnames(count_table_female))
count_table_female_72h=count_table_female[,-rm]
rm=grep("24h",colnames(count_table_female_72h))
count_table_female_72h=count_table_female_72h[,-rm]
rm=grep("48h",colnames(count_table_female_72h))
count_table_female_72h=count_table_female_72h[,-rm]


sampleTable=read.csv("Sample_female_72h.csv",header=T,sep = ",",check.names = F)
rownames(sampleTable)=sampleTable$Sample
colnames(count_table_female_72h)==rownames(sampleTable) #check if the count table has the same order as sample table
str(sampleTable)
id=colnames(count_table_female_72h)
id=data.frame(id)
sampleTable=sampleTable[id$id,]
colnames(count_table_female_72h)==rownames(sampleTable)
#-----------------#

dds_female_72h_count=DESeqDataSetFromMatrix(countData=count_table_female_72h, colData=sampleTable, design = ~Salinity)
dds_female_72h_count=DESeq(dds_female_72h_count)
female_72h_result=results(dds_female_72h_count, contrast = c("Salinity", "12ppt", "0ppt"))
DEgenes_female_72h=as.data.frame(female_72h_result)


write.csv(DEgenes_female_72h,file="DEgenes_female_72h.csv")


resSig_female_72h=female_72h_result[which(female_72h_result$padj<0.05),]
head(resSig_female_72h[order(resSig_female_72h$log2FoldChange, decreasing = TRUE), ])
head(resSig_female_72h[order(resSig_female_72h$log2FoldChange, decreasing = FALSE), ])
table(resSig_female_72h$padj < 0.05)


write.csv(resSig_female_72h,file="DEgenes_female_0.05_72h.csv")
#---------------#

nt <- normTransform(dds_female_72h_count, f = log2, pc = 1) # defaults to log2(x+1)
log2.norm.counts <- assay(nt)[rownames(as.data.frame(resSig_female_72h)),]
coldata_female_72h=sampleTable[,c("Sample","Salinity")]
coldata_female_72h=as.data.frame(coldata_female_72h)
#---------------#

ann_colors=list(
  Salinity=c('0ppt'='#00FFFF','12ppt'='#FF0000')
)
rownames(coldata_female_72h)=coldata_female_72h$Sample
coldata_female_72h$Sample=NULL

png("./female_72h_heatmap.png", res=600, height =5, width = 5,units="in")
pheatmap(log2.norm.counts, 
         cluster_rows=TRUE, 
         show_rownames=FALSE,
         cluster_cols=T,
         scale = "row",
         clustering_distance_cols="euclidean",
         clustering_distance_rows="euclidean",
         clustering_method = "complete",
         annotation_col=coldata_female_72h,
         annotation_colors = ann_colors,
         border_color=NA)
dev.off()


#combine female data
setwd("E:/filter")
files <- list.files(pattern = "DEgenes_female_0.05_.*\\.csv$")


if (length(files) >= 5) {
  
  data_list_female_deseq <- lapply(files[1:5], read.csv)
  
  
  combined_female_data_deseq <- do.call(rbind, data_list_female_deseq)
  
  head(combined_female_data_deseq)
  
  
  write.csv(combined_female_data_deseq, "combined_female_data_deseq.csv", row.names = FALSE)
} else {
  stop("Less than 5 CSV files found.")
}


deseq2_genes_female <- subset(combined_female_data_deseq, select = X)


duplicated_rows <- duplicated(deseq2_genes_female$X) | duplicated(deseq2_genes_female$X, fromLast = TRUE)

deseq2_genes_female_unique <- deseq2_genes_female[!duplicated_rows, ]

write.table(deseq2_genes_female_unique, file = "deseq2_genes_female_unique.txt", sep = "\t", row.names = FALSE, quote = FALSE)


#Male across time DEseq2
rm=grep("F", colnames(count_table))
count_table_male=count_table[,-rm]


sampleTable=read.csv("Sample_1.csv",header=T,sep = ",",check.names = F)
rownames(sampleTable)=sampleTable$Sample
colnames(count_table_male)==rownames(sampleTable) #check if the count table has the same order as sample table
str(sampleTable)
id=colnames(count_table_male)
id=data.frame(id)
sampleTable=sampleTable[id$id,]
colnames(count_table_male)==rownames(sampleTable)
#-----------------#

dds_male_count=DESeqDataSetFromMatrix(countData=count_table_male, colData=sampleTable, design = ~Time+Salinity)
dds_male_count=DESeq(dds_male_count)
male_result=results(dds_male_count, contrast = c("Salinity", "12ppt", "0ppt"))
DEgenes_male=as.data.frame(male_result)


write.csv(DEgenes_male,file="DEgenes_male_all.csv")

resSig_male=male_result[which(male_result$padj<0.05),]
head(resSig_male[order(resSig_male$log2FoldChange, decreasing = TRUE), ])
head(resSig_male[order(resSig_male$log2FoldChange, decreasing = FALSE), ])
table(resSig_male$padj < 0.05)


write.csv(resSig_male,file="DEgenes_male_0.05_all.csv")
#---------------#

nt <- normTransform(dds_male_count, f = log2, pc = 1) # defaults to log2(x+1)
log2.norm.counts <- assay(nt)[rownames(as.data.frame(resSig_male)),]
coldata_male=sampleTable[,c("Sample","Salinity")]
coldata_male=as.data.frame(coldata_male)
#---------------#

ann_colors=list(
  Salinity=c('0ppt'='#00FFFF','12ppt'='#FF0000')
)
rownames(coldata_male)=coldata_male$Sample
coldata_male$Sample=NULL

png("./male_heatmap.png", res=600, height =5, width = 5,units="in")
pheatmap(log2.norm.counts, 
         cluster_rows=TRUE, 
         show_rownames=FALSE,
         cluster_cols=T,
         scale = "row",
         clustering_distance_cols="euclidean",
         clustering_distance_rows="euclidean",
         clustering_method = "complete",
         annotation_col=coldata_male,
         annotation_colors = ann_colors,
         border_color=NA)
dev.off()


#male 1h DEseq2
rm=grep("24h",colnames(count_table_male))
count_table_male_1h=count_table_male[,-rm]
rm=grep("48h",colnames(count_table_male_1h))
count_table_male_1h=count_table_male_1h[,-rm]
rm=grep("72h",colnames(count_table_male_1h))
count_table_male_1h=count_table_male_1h[,-rm]


sampleTable=read.csv("Sample_male_1h.csv",header=T,sep = ",",check.names = F)
rownames(sampleTable)=sampleTable$Sample
colnames(count_table_male_1h)==rownames(sampleTable) #check if the count table has the same order as sample table
str(sampleTable)
id=colnames(count_table_male_1h)
id=data.frame(id)
sampleTable=sampleTable[id$id,]
colnames(count_table_male_1h)==rownames(sampleTable)
#-----------------#

dds_male_1h_count=DESeqDataSetFromMatrix(countData=count_table_male_1h, colData=sampleTable, design = ~Salinity)
dds_male_1h_count=DESeq(dds_male_1h_count)
male_1h_result=results(dds_male_1h_count, contrast = c("Salinity", "12ppt", "0ppt"))
DEgenes_male_1h=as.data.frame(male_1h_result)


write.csv(DEgenes_male_1h,file="DEgenes_male_1h.csv")


resSig_male_1h=male_1h_result[which(male_1h_result$padj<0.05),]
head(resSig_male_1h[order(resSig_male_1h$log2FoldChange, decreasing = TRUE), ])
head(resSig_male_1h[order(resSig_male_1h$log2FoldChange, decreasing = FALSE), ])
table(resSig_male_1h$padj < 0.05)


write.csv(resSig_male_1h,file="DEgenes_male_0.05_1h.csv")
#---------------#

nt <- normTransform(dds_male_1h_count, f = log2, pc = 1) # defaults to log2(x+1)
log2.norm.counts <- assay(nt)[rownames(as.data.frame(resSig_male_1h)),]
coldata_male_1h=sampleTable[,c("Sample","Salinity")]
coldata_male_1h=as.data.frame(coldata_male_1h)
#---------------#

ann_colors=list(
  Salinity=c('0ppt'='#00FFFF','12ppt'='#FF0000')
)
rownames(coldata_male_1h)=coldata_male_1h$Sample
coldata_male_1h$Sample=NULL

png("./male_1h_heatmap.png", res=600, height =5, width = 5,units="in")
pheatmap(log2.norm.counts, 
         cluster_rows=TRUE, 
         show_rownames=FALSE,
         cluster_cols=T,
         scale = "row",
         clustering_distance_cols="euclidean",
         clustering_distance_rows="euclidean",
         clustering_method = "complete",
         annotation_col=coldata_male_1h,
         annotation_colors = ann_colors,
         border_color=NA)
dev.off()


#male 24h DEseq2
rm=grep("1h",colnames(count_table_male))
count_table_male_24h=count_table_male[,-rm]
rm=grep("48h",colnames(count_table_male_24h))
count_table_male_24h=count_table_male_24h[,-rm]
rm=grep("72h",colnames(count_table_male_24h))
count_table_male_24h=count_table_male_24h[,-rm]


sampleTable=read.csv("Sample_male_24h.csv",header=T,sep = ",",check.names = F)
rownames(sampleTable)=sampleTable$Sample
colnames(count_table_male_24h)==rownames(sampleTable) #check if the count table has the same order as sample table
str(sampleTable)
id=colnames(count_table_male_24h)
id=data.frame(id)
sampleTable=sampleTable[id$id,]
colnames(count_table_male_24h)==rownames(sampleTable)
#-----------------#

dds_male_24h_count=DESeqDataSetFromMatrix(countData=count_table_male_24h, colData=sampleTable, design = ~Salinity)
dds_male_24h_count=DESeq(dds_male_24h_count)
male_24h_result=results(dds_male_24h_count, contrast = c("Salinity", "12ppt", "0ppt"))
DEgenes_male_24h=as.data.frame(male_24h_result)


write.csv(DEgenes_male_24h,file="DEgenes_male_24h.csv")


resSig_male_24h=male_24h_result[which(male_24h_result$padj<0.05),]
head(resSig_male_24h[order(resSig_male_24h$log2FoldChange, decreasing = TRUE), ])
head(resSig_male_24h[order(resSig_male_24h$log2FoldChange, decreasing = FALSE), ])
table(resSig_male_24h$padj < 0.05)


write.csv(resSig_male_24h,file="DEgenes_male_0.05_24h.csv")
#---------------#

nt <- normTransform(dds_male_24h_count, f = log2, pc = 1) # defaults to log2(x+1)
log2.norm.counts <- assay(nt)[rownames(as.data.frame(resSig_male_24h)),]
coldata_male_24h=sampleTable[,c("Sample","Salinity")]
coldata_male_24h=as.data.frame(coldata_male_24h)
#---------------#

ann_colors=list(
  Salinity=c('0ppt'='#00FFFF','12ppt'='#FF0000')
)
rownames(coldata_male_24h)=coldata_male_24h$Sample
coldata_male_24h$Sample=NULL

png("./male_24h_heatmap.png", res=600, height =5, width = 5,units="in")
pheatmap(log2.norm.counts, 
         cluster_rows=TRUE, 
         show_rownames=FALSE,
         cluster_cols=T,
         scale = "row",
         clustering_distance_cols="euclidean",
         clustering_distance_rows="euclidean",
         clustering_method = "complete",
         annotation_col=coldata_male_24h,
         annotation_colors = ann_colors,
         border_color=NA)
dev.off()


#male 48h DEseq2
rm=grep("1h",colnames(count_table_male))
count_table_male_48h=count_table_male[,-rm]
rm=grep("24h",colnames(count_table_male_48h))
count_table_male_48h=count_table_male_48h[,-rm]
rm=grep("72h",colnames(count_table_male_48h))
count_table_male_48h=count_table_male_48h[,-rm]


sampleTable=read.csv("Sample_male_48h.csv",header=T,sep = ",",check.names = F)
rownames(sampleTable)=sampleTable$Sample
colnames(count_table_male_48h)==rownames(sampleTable) #check if the count table has the same order as sample table
str(sampleTable)
id=colnames(count_table_male_48h)
id=data.frame(id)
sampleTable=sampleTable[id$id,]
colnames(count_table_male_48h)==rownames(sampleTable)
#-----------------#

dds_male_48h_count=DESeqDataSetFromMatrix(countData=count_table_male_48h, colData=sampleTable, design = ~Salinity)
dds_male_48h_count=DESeq(dds_male_48h_count)
male_48h_result=results(dds_male_48h_count, contrast = c("Salinity", "12ppt", "0ppt"))
DEgenes_male_48h=as.data.frame(male_48h_result)


write.csv(DEgenes_male_48h,file="DEgenes_male_48h.csv")


resSig_male_48h=male_48h_result[which(male_48h_result$padj<0.05),]
head(resSig_male_48h[order(resSig_male_48h$log2FoldChange, decreasing = TRUE), ])
head(resSig_male_48h[order(resSig_male_48h$log2FoldChange, decreasing = FALSE), ])
table(resSig_male_48h$padj < 0.05)


write.csv(resSig_male_48h,file="DEgenes_male_0.05_48h.csv")
#---------------#

nt <- normTransform(dds_male_48h_count, f = log2, pc = 1) # defaults to log2(x+1)
log2.norm.counts <- assay(nt)[rownames(as.data.frame(resSig_male_48h)),]
coldata_male_48h=sampleTable[,c("Sample","Salinity")]
coldata_male_48h=as.data.frame(coldata_male_48h)
#---------------#

ann_colors=list(
  Salinity=c('0ppt'='#00FFFF','12ppt'='#FF0000')
)
rownames(coldata_male_48h)=coldata_male_48h$Sample
coldata_male_48h$Sample=NULL

png("./male_48h_heatmap.png", res=600, height =5, width = 5,units="in")
pheatmap(log2.norm.counts, 
         cluster_rows=TRUE, 
         show_rownames=FALSE,
         cluster_cols=T,
         scale = "row",
         clustering_distance_cols="euclidean",
         clustering_distance_rows="euclidean",
         clustering_method = "complete",
         annotation_col=coldata_male_48h,
         annotation_colors = ann_colors,
         border_color=NA)
dev.off()


#male 72h DEseq2
rm=grep("1h",colnames(count_table_male))
count_table_male_72h=count_table_male[,-rm]
rm=grep("24h",colnames(count_table_male_72h))
count_table_male_72h=count_table_male_72h[,-rm]
rm=grep("48h",colnames(count_table_male_72h))
count_table_male_72h=count_table_male_72h[,-rm]


sampleTable=read.csv("Sample_male_72h.csv",header=T,sep = ",",check.names = F)
rownames(sampleTable)=sampleTable$Sample
colnames(count_table_male_72h)==rownames(sampleTable) #check if the count table has the same order as sample table
str(sampleTable)
id=colnames(count_table_male_72h)
id=data.frame(id)
sampleTable=sampleTable[id$id,]
colnames(count_table_male_72h)==rownames(sampleTable)
#-----------------#

dds_male_72h_count=DESeqDataSetFromMatrix(countData=count_table_male_72h, colData=sampleTable, design = ~Salinity)
dds_male_72h_count=DESeq(dds_male_72h_count)
male_72h_result=results(dds_male_72h_count, contrast = c("Salinity", "12ppt", "0ppt"))
DEgenes_male_72h=as.data.frame(male_72h_result)


write.csv(DEgenes_male_72h,file="DEgenes_male_72h.csv")


resSig_male_72h=male_72h_result[which(male_72h_result$padj<0.05),]
head(resSig_male_72h[order(resSig_male_72h$log2FoldChange, decreasing = TRUE), ])
head(resSig_male_72h[order(resSig_male_72h$log2FoldChange, decreasing = FALSE), ])
table(resSig_male_72h$padj < 0.05)


write.csv(resSig_male_72h,file="DEgenes_male_0.05_72h.csv")
#---------------#

nt <- normTransform(dds_male_72h_count, f = log2, pc = 1) # defaults to log2(x+1)
log2.norm.counts <- assay(nt)[rownames(as.data.frame(resSig_male_72h)),]
coldata_male_72h=sampleTable[,c("Sample","Salinity")]
coldata_male_72h=as.data.frame(coldata_male_72h)
#---------------#

ann_colors=list(
  Salinity=c('0ppt'='#00FFFF','12ppt'='#FF0000')
)
rownames(coldata_male_72h)=coldata_male_72h$Sample
coldata_male_72h$Sample=NULL

png("./male_72h_heatmap.png", res=600, height =5, width = 5,units="in")
pheatmap(log2.norm.counts, 
         cluster_rows=TRUE, 
         show_rownames=FALSE,
         cluster_cols=T,
         scale = "row",
         clustering_distance_cols="euclidean",
         clustering_distance_rows="euclidean",
         clustering_method = "complete",
         annotation_col=coldata_male_72h,
         annotation_colors = ann_colors,
         border_color=NA)
dev.off()


#combine male data
setwd("E:/filter")
files <- list.files(pattern = "DEgenes_male_0.05_.*\\.csv$")


if (length(files) >= 5) {
  
  data_list_male_deseq <- lapply(files[1:5], read.csv)
  
  
  combined_male_data_deseq <- do.call(rbind, data_list_male_deseq)
  
  head(combined_male_data_deseq)
  
  
  write.csv(combined_male_data_deseq, "combined_male_data_deseq.csv", row.names = FALSE)
} else {
  stop("Less than 5 CSV files found.")
}


deseq2_genes_male <- subset(combined_male_data_deseq, select = X)


duplicated_rows <- duplicated(deseq2_genes_male$X) | duplicated(deseq2_genes_male$X, fromLast = TRUE)

deseq2_genes_male_unique <- deseq2_genes_male[!duplicated_rows, ]

write.table(deseq2_genes_male_unique, file = "deseq2_genes_male_unique.txt", sep = "\t", row.names = FALSE, quote = FALSE)
