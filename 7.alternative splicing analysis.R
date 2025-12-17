suppressPackageStartupMessages(library("DEXSeq"))
library(BiocParallel)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(VennDiagram)
library(rstatix)
library(car)
library(pheatmap)
library(vegan)
library(ggrepel)


#####1 data input#####
#all locations
setwd("E:/splice")
countFiles=list.files(pattern = ".txt")#需要用DEXSEQ自带的脚本获取count文件
countFiles
setwd("E:/exon")
sampleTable=read.table(file ="sampleTable.csv",header=T,sep=",")
gffFile=list.files(pattern = ".gff")


#####2 DEXSeq dataset ####
dxd=DEXSeqDataSetFromHTSeq(
  countFiles,
  sampleData=sampleTable,
  design= ~sample + exon + Salinity:exon,
  flattenedfile=gffFile)
dxd#共有243197个exon


#####3 filtering(dxd之后进行filter)#####
head(dxd)
head(featureCounts(dxd))# check countdata
#each exon >= 1 across all samples(32)
filter_dxd=subset(dxd,apply(featureCounts(dxd),1, function(x) sum(x >=1))>= 32)
filter_dxd#121996 exons
head(featureCounts(filter_dxd),3)
filtered_count = featureCounts(filter_dxd)
colnames(filtered_count)=sampleTable$Sample
write.csv(filtered_count,"filtered_count.csv",row.names = TRUE, quote = FALSE)


#####4 PCA/RDA based on the prcomp() function#####

####4.1 all PCA ####
all_dds = DESeqDataSetFromMatrix(countData = filtered_count, colData = sampleTable, design = ~Salinity+Time+Sex )
head(all_dds)
all_dds$Salinity=relevel(all_dds$Salinity, ref ="0ppt")
all_dds$Salinity
all_vst=vst(all_dds,blind=F)
head(all_vst)
head(assay(all_vst), 3)

#all PCA
pcadata=all_vst@assays@data@listData[[1]]
pcadata=t(pcadata)
id=rownames(pcadata)
id <- as.vector(id)
sample_ids <- data.frame(Sample = sampleTable$Sample)
sample_indices <- match(sample_ids$Sample, id)
pcadata <- pcadata[sample_indices, ]
pcadata=prcomp(pcadata)
summary(pcadata)#pc1~5
write.table(pcadata[["x"]], file="pcadata_gene.txt",sep="\t",quote=FALSE,row.names=T)

rda=pcadata$x[,1:5]
write.csv(rda,file="rda_all.csv")
genet=read.table("rda_all.csv",header=TRUE,check.names=FALSE,row.names = 1,sep=",")
design_test=read.table("sampleTable.csv",sep=",",header = TRUE,check.names = T,row.names = 1)
genet=genet[rownames(design_test),]

rda.vpa=varpart(genet[,1:5],~Salinity,~Time,~Sex,data=design_test)
str(design_test)
rda.vpa
str(rda.vpa)
summary(rda.vpa)


#rda
rda=rda(genet[,1:5]~+Salinity+Time+Sex,data=design_test)
summary(rda)
RsquareAdj(rda)#0.1305833
anova(rda)#1.9312  0.012 *
plot(rda)

### look the sig###
#salinity
rda_sal=rda(genet[,1:5]~Salinity,data=design_test)
summary(rda_sal)
RsquareAdj(rda_sal)#0.006486492
anova(rda_sal)#1.2024  0.317
plot(rda_sal)
rda_sal_partial=rda(genet[,1:5]~Salinity+Condition(Time+Sex),data=design_test)
summary(rda_sal_partial)
RsquareAdj(rda_sal_partial)#0.01204363
anova(rda_sal_partial)#1.374  0.243
plot(rda_sal_partial)

#time
rda_time=rda(genet[,1:5]~Time,data=design_test)
summary(rda_time)
RsquareAdj(rda_time)#0.1162975
anova(rda_time)#2.3599  0.007 **
plot(rda_time)
rda_time_partial=rda(genet[,1:5]~Time+Condition(Salinity+Sex),data=design_test)
summary(rda_time_partial)
RsquareAdj(rda_time_partial)#0.1257959
anova(rda_time_partial)#2.3987  0.005 **
plot(rda_time_partial)

#sex
rda_sex=rda(genet[,1:5]~Sex,data=design_test)
summary(rda_sex)
RsquareAdj(rda_sex)#-0.001858636
anova(rda_sex)#0.9425   0.44
plot(rda_sex)
rda_sex_partial=rda(genet[,1:5]~Sex+Condition(Salinity+Time),data=design_test)
summary(rda_sex_partial)
RsquareAdj(rda_sex_partial)#0.002771268
anova(rda_sex_partial)#1.0861  0.354
plot(rda_sex_partial)

#
#
sommaire = summary(rda)
sommaire
df1  <- data.frame(sommaire$sites[,1:2])      
df2_plot<- data.frame(sommaire$biplot[,1:2])
df2  <- data.frame(sommaire$species[,1:2])


# prepare df1 with group info
df1<-merge(df1, design_test, by=0, all=TRUE)
rownames(df1)=df1$Row.names
df1$Row.names <- NULL
colnames(df1)

df1$Time<-as.factor(df1$Time)
df1$Salinity<-as.factor(df1$Salinity)
df1$Sex<-as.factor(df1$Sex)

rownames(df2)

df2subset_plot<-df2_plot[c("Salinity12ppt","Time48h","SexMale"),]
#title_lab=expression(bold(paste(" adj. ",R^{2}," =13.06%;",italic(p),"-value"," = 0.012"))) # fill with you own stats in sommaire

rda1 =round(rda$CCA$eig[1]/sum(rda$CCA$eig)*100,2) # 52.4
rda2 =round(rda$CCA$eig[2]/sum(rda$CCA$eig)*100,2) # 21.0
theme = theme(
  panel.grid.major = element_blank(),
  panel.background = element_blank(),
  panel.grid.minor = element_blank(),
  strip.background = element_blank(),
  axis.text.x = element_text(colour = "black", size = 14),  # 调大坐标轴字体
  axis.text.y = element_text(colour = "black", size = 14),  # 调大坐标轴字体
  axis.title.x = element_text(colour = "black", size = 12, face = "bold"),
  axis.title.y = element_text(colour = "black", size = 12, face = "bold"),
  axis.ticks = element_line(colour = "black"),
  plot.margin = unit(c(1,1,1,1), "line"),
  plot.title = element_text(size = 12, hjust = 0),
  aspect.ratio = 1,
  legend.background = element_rect(fill = "NA"),
  panel.border = element_rect(colour = "black", fill = NA, size = 1),
  legend.text = element_text(size = 10),
  legend.title = element_text(size = 10, face = "bold"),
  legend.key = element_rect(fill = NA)
)

b=ggplot() +
  # 增大点的大小（size=4）
  geom_point(size = 4, data = df1, aes(x = RDA1, y = RDA2, fill = Time, shape = Sex, color = Salinity)) +  
  # 使用更鲜明的四色方案
  scale_fill_manual(values = c("#FF6B6B", "#4ECDC4", "#556270", "#C7F464"),
                    name = "Time", labels = c("1h", "24h", "48h", "72h")) +
  scale_shape_manual(values = c(21, 24), name = "Sex", labels = c("Female", "Male")) +
  # 设置Salinity图例为空心
  scale_color_manual(values = c('blue', "red"), 
                     name = "Salinity", labels = c("0ppt", "12ppt"),
                     guide = guide_legend(
                       override.aes = list(shape = 21, fill = NA, size = 4))) +  # 空心设置
  
  guides(
    fill = guide_legend(override.aes = list(shape = 21, size = 4)),  # 调大图例点大小
    shape = guide_legend(override.aes = list(size = 4))  # 调大图例点大小
  ) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  coord_fixed() +
  theme +
  geom_segment(data = df2subset_plot, 
               aes(x = 0, xend = RDA1, y = 0, yend = RDA2),
               color = "black", arrow = arrow(length = unit(0.01, "npc"))) +
  labs(x = "RDA1 (52.4%)", y = "RDA2 (21.0%)") +
  geom_text_repel(data = df2subset_plot,
                  aes(x = RDA1, y = RDA2, label = c("Sex", "Time**", "Salinity")),
                  size = 6) +
  theme(legend.spacing.x = unit(0.5, 'cm'), 
        legend.spacing.y = unit(0, 'cm'))

ggsave("./RDA_DS_rename_pca2.png", device="png",
       #path = "xxx",
       width = 6.2, height = 6.2, units="in",#scaling = 0.7,
       dpi=300)


####5 female DEU####
filter_dxd_female = filter_dxd[,colData(filter_dxd)$Sex %in% "Female"]
head(featureCounts(filter_dxd_female))

filter_dxd_female = estimateSizeFactors(filter_dxd_female)

####5.1 Female across time: DEU for salinity using ReducedModel####
sampleAnnotation(filter_dxd_female)
#We specify two design formula, which indicate that the time factor should be treated as a blocking factor
formulaFullModel    =  ~ sample + exon + Time:exon + Salinity:exon
formulaReducedModel =  ~ sample + exon + Time:exon

dxd_Fe_sal = estimateDispersions(filter_dxd_female, formula = formulaFullModel )
dxd_Fe_sal_2 = testForDEU(dxd_Fe_sal, 
                            reducedModel = formulaReducedModel, 
                            fullModel = formulaFullModel )
dxr2 = DEXSeqResults(dxd_Fe_sal_2)
table(dxr2$padj < 0.01)#0
table(dxr2$padj < 0.05 )#1

FeSig_sal_0.05=dxr2[which(dxr2$padj<0.05),]
FeSig_sal_0.05=as.data.frame(FeSig_sal_0.05)
write.csv(FeSig_sal_0.05, file="./FeSig_sal_0.05.csv")
write.csv(dxr2, file="./Fe_all.csv")

hist(FeSig_sal_0.05$pvalue)# data structure
hist(FeSig_sal_0.05$stat)# data structure


#5.2 female 1h DEU
BPPARAM = SnowParam(8)

filter_dxd_female_1h = filter_dxd_female[,colData(filter_dxd_female)$Time %in% "1h"]
head(featureCounts(filter_dxd_female_1h))


filter_dxd_female_1h = estimateSizeFactors(filter_dxd_female_1h)
filter_dxd_female_1h = estimateDispersions(filter_dxd_female_1h,BPPARAM=BPPARAM)
filter_dxd_female_1h = testForDEU(filter_dxd_female_1h,BPPARAM=BPPARAM) 

filter_dxd_female_1h = estimateExonFoldChanges(filter_dxd_female_1h, fitExpToVar="Salinity",BPPARAM=BPPARAM)
filter_dxr_female_1h = DEXSeqResults(filter_dxd_female_1h)
filter_dxr_female_1h
table(filter_dxr_female_1h$padj < 0.01)#25
table(filter_dxr_female_1h$padj < 0.05 )#37


Fe1hSig_Sal_0.05=filter_dxr_female_1h[which(filter_dxr_female_1h$padj<0.05),]
Fe1hSig_Sal_0.05=as.data.frame(Fe1hSig_Sal_0.05)
write.csv(Fe1hSig_Sal_0.05, file="./Fe1hSig_Sal_0.05.csv")
write.csv(filter_dxr_female_1h, file="./Fe1h_all.csv")

hist(Fe1hSig_Sal_0.05$pvalue)# data structure
hist(Fe1hSig_Sal_0.05$stat)# data structure


#5.3 female 24h DEU
BPPARAM = SnowParam(8)

filter_dxd_female_24h = filter_dxd_female[,colData(filter_dxd_female)$Time %in% "24h"]
head(featureCounts(filter_dxd_female_24h))


filter_dxd_female_24h = estimateSizeFactors(filter_dxd_female_24h)
filter_dxd_female_24h = estimateDispersions(filter_dxd_female_24h,BPPARAM=BPPARAM)
filter_dxd_female_24h = testForDEU(filter_dxd_female_24h,BPPARAM=BPPARAM) 

filter_dxd_female_24h = estimateExonFoldChanges(filter_dxd_female_24h, fitExpToVar="Salinity",BPPARAM=BPPARAM)
filter_dxr_female_24h = DEXSeqResults(filter_dxd_female_24h) 
filter_dxr_female_24h
table(filter_dxr_female_24h$padj < 0.01)#97
table(filter_dxr_female_24h$padj < 0.05 )#186


Fe24hSig_Sal_0.05=filter_dxr_female_24h[which(filter_dxr_female_24h$padj<0.05),]
Fe24hSig_Sal_0.05=as.data.frame(Fe24hSig_Sal_0.05)
write.csv(Fe24hSig_Sal_0.05, file="./Fe24hSig_Sal_0.05.csv")
write.csv(filter_dxr_female_24h, file="./Fe24h_all.csv")

hist(Fe24hSig_Sal_0.05$pvalue)# data structure
hist(Fe24hSig_Sal_0.05$stat)# data structure


#5.4 female 48h DEU
BPPARAM = SnowParam(8)

filter_dxd_female_48h = filter_dxd_female[,colData(filter_dxd_female)$Time %in% "48h"]
head(featureCounts(filter_dxd_female_48h))


filter_dxd_female_48h = estimateSizeFactors(filter_dxd_female_48h)
filter_dxd_female_48h = estimateDispersions(filter_dxd_female_48h,BPPARAM=BPPARAM)
filter_dxd_female_48h = testForDEU(filter_dxd_female_48h,BPPARAM=BPPARAM) 

filter_dxd_female_48h = estimateExonFoldChanges(filter_dxd_female_48h, fitExpToVar="Salinity",BPPARAM=BPPARAM)
filter_dxr_female_48h = DEXSeqResults(filter_dxd_female_48h) 
filter_dxr_female_48h
table(filter_dxr_female_48h$padj < 0.01)#26
table(filter_dxr_female_48h$padj < 0.05 )#49


Fe48hSig_Sal_0.05=filter_dxr_female_48h[which(filter_dxr_female_48h$padj<0.05),]
Fe48hSig_Sal_0.05=as.data.frame(Fe48hSig_Sal_0.05)
write.csv(Fe48hSig_Sal_0.05, file="./Fe48hSig_Sal_0.05.csv")
write.csv(filter_dxr_female_48h, file="./Fe48h_all.csv")

hist(Fe48hSig_Sal_0.05$pvalue)# data structure
hist(Fe48hSig_Sal_0.05$stat)# data structure


#5.5 female 72h DEU
BPPARAM = SnowParam(8)

filter_dxd_female_72h = filter_dxd_female[,colData(filter_dxd_female)$Time %in% "72h"]
head(featureCounts(filter_dxd_female_72h))


filter_dxd_female_72h = estimateSizeFactors(filter_dxd_female_72h)
filter_dxd_female_72h = estimateDispersions(filter_dxd_female_72h,BPPARAM=BPPARAM)
filter_dxd_female_72h = testForDEU(filter_dxd_female_72h,BPPARAM=BPPARAM) 

filter_dxd_female_72h = estimateExonFoldChanges(filter_dxd_female_72h, fitExpToVar="Salinity",BPPARAM=BPPARAM)
filter_dxr_female_72h = DEXSeqResults(filter_dxd_female_72h) 
filter_dxr_female_72h
table(filter_dxr_female_72h$padj < 0.01)#171
table(filter_dxr_female_72h$padj < 0.05 )#370


Fe72hSig_Sal_0.05=filter_dxr_female_72h[which(filter_dxr_female_72h$padj<0.05),]
Fe72hSig_Sal_0.05=as.data.frame(Fe72hSig_Sal_0.05)
write.csv(Fe72hSig_Sal_0.05, file="./Fe72hSig_Sal_0.05.csv")
write.csv(filter_dxr_female_72h, file="./Fe72h_all.csv")

hist(Fe72hSig_Sal_0.05$pvalue)# data structure
hist(Fe72hSig_Sal_0.05$stat)# data structure


####6 male DEU####
filter_dxd_male = filter_dxd[,colData(filter_dxd)$Sex %in% "Male"]
head(featureCounts(filter_dxd_male))

filter_dxd_male = estimateSizeFactors(filter_dxd_male)

####6.1 male across time: DEU for salinity using ReducedModel####
sampleAnnotation(filter_dxd_male)
#We specify two design formula, which indicate that the time factor should be treated as a blocking factor
formulaFullModel    =  ~ sample + exon + Time:exon + Salinity:exon
formulaReducedModel =  ~ sample + exon + Time:exon

dxd_Male_sal = estimateDispersions(filter_dxd_male, formula = formulaFullModel )
dxd_Male_sal_2 = testForDEU(dxd_Male_sal, 
                          reducedModel = formulaReducedModel, 
                          fullModel = formulaFullModel )
dxr2 = DEXSeqResults(dxd_Male_sal_2)
table(dxr2$padj < 0.01)#1
table(dxr2$padj < 0.05 )#1

MaleSig_sal_0.05=dxr2[which(dxr2$padj<0.05),]
MaleSig_sal_0.05=as.data.frame(MaleSig_sal_0.05)
write.csv(MaleSig_sal_0.05, file="./MaleSig_sal_0.05.csv")
write.csv(dxr2, file="./Male_all.csv")

hist(MaleSig_sal_0.05$pvalue)# data structure
hist(MaleSig_sal_0.05$stat)# data structure


#6.2 male 1h DEU
BPPARAM = SnowParam(8)

filter_dxd_male_1h = filter_dxd_male[,colData(filter_dxd_male)$Time %in% "1h"]
head(featureCounts(filter_dxd_male_1h))


filter_dxd_male_1h = estimateSizeFactors(filter_dxd_male_1h)
filter_dxd_male_1h = estimateDispersions(filter_dxd_male_1h,BPPARAM=BPPARAM)
filter_dxd_male_1h = testForDEU(filter_dxd_male_1h,BPPARAM=BPPARAM) 

filter_dxd_male_1h = estimateExonFoldChanges(filter_dxd_male_1h, fitExpToVar="Salinity",BPPARAM=BPPARAM)
filter_dxr_male_1h = DEXSeqResults(filter_dxd_male_1h) 
filter_dxr_male_1h
table(filter_dxr_male_1h$padj < 0.01)#44
table(filter_dxr_male_1h$padj < 0.05 )#75


Male1hSig_Sal_0.05=filter_dxr_male_1h[which(filter_dxr_male_1h$padj<0.05),]
Male1hSig_Sal_0.05=as.data.frame(Male1hSig_Sal_0.05)
write.csv(Male1hSig_Sal_0.05, file="./Male1hSig_Sal_0.05.csv")
write.csv(filter_dxr_male_1h, file="./Male1h_all.csv")

hist(Male1hSig_Sal_0.05$pvalue)# data structure
hist(Male1hSig_Sal_0.05$stat)# data structure


#6.3 male 24h DEU
BPPARAM = SnowParam(8)

filter_dxd_male_24h = filter_dxd_male[,colData(filter_dxd_male)$Time %in% "24h"]
head(featureCounts(filter_dxd_male_24h))


filter_dxd_male_24h = estimateSizeFactors(filter_dxd_male_24h)
filter_dxd_male_24h = estimateDispersions(filter_dxd_male_24h,BPPARAM=BPPARAM)
filter_dxd_male_24h = testForDEU(filter_dxd_male_24h,BPPARAM=BPPARAM) 

filter_dxd_male_24h = estimateExonFoldChanges(filter_dxd_male_24h, fitExpToVar="Salinity",BPPARAM=BPPARAM)
filter_dxr_male_24h = DEXSeqResults(filter_dxd_male_24h) 
filter_dxr_male_24h
table(filter_dxr_male_24h$padj < 0.01)#131
table(filter_dxr_male_24h$padj < 0.05 )#205


Male24hSig_Sal_0.05=filter_dxr_male_24h[which(filter_dxr_male_24h$padj<0.05),]
Male24hSig_Sal_0.05=as.data.frame(Male24hSig_Sal_0.05)
write.csv(Male24hSig_Sal_0.05, file="./Male24hSig_Sal_0.05.csv")
write.csv(filter_dxr_male_24h, file="./Male24h_all.csv")

hist(Male24hSig_Sal_0.05$pvalue)# data structure
hist(Male24hSig_Sal_0.05$stat)# data structure


#6.4 male 48h DEU
BPPARAM = SnowParam(8)

filter_dxd_male_48h = filter_dxd_male[,colData(filter_dxd_male)$Time %in% "48h"]
head(featureCounts(filter_dxd_male_48h))


filter_dxd_male_48h = estimateSizeFactors(filter_dxd_male_48h)
filter_dxd_male_48h = estimateDispersions(filter_dxd_male_48h,BPPARAM=BPPARAM)
filter_dxd_male_48h = testForDEU(filter_dxd_male_48h,BPPARAM=BPPARAM) 

filter_dxd_male_48h = estimateExonFoldChanges(filter_dxd_male_48h, fitExpToVar="Salinity",BPPARAM=BPPARAM)
filter_dxr_male_48h = DEXSeqResults(filter_dxd_male_48h) 
filter_dxr_male_48h
table(filter_dxr_male_48h$padj < 0.01)#34
table(filter_dxr_male_48h$padj < 0.05 )#59


Male48hSig_Sal_0.05=filter_dxr_male_48h[which(filter_dxr_male_48h$padj<0.05),]
Male48hSig_Sal_0.05=as.data.frame(Male48hSig_Sal_0.05)
write.csv(Male48hSig_Sal_0.05, file="./Male48hSig_Sal_0.05.csv")
write.csv(filter_dxr_male_48h, file="./Male48h_all.csv")

hist(Male48hSig_Sal_0.05$pvalue)# data structure
hist(Male48hSig_Sal_0.05$stat)# data structure


#6.5 male 72h DEU
BPPARAM = SnowParam(8)

filter_dxd_male_72h = filter_dxd_male[,colData(filter_dxd_male)$Time %in% "72h"]
head(featureCounts(filter_dxd_male_72h))


filter_dxd_male_72h = estimateSizeFactors(filter_dxd_male_72h)
filter_dxd_male_72h = estimateDispersions(filter_dxd_male_72h,BPPARAM=BPPARAM)
filter_dxd_male_72h = testForDEU(filter_dxd_male_72h,BPPARAM=BPPARAM) 

filter_dxd_male_72h = estimateExonFoldChanges(filter_dxd_male_72h, fitExpToVar="Salinity",BPPARAM=BPPARAM)
filter_dxr_male_72h = DEXSeqResults(filter_dxd_male_72h) 
filter_dxr_male_72h
table(filter_dxr_male_72h$padj < 0.01)#2070
table(filter_dxr_male_72h$padj < 0.05 )#4152


Male72hSig_Sal_0.05=filter_dxr_male_72h[which(filter_dxr_male_72h$padj<0.05),]
Male72hSig_Sal_0.05=as.data.frame(Male72hSig_Sal_0.05)
write.csv(Male72hSig_Sal_0.05, file="./Male72hSig_Sal_0.05.csv")
write.csv(filter_dxr_male_72h, file="./Male72h_all.csv")

hist(Male72hSig_Sal_0.05$pvalue)# data structure
hist(Male72hSig_Sal_0.05$stat)# data structure


#7 combine DEU data
#7.1 combine female DEU data
setwd("E:/exon")

deu_across <- read.csv("Fe_Sig_Sal_0.05_acrosstime.csv", header = TRUE, sep = ",")
deu_1h <- read.csv("Fe_Sig_Sal_0.05_1h.csv", header = TRUE, sep = ",")
deu_24h <- read.csv("Fe_Sig_Sal_0.05_24h.csv", header = TRUE, sep = ",")
deu_48h <- read.csv("Fe_Sig_Sal_0.05_48h.csv", header = TRUE, sep = ",")
deu_72h <- read.csv("Fe_Sig_Sal_0.05_72h.csv", header = TRUE, sep = ",")


exon_across <- deu_across$Exon
exon_1h <- deu_1h$Exon
exon_24h <- deu_24h$Exon
exon_48h <- deu_48h$Exon
exon_72h <- deu_72h$Exon

exon_combined <- c(exon_across, exon_1h, exon_24h, exon_48h, exon_72h)
exon_combined


duplicates_in_exon <- duplicated(exon_combined)

duplicate_positions <- which(duplicates_in_exon)
print(duplicate_positions)


if (length(duplicate_positions) > 0) {
  print(exon_combined[duplicate_positions])
} else {
  print("No duplicates found.")
}

exon_combined_unique <- unique(exon_combined)
print(exon_combined_unique)

exon_combined_unique <- exon_combined_unique[!grepl("\\+", exon_combined_unique)]
print(exon_combined_unique)


write.table(exon_combined_unique, file = "exon_combined_female_unique.txt", row.names = FALSE, col.names = FALSE, sep = "\t")

#7.2 combine male DEU data
setwd("E:/exon")

deu_across <- read.csv("Male_Sig_Sal_0.05_acrosstime.csv", header = TRUE, sep = ",")
deu_1h <- read.csv("Male_Sig_Sal_0.05_1h.csv", header = TRUE, sep = ",")
deu_24h <- read.csv("Male_Sig_Sal_0.05_24h.csv", header = TRUE, sep = ",")
deu_48h <- read.csv("Male_Sig_Sal_0.05_48h.csv", header = TRUE, sep = ",")
deu_72h <- read.csv("Male_Sig_Sal_0.05_72h.csv", header = TRUE, sep = ",")


exon_across <- deu_across$Exon
exon_1h <- deu_1h$Exon
exon_24h <- deu_24h$Exon
exon_48h <- deu_48h$Exon
exon_72h <- deu_72h$Exon

exon_combined <- c(exon_across, exon_1h, exon_24h, exon_48h, exon_72h)
exon_combined


duplicates_in_exon <- duplicated(exon_combined)

duplicate_positions <- which(duplicates_in_exon)
print(duplicate_positions)


if (length(duplicate_positions) > 0) {
  print(exon_combined[duplicate_positions])
} else {
  print("No duplicates found.")
}

exon_combined_unique <- unique(exon_combined)
print(exon_combined_unique)

exon_combined_unique <- exon_combined_unique[!grepl("\\+", exon_combined_unique)]
print(exon_combined_unique)


write.table(exon_combined_unique, file = "exon_combined_male_unique.txt", row.names = FALSE, col.names = FALSE, sep = "\t")


#8 combine DEU & WGCNA exons
#8.1 Female DEU & WGCNA exons
setwd("E:/exon")
DEU_exons_female <- read.table("exon_combined_female_unique.txt", header = F, sep = "\t", col.names = c("Exon"))
setwd("E:/exon/female")
combined_wgcna_female_file <- read.table("combined_wgcna_female_file.txt", header = F, sep = "\t", col.names = c("Exon"))


merged_female_exons <- rbind(DEU_exons_female, combined_wgcna_female_file)
merged_female_exons_unique <- merged_female_exons[!duplicated(merged_female_exons$Exon), ]
write.table(merged_female_exons_unique, file = "merged_female_exons_unique.txt", sep = "\t", row.names = FALSE, quote = FALSE)
merged_female_exons_unique <- read.table("merged_female_exons_unique.txt", header = T, sep = "\t")


#8.2 Male DEU & WGCNA exons
setwd("E:/exon")
DEU_exons_male <- read.table("exon_combined_male_unique.txt", header = F, sep = "\t", col.names = c("Exon"))
setwd("E:/exon/male")
combined_wgcna_male_file <- read.table("combined_wgcna_male_file.txt", header = F, sep = "\t", col.names = c("Exon"))


merged_male_exons <- rbind(DEU_exons_male, combined_wgcna_male_file)
merged_male_exons_unique <- merged_male_exons[!duplicated(merged_male_exons$Exon), ]
write.table(merged_male_exons_unique, file = "merged_male_exons_unique.txt", sep = "\t", row.names = FALSE, quote = FALSE)
merged_male_exons_unique <- read.table("merged_male_exons_unique.txt", header = T, sep = "\t")


#
female_exons <- merged_female_exons_unique$x
male_exons <- merged_male_exons_unique$x


#choose shared exons in female and male
common_exons <- intersect(female_exons, male_exons)
print(common_exons)

setwd("E:/exon")
write.table(common_exons, file = "shared_exons.txt", sep = "\t", row.names = FALSE, quote = FALSE)
shared_exons <- read.table("shared_exons.txt", header = T, sep = "\t")

#
filtered_count=read.csv("filtered_count.csv",header = T)#121996
final_filtered_count <- filtered_count[!grepl("\\+", filtered_count[, 1]), ]#118223
write.csv(final_filtered_count,"final_filtered_count.csv",row.names = TRUE, quote = FALSE)

#9 direction analysis ####
#1h exon direction
library(ggplot2)
setwd("E:/exon")
shared_exons <- read.table("shared_exons.txt", header = T, sep = "\t")

Sal1=read.csv("Fe1h_all.csv",header = T)
Sal2=read.csv("Male1h_all.csv",header = T)

names(Sal1)[names(Sal1) == 'X'] <- 'exonname.Sal1'
names(Sal2)[names(Sal2) == 'X'] <- 'exonname.Sal2'
names(Sal1)[names(Sal1) == 'log2fold_12ppt_0ppt'] <- 'log2FoldChange.Sal1'
names(Sal2)[names(Sal2) == 'log2fold_12ppt_0ppt'] <- 'log2FoldChange.Sal2'
Sal1 <- Sal1[, c("exonname.Sal1", "log2FoldChange.Sal1")]
Sal2 <- Sal2[, c("exonname.Sal2", "log2FoldChange.Sal2")]
Sal1 <- subset(Sal1, !grepl("\\+", exonname.Sal1))
Sal2 <- subset(Sal2, !grepl("\\+", exonname.Sal2))
Sal1 <- na.omit(Sal1, cols = "log2FoldChange.Sal1")
Sal2 <- na.omit(Sal2, cols = "log2FoldChange.Sal2")


Sal1_Sal2_Merged <- cbind(Sal1, Sal2)
Sal1_Sal2_Merged <- Sal1_Sal2_Merged[Sal1_Sal2_Merged$exonname.Sal1 %in% shared_exons$x, ]
Sal1_Sal2_Merged <- Sal1_Sal2_Merged[Sal1_Sal2_Merged$exonname.Sal2 %in% shared_exons$x, ]

#shared_exons_final
shared_exons_final <- merge(shared_exons, Sal1_Sal2_Merged, by.x = "x", by.y = "exonname.Sal1")
shared_exons_final <- shared_exons_final$x
setwd("E:/exon")
write.table(shared_exons_final, file = "shared_exons_final.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#The non-logFC exon
unique_to_shared_exons <- setdiff(shared_exons$x, Sal1_Sal2_Merged$exonname.Sal1)
write.table(unique_to_shared_exons, file = "non_logFC_exon.txt", sep = "\t", row.names = FALSE, quote = FALSE)
setwd("E:/exon")
unique_to_shared_exons <- read.table("non_logFC_exon.txt", header = T, sep = "\t")

#Female salinity responsive exons
setwd("E:/exon/female")
merged_female_exons_unique <- read.table("merged_female_exons_unique.txt", header = T, sep = "\t")
merged_female_exons_unique <- merged_female_exons_unique[!merged_female_exons_unique$x %in% unique_to_shared_exons$x, ]#40669
write.table(merged_female_exons_unique, file = "merged_female_exons_unique.txt", sep = "\t", row.names = FALSE, quote = FALSE)
merged_female_exons_unique <- read.table("merged_female_exons_unique.txt", header = T, sep = "\t")
#Male salinity responsive exons
setwd("E:/exon/male")
merged_male_exons_unique <- read.table("merged_male_exons_unique.txt", header = T, sep = "\t")
merged_male_exons_unique <- merged_male_exons_unique[!merged_male_exons_unique$x %in% unique_to_shared_exons$x, ]#43378
write.table(merged_male_exons_unique, file = "merged_male_exons_unique.txt", sep = "\t", row.names = FALSE, quote = FALSE)
merged_male_exons_unique <- read.table("merged_male_exons_unique.txt", header = T, sep = "\t")

#The Hypergeometric Distribution: phyper: shared exon
phyper(14868, 43378, 74845, 40669, lower.tail = FALSE)#0.7521463

#venn for shared exons
library(VennDiagram)
list2=list(merged_female_exons_unique$x,merged_male_exons_unique$x)
names(list2) <- c("female","male")
setwd("E:/exon")
venn = venn.diagram(
  list2,
  filename = "venn_responsive_exon.tiff",
  fontface="bold",
  cat.fontface="bold",
  main.fontface="bold",
  fontfamily = "Arial",
  cat.fontfamily="Arial",
  main.fontfamily = "Arial",
  fill=c("white","white"),
  cat.cex = 0.5,
  ext.text = F
)



Sal1_Sal2_Merged$direction[which(Sal1_Sal2_Merged$log2FoldChange.Sal1*Sal1_Sal2_Merged$log2FoldChange.Sal2 >0)] <- "same"
Sal1_Sal2_Merged$direction[which(Sal1_Sal2_Merged$log2FoldChange.Sal1*Sal1_Sal2_Merged$log2FoldChange.Sal2 <0)] <- "opposite"


###count the number of exon in different quadrant
Sal1_Sal2_Merged$quarant=ifelse(Sal1_Sal2_Merged$log2FoldChange.Sal1 >0 & Sal1_Sal2_Merged$log2FoldChange.Sal2 >0, "first",
                                ifelse(Sal1_Sal2_Merged$log2FoldChange.Sal1 <0 & Sal1_Sal2_Merged$log2FoldChange.Sal2 >0, "second",
                                       ifelse(Sal1_Sal2_Merged$log2FoldChange.Sal1 <0 & Sal1_Sal2_Merged$log2FoldChange.Sal2 <0, "third", "fourth")))
nrow(Sal1_Sal2_Merged[Sal1_Sal2_Merged$quarant=="first",])#3362
nrow(Sal1_Sal2_Merged[Sal1_Sal2_Merged$quarant=="second",])#4153
nrow(Sal1_Sal2_Merged[Sal1_Sal2_Merged$quarant=="third",])#4009
nrow(Sal1_Sal2_Merged[Sal1_Sal2_Merged$quarant=="fourth",])#3345

p <- ggplot(Sal1_Sal2_Merged, aes(x=log2FoldChange.Sal1, y=log2FoldChange.Sal2, fill=direction)) +
  geom_hline(yintercept=0, col="grey50") +
  geom_vline(xintercept=0, col="grey50") +
  geom_point(colour = "black",shape = 21,size=3, alpha=1)+
  scale_fill_manual(values=c("white","grey"))+
  xlab(expression("log"["2"]~"FC"~"for salinity-responsive spliced exons in female"))+
  ylab(expression("log"["2"]~"FC"~"for salinity-responsive spliced exons in male"))+
  #geom_abline(intercept=0, slope=1, show.legend = FALSE, lty=2) +
  #geom_smooth(method="lm", se=FALSE,aes(group=1), color="black", show.legend = F,formula = y ~ x) +
  #stat_poly_eq(aes(group=1,label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~~')),formula = y ~ x, parse = T) +
  theme_bw() +
  #guides(shape = guide_legend(override.aes = list(size = 5, alpha=1),
  #fill = guide_legend(override.aes = list(size = 5, alpha=1))))+
  theme(
    strip.background = element_blank(),
    strip.text = element_blank(),
    legend.position = "none",
    aspect.ratio = 1,  # 设置宽高比为1
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    plot.title = element_text(
      size = 22,            # 标题字体大小
      hjust = 0.5,          # 水平居中
      margin = margin(b = 15),# 下边距15点（增加与图的间距）
      face = "bold"
    )
  ) +
  guides(fill=guide_legend(title=NULL)) +
  # 添加"1 h"标题（顶部居中）
  ggtitle("1 h") +  # 关键添加
  geom_text(x=4.3,y=3.3,label="3362",size=6)+
  geom_text(x=-3.1,y=3.3,label="4153",size=6)+
  geom_text(x=-3.1,y=-4.4,label="4009",size=6)+
  geom_text(x=4.3,y=-4.4,label="3345",size=6)

setwd("E:/paper_plot")
ggsave("plot_1h_exon.png", p, width = 10, height = 10, dpi = 300)

#24h exon direction
setwd("E:/exon")
shared_exons <- read.table("shared_exons.txt", header = T, sep = "\t")
Sal1=read.csv("Fe24h_all.csv",header = T)
Sal2=read.csv("Male24h_all.csv",header = T)

names(Sal1)[names(Sal1) == 'X'] <- 'exonname.Sal1'
names(Sal2)[names(Sal2) == 'X'] <- 'exonname.Sal2'
names(Sal1)[names(Sal1) == 'log2fold_12ppt_0ppt'] <- 'log2FoldChange.Sal1'
names(Sal2)[names(Sal2) == 'log2fold_12ppt_0ppt'] <- 'log2FoldChange.Sal2'
Sal1 <- Sal1[, c("exonname.Sal1", "log2FoldChange.Sal1")]
Sal2 <- Sal2[, c("exonname.Sal2", "log2FoldChange.Sal2")]
Sal1 <- subset(Sal1, !grepl("\\+", exonname.Sal1))
Sal2 <- subset(Sal2, !grepl("\\+", exonname.Sal2))
Sal1 <- na.omit(Sal1, cols = "log2FoldChange.Sal1")
Sal2 <- na.omit(Sal2, cols = "log2FoldChange.Sal2")


Sal1_Sal2_Merged <- cbind(Sal1, Sal2)
Sal1_Sal2_Merged <- Sal1_Sal2_Merged[Sal1_Sal2_Merged$exonname.Sal1 %in% shared_exons$x, ]
Sal1_Sal2_Merged <- Sal1_Sal2_Merged[Sal1_Sal2_Merged$exonname.Sal2 %in% shared_exons$x, ]

Sal1_Sal2_Merged$direction[which(Sal1_Sal2_Merged$log2FoldChange.Sal1*Sal1_Sal2_Merged$log2FoldChange.Sal2 >0)] <- "same"
Sal1_Sal2_Merged$direction[which(Sal1_Sal2_Merged$log2FoldChange.Sal1*Sal1_Sal2_Merged$log2FoldChange.Sal2 <0)] <- "opposite"


###count the number of exon in different quadrant
Sal1_Sal2_Merged$quarant=ifelse(Sal1_Sal2_Merged$log2FoldChange.Sal1 >0 & Sal1_Sal2_Merged$log2FoldChange.Sal2 >0, "first",
                                ifelse(Sal1_Sal2_Merged$log2FoldChange.Sal1 <0 & Sal1_Sal2_Merged$log2FoldChange.Sal2 >0, "second",
                                       ifelse(Sal1_Sal2_Merged$log2FoldChange.Sal1 <0 & Sal1_Sal2_Merged$log2FoldChange.Sal2 <0, "third", "fourth")))
nrow(Sal1_Sal2_Merged[Sal1_Sal2_Merged$quarant=="first",])#3705
nrow(Sal1_Sal2_Merged[Sal1_Sal2_Merged$quarant=="second",])#3447
nrow(Sal1_Sal2_Merged[Sal1_Sal2_Merged$quarant=="third",])#3680
nrow(Sal1_Sal2_Merged[Sal1_Sal2_Merged$quarant=="fourth",])#4037

p <- ggplot(Sal1_Sal2_Merged, aes(x=log2FoldChange.Sal1, y=log2FoldChange.Sal2, fill=direction)) +
  geom_hline(yintercept=0, col="grey50") +
  geom_vline(xintercept=0, col="grey50") +
  geom_point(colour = "black",shape = 21,size=3, alpha=1)+
  scale_fill_manual(values=c("white","grey"))+
  xlab(expression("log"["2"]~"FC"~"for salinity-responsive spliced exons in female"))+
  ylab(expression("log"["2"]~"FC"~"for salinity-responsive spliced exons in male"))+
  #geom_abline(intercept=0, slope=1, show.legend = FALSE, lty=2) +
  #geom_smooth(method="lm", se=FALSE,aes(group=1), color="black", show.legend = F,formula = y ~ x) +
  #stat_poly_eq(aes(group=1,label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~~')),formula = y ~ x, parse = T) +
  theme_bw() +
  #guides(shape = guide_legend(override.aes = list(size = 5, alpha=1),
  #fill = guide_legend(override.aes = list(size = 5, alpha=1))))+
  theme(
    strip.background = element_blank(),
    strip.text = element_blank(),
    legend.position = "none",
    aspect.ratio = 1,  # 设置宽高比为1
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    plot.title = element_text(
      size = 22,            # 标题字体大小
      hjust = 0.5,          # 水平居中
      margin = margin(b = 15),# 下边距15点（增加与图的间距）
      face = "bold"
    )
  ) +
  guides(fill=guide_legend(title=NULL)) +
  # 添加"24 h"标题（顶部居中）
  ggtitle("24 h") +  # 关键添加
  geom_text(x=3.3,y=3.1,label="3705",size=6)+
  geom_text(x=-2.8,y=3.1,label="3447",size=6)+
  geom_text(x=-2.8,y=-4,label="3680",size=6)+
  geom_text(x=3.3,y=-4,label="4037",size=6)

setwd("E:/paper_plot")
ggsave("plot_24h_exon.png", p, width = 10, height = 10, dpi = 300)

#48h exon direction
setwd("E:/exon")
shared_exons <- read.table("shared_exons.txt", header = T, sep = "\t")
Sal1=read.csv("Fe48h_all.csv",header = T)
Sal2=read.csv("Male48h_all.csv",header = T)

names(Sal1)[names(Sal1) == 'X'] <- 'exonname.Sal1'
names(Sal2)[names(Sal2) == 'X'] <- 'exonname.Sal2'
names(Sal1)[names(Sal1) == 'log2fold_12ppt_0ppt'] <- 'log2FoldChange.Sal1'
names(Sal2)[names(Sal2) == 'log2fold_12ppt_0ppt'] <- 'log2FoldChange.Sal2'
Sal1 <- Sal1[, c("exonname.Sal1", "log2FoldChange.Sal1")]
Sal2 <- Sal2[, c("exonname.Sal2", "log2FoldChange.Sal2")]
Sal1 <- subset(Sal1, !grepl("\\+", exonname.Sal1))
Sal2 <- subset(Sal2, !grepl("\\+", exonname.Sal2))
Sal1 <- na.omit(Sal1, cols = "log2FoldChange.Sal1")
Sal2 <- na.omit(Sal2, cols = "log2FoldChange.Sal2")


Sal1_Sal2_Merged <- cbind(Sal1, Sal2)
Sal1_Sal2_Merged <- Sal1_Sal2_Merged[Sal1_Sal2_Merged$exonname.Sal1 %in% shared_exons$x, ]
Sal1_Sal2_Merged <- Sal1_Sal2_Merged[Sal1_Sal2_Merged$exonname.Sal2 %in% shared_exons$x, ]

Sal1_Sal2_Merged$direction[which(Sal1_Sal2_Merged$log2FoldChange.Sal1*Sal1_Sal2_Merged$log2FoldChange.Sal2 >0)] <- "same"
Sal1_Sal2_Merged$direction[which(Sal1_Sal2_Merged$log2FoldChange.Sal1*Sal1_Sal2_Merged$log2FoldChange.Sal2 <0)] <- "opposite"


###count the number of exon in different quadrant
Sal1_Sal2_Merged$quarant=ifelse(Sal1_Sal2_Merged$log2FoldChange.Sal1 >0 & Sal1_Sal2_Merged$log2FoldChange.Sal2 >0, "first",
                                ifelse(Sal1_Sal2_Merged$log2FoldChange.Sal1 <0 & Sal1_Sal2_Merged$log2FoldChange.Sal2 >0, "second",
                                       ifelse(Sal1_Sal2_Merged$log2FoldChange.Sal1 <0 & Sal1_Sal2_Merged$log2FoldChange.Sal2 <0, "third", "fourth")))
nrow(Sal1_Sal2_Merged[Sal1_Sal2_Merged$quarant=="first",])#4033
nrow(Sal1_Sal2_Merged[Sal1_Sal2_Merged$quarant=="second",])#3603
nrow(Sal1_Sal2_Merged[Sal1_Sal2_Merged$quarant=="third",])#3363
nrow(Sal1_Sal2_Merged[Sal1_Sal2_Merged$quarant=="fourth",])#3870

p <- ggplot(Sal1_Sal2_Merged, aes(x=log2FoldChange.Sal1, y=log2FoldChange.Sal2, fill=direction)) +
  geom_hline(yintercept=0, col="grey50") +
  geom_vline(xintercept=0, col="grey50") +
  geom_point(colour = "black",shape = 21,size=3, alpha=1)+
  scale_fill_manual(values=c("white","grey"))+
  xlab(expression("log"["2"]~"FC"~"for salinity-responsive spliced exons in female"))+
  ylab(expression("log"["2"]~"FC"~"for salinity-responsive spliced exons in male"))+
  #geom_abline(intercept=0, slope=1, show.legend = FALSE, lty=2) +
  #geom_smooth(method="lm", se=FALSE,aes(group=1), color="black", show.legend = F,formula = y ~ x) +
  #stat_poly_eq(aes(group=1,label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~~')),formula = y ~ x, parse = T) +
  theme_bw() +
  #guides(shape = guide_legend(override.aes = list(size = 5, alpha=1),
  #fill = guide_legend(override.aes = list(size = 5, alpha=1))))+
  theme(
    strip.background = element_blank(),
    strip.text = element_blank(),
    legend.position = "none",
    aspect.ratio = 1,  # 设置宽高比为1
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    plot.title = element_text(
      size = 22,            # 标题字体大小
      hjust = 0.5,          # 水平居中
      margin = margin(b = 15),# 下边距15点（增加与图的间距）
      face = "bold"
    )
  ) +
  guides(fill=guide_legend(title=NULL)) +
  # 添加"48 h"标题（顶部居中）
  ggtitle("48 h") +  # 关键添加
  geom_text(x=4.15,y=4.9,label="4033",size=6)+
  geom_text(x=-3.3,y=4.9,label="3603",size=6)+
  geom_text(x=-3.3,y=-2.9,label="3363",size=6)+
  geom_text(x=4.15,y=-2.9,label="3870",size=6)

setwd("E:/paper_plot")
ggsave("plot_48h_exon.png", p, width = 10, height = 10, dpi = 300)

#72h exon direction
setwd("E:/exon")
shared_exons <- read.table("shared_exons.txt", header = T, sep = "\t")
Sal1=read.csv("Fe72h_all.csv",header = T)
Sal2=read.csv("Male72h_all.csv",header = T)

names(Sal1)[names(Sal1) == 'X'] <- 'exonname.Sal1'
names(Sal2)[names(Sal2) == 'X'] <- 'exonname.Sal2'
names(Sal1)[names(Sal1) == 'log2fold_12ppt_0ppt'] <- 'log2FoldChange.Sal1'
names(Sal2)[names(Sal2) == 'log2fold_12ppt_0ppt'] <- 'log2FoldChange.Sal2'
Sal1 <- Sal1[, c("exonname.Sal1", "log2FoldChange.Sal1")]
Sal2 <- Sal2[, c("exonname.Sal2", "log2FoldChange.Sal2")]
Sal1 <- subset(Sal1, !grepl("\\+", exonname.Sal1))
Sal2 <- subset(Sal2, !grepl("\\+", exonname.Sal2))
Sal1 <- na.omit(Sal1, cols = "log2FoldChange.Sal1")
Sal2 <- na.omit(Sal2, cols = "log2FoldChange.Sal2")

###
rownames(Sal1) <- Sal1$exonname.Sal1
rownames(Sal2) <- Sal2$exonname.Sal2
common_rownames <- intersect(rownames(Sal1), rownames(Sal2))
Sal1_common <- Sal1[rownames(Sal1) %in% common_rownames, ]
Sal2_common <- Sal2[rownames(Sal2) %in% common_rownames, ]
rownames(Sal1_common) <- NULL
rownames(Sal2_common) <- NULL
Sal1 <- Sal1_common
Sal2 <- Sal2_common
###

Sal1_Sal2_Merged <- cbind(Sal1, Sal2)
Sal1_Sal2_Merged <- Sal1_Sal2_Merged[Sal1_Sal2_Merged$exonname.Sal1 %in% shared_exons$x, ]
Sal1_Sal2_Merged <- Sal1_Sal2_Merged[Sal1_Sal2_Merged$exonname.Sal2 %in% shared_exons$x, ]

Sal1_Sal2_Merged$direction[which(Sal1_Sal2_Merged$log2FoldChange.Sal1*Sal1_Sal2_Merged$log2FoldChange.Sal2 >0)] <- "same"
Sal1_Sal2_Merged$direction[which(Sal1_Sal2_Merged$log2FoldChange.Sal1*Sal1_Sal2_Merged$log2FoldChange.Sal2 <0)] <- "opposite"


###count the number of exon in different quadrant
Sal1_Sal2_Merged$quarant=ifelse(Sal1_Sal2_Merged$log2FoldChange.Sal1 >0 & Sal1_Sal2_Merged$log2FoldChange.Sal2 >0, "first",
                                ifelse(Sal1_Sal2_Merged$log2FoldChange.Sal1 <0 & Sal1_Sal2_Merged$log2FoldChange.Sal2 >0, "second",
                                       ifelse(Sal1_Sal2_Merged$log2FoldChange.Sal1 <0 & Sal1_Sal2_Merged$log2FoldChange.Sal2 <0, "third", "fourth")))
nrow(Sal1_Sal2_Merged[Sal1_Sal2_Merged$quarant=="first",])#4391
nrow(Sal1_Sal2_Merged[Sal1_Sal2_Merged$quarant=="second",])#4525
nrow(Sal1_Sal2_Merged[Sal1_Sal2_Merged$quarant=="third",])#2872
nrow(Sal1_Sal2_Merged[Sal1_Sal2_Merged$quarant=="fourth",])#3081

p <- ggplot(Sal1_Sal2_Merged, aes(x=log2FoldChange.Sal1, y=log2FoldChange.Sal2, fill=direction)) +
  geom_hline(yintercept=0, col="grey50") +
  geom_vline(xintercept=0, col="grey50") +
  geom_point(colour = "black",shape = 21,size=3, alpha=1)+
  scale_fill_manual(values=c("white","grey"))+
  xlab(expression("log"["2"]~"FC"~"for salinity-responsive spliced exons in female"))+
  ylab(expression("log"["2"]~"FC"~"for salinity-responsive spliced exons in male"))+
  #geom_abline(intercept=0, slope=1, show.legend = FALSE, lty=2) +
  #geom_smooth(method="lm", se=FALSE,aes(group=1), color="black", show.legend = F,formula = y ~ x) +
  #stat_poly_eq(aes(group=1,label = paste(..eq.label.., ..adj.rr.label.., sep = '~~~~')),formula = y ~ x, parse = T) +
  theme_bw() +
  #guides(shape = guide_legend(override.aes = list(size = 5, alpha=1),
  #fill = guide_legend(override.aes = list(size = 5, alpha=1))))+
  theme(
    strip.background = element_blank(),
    strip.text = element_blank(),
    legend.position = "none",
    aspect.ratio = 1,  # 设置宽高比为1
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    plot.title = element_text(
      size = 22,            # 标题字体大小
      hjust = 0.5,          # 水平居中
      margin = margin(b = 15),# 下边距15点（增加与图的间距）
      face = "bold"
    )
  ) +
  guides(fill=guide_legend(title=NULL)) +
  # 添加"72 h"标题（顶部居中）
  ggtitle("72 h") +  # 关键添加
  geom_text(x=4.1,y=5.05,label="4391",size=6)+
  geom_text(x=-2.15,y=5.05,label="4525",size=6)+
  geom_text(x=-2.15,y=-3.05,label="2872",size=6)+
  geom_text(x=4.1,y=-3.05,label="3081",size=6)

setwd("E:/paper_plot")
ggsave("plot_72h_exon.png", p, width = 10, height = 10, dpi = 300)


#### 10 G-test ####
#1h
library(RVAideMemoire)

observed = c(7371, 7498)    # observed frequencies
expected = c(0.5, 0.5)      # expected proportions

G.test(x=observed,
       p=expected)
#G = 1.0848, df = 1, p-value = 0.2976


#24h
observed = c(7385, 7484)    # observed frequencies
expected = c(0.5, 0.5)      # expected proportions

G.test(x=observed,
       p=expected)
#G = 0.65916, df = 1, p-value = 0.4169


#48h
observed = c(7396, 7473)    # observed frequencies
expected = c(0.5, 0.5)      # expected proportions

G.test(x=observed,
       p=expected)
#G = 0.39875, df = 1, p-value = 0.5277


#72h
observed = c(7263, 7606)    # observed frequencies
expected = c(0.5, 0.5)      # expected proportions

G.test(x=observed,
       p=expected)
#G = 7.9131, df = 1, p-value = 0.004908


### 11 boxplot of exon ###
#boxplot_1h
library(ggplot2)
library(ggpubr)
library(Rmisc)

Logfc_female_1h <- as.data.frame(cbind(Sal1_Sal2_Merged[["exonname.Sal1"]], Sal1_Sal2_Merged[["log2FoldChange.Sal1"]]))
Logfc_male_1h <- as.data.frame(cbind(Sal1_Sal2_Merged[["exonname.Sal2"]], Sal1_Sal2_Merged[["log2FoldChange.Sal2"]]))

names(Logfc_female_1h) <- c("exonname", "log2FoldChange")
names(Logfc_male_1h) <- c("exonname", "log2FoldChange")

#
library(gridExtra)
library(gapminder)
library(dplyr)
library(tidyr)

Logfc_female_1h$Sex="Female"
Logfc_male_1h$Sex="Male"

dat_my=rbind(Logfc_female_1h,Logfc_male_1h)

#logFC取绝对值：变化的倍数
dat_test=dat_my
dat_test$log2FoldChange <- as.numeric(dat_test$log2FoldChange)
dat_test$log2FoldChange=abs(dat_test$log2FoldChange)

df_summary <-  dat_test %>%
  group_by(Sex) %>%
  summarise(mean_value = mean(log2FoldChange, na.rm = TRUE))#平均值 F:0.2548421, M:0.2481658


scaleFUN <- function(x) sprintf("%.1f", x)
plas_1=ggplot(dat_test, aes(x = Sex, y = log2FoldChange)) +
  #geom_boxplot() +
  stat_boxplot(geom="errorbar",position=position_dodge(width=0.2),width=0.1)+
  geom_boxplot(position=position_dodge(width =0.2),width=0.4)+#,outlier.shape=NA
  #facet_wrap(~variable, scale = "free") +
  #scale_fill_manual(values = c("#d6604d","#217db4"))+
  #scale_y_continuous(name = "Log2Fold change") +
  #scale_y_continuous(labels=scaleFUN,name="Absolute log2 fold change")+
  ylab(expression("|Log"["2"]~"Fold Change|"))+
  scale_y_continuous(labels = scaleFUN)+
  #scale_x_discrete(labels = abbreviate, name = "Sex")+
  stat_compare_means(aes(group = Sex),
                     size = 3,#method = "t.test"，
                     #digits=2,
                     method = "wilcox.test",
                     method.args = list(paired = FALSE),
                     hide.ns = F,
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), 
                                      symbols = c("****", "***", "**", "*", "NS")))+
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x=element_text(vjust=1,size=13),
        axis.text.y=element_text(vjust=1,size=13),
        axis.title.y=element_text(vjust=1,size=14),
        axis.title.x=element_text(vjust=1,size=13),
        plot.title = element_text(
          size = 16,
          hjust = 0.5,
          face = "bold",
          margin = margin(b = 10)
        )) +  # 结束theme()
  xlab(NULL)+
  # 添加"1 h"标题
  ggtitle("1 h")

wilcox_result <- wilcox.test(dat_test$log2FoldChange ~ dat_test$Sex, paired = FALSE)
print(wilcox_result)#W = 112716857, p-value = 0.003324

ggsave("E:/paper_plot/boxplot_1h_exon.png", device="png",
       #path = "xxx",
       width = 8, height = 4, units="in",#scaling = 0.7,
       dpi=600)
#------------------------------#
#boxplot_24h
Logfc_female_24h <- as.data.frame(cbind(Sal1_Sal2_Merged[["exonname.Sal1"]], Sal1_Sal2_Merged[["log2FoldChange.Sal1"]]))
Logfc_male_24h <- as.data.frame(cbind(Sal1_Sal2_Merged[["exonname.Sal2"]], Sal1_Sal2_Merged[["log2FoldChange.Sal2"]]))

names(Logfc_female_24h) <- c("exonname", "log2FoldChange")
names(Logfc_male_24h) <- c("exonname", "log2FoldChange")

#
Logfc_female_24h$Sex="Female"
Logfc_male_24h$Sex="Male"

dat_my=rbind(Logfc_female_24h,Logfc_male_24h)

#logFC取绝对值：变化的倍数
dat_test=dat_my
dat_test$log2FoldChange <- as.numeric(dat_test$log2FoldChange)
dat_test$log2FoldChange=abs(dat_test$log2FoldChange)

df_summary <-  dat_test %>%
  group_by(Sex) %>%
  summarise(mean_value = mean(log2FoldChange, na.rm = TRUE))#平均值 F:0.2608730, M:0.2642746


scaleFUN <- function(x) sprintf("%.1f", x)
plas_1=ggplot(dat_test, aes(x = Sex, y = log2FoldChange)) +
  #geom_boxplot() +
  stat_boxplot(geom="errorbar",position=position_dodge(width=0.2),width=0.1)+
  geom_boxplot(position=position_dodge(width =0.2),width=0.4)+#,outlier.shape=NA
  #facet_wrap(~variable, scale = "free") +
  #scale_fill_manual(values = c("#d6604d","#217db4"))+
  #scale_y_continuous(name = "Log2Fold change") +
  #scale_y_continuous(labels=scaleFUN,name="Absolute log2 fold change")+
  ylab(expression("|Log"["2"]~"Fold Change|"))+
  scale_y_continuous(labels = scaleFUN)+
  #scale_x_discrete(labels = abbreviate, name = "Sex")+
  stat_compare_means(aes(group = Sex),
                     size = 3,#method = "t.test"，
                     #digits=2,
                     method = "wilcox.test",
                     method.args = list(paired = FALSE),
                     hide.ns = F,
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), 
                                      symbols = c("****", "***", "**", "*", "NS")))+
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x=element_text(vjust=1,size=13),
        axis.text.y=element_text(vjust=1,size=13),
        axis.title.y=element_text(vjust=1,size=14),
        axis.title.x=element_text(vjust=1,size=13),
        plot.title = element_text(
          size = 16,
          hjust = 0.5,
          face = "bold",
          margin = margin(b = 10)
        )) +  # 结束theme()
  xlab(NULL)+
  # 添加"24 h"标题
  ggtitle("24 h")

wilcox_result <- wilcox.test(dat_test$log2FoldChange ~ dat_test$Sex, paired = FALSE)
print(wilcox_result)#W = 109938500, p-value = 0.4137

ggsave("E:/paper_plot/boxplot_24h_exon.png", device="png",
       #path = "xxx",
       width = 8, height = 4, units="in",#scaling = 0.7,
       dpi=600)
#------------------------------#
#boxplot_48h
Logfc_female_48h <- as.data.frame(cbind(Sal1_Sal2_Merged[["exonname.Sal1"]], Sal1_Sal2_Merged[["log2FoldChange.Sal1"]]))
Logfc_male_48h <- as.data.frame(cbind(Sal1_Sal2_Merged[["exonname.Sal2"]], Sal1_Sal2_Merged[["log2FoldChange.Sal2"]]))

names(Logfc_female_48h) <- c("exonname", "log2FoldChange")
names(Logfc_male_48h) <- c("exonname", "log2FoldChange")

#
Logfc_female_48h$Sex="Female"
Logfc_male_48h$Sex="Male"

dat_my=rbind(Logfc_female_48h,Logfc_male_48h)

#logFC取绝对值：变化的倍数
dat_test=dat_my
dat_test$log2FoldChange <- as.numeric(dat_test$log2FoldChange)
dat_test$log2FoldChange=abs(dat_test$log2FoldChange)

df_summary <-  dat_test %>%
  group_by(Sex) %>%
  summarise(mean_value = mean(log2FoldChange, na.rm = TRUE))#平均值 F:0.1998416, M:0.2684339


scaleFUN <- function(x) sprintf("%.1f", x)
plas_1=ggplot(dat_test, aes(x = Sex, y = log2FoldChange)) +
  #geom_boxplot() +
  stat_boxplot(geom="errorbar",position=position_dodge(width=0.2),width=0.1)+
  geom_boxplot(position=position_dodge(width =0.2),width=0.4)+#,outlier.shape=NA
  #facet_wrap(~variable, scale = "free") +
  #scale_fill_manual(values = c("#d6604d","#217db4"))+
  #scale_y_continuous(name = "Log2Fold change") +
  #scale_y_continuous(labels=scaleFUN,name="Absolute log2 fold change")+
  ylab(expression("|Log"["2"]~"Fold Change|"))+
  scale_y_continuous(labels = scaleFUN)+
  #scale_x_discrete(labels = abbreviate, name = "Sex")+
  stat_compare_means(aes(group = Sex),
                     size = 3,#method = "t.test"，
                     #digits=2,
                     method = "wilcox.test",
                     method.args = list(paired = FALSE),
                     hide.ns = F,
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), 
                                      symbols = c("****", "***", "**", "*", "NS")))+
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x=element_text(vjust=1,size=13),
        axis.text.y=element_text(vjust=1,size=13),
        axis.title.y=element_text(vjust=1,size=14),
        axis.title.x=element_text(vjust=1,size=13),
        plot.title = element_text(
          size = 16,
          hjust = 0.5,
          face = "bold",
          margin = margin(b = 10)
        )) +  # 结束theme()
  xlab(NULL)+
  # 添加"48 h"标题
  ggtitle("48 h")


wilcox_result <- wilcox.test(dat_test$log2FoldChange ~ dat_test$Sex, paired = FALSE)
print(wilcox_result)#W = 91984315, p-value < 2.2e-16

ggsave("E:/paper_plot/boxplot_48h_exon.png", device="png",
       #path = "xxx",
       width = 8, height = 4, units="in",#scaling = 0.7,
       dpi=600)
#------------------------------#
#boxplot_72h
Logfc_female_72h <- as.data.frame(cbind(Sal1_Sal2_Merged[["exonname.Sal1"]], Sal1_Sal2_Merged[["log2FoldChange.Sal1"]]))
Logfc_male_72h <- as.data.frame(cbind(Sal1_Sal2_Merged[["exonname.Sal2"]], Sal1_Sal2_Merged[["log2FoldChange.Sal2"]]))

names(Logfc_female_72h) <- c("exonname", "log2FoldChange")
names(Logfc_male_72h) <- c("exonname", "log2FoldChange")

#
Logfc_female_72h$Sex="Female"
Logfc_male_72h$Sex="Male"

dat_my=rbind(Logfc_female_72h,Logfc_male_72h)

#logFC取绝对值：变化的倍数
dat_test=dat_my
dat_test$log2FoldChange <- as.numeric(dat_test$log2FoldChange)
dat_test$log2FoldChange=abs(dat_test$log2FoldChange)

df_summary <-  dat_test %>%
  group_by(Sex) %>%
  summarise(mean_value = mean(log2FoldChange, na.rm = TRUE))#平均值 F:0.2226624, M:0.2992701


scaleFUN <- function(x) sprintf("%.1f", x)
plas_1=ggplot(dat_test, aes(x = Sex, y = log2FoldChange)) +
  #geom_boxplot() +
  stat_boxplot(geom="errorbar",position=position_dodge(width=0.2),width=0.1)+
  geom_boxplot(position=position_dodge(width =0.2),width=0.4)+#,outlier.shape=NA
  #facet_wrap(~variable, scale = "free") +
  #scale_fill_manual(values = c("#d6604d","#217db4"))+
  #scale_y_continuous(name = "Log2Fold change") +
  #scale_y_continuous(labels=scaleFUN,name="Absolute log2 fold change")+
  ylab(expression("|Log"["2"]~"Fold Change|"))+
  scale_y_continuous(labels = scaleFUN)+
  #scale_x_discrete(labels = abbreviate, name = "Sex")+
  stat_compare_means(aes(group = Sex),
                     size = 3,#method = "t.test"，
                     #digits=2,
                     method = "wilcox.test",
                     method.args = list(paired = FALSE),
                     hide.ns = F,
                     symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), 
                                      symbols = c("****", "***", "**", "*", "NS")))+
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x=element_text(vjust=1,size=13),
        axis.text.y=element_text(vjust=1,size=13),
        axis.title.y=element_text(vjust=1,size=14),
        axis.title.x=element_text(vjust=1,size=13),
        plot.title = element_text(
          size = 16,
          hjust = 0.5,
          face = "bold",
          margin = margin(b = 10)
        )) +  # 结束theme()
  xlab(NULL)+
  # 添加"72 h"标题
  ggtitle("72 h")

wilcox_result <- wilcox.test(dat_test$log2FoldChange ~ dat_test$Sex, paired = FALSE)
print(wilcox_result)#W = 92938244, p-value < 2.2e-16

ggsave("E:/paper_plot/boxplot_72h_exon.png", device="png",
       #path = "xxx",
       width = 8, height = 4, units="in",#scaling = 0.7,
       dpi=600)


# 12 parallel DEUs and parallel DE genes####
# 12.1 parallel DEUs
setwd("E:/exon")
#1h
parallelDEU_1h=Sal1_Sal2_Merged[which(Sal1_Sal2_Merged$direction=="same"),]
#24h
parallelDEU_24h=Sal1_Sal2_Merged[which(Sal1_Sal2_Merged$direction=="same"),]
#48h
parallelDEU_48h=Sal1_Sal2_Merged[which(Sal1_Sal2_Merged$direction=="same"),]
#72h
parallelDEU_72h=Sal1_Sal2_Merged[which(Sal1_Sal2_Merged$direction=="same"),]

rownames(parallelDEU_1h) <- parallelDEU_1h$exonname.Sal1
rownames(parallelDEU_24h) <- parallelDEU_24h$exonname.Sal1
rownames(parallelDEU_48h) <- parallelDEU_48h$exonname.Sal1
rownames(parallelDEU_72h) <- parallelDEU_72h$exonname.Sal1

common_exonnames <- Reduce(intersect, list(rownames(parallelDEU_1h), rownames(parallelDEU_24h), rownames(parallelDEU_48h), rownames(parallelDEU_72h)))
write.csv(common_exonnames, "common_exonnames.csv", row.names = F)

common_exonnames=read.csv("common_exonnames.csv",header = T)
common_DEU_genenames <- common_exonnames %>%
  separate(x, into = c("genename", "exonpart"), sep = ":", extra = "drop")

genename_DEU <- common_DEU_genenames$genename
genename_DEU
genename_DEU_unique <- unique(genename_DEU)
write.csv(genename_DEU_unique, "common_parallel_splicing_genenames.csv", row.names = F)

# 12.2 parallel DEGs
setwd("E:/filter")
#1h
parallelDEG_1h=Sal1_Sal2_Merged[which(Sal1_Sal2_Merged$direction=="same"),]
#24h
parallelDEG_24h=Sal1_Sal2_Merged[which(Sal1_Sal2_Merged$direction=="same"),]
#48h
parallelDEG_48h=Sal1_Sal2_Merged[which(Sal1_Sal2_Merged$direction=="same"),]
#72h
parallelDEG_72h=Sal1_Sal2_Merged[which(Sal1_Sal2_Merged$direction=="same"),]

rownames(parallelDEG_1h) <- parallelDEG_1h$genename.Sal1
rownames(parallelDEG_24h) <- parallelDEG_24h$genename.Sal1
rownames(parallelDEG_48h) <- parallelDEG_48h$genename.Sal1
rownames(parallelDEG_72h) <- parallelDEG_72h$genename.Sal1

common_genenames <- Reduce(intersect, list(rownames(parallelDEG_1h), rownames(parallelDEG_24h), rownames(parallelDEG_48h), rownames(parallelDEG_72h)))
write.csv(common_genenames, "common_genenames.csv", row.names = F)

common_genenames=read.csv("common_genenames.csv",header = T)
genename_DEG <- common_genenames$x
genename_DEG
genename_DEG_unique <- unique(genename_DEG)
genename_DEG_unique
write.csv(genename_DEG_unique, "common_parallel_express_genenames.csv", row.names = F)
