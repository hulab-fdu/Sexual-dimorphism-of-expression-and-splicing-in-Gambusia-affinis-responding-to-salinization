library(reshape2)
library(DESeq2)
library(grid)
library(genefilter)
library(pheatmap)
library(VennDiagram)
library(ggplot2)
library("RColorBrewer")
library(EnhancedVolcano)
library(tidyverse)
library(ggpubr)
library(gridExtra)
library(yarrr)
library(plyr)
library(goseq)
library(clusterProfiler)
library(tximport)
library(GenomicFeatures)
library(ggfortify)
library(rgl)
library(car)
library(plot3D)
library(vegan)
library(ggrepel)
library(dcast)
library(reshape2)
library(vegan)

####1 package:reshape2:count_table and filter ####
setwd("E:/filter")
allcount=read.csv('all_sample.txt',header=F,sep="\t")
colnames(allcount)=c('sample','gene.ID','count')
allcount=dcast(allcount,formula=gene.ID~sample,value.var="count") #gene.ID as row-,sample as col-,count as the variate
rownames(allcount)=allcount$gene.ID
noint = rownames(allcount) %in% c("__alignment_not_unique","__ambiguous","__no_feature","__not_aligned","__too_low_aQual")
head(noint,6)
count_table=allcount[!noint,]
count_table=count_table[,-1]
rm=grep("G",colnames(count_table))
count_table=count_table[,-rm]
write.table(count_table, file="counts.txt",sep="\t",quote=FALSE,row.names=T)


sampleTable=read.csv("Sample.csv",header=T,sep = ",",check.names = F)
rownames(sampleTable)=sampleTable$Sample
colnames(count_table)==rownames(sampleTable) #check if the count table has the same order as sample table
str(sampleTable)
id=colnames(count_table)
id=data.frame(id)
sampleTable=sampleTable[id$id,]
colnames(count_table)==rownames(sampleTable)


#count_filter: gene >1 at all samples
count_table=subset(count_table,apply(count_table,1, function(x) sum(x >= 1))>= 32)
write.table(count_table, file="filter_table.txt",sep="\t",quote=FALSE,row.names=T)


####1.1 all_PCA####
all=DESeqDataSetFromMatrix(countData=count_table, colData=sampleTable, design = ~Salinity+Time+Sex)
all$Salinity=relevel(all$Salinity, ref ="0ppt")
all$Salinity
vst_all=vst(all,blind=F)
write.table(vst_all@assays@data@listData[[1]], file="vst_count_table.txt",sep="\t",quote=FALSE,row.names=T)

#PCA based on the prcomp() function
pcadata <- vst_all@assays@data@listData[[1]]
pcadata <- t(pcadata)
id <- rownames(pcadata)
id <- as.vector(id)
sample_ids <- data.frame(Sample = sampleTable$Sample)
sample_indices <- match(sample_ids$Sample, id)
pcadata <- pcadata[sample_indices, ]
pcadata=prcomp(pcadata)
summary(pcadata)#pc1~3
write.table(pcadata[["x"]], file="pcadata_gene.txt",sep="\t",quote=FALSE,row.names=T)


rda=pcadata$x[,1:3]
write.csv(rda,file="rda_all.csv")
genet=read.table("rda_all.csv",header=TRUE,check.names=FALSE,row.names = 1,sep=",")
design_test=read.table("Sample.csv",sep=",",header = TRUE,check.names = T,row.names = 1)
genet=genet[rownames(design_test),]

rda.vpa=varpart(genet[,1:3],~Salinity,~Time,~Sex,data=design_test)
str(design_test)
rda.vpa
str(rda.vpa)
summary(rda.vpa)


#rda
rda=rda(genet[,1:3]~+Salinity+Time+Sex,data=design_test)
summary(rda)
RsquareAdj(rda)#0.1445525
anova(rda)#2.0477  0.032 *
plot(rda)

### look the sig###
#salinity
rda_sal=rda(genet[,1:3]~Salinity,data=design_test)
summary(rda_sal)
RsquareAdj(rda_sal)#0.005785041
anova(rda_sal)#1.1804  0.276
plot(rda_sal)
rda_sal_partial=rda(genet[,1:3]~Salinity+Condition(Time+Sex),data=design_test)
summary(rda_sal_partial)
RsquareAdj(rda_sal_partial)#0.01178162
anova(rda_sal_partial)#1.3719  0.246
plot(rda_sal_partial)

#time
rda_time=rda(genet[,1:3]~Time,data=design_test)
summary(rda_time)
RsquareAdj(rda_time)#0.1194538
anova(rda_time)#2.4018  0.026 *
plot(rda_time)
rda_time_partial=rda(genet[,1:3]~Time+Condition(Salinity+Sex),data=design_test)
summary(rda_time_partial)
RsquareAdj(rda_time_partial)#0.1302884
anova(rda_time_partial)#2.4723  0.025 *
plot(rda_time_partial)

#sex
rda_sex=rda(genet[,1:3]~Sex,data=design_test)
summary(rda_sex)
RsquareAdj(rda_sex)#0.008003572
anova(rda_sex)#1.2501  0.297
plot(rda_sex)
rda_sex_partial=rda(genet[,1:3]~Sex+Condition(Salinity+Time),data=design_test)
summary(rda_sex_partial)
RsquareAdj(rda_sex_partial)#0.01424666
anova(rda_sex_partial)#1.4497  0.229
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
#title_lab=expression(bold(paste(" adj. ",R^{2}," =14.46%;",italic(p),"-value"," = 0.032"))) # fill with you own stats in sommaire

rda1 =round(rda$CCA$eig[1]/sum(rda$CCA$eig)*100,2) # 65.5
rda2 =round(rda$CCA$eig[2]/sum(rda$CCA$eig)*100,2) # 25.4
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

a=ggplot() +
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
  labs(x="RDA1 (65.5%)",y="RDA2 (25.4%)") +
  geom_text_repel(data=df2subset_plot,aes(x=RDA1, y=RDA2,label=c("Sex","Time*","Salinity")),size=6)+
  theme(legend.spacing.x = unit(0.5, 'cm'),legend.spacing.y = unit(0, 'cm'))


ggsave("./RDA_DE_rename_pca_1.png", device="png",
       #path = "xxx",
       width = 6.2, height = 6.2, units="in",#scaling = 0.7,
       dpi=300)
