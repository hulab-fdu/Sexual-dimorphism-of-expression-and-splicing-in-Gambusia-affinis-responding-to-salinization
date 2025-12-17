#WGCNA_exon_female
setwd("E:/exon")
library(DEXSeq)
library(reshape2)
library(DESeq2)
library(WGCNA)
library(dbplyr)


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
colnames(filtered_count)
rm=grep("M", colnames(filtered_count))
filtered_count_female=filtered_count[,-rm]


sampleTable=read.csv("sampleTable_F.csv",header=T,sep = ",",check.names = F)
rownames(sampleTable)=sampleTable$sampleID
colnames(filtered_count_female)==rownames(sampleTable)#check if the count table has the same order as sample table
str(sampleTable)
id=colnames(filtered_count_female)
id=data.frame(id)
sampleTable=sampleTable[id$id,]
colnames(filtered_count_female)==rownames(sampleTable)

sampleTable$salinity <- factor(sampleTable$salinity)
rlog_count=rlogTransformation(filtered_count_female,blind=F)
head(rlog_count,3)
rlog_count_data=as.data.frame(rlog_count)
#-----------------#


####1 data loading ####
options(stringsAsFactors = FALSE)
enableWGCNAThreads()#multiple thread
disableWGCNAThreads()
datExpr_count= as.data.frame(t(rlog_count_data))#row:sample,col:exon
gsg = goodSamplesGenes(datExpr_count, verbose = 3)
gsg$allOK#TRUE


sampleTree = hclust(dist(datExpr_count,method = "euclidean"), method = "ward.D2") 
sizeGrWindow(12,9) 
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)


#trait data loading
rownames(datExpr_count)==rownames(sampleTable)#TRUE


sampleTree2 = hclust(dist(datExpr_count), method = "average")
traitColors = numbers2colors(as.numeric(sampleTable$salinity), signed = TRUE) 
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(sampleTable),
                    main = "Sample dendrogram and trait heatmap")

#data save
save(datExpr_count, sampleTable, file = "wgcna_female_dataInput.RData")


####network construct####
#soft-power selection
gc()
memory.limit(999999999)
powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft = pickSoftThreshold(datExpr_count, powerVector = powers, verbose = 5, networkType="signed", allowWGCNAThreads(nThreads = 6))
save(sft, file = "sft_female.RData")
sizeGrWindow(12,12)
par(mfrow = c(1,2))
cex1 = 0.9
#Soft Threshold
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
      xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
      main = paste("Scale independence"));
 text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
      labels=powers,cex=cex1,col="red");
 abline(h=0.90,col="red")
#Mean connectivity
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

sft$powerEstimate # 16
gc()

#one_step
softPower = sft$powerEstimate
net=blockwiseModules(datExpr_count, power = softPower, networkType = "signed",
                     TOMType = "signed", minModuleSize = 30, maxBlockSize = 5000,
                     reassignThreshold = 0, mergeCutHeight = 0.25, numericLabels = TRUE,
                     pamRespectsDendro = FALSE, saveTOMs = TRUE, saveTOMFileBase = "methAfterFilterTOM0.15",
                     verbose = 3)
## time consuming; save result to avoid running it again!
save(net, file = "net_female_exon.RData")


# get module label for each gene ('0' corresponding to no
# module) and convert label to a module color ('grey'
# corresponding to no module)
moduleLabels=net$colors
moduleColors=labels2colors(net$colors)
## number of sites without module (assigned to grey module):
length(moduleLabels[moduleLabels == 0])#763
# get genes per module #table(net$colors)
ExonsperModule <- NULL
for (c in seq(1, length(levels(as.factor(moduleColors))), by = 1)) {
  color <- levels(as.factor(moduleColors))[c]
  length <- length(moduleColors[moduleColors == color])
  module <- data.frame(Module = color, NoGenes = length)
  ExonsperModule <- rbind(ExonsperModule, module)
}


# Calculate eigengenes
MEList = moduleEigengenes(datExpr_count, colors = moduleColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")


#Quantify module similarity by eigengene correlation. Eigengenes: Module representatives
sizeGrWindow(6, 6)
par(cex = 1)
png("Eigengene adjacency heatmap.png",units="in", width=10, height=10,res=400)
plotEigengeneNetworks(orderMEs(MEs), "", plotAdjacency = T, colorLabels = F,
                      marDendro = c(0, 4, 2, 4), marHeatmap = c(3, 4, 1, 1.5),
                      plotDendrograms = T, xLabelsAngle = 90)
dev.off()
# get correlation between modules
cor(orderMEs(MEs))


# quantify module-trait associations
# Define numbers of exons and samples
nExons = ncol(datExpr_count)
nSamples = nrow(datExpr_count)

# Recalculate MEs with color labels
MEs = orderMEs(MEs)

sampleTable=read.csv("sampleTable_F.csv",header=T,sep = ",",check.names = F)
rownames(sampleTable)=sampleTable$sampleID
str(sampleTable)
design_female=sampleTable

rownames(design_female)=design_female$sampleID
design_female$sampleID=NULL

design_female

colnames(MEs) <- sub("^ME", "", colnames(MEs))
moduleTraitCor = cor(MEs, design_female, use = "p");

moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

Sig.module_female=as.data.frame(moduleTraitPvalue)[which(as.data.frame(moduleTraitPvalue)$'salinity'<0.05|
                                                         as.data.frame(moduleTraitPvalue)$'salinity:1 h'<0.05|
                                                         as.data.frame(moduleTraitPvalue)$'salinity:24 h'<0.05|
                                                         as.data.frame(moduleTraitPvalue)$'salinity:48 h'<0.05|
                                                         as.data.frame(moduleTraitPvalue)$'salinity:72 h'<0.05),]
rownames(Sig.module_female)#140

Sig.moduleTraitCor <- moduleTraitCor[row.names(moduleTraitCor) %in% row.names(Sig.module_female), ]
Sig.moduleTraitPvalue <- moduleTraitPvalue[row.names(moduleTraitPvalue) %in% row.names(Sig.module_female), ]


# Will display correlations and their p-values
# 创建一个矩阵，用于存储最终的注释文本
textMatrix <- matrix("-", nrow = nrow(Sig.moduleTraitPvalue), ncol = ncol(Sig.moduleTraitPvalue))

# 找出 Sig.moduleTraitPvalue 中小于 0.05 的位置
significant_positions <- which(Sig.moduleTraitPvalue < 0.05, arr.ind = TRUE)

# 对于小于 0.05 的值，保留两位有效数字并存储到 textMatrix 中
textMatrix[significant_positions] <- signif(Sig.moduleTraitPvalue[significant_positions], 2)

# 确保 textMatrix 的维度与 Sig.moduleTraitPvalue 一致
dim(textMatrix) <- dim(Sig.moduleTraitPvalue)

# Display the correlation values within a heatmap plot
png("WGCNA_female_heatmap1.png", width = 3000, height = 6000, res = 300)
par(mar = c(6, 15, 3, 3));
labeledHeatmap(Matrix = Sig.moduleTraitCor,
               xLabels = names(as.data.frame(design_female)),
               yLabels = rownames(Sig.moduleTraitCor),
               ySymbols = rownames(Sig.moduleTraitCor),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()


# $ Choose interesting modules for downstream GO analysis
setwd("E:/exon/female")
allLLIDs=names(datExpr_count)

rownames(Sig.module_female) <- gsub("^ME", "", rownames(Sig.module_female))

intModules=rownames(Sig.module_female)

intModules

for (module in intModules)
{
  # Select module probes
  modExons = (moduleColors==module)
  # Get their entrez ID codes
  modLLIDs = allLLIDs[modExons];
  # Write them into a file
  fileName = paste("LocusLinkIDsfemale-", module, ".txt", sep="");
  write.table(as.data.frame(modLLIDs), file = fileName,
              row.names = FALSE, col.names = FALSE,quote = F)
}


# As background in the enrichment analysis, we will use all filtered genes in the analysis.
fileName = paste("LocusLinkIDsFemale-all.txt", sep="");
write.table(as.data.frame(allLLIDs), file = fileName,
            row.names = FALSE, col.names = FALSE,quote = F)


# Combine the exons of female WGCNA
files <- list.files(pattern = "LocusLinkIDsfemale-.*\\.txt")
data_list_female <- lapply(files, read.table, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
combined_data_wgcna_female <- do.call(rbind, data_list_female)


duplicates <- duplicated(combined_data_wgcna_female) | duplicated(combined_data_wgcna_female, fromLast = TRUE)
number_of_duplicates <- sum(duplicates)
number_of_duplicates #0

#combined_data_wgcna_female_unique <- unique(combined_data_wgcna_female)


combined_data_wgcna_female <- combined_data_wgcna_female[!grepl("\\+", combined_data_wgcna_female[, 1]), ]
combined_data_wgcna_female


write.table(combined_data_wgcna_female, "combined_wgcna_female_file.txt", col.names = FALSE, row.names = FALSE, sep = "\t")


#######################################################
#WGCNA_exon_male
setwd("E:/exon")
library(DEXSeq)
library(reshape2)
library(DESeq2)
library(WGCNA)
library(dbplyr)


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
colnames(filtered_count)
rm=grep("F", colnames(filtered_count))
filtered_count_male=filtered_count[,-rm]


sampleTable=read.csv("sampleTable_M.csv",header=T,sep = ",",check.names = F)
rownames(sampleTable)=sampleTable$sampleID
colnames(filtered_count_male)==rownames(sampleTable)#check if the count table has the same order as sample table
str(sampleTable)
id=colnames(filtered_count_male)
id=data.frame(id)
sampleTable=sampleTable[id$id,]
colnames(filtered_count_male)==rownames(sampleTable)

sampleTable$salinity <- factor(sampleTable$salinity)
rlog_count=rlogTransformation(filtered_count_male,blind=F)
head(rlog_count,3)
rlog_count_data=as.data.frame(rlog_count)
#-----------------#


####1 data loading ####
options(stringsAsFactors = FALSE)
enableWGCNAThreads()#multiple thread
disableWGCNAThreads()
datExpr_count= as.data.frame(t(rlog_count_data))#row:sample,col:exon

gsg = goodSamplesGenes(datExpr_count, verbose = 3)
gsg$allOK#TRUE


sampleTree = hclust(dist(datExpr_count,method = "euclidean"), method = "ward.D2") 
sizeGrWindow(12,9) 
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)


#trait data loading
rownames(datExpr_count)==rownames(sampleTable)#TRUE


sampleTree2 = hclust(dist(datExpr_count), method = "average")
traitColors = numbers2colors(as.numeric(sampleTable$salinity), signed = TRUE) 
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(sampleTable),
                    main = "Sample dendrogram and trait heatmap")

#data save
save(datExpr_count, sampleTable, file = "wgcna_male_dataInput.RData")


####network construct####
#soft-power selection
gc()
memory.limit(999999999)
powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft = pickSoftThreshold(datExpr_count, powerVector = powers, verbose = 5, networkType="signed", allowWGCNAThreads(nThreads = 6))
save(sft, file = "sft_male.RData")
sizeGrWindow(12,12)
par(mfrow = c(1,2))
cex1 = 0.9
#Soft Threshold
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red")
#Mean connectivity
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

sft$powerEstimate # 18
gc()

#one_step
softPower = sft$powerEstimate

#
net=blockwiseModules(datExpr_count, power = softPower, networkType = "signed",
                     TOMType = "signed", minModuleSize = 30, maxBlockSize = 5000,
                     reassignThreshold = 0, mergeCutHeight = 0.25, numericLabels = TRUE,
                     pamRespectsDendro = FALSE, saveTOMs = TRUE, saveTOMFileBase = "methAfterFilterTOM0.15",
                     verbose = 3)
## time consuming; save result to avoid running it again!
save(net, file = "net_male_exon.RData")


# get module label for each gene ('0' corresponding to no
# module) and convert label to a module color ('grey'
# corresponding to no module)
moduleLabels=net$colors
moduleColors=labels2colors(net$colors)
## number of sites without module (assigned to grey module):
length(moduleLabels[moduleLabels == 0])#472
# get genes per module #table(net$colors)
ExonsperModule <- NULL
for (c in seq(1, length(levels(as.factor(moduleColors))), by = 1)) {
  color <- levels(as.factor(moduleColors))[c]
  length <- length(moduleColors[moduleColors == color])
  module <- data.frame(Module = color, NoGenes = length)
  ExonsperModule <- rbind(ExonsperModule, module)
}


# Calculate eigengenes
MEList = moduleEigengenes(datExpr_count, colors = moduleColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")


#Quantify module similarity by eigengene correlation. Eigengenes: Module representatives
sizeGrWindow(6, 6)
par(cex = 1)
png("Eigengene adjacency heatmap.png",units="in", width=10, height=10,res=400)
plotEigengeneNetworks(orderMEs(MEs), "", plotAdjacency = T, colorLabels = F,
                      marDendro = c(0, 4, 2, 4), marHeatmap = c(3, 4, 1, 1.5),
                      plotDendrograms = T, xLabelsAngle = 90)
dev.off()
# get correlation between modules
cor(orderMEs(MEs))


# quantify module-trait associations
# Define numbers of exons and samples
nExons = ncol(datExpr_count)
nSamples = nrow(datExpr_count)

# Recalculate MEs with color labels
MEs = orderMEs(MEs)

sampleTable=read.csv("sampleTable_M.csv",header=T,sep = ",",check.names = F)
rownames(sampleTable)=sampleTable$sampleID
str(sampleTable)
design_male=sampleTable

rownames(design_male)=design_male$sampleID
design_male$sampleID=NULL

design_male

colnames(MEs) <- sub("^ME", "", colnames(MEs))
moduleTraitCor = cor(MEs, design_male, use = "p");

moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

Sig.module_male=as.data.frame(moduleTraitPvalue)[which(as.data.frame(moduleTraitPvalue)$'salinity'<0.05|
                                                           as.data.frame(moduleTraitPvalue)$'salinity:1 h'<0.05|
                                                           as.data.frame(moduleTraitPvalue)$'salinity:24 h'<0.05|
                                                           as.data.frame(moduleTraitPvalue)$'salinity:48 h'<0.05|
                                                           as.data.frame(moduleTraitPvalue)$'salinity:72 h'<0.05),]
rownames(Sig.module_male)#72

Sig.moduleTraitCor <- moduleTraitCor[row.names(moduleTraitCor) %in% row.names(Sig.module_male), ]
Sig.moduleTraitPvalue <- moduleTraitPvalue[row.names(moduleTraitPvalue) %in% row.names(Sig.module_male), ]


# Will display correlations and their p-values
# 创建一个矩阵，用于存储最终的注释文本
textMatrix <- matrix("-", nrow = nrow(Sig.moduleTraitPvalue), ncol = ncol(Sig.moduleTraitPvalue))

# 找出 Sig.moduleTraitPvalue 中小于 0.05 的位置
significant_positions <- which(Sig.moduleTraitPvalue < 0.05, arr.ind = TRUE)

# 对于小于 0.05 的值，保留两位有效数字并存储到 textMatrix 中
textMatrix[significant_positions] <- signif(Sig.moduleTraitPvalue[significant_positions], 2)

# 确保 textMatrix 的维度与 Sig.moduleTraitPvalue 一致
dim(textMatrix) <- dim(Sig.moduleTraitPvalue)


# Display the correlation values within a heatmap plot
png("WGCNA_male_heatmap1.png", width = 3000, height = 6000, res = 300)
par(mar = c(6, 15, 3, 3));
labeledHeatmap(Matrix = Sig.moduleTraitCor,
               xLabels = names(as.data.frame(design_male)),
               yLabels = rownames(Sig.moduleTraitCor),
               ySymbols = rownames(Sig.moduleTraitCor),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()


# $ Choose interesting modules for downstream GO analysis
setwd("E:/exon/male")
allLLIDs=names(datExpr_count)

rownames(Sig.module_male) <- gsub("^ME", "", rownames(Sig.module_male))

intModules=rownames(Sig.module_male)

intModules

for (module in intModules)
{
  # Select module probes
  modExons = (moduleColors==module)
  # Get their entrez ID codes
  modLLIDs = allLLIDs[modExons];
  # Write them into a file
  fileName = paste("LocusLinkIDsmale-", module, ".txt", sep="");
  write.table(as.data.frame(modLLIDs), file = fileName,
              row.names = FALSE, col.names = FALSE,quote = F)
}


# As background in the enrichment analysis, we will use all filtered genes in the analysis.
fileName = paste("LocusLinkIDsMale-all.txt", sep="");
write.table(as.data.frame(allLLIDs), file = fileName,
            row.names = FALSE, col.names = FALSE,quote = F)


# Combine the exons of male WGCNA
files <- list.files(pattern = "LocusLinkIDsmale-.*\\.txt")
data_list_male <- lapply(files, read.table, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
combined_data_wgcna_male <- do.call(rbind, data_list_male)


duplicates <- duplicated(combined_data_wgcna_male) | duplicated(combined_data_wgcna_male, fromLast = TRUE)
number_of_duplicates <- sum(duplicates)
number_of_duplicates #0

#combined_data_wgcna_male_unique <- unique(combined_data_wgcna_male)


combined_data_wgcna_male <- combined_data_wgcna_male[!grepl("\\+", combined_data_wgcna_male[, 1]), ]
combined_data_wgcna_male


write.table(combined_data_wgcna_male, "combined_wgcna_male_file.txt", col.names = FALSE, row.names = FALSE, sep = "\t")

