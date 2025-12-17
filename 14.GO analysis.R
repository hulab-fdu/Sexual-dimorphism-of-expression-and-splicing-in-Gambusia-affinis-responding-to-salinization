library(tidyr)
library(rtracklayer)
library(data.table)
library(topGO)
library(ggpubr)
library(GenomicRanges)
library(biomaRt)
library(topGO)
library(ggsci)


# 1 GO analysis: all parallel express gene #
setwd("E:/filter")
counts=read.csv('filter_table.txt',header=T,sep="\t")
# import Gambusia annotation, extract filtered genes info
setwd("E:/goterm")
gtf = rtracklayer::import("Gambusia_affinis.ASM309773v1.109.gtf")
gtf_df=as.data.frame(gtf)
gtf_df_gene=gtf_df[gtf$type=="gene",] # only need gene anno
gene_pool=gtf_df_gene[gtf_df_gene$gene_id %in% rownames(counts),]
rownames(gene_pool)=c(1:nrow(gene_pool))
gene_pool_grange=as(gene_pool[,1:3], "GRanges")

# import liftover annotation from platyfish, a.k.a. xm
liftover.gtf=rtracklayer::import("mosquito.ensembl.liftover.gtf")
liftover.gtf_df=as.data.frame(liftover.gtf)
liftover.gtf_df_gene=liftover.gtf_df[liftover.gtf_df$type=="gene",]
rownames(liftover.gtf_df_gene)=c(1:nrow(liftover.gtf_df_gene))
liftover.gtf_df_gene_grange=as(liftover.gtf_df_gene[,1:3], "GRanges")

# find overlaps between filtered genes info and liftover annotation
library(GenomicRanges)
overlap_GO=findOverlaps(gene_pool_grange, liftover.gtf_df_gene_grange, type = c("within"), ignore.strand=T)
overlap_GO_df=data.frame(overlap_GO)

# individually extract ensembl gene names from gene pool and xm
gene_pool_filtered=gene_pool[overlap_GO_df$queryHits,]
liftover.gtf_df_gene_filtered=liftover.gtf_df_gene[overlap_GO_df$subjectHits,]

# this is the key for GO term analysis
combine=data.frame(ga=gene_pool_filtered$gene_id, xm=liftover.gtf_df_gene_filtered$gene_id)

#get GO id for each gene
library(biomaRt)
mart=useDataset("xmaculatus_gene_ensembl", useMart("ensembl"))

#get all go id for universe
ensembl_genes_universe=combine$xm

#用getBM()函数获取注释（如果报错，可能是网络问题，通过连接vpn解决，或者多试几次）
Gene_pool=getBM(filters = "ensembl_gene_id", #已知的注释类型
                attributes = c("ensembl_gene_id","go_id"), #要获取的注释类型
                values = ensembl_genes_universe,#已知的注释类型的数据，把上面我们通过数据处理得到的ensembl基因序号作为ensembl_gene_id 的值 
                mart = mart)#选定的数据库的基因组

# GO: analysis ####

# overlap_GO genes
setwd("E:/liftoff")
parallel_express_gene=read.csv('common_parallel_express_genenames.csv',header=T)
setwd("E:/goterm")
temp1=list(parallel_express_gene)

for (i in 1:length(temp1)) {
  temp_3=temp1[[i]]
  name = paste("group", i, sep="")
  #get all go id for DE genes
  group1_xm=merge(as.data.frame(temp_3),combine,by.x="x",by.y="ga")###这里会失去部分genes
  if(length(is.na(group1_xm))==0){
    next
  }
  ensembl_genes_group1=group1_xm$xm
  group1.gene=getBM(filters = "ensembl_gene_id", 
                    attributes = c("ensembl_gene_id","go_id"), 
                    values = ensembl_genes_group1, 
                    mart = mart)
  #err = try(group1.gene,TRUE)
  
  # Remove blank entries
  Gene_pool <- Gene_pool[Gene_pool$go_id != '',]
  group1.gene <- group1.gene[group1.gene$go_id != '',]
  
  # convert from table format to list format
  geneID2GO <- by(Gene_pool$go_id,
                  Gene_pool$ensembl_gene_id,
                  function(x) as.character(x))
  myInterestingGenes=by(group1.gene$go_id,
                        group1.gene$ensembl_gene_id,
                        function(x) as.character(x))
  # examine result
  head(geneID2GO)
  head(myInterestingGenes)
  
  correction<-"fdr" #多重检验校正
  geneNames = names(geneID2GO)
  myInterestingGenesNames=names(myInterestingGenes)
  geneList = factor(as.integer(geneNames %in% myInterestingGenesNames))#定义背景基因和感兴趣基因
  #有两个levels：0和1，其中1表示感兴趣的基因
  names(geneList) <- geneNames
  head(geneList,3)
  
  ontology=c("MF","BP","CC")
  for (i in 1:length(ontology)) {
    tgData = new("topGOdata", 
                 ontology = ontology[i], 
                 allGenes = geneList, 
                 annot = annFUN.gene2GO, #基因对应的GO注释读取（从geneID2GO中读取）
                 gene2GO = geneID2GO)
    fisherRes = runTest(tgData, algorithm="classic", statistic="fisher")
    fisherResCor = p.adjust(score(fisherRes), method=correction)
    weightRes = runTest(tgData, algorithm="weight01", statistic="fisher")
    weightResCor = p.adjust(score(weightRes), method=correction)
    allRes = GenTable(tgData, 
                      classic=fisherRes, 
                      weight=weightRes, 
                      orderBy="weight", 
                      ranksOf="classic", 
                      topNodes=200)#提取top 200个GO.如何判断显著富集？fisher.COR<0.1
    allRes$fisher.COR = fisherResCor[allRes$GO.ID]
    allRes$weight.COR = weightResCor[allRes$GO.ID]
    write.csv(allRes, paste(name,ontology[i],"csv",sep="."))}
}

#MF
express_MF=read.table("./express.MF.csv",header = T,row.names = 1,sep = ",")
express_MF_gene2GO=merge(express_MF,Gene_pool,by.x="GO.ID",by.y="go_id")
express_MF_gene2GO_final=merge(express_MF_gene2GO,combine,by.x="ensembl_gene_id",by.y="xm")
express_MF_gene2GO_final=express_MF_gene2GO_final[,c("ensembl_gene_id","ga","GO.ID","Term","Annotated","Significant","Expected",
                                                     "classic","weight","fisher.COR","weight.COR" )]

express_MF_gene2GO_final=merge(express_MF_gene2GO_final,as.data.frame(parallel_express_gene),by.x = "ga",by.y = "x")
express_MF_gene2GO_final=express_MF_gene2GO_final[order(express_MF_gene2GO_final$weight),]
express_MF_gene2GO_final=express_MF_gene2GO_final[which(express_MF_gene2GO_final$weight<0.05),]
length(unique(express_MF_gene2GO_final$ga))#40
length(unique(express_MF_gene2GO_final$GO.ID))#31
write.csv(express_MF_gene2GO_final,"./express_MF_gene2GO_final.csv",row.names = F)

#BP
express_BP=read.table("./express.BP.csv",header = T,row.names = 1,sep = ",")
express_BP_gene2GO=merge(express_BP,Gene_pool,by.x="GO.ID",by.y="go_id")
express_BP_gene2GO_final=merge(express_BP_gene2GO,combine,by.x="ensembl_gene_id",by.y="xm")
express_BP_gene2GO_final=express_BP_gene2GO_final[,c("ensembl_gene_id","ga","GO.ID","Term","Annotated","Significant","Expected",
                                                       "classic","weight","fisher.COR","weight.COR" )]

express_BP_gene2GO_final=merge(express_BP_gene2GO_final,as.data.frame(parallel_express_gene),by.x = "ga",by.y = "x")
express_BP_gene2GO_final=express_BP_gene2GO_final[order(express_BP_gene2GO_final$weight),]
express_BP_gene2GO_final=express_BP_gene2GO_final[which(express_BP_gene2GO_final$weight<0.05),]
length(unique(express_BP_gene2GO_final$ga))#43
length(unique(express_BP_gene2GO_final$GO.ID))#54
write.csv(express_BP_gene2GO_final,"./express_BP_gene2GO_final.csv",row.names = F)

#CC
express_CC=read.table("./express.CC.csv",header = T,row.names = 1,sep = ",")
express_CC_gene2GO=merge(express_CC,Gene_pool,by.x="GO.ID",by.y="go_id")
express_CC_gene2GO_final=merge(express_CC_gene2GO,combine,by.x="ensembl_gene_id",by.y="xm")
express_CC_gene2GO_final=express_CC_gene2GO_final[,c("ensembl_gene_id","ga","GO.ID","Term","Annotated","Significant","Expected",
                                                     "classic","weight","fisher.COR","weight.COR" )]

express_CC_gene2GO_final=merge(express_CC_gene2GO_final,as.data.frame(parallel_express_gene),by.x = "ga",by.y = "x")
express_CC_gene2GO_final=express_CC_gene2GO_final[order(express_CC_gene2GO_final$weight),]
express_CC_gene2GO_final=express_CC_gene2GO_final[which(express_CC_gene2GO_final$weight<0.05),]
length(unique(express_CC_gene2GO_final$ga))#51
length(unique(express_CC_gene2GO_final$GO.ID))#10
write.csv(express_CC_gene2GO_final,"./express_CC_gene2GO_final.csv",row.names = F)

# 2 GO analysis: all nonparallel express gene #
setwd("E:/filter")
counts=read.csv('filter_table.txt',header=T,sep="\t")
# import Gambusia annotation, extract filtered genes info
setwd("E:/goterm")
gtf = rtracklayer::import("Gambusia_affinis.ASM309773v1.109.gtf")
gtf_df=as.data.frame(gtf)
gtf_df_gene=gtf_df[gtf$type=="gene",] # only need gene anno
gene_pool=gtf_df_gene[gtf_df_gene$gene_id %in% rownames(counts),]
rownames(gene_pool)=c(1:nrow(gene_pool))
gene_pool_grange=as(gene_pool[,1:3], "GRanges")

# import liftover annotation from platyfish, a.k.a. xm
liftover.gtf=rtracklayer::import("mosquito.ensembl.liftover.gtf")
liftover.gtf_df=as.data.frame(liftover.gtf)
liftover.gtf_df_gene=liftover.gtf_df[liftover.gtf_df$type=="gene",]
rownames(liftover.gtf_df_gene)=c(1:nrow(liftover.gtf_df_gene))
liftover.gtf_df_gene_grange=as(liftover.gtf_df_gene[,1:3], "GRanges")

# find overlaps between filtered genes info and liftover annotation
library(GenomicRanges)
overlap_GO=findOverlaps(gene_pool_grange, liftover.gtf_df_gene_grange, type = c("within"), ignore.strand=T)
overlap_GO_df=data.frame(overlap_GO)

# individually extract ensembl gene names from gene pool and xm
gene_pool_filtered=gene_pool[overlap_GO_df$queryHits,]
liftover.gtf_df_gene_filtered=liftover.gtf_df_gene[overlap_GO_df$subjectHits,]

# this is the key for GO term analysis
combine=data.frame(ga=gene_pool_filtered$gene_id, xm=liftover.gtf_df_gene_filtered$gene_id)

#get GO id for each gene
library(biomaRt)
mart=useDataset("xmaculatus_gene_ensembl", useMart("ensembl"))

#get all go id for universe
ensembl_genes_universe=combine$xm

#用getBM()函数获取注释（如果报错，可能是网络问题，通过连接vpn解决，或者多试几次）
Gene_pool=getBM(filters = "ensembl_gene_id", #已知的注释类型
                attributes = c("ensembl_gene_id","go_id"), #要获取的注释类型
                values = ensembl_genes_universe,#已知的注释类型的数据，把上面我们通过数据处理得到的ensembl基因序号作为ensembl_gene_id 的值 
                mart = mart)#选定的数据库的基因组

# GO: analysis ####

# overlap_GO genes
setwd("E:/liftoff")
nonparallel_express_gene=read.csv('non_pa_express_genenames.csv',header=T)
setwd("E:/goterm")
temp1=list(nonparallel_express_gene)

for (i in 1:length(temp1)) {
  temp_3=temp1[[i]]
  name = paste("group", i, sep="")
  #get all go id for DE genes
  group1_xm=merge(as.data.frame(temp_3),combine,by.x="x",by.y="ga")###这里会失去部分genes
  if(length(is.na(group1_xm))==0){
    next
  }
  ensembl_genes_group1=group1_xm$xm
  group1.gene=getBM(filters = "ensembl_gene_id", 
                    attributes = c("ensembl_gene_id","go_id"), 
                    values = ensembl_genes_group1, 
                    mart = mart)
  #err = try(group1.gene,TRUE)
  
  # Remove blank entries
  Gene_pool <- Gene_pool[Gene_pool$go_id != '',]
  group1.gene <- group1.gene[group1.gene$go_id != '',]
  
  # convert from table format to list format
  geneID2GO <- by(Gene_pool$go_id,
                  Gene_pool$ensembl_gene_id,
                  function(x) as.character(x))
  myInterestingGenes=by(group1.gene$go_id,
                        group1.gene$ensembl_gene_id,
                        function(x) as.character(x))
  # examine result
  head(geneID2GO)
  head(myInterestingGenes)
  
  correction<-"fdr" #多重检验校正
  geneNames = names(geneID2GO)
  myInterestingGenesNames=names(myInterestingGenes)
  geneList = factor(as.integer(geneNames %in% myInterestingGenesNames))#定义背景基因和感兴趣基因
  #有两个levels：0和1，其中1表示感兴趣的基因
  names(geneList) <- geneNames
  head(geneList,3)
  
  ontology=c("MF","BP","CC")
  for (i in 1:length(ontology)) {
    tgData = new("topGOdata", 
                 ontology = ontology[i], 
                 allGenes = geneList, 
                 annot = annFUN.gene2GO, #基因对应的GO注释读取（从geneID2GO中读取）
                 gene2GO = geneID2GO)
    fisherRes = runTest(tgData, algorithm="classic", statistic="fisher")
    fisherResCor = p.adjust(score(fisherRes), method=correction)
    weightRes = runTest(tgData, algorithm="weight01", statistic="fisher")
    weightResCor = p.adjust(score(weightRes), method=correction)
    allRes = GenTable(tgData, 
                      classic=fisherRes, 
                      weight=weightRes, 
                      orderBy="weight", 
                      ranksOf="classic", 
                      topNodes=200)#提取top 200个GO.如何判断显著富集？fisher.COR<0.1
    allRes$fisher.COR = fisherResCor[allRes$GO.ID]
    allRes$weight.COR = weightResCor[allRes$GO.ID]
    write.csv(allRes, paste(name,ontology[i],"csv",sep="."))}
}

#MF
express_MF=read.table("./nonexpress.MF.csv",header = T,row.names = 1,sep = ",")
express_MF_gene2GO=merge(express_MF,Gene_pool,by.x="GO.ID",by.y="go_id")
express_MF_gene2GO_final=merge(express_MF_gene2GO,combine,by.x="ensembl_gene_id",by.y="xm")
express_MF_gene2GO_final=express_MF_gene2GO_final[,c("ensembl_gene_id","ga","GO.ID","Term","Annotated","Significant","Expected",
                                                     "classic","weight","fisher.COR","weight.COR" )]

express_MF_gene2GO_final=merge(express_MF_gene2GO_final,as.data.frame(nonparallel_express_gene),by.x = "ga",by.y = "x")
express_MF_gene2GO_final=express_MF_gene2GO_final[order(express_MF_gene2GO_final$weight),]
express_MF_gene2GO_final=express_MF_gene2GO_final[which(express_MF_gene2GO_final$weight<0.05),]
length(unique(express_MF_gene2GO_final$ga))#163
length(unique(express_MF_gene2GO_final$GO.ID))#16
write.csv(express_MF_gene2GO_final,"./nonexpress_MF_gene2GO_final.csv",row.names = F)

#BP
express_BP=read.table("./nonexpress.BP.csv",header = T,row.names = 1,sep = ",")
express_BP_gene2GO=merge(express_BP,Gene_pool,by.x="GO.ID",by.y="go_id")
express_BP_gene2GO_final=merge(express_BP_gene2GO,combine,by.x="ensembl_gene_id",by.y="xm")
express_BP_gene2GO_final=express_BP_gene2GO_final[,c("ensembl_gene_id","ga","GO.ID","Term","Annotated","Significant","Expected",
                                                     "classic","weight","fisher.COR","weight.COR" )]

express_BP_gene2GO_final=merge(express_BP_gene2GO_final,as.data.frame(nonparallel_express_gene),by.x = "ga",by.y = "x")
express_BP_gene2GO_final=express_BP_gene2GO_final[order(express_BP_gene2GO_final$weight),]
express_BP_gene2GO_final=express_BP_gene2GO_final[which(express_BP_gene2GO_final$weight<0.05),]
length(unique(express_BP_gene2GO_final$ga))#222
length(unique(express_BP_gene2GO_final$GO.ID))#36
write.csv(express_BP_gene2GO_final,"./nonexpress_BP_gene2GO_final.csv",row.names = F)

#CC
express_CC=read.table("./nonexpress.CC.csv",header = T,row.names = 1,sep = ",")
express_CC_gene2GO=merge(express_CC,Gene_pool,by.x="GO.ID",by.y="go_id")
express_CC_gene2GO_final=merge(express_CC_gene2GO,combine,by.x="ensembl_gene_id",by.y="xm")
express_CC_gene2GO_final=express_CC_gene2GO_final[,c("ensembl_gene_id","ga","GO.ID","Term","Annotated","Significant","Expected",
                                                     "classic","weight","fisher.COR","weight.COR" )]

express_CC_gene2GO_final=merge(express_CC_gene2GO_final,as.data.frame(nonparallel_express_gene),by.x = "ga",by.y = "x")
express_CC_gene2GO_final=express_CC_gene2GO_final[order(express_CC_gene2GO_final$weight),]
express_CC_gene2GO_final=express_CC_gene2GO_final[which(express_CC_gene2GO_final$weight<0.05),]
length(unique(express_CC_gene2GO_final$ga))#596
length(unique(express_CC_gene2GO_final$GO.ID))#20
write.csv(express_CC_gene2GO_final,"./nonexpress_CC_gene2GO_final.csv",row.names = F)


# 3 GO analysis: all parallel splice gene #
setwd("E:/exon")
counts=read.csv('filtered_splice_genes.csv',header=T)
rownames(counts) <- counts[, 1]
# import Gambusia annotation, extract filtered genes info
setwd("E:/goterm")
gtf = rtracklayer::import("Gambusia_affinis.ASM309773v1.109.gtf")
gtf_df=as.data.frame(gtf)
gtf_df_gene=gtf_df[gtf$type=="gene",] # only need gene anno
gene_pool=gtf_df_gene[gtf_df_gene$gene_id %in% rownames(counts),]
rownames(gene_pool)=c(1:nrow(gene_pool))
gene_pool_grange=as(gene_pool[,1:3], "GRanges")

# import liftover annotation from platyfish, a.k.a. xm
liftover.gtf=rtracklayer::import("mosquito.ensembl.liftover.gtf")
liftover.gtf_df=as.data.frame(liftover.gtf)
liftover.gtf_df_gene=liftover.gtf_df[liftover.gtf_df$type=="gene",]
rownames(liftover.gtf_df_gene)=c(1:nrow(liftover.gtf_df_gene))
liftover.gtf_df_gene_grange=as(liftover.gtf_df_gene[,1:3], "GRanges")

# find overlaps between filtered genes info and liftover annotation
library(GenomicRanges)
overlap_GO=findOverlaps(gene_pool_grange, liftover.gtf_df_gene_grange, type = c("within"), ignore.strand=T)
overlap_GO_df=data.frame(overlap_GO)

# individually extract ensembl gene names from gene pool and xm
gene_pool_filtered=gene_pool[overlap_GO_df$queryHits,]
liftover.gtf_df_gene_filtered=liftover.gtf_df_gene[overlap_GO_df$subjectHits,]

# this is the key for GO term analysis
combine=data.frame(ga=gene_pool_filtered$gene_id, xm=liftover.gtf_df_gene_filtered$gene_id)

#get GO id for each gene
library(biomaRt)
mart=useDataset("xmaculatus_gene_ensembl", useMart("ensembl"))

#get all go id for universe
ensembl_genes_universe=combine$xm

#用getBM()函数获取注释（如果报错，可能是网络问题，通过连接vpn解决，或者多试几次）
Gene_pool=getBM(filters = "ensembl_gene_id", #已知的注释类型
                attributes = c("ensembl_gene_id","go_id"), #要获取的注释类型
                values = ensembl_genes_universe,#已知的注释类型的数据，把上面我们通过数据处理得到的ensembl基因序号作为ensembl_gene_id 的值 
                mart = mart)#选定的数据库的基因组

# GO: analysis ####

# overlap_GO genes
setwd("E:/liftoff")
parallel_splice_gene=read.csv('common_parallel_splicing_genenames.csv',header=T)
setwd("E:/goterm")
temp1=list(parallel_splice_gene)

for (i in 1:length(temp1)) {
  temp_3=temp1[[i]]
  name = paste("group", i, sep="")
  #get all go id for DE genes
  group1_xm=merge(as.data.frame(temp_3),combine,by.x="x",by.y="ga")###这里会失去部分genes
  if(length(is.na(group1_xm))==0){
    next
  }
  ensembl_genes_group1=group1_xm$xm
  group1.gene=getBM(filters = "ensembl_gene_id", 
                    attributes = c("ensembl_gene_id","go_id"), 
                    values = ensembl_genes_group1, 
                    mart = mart)
  #err = try(group1.gene,TRUE)
  
  # Remove blank entries
  Gene_pool <- Gene_pool[Gene_pool$go_id != '',]
  group1.gene <- group1.gene[group1.gene$go_id != '',]
  
  # convert from table format to list format
  geneID2GO <- by(Gene_pool$go_id,
                  Gene_pool$ensembl_gene_id,
                  function(x) as.character(x))
  myInterestingGenes=by(group1.gene$go_id,
                        group1.gene$ensembl_gene_id,
                        function(x) as.character(x))
  # examine result
  head(geneID2GO)
  head(myInterestingGenes)
  
  correction<-"fdr" #多重检验校正
  geneNames = names(geneID2GO)
  myInterestingGenesNames=names(myInterestingGenes)
  geneList = factor(as.integer(geneNames %in% myInterestingGenesNames))#定义背景基因和感兴趣基因
  #有两个levels：0和1，其中1表示感兴趣的基因
  names(geneList) <- geneNames
  head(geneList,3)
  
  ontology=c("MF","BP","CC")
  for (i in 1:length(ontology)) {
    tgData = new("topGOdata", 
                 ontology = ontology[i], 
                 allGenes = geneList, 
                 annot = annFUN.gene2GO, #基因对应的GO注释读取（从geneID2GO中读取）
                 gene2GO = geneID2GO)
    fisherRes = runTest(tgData, algorithm="classic", statistic="fisher")
    fisherResCor = p.adjust(score(fisherRes), method=correction)
    weightRes = runTest(tgData, algorithm="weight01", statistic="fisher")
    weightResCor = p.adjust(score(weightRes), method=correction)
    allRes = GenTable(tgData, 
                      classic=fisherRes, 
                      weight=weightRes, 
                      orderBy="weight", 
                      ranksOf="classic", 
                      topNodes=200)#提取top 200个GO.如何判断显著富集？fisher.COR<0.1
    allRes$fisher.COR = fisherResCor[allRes$GO.ID]
    allRes$weight.COR = weightResCor[allRes$GO.ID]
    write.csv(allRes, paste(name,ontology[i],"csv",sep="."))}
}

#MF
splice_MF=read.table("./splice.MF.csv",header = T,row.names = 1,sep = ",")
splice_MF_gene2GO=merge(splice_MF,Gene_pool,by.x="GO.ID",by.y="go_id")
splice_MF_gene2GO_final=merge(splice_MF_gene2GO,combine,by.x="ensembl_gene_id",by.y="xm")
splice_MF_gene2GO_final=splice_MF_gene2GO_final[,c("ensembl_gene_id","ga","GO.ID","Term","Annotated","Significant","Expected",
                                                   "classic","weight","fisher.COR","weight.COR" )]

splice_MF_gene2GO_final=merge(splice_MF_gene2GO_final,as.data.frame(parallel_splice_gene),by.x = "ga",by.y = "x")
splice_MF_gene2GO_final=splice_MF_gene2GO_final[order(splice_MF_gene2GO_final$weight),]
splice_MF_gene2GO_final=splice_MF_gene2GO_final[which(splice_MF_gene2GO_final$weight<0.05),]
length(unique(splice_MF_gene2GO_final$ga))#52
length(unique(splice_MF_gene2GO_final$GO.ID))#30
write.csv(splice_MF_gene2GO_final,"./splice_MF_gene2GO_final.csv",row.names = F)

#BP
splice_BP=read.table("./splice.BP.csv",header = T,row.names = 1,sep = ",")
splice_BP_gene2GO=merge(splice_BP,Gene_pool,by.x="GO.ID",by.y="go_id")
splice_BP_gene2GO_final=merge(splice_BP_gene2GO,combine,by.x="ensembl_gene_id",by.y="xm")
splice_BP_gene2GO_final=splice_BP_gene2GO_final[,c("ensembl_gene_id","ga","GO.ID","Term","Annotated","Significant","Expected",
                                                     "classic","weight","fisher.COR","weight.COR" )]

splice_BP_gene2GO_final=merge(splice_BP_gene2GO_final,as.data.frame(parallel_splice_gene),by.x = "ga",by.y = "x")
splice_BP_gene2GO_final=splice_BP_gene2GO_final[order(splice_BP_gene2GO_final$weight),]
splice_BP_gene2GO_final=splice_BP_gene2GO_final[which(splice_BP_gene2GO_final$weight<0.05),]
length(unique(splice_BP_gene2GO_final$ga))#16
length(unique(splice_BP_gene2GO_final$GO.ID))#8
write.csv(splice_BP_gene2GO_final,"./splice_BP_gene2GO_final.csv",row.names = F)

#CC
splice_CC=read.table("./splice.CC.csv",header = T,row.names = 1,sep = ",")
splice_CC_gene2GO=merge(splice_CC,Gene_pool,by.x="GO.ID",by.y="go_id")
splice_CC_gene2GO_final=merge(splice_CC_gene2GO,combine,by.x="ensembl_gene_id",by.y="xm")
splice_CC_gene2GO_final=splice_CC_gene2GO_final[,c("ensembl_gene_id","ga","GO.ID","Term","Annotated","Significant","Expected",
                                                   "classic","weight","fisher.COR","weight.COR" )]

splice_CC_gene2GO_final=merge(splice_CC_gene2GO_final,as.data.frame(parallel_splice_gene),by.x = "ga",by.y = "x")
splice_CC_gene2GO_final=splice_CC_gene2GO_final[order(splice_CC_gene2GO_final$weight),]
splice_CC_gene2GO_final=splice_CC_gene2GO_final[which(splice_CC_gene2GO_final$weight<0.05),]
length(unique(splice_CC_gene2GO_final$ga))#19
length(unique(splice_CC_gene2GO_final$GO.ID))#13
write.csv(splice_CC_gene2GO_final,"./splice_CC_gene2GO_final.csv",row.names = F)

# 4 GO analysis: all nonparallel splice gene #
setwd("E:/exon")
counts=read.csv('filtered_splice_genes.csv',header=T)
rownames(counts) <- counts[, 1]
# import Gambusia annotation, extract filtered genes info
setwd("E:/goterm")
gtf = rtracklayer::import("Gambusia_affinis.ASM309773v1.109.gtf")
gtf_df=as.data.frame(gtf)
gtf_df_gene=gtf_df[gtf$type=="gene",] # only need gene anno
gene_pool=gtf_df_gene[gtf_df_gene$gene_id %in% rownames(counts),]
rownames(gene_pool)=c(1:nrow(gene_pool))
gene_pool_grange=as(gene_pool[,1:3], "GRanges")

# import liftover annotation from platyfish, a.k.a. xm
liftover.gtf=rtracklayer::import("mosquito.ensembl.liftover.gtf")
liftover.gtf_df=as.data.frame(liftover.gtf)
liftover.gtf_df_gene=liftover.gtf_df[liftover.gtf_df$type=="gene",]
rownames(liftover.gtf_df_gene)=c(1:nrow(liftover.gtf_df_gene))
liftover.gtf_df_gene_grange=as(liftover.gtf_df_gene[,1:3], "GRanges")

# find overlaps between filtered genes info and liftover annotation
library(GenomicRanges)
overlap_GO=findOverlaps(gene_pool_grange, liftover.gtf_df_gene_grange, type = c("within"), ignore.strand=T)
overlap_GO_df=data.frame(overlap_GO)

# individually extract ensembl gene names from gene pool and xm
gene_pool_filtered=gene_pool[overlap_GO_df$queryHits,]
liftover.gtf_df_gene_filtered=liftover.gtf_df_gene[overlap_GO_df$subjectHits,]

# this is the key for GO term analysis
combine=data.frame(ga=gene_pool_filtered$gene_id, xm=liftover.gtf_df_gene_filtered$gene_id)

#get GO id for each gene
library(biomaRt)
mart=useDataset("xmaculatus_gene_ensembl", useMart("ensembl"))

#get all go id for universe
ensembl_genes_universe=combine$xm

#用getBM()函数获取注释（如果报错，可能是网络问题，通过连接vpn解决，或者多试几次）
Gene_pool=getBM(filters = "ensembl_gene_id", #已知的注释类型
                attributes = c("ensembl_gene_id","go_id"), #要获取的注释类型
                values = ensembl_genes_universe,#已知的注释类型的数据，把上面我们通过数据处理得到的ensembl基因序号作为ensembl_gene_id 的值 
                mart = mart)#选定的数据库的基因组

# GO: analysis ####

# overlap_GO genes
setwd("E:/liftoff")
nonparallel_splice_gene=read.csv('non_pa_splice_genenames.csv',header=T)
setwd("E:/goterm")
temp1=list(nonparallel_splice_gene)

for (i in 1:length(temp1)) {
  temp_3=temp1[[i]]
  name = paste("group", i, sep="")
  #get all go id for DE genes
  group1_xm=merge(as.data.frame(temp_3),combine,by.x="x",by.y="ga")###这里会失去部分genes
  if(length(is.na(group1_xm))==0){
    next
  }
  ensembl_genes_group1=group1_xm$xm
  group1.gene=getBM(filters = "ensembl_gene_id", 
                    attributes = c("ensembl_gene_id","go_id"), 
                    values = ensembl_genes_group1, 
                    mart = mart)
  #err = try(group1.gene,TRUE)
  
  # Remove blank entries
  Gene_pool <- Gene_pool[Gene_pool$go_id != '',]
  group1.gene <- group1.gene[group1.gene$go_id != '',]
  
  # convert from table format to list format
  geneID2GO <- by(Gene_pool$go_id,
                  Gene_pool$ensembl_gene_id,
                  function(x) as.character(x))
  myInterestingGenes=by(group1.gene$go_id,
                        group1.gene$ensembl_gene_id,
                        function(x) as.character(x))
  # examine result
  head(geneID2GO)
  head(myInterestingGenes)
  
  correction<-"fdr" #多重检验校正
  geneNames = names(geneID2GO)
  myInterestingGenesNames=names(myInterestingGenes)
  geneList = factor(as.integer(geneNames %in% myInterestingGenesNames))#定义背景基因和感兴趣基因
  #有两个levels：0和1，其中1表示感兴趣的基因
  names(geneList) <- geneNames
  head(geneList,3)
  
  ontology=c("MF","BP","CC")
  for (i in 1:length(ontology)) {
    tgData = new("topGOdata", 
                 ontology = ontology[i], 
                 allGenes = geneList, 
                 annot = annFUN.gene2GO, #基因对应的GO注释读取（从geneID2GO中读取）
                 gene2GO = geneID2GO)
    fisherRes = runTest(tgData, algorithm="classic", statistic="fisher")
    fisherResCor = p.adjust(score(fisherRes), method=correction)
    weightRes = runTest(tgData, algorithm="weight01", statistic="fisher")
    weightResCor = p.adjust(score(weightRes), method=correction)
    allRes = GenTable(tgData, 
                      classic=fisherRes, 
                      weight=weightRes, 
                      orderBy="weight", 
                      ranksOf="classic", 
                      topNodes=200)#提取top 200个GO.如何判断显著富集？fisher.COR<0.1
    allRes$fisher.COR = fisherResCor[allRes$GO.ID]
    allRes$weight.COR = weightResCor[allRes$GO.ID]
    write.csv(allRes, paste(name,ontology[i],"csv",sep="."))}
}

#MF
splice_MF=read.table("./nonsplice.MF.csv",header = T,row.names = 1,sep = ",")
splice_MF_gene2GO=merge(splice_MF,Gene_pool,by.x="GO.ID",by.y="go_id")
splice_MF_gene2GO_final=merge(splice_MF_gene2GO,combine,by.x="ensembl_gene_id",by.y="xm")
splice_MF_gene2GO_final=splice_MF_gene2GO_final[,c("ensembl_gene_id","ga","GO.ID","Term","Annotated","Significant","Expected",
                                                   "classic","weight","fisher.COR","weight.COR" )]

splice_MF_gene2GO_final=merge(splice_MF_gene2GO_final,as.data.frame(nonparallel_splice_gene),by.x = "ga",by.y = "x")
splice_MF_gene2GO_final=splice_MF_gene2GO_final[order(splice_MF_gene2GO_final$weight),]
splice_MF_gene2GO_final=splice_MF_gene2GO_final[which(splice_MF_gene2GO_final$weight<0.05),]
length(unique(splice_MF_gene2GO_final$ga))#706
length(unique(splice_MF_gene2GO_final$GO.ID))#19
write.csv(splice_MF_gene2GO_final,"./nonsplice_MF_gene2GO_final.csv",row.names = F)

#BP
splice_BP=read.table("./nonsplice.BP.csv",header = T,row.names = 1,sep = ",")
splice_BP_gene2GO=merge(splice_BP,Gene_pool,by.x="GO.ID",by.y="go_id")
splice_BP_gene2GO_final=merge(splice_BP_gene2GO,combine,by.x="ensembl_gene_id",by.y="xm")
splice_BP_gene2GO_final=splice_BP_gene2GO_final[,c("ensembl_gene_id","ga","GO.ID","Term","Annotated","Significant","Expected",
                                                   "classic","weight","fisher.COR","weight.COR" )]

splice_BP_gene2GO_final=merge(splice_BP_gene2GO_final,as.data.frame(nonparallel_splice_gene),by.x = "ga",by.y = "x")
splice_BP_gene2GO_final=splice_BP_gene2GO_final[order(splice_BP_gene2GO_final$weight),]
splice_BP_gene2GO_final=splice_BP_gene2GO_final[which(splice_BP_gene2GO_final$weight<0.05),]
length(unique(splice_BP_gene2GO_final$ga))#279
length(unique(splice_BP_gene2GO_final$GO.ID))#35
write.csv(splice_BP_gene2GO_final,"./nonsplice_BP_gene2GO_final.csv",row.names = F)

#CC
splice_CC=read.table("./nonsplice.CC.csv",header = T,row.names = 1,sep = ",")
splice_CC_gene2GO=merge(splice_CC,Gene_pool,by.x="GO.ID",by.y="go_id")
splice_CC_gene2GO_final=merge(splice_CC_gene2GO,combine,by.x="ensembl_gene_id",by.y="xm")
splice_CC_gene2GO_final=splice_CC_gene2GO_final[,c("ensembl_gene_id","ga","GO.ID","Term","Annotated","Significant","Expected",
                                                   "classic","weight","fisher.COR","weight.COR" )]

splice_CC_gene2GO_final=merge(splice_CC_gene2GO_final,as.data.frame(nonparallel_splice_gene),by.x = "ga",by.y = "x")
splice_CC_gene2GO_final=splice_CC_gene2GO_final[order(splice_CC_gene2GO_final$weight),]
splice_CC_gene2GO_final=splice_CC_gene2GO_final[which(splice_CC_gene2GO_final$weight<0.05),]
length(unique(splice_CC_gene2GO_final$ga))#544
length(unique(splice_CC_gene2GO_final$GO.ID))#13
write.csv(splice_CC_gene2GO_final,"./nonsplice_CC_gene2GO_final.csv",row.names = F)


# 5 GO analysis: female parallel express gene (regulated by miRNA) #
setwd("E:/filter")
counts=read.csv('filter_table.txt',header=T,sep="\t")
# import Gambusia annotation, extract filtered genes info
setwd("E:/goterm")
gtf = rtracklayer::import("Gambusia_affinis.ASM309773v1.109.gtf")
gtf_df=as.data.frame(gtf)
gtf_df_gene=gtf_df[gtf$type=="gene",] # only need gene anno
gene_pool=gtf_df_gene[gtf_df_gene$gene_id %in% rownames(counts),]
rownames(gene_pool)=c(1:nrow(gene_pool))
gene_pool_grange=as(gene_pool[,1:3], "GRanges")

# import liftover annotation from platyfish, a.k.a. xm
liftover.gtf=rtracklayer::import("mosquito.ensembl.liftover.gtf")
liftover.gtf_df=as.data.frame(liftover.gtf)
liftover.gtf_df_gene=liftover.gtf_df[liftover.gtf_df$type=="gene",]
rownames(liftover.gtf_df_gene)=c(1:nrow(liftover.gtf_df_gene))
liftover.gtf_df_gene_grange=as(liftover.gtf_df_gene[,1:3], "GRanges")

# find overlaps between filtered genes info and liftover annotation
library(GenomicRanges)
overlap_GO=findOverlaps(gene_pool_grange, liftover.gtf_df_gene_grange, type = c("within"), ignore.strand=T)
overlap_GO_df=data.frame(overlap_GO)

# individually extract ensembl gene names from gene pool and xm
gene_pool_filtered=gene_pool[overlap_GO_df$queryHits,]
liftover.gtf_df_gene_filtered=liftover.gtf_df_gene[overlap_GO_df$subjectHits,]

# this is the key for GO term analysis
combine=data.frame(ga=gene_pool_filtered$gene_id, xm=liftover.gtf_df_gene_filtered$gene_id)

#get GO id for female
library(biomaRt)
mart=useDataset("xmaculatus_gene_ensembl", useMart("ensembl"))

#get all go id for universe
ensembl_genes_universe=combine$xm

#用getBM()函数获取注释（如果报错，可能是网络问题，通过连接vpn解决，或者多试几次）
Gene_pool=getBM(filters = "ensembl_gene_id", #已知的注释类型
                attributes = c("ensembl_gene_id","go_id","name_1006","namespace_1003"), #要获取的注释类型
                values = ensembl_genes_universe,
                mart = mart,
                curl = curl::new_handle(timeout = 10000))#选定的数据库的基因组


setwd("E:/goterm")
cor_female_gene=read.csv('goterm_female.csv',header=T)

# 加载必要的库
library(dplyr)
# 使用 merge() 函数进行内连接，找出匹配的行
matched_rows <- merge(cor_female_gene, combine, by.x = "mR", by.y = "ga")#6
# 查看结果
print(matched_rows)

female_GO_final=merge(matched_rows,Gene_pool,by.x="xm",by.y="ensembl_gene_id")#21
female_GO_final <- female_GO_final[female_GO_final$go_id != "", ]#19

setwd("E:/goterm")
write.csv(female_GO_final,"./female_GO_final.csv",row.names = F)


# 6 GO analysis: male parallel express gene (regulated by miRNA) #
setwd("E:/filter")
counts=read.csv('filter_table.txt',header=T,sep="\t")
# import Gambusia annotation, extract filtered genes info
setwd("E:/goterm")
gtf = rtracklayer::import("Gambusia_affinis.ASM309773v1.109.gtf")
gtf_df=as.data.frame(gtf)
gtf_df_gene=gtf_df[gtf$type=="gene",] # only need gene anno
gene_pool=gtf_df_gene[gtf_df_gene$gene_id %in% rownames(counts),]
rownames(gene_pool)=c(1:nrow(gene_pool))
gene_pool_grange=as(gene_pool[,1:3], "GRanges")

# import liftover annotation from platyfish, a.k.a. xm
liftover.gtf=rtracklayer::import("mosquito.ensembl.liftover.gtf")
liftover.gtf_df=as.data.frame(liftover.gtf)
liftover.gtf_df_gene=liftover.gtf_df[liftover.gtf_df$type=="gene",]
rownames(liftover.gtf_df_gene)=c(1:nrow(liftover.gtf_df_gene))
liftover.gtf_df_gene_grange=as(liftover.gtf_df_gene[,1:3], "GRanges")

# find overlaps between filtered genes info and liftover annotation
library(GenomicRanges)
overlap_GO=findOverlaps(gene_pool_grange, liftover.gtf_df_gene_grange, type = c("within"), ignore.strand=T)
overlap_GO_df=data.frame(overlap_GO)

# individually extract ensembl gene names from gene pool and xm
gene_pool_filtered=gene_pool[overlap_GO_df$queryHits,]
liftover.gtf_df_gene_filtered=liftover.gtf_df_gene[overlap_GO_df$subjectHits,]

# this is the key for GO term analysis
combine=data.frame(ga=gene_pool_filtered$gene_id, xm=liftover.gtf_df_gene_filtered$gene_id)

#get GO id for male
library(biomaRt)
mart=useDataset("xmaculatus_gene_ensembl", useMart("ensembl"))

#get all go id for universe
ensembl_genes_universe=combine$xm

#用getBM()函数获取注释（如果报错，可能是网络问题，通过连接vpn解决，或者多试几次）
Gene_pool=getBM(filters = "ensembl_gene_id", #已知的注释类型
                attributes = c("ensembl_gene_id","go_id","name_1006","namespace_1003"), #要获取的注释类型
                values = ensembl_genes_universe,
                mart = mart,
                curl = curl::new_handle(timeout = 10000))#选定的数据库的基因组

setwd("E:/goterm")
cor_male_gene=read.csv('goterm_male.csv',header=T)

# 加载必要的库
library(dplyr)
# 使用 merge() 函数进行内连接，找出匹配的行
matched_rows <- merge(cor_male_gene, combine, by.x = "mR", by.y = "ga")#11
# 查看结果
print(matched_rows)

male_GO_final=merge(matched_rows,Gene_pool,by.x="xm",by.y="ensembl_gene_id")#39
male_GO_final <- male_GO_final[male_GO_final$go_id != "", ]#38

setwd("E:/goterm")
write.csv(male_GO_final,"./male_GO_final.csv",row.names = F)


##female:富集分析(运行完5的前一部分code，再运行下列code)
#get all go id for universe
ensembl_genes_universe=combine$xm

#用getBM()函数获取注释（如果报错，可能是网络问题，通过连接vpn解决，或者多试几次）
Gene_pool=getBM(filters = "ensembl_gene_id", #已知的注释类型
                attributes = c("ensembl_gene_id","go_id"), #要获取的注释类型
                values = ensembl_genes_universe,#已知的注释类型的数据，把上面我们通过数据处理得到的ensembl基因序号作为ensembl_gene_id 的值 
                mart = mart)#选定的数据库的基因组

# GO: analysis ####

# overlap_GO genes
setwd("E:/goterm")
parallel_female_gene=read.csv('goterm_female.csv',header=T)
temp1=list(female_parallel_express_gene)

for (i in 1:length(temp1)) {
  temp_3=temp1[[i]]
  name = paste("group", i, sep="")
  #get all go id for DE genes
  group1_xm=merge(as.data.frame(temp_3),combine,by.x="mR",by.y="ga")###这里会失去部分genes
  if(length(is.na(group1_xm))==0){
    next
  }
  ensembl_genes_group1=group1_xm$xm
  group1.gene=getBM(filters = "ensembl_gene_id", 
                    attributes = c("ensembl_gene_id","go_id"), 
                    values = ensembl_genes_group1, 
                    mart = mart)
  #err = try(group1.gene,TRUE)
  
  # Remove blank entries
  Gene_pool <- Gene_pool[Gene_pool$go_id != '',]
  group1.gene <- group1.gene[group1.gene$go_id != '',]
  
  # convert from table format to list format
  geneID2GO <- by(Gene_pool$go_id,
                  Gene_pool$ensembl_gene_id,
                  function(x) as.character(x))
  myInterestingGenes=by(group1.gene$go_id,
                        group1.gene$ensembl_gene_id,
                        function(x) as.character(x))
  # examine result
  head(geneID2GO)
  head(myInterestingGenes)
  
  correction<-"fdr" #多重检验校正
  geneNames = names(geneID2GO)
  myInterestingGenesNames=names(myInterestingGenes)
  geneList = factor(as.integer(geneNames %in% myInterestingGenesNames))#定义背景基因和感兴趣基因
  #有两个levels：0和1，其中1表示感兴趣的基因
  names(geneList) <- geneNames
  head(geneList,3)
  
  ontology=c("MF","BP","CC")
  for (i in 1:length(ontology)) {
    tgData = new("topGOdata", 
                 ontology = ontology[i], 
                 allGenes = geneList, 
                 annot = annFUN.gene2GO, #基因对应的GO注释读取（从geneID2GO中读取）
                 gene2GO = geneID2GO)
    fisherRes = runTest(tgData, algorithm="classic", statistic="fisher")
    fisherResCor = p.adjust(score(fisherRes), method=correction)
    weightRes = runTest(tgData, algorithm="weight01", statistic="fisher")
    weightResCor = p.adjust(score(weightRes), method=correction)
    allRes = GenTable(tgData, 
                      classic=fisherRes, 
                      weight=weightRes, 
                      orderBy="weight", 
                      ranksOf="classic", 
                      topNodes=200)#提取top 200个GO.如何判断显著富集？fisher.COR<0.1
    allRes$fisher.COR = fisherResCor[allRes$GO.ID]
    allRes$weight.COR = weightResCor[allRes$GO.ID]
    write.csv(allRes, paste(name,ontology[i],"csv",sep="."))}
}

#MF
female_MF=read.table("./female.MF.csv",header = T,row.names = 1,sep = ",")
female_MF_gene2GO=merge(female_MF,Gene_pool,by.x="GO.ID",by.y="go_id")
female_MF_gene2GO_final=merge(female_MF_gene2GO,combine,by.x="ensembl_gene_id",by.y="xm")
female_MF_gene2GO_final=female_MF_gene2GO_final[,c("ensembl_gene_id","ga","GO.ID","Term","Annotated","Significant","Expected",
                                                     "classic","weight","fisher.COR","weight.COR" )]

female_MF_gene2GO_final=merge(female_MF_gene2GO_final,as.data.frame(parallel_female_gene),by.x = "ga",by.y = "mR")
female_MF_gene2GO_final=female_MF_gene2GO_final[order(female_MF_gene2GO_final$weight),]
female_MF_gene2GO_final=female_MF_gene2GO_final[which(female_MF_gene2GO_final$weight<0.05),]
length(unique(female_MF_gene2GO_final$ga))#2
length(unique(female_MF_gene2GO_final$GO.ID))#2
write.csv(female_MF_gene2GO_final,"./female_MF_gene2GO_final.csv",row.names = F)

#BP
female_BP=read.table("./female.BP.csv",header = T,row.names = 1,sep = ",")
female_BP_gene2GO=merge(female_BP,Gene_pool,by.x="GO.ID",by.y="go_id")
female_BP_gene2GO_final=merge(female_BP_gene2GO,combine,by.x="ensembl_gene_id",by.y="xm")
female_BP_gene2GO_final=female_BP_gene2GO_final[,c("ensembl_gene_id","ga","GO.ID","Term","Annotated","Significant","Expected",
                                                     "classic","weight","fisher.COR","weight.COR" )]

female_BP_gene2GO_final=merge(female_BP_gene2GO_final,as.data.frame(parallel_female_gene),by.x = "ga",by.y = "mR")
female_BP_gene2GO_final=female_BP_gene2GO_final[order(female_BP_gene2GO_final$weight),]
female_BP_gene2GO_final=female_BP_gene2GO_final[which(female_BP_gene2GO_final$weight<0.05),]
length(unique(female_BP_gene2GO_final$ga))#1
length(unique(female_BP_gene2GO_final$GO.ID))#1
write.csv(female_BP_gene2GO_final,"./female_BP_gene2GO_final.csv",row.names = F)

#CC
female_CC=read.table("./female.CC.csv",header = T,row.names = 1,sep = ",")
female_CC_gene2GO=merge(female_CC,Gene_pool,by.x="GO.ID",by.y="go_id")
female_CC_gene2GO_final=merge(female_CC_gene2GO,combine,by.x="ensembl_gene_id",by.y="xm")
female_CC_gene2GO_final=female_CC_gene2GO_final[,c("ensembl_gene_id","ga","GO.ID","Term","Annotated","Significant","Expected",
                                                     "classic","weight","fisher.COR","weight.COR" )]

female_CC_gene2GO_final=merge(female_CC_gene2GO_final,as.data.frame(parallel_female_gene),by.x = "ga",by.y = "mR")
female_CC_gene2GO_final=female_CC_gene2GO_final[order(female_CC_gene2GO_final$weight),]
female_CC_gene2GO_final=female_CC_gene2GO_final[which(female_CC_gene2GO_final$weight<0.05),]
length(unique(female_CC_gene2GO_final$ga))#1
length(unique(female_CC_gene2GO_final$GO.ID))#1
write.csv(female_CC_gene2GO_final,"./female_CC_gene2GO_final.csv",row.names = F)


##male:富集分析(运行完6的前一部分code，再运行下列code)
#get all go id for universe
ensembl_genes_universe=combine$xm

#用getBM()函数获取注释（如果报错，可能是网络问题，通过连接vpn解决，或者多试几次）
Gene_pool=getBM(filters = "ensembl_gene_id", #已知的注释类型
                attributes = c("ensembl_gene_id","go_id"), #要获取的注释类型
                values = ensembl_genes_universe,#已知的注释类型的数据，把上面我们通过数据处理得到的ensembl基因序号作为ensembl_gene_id 的值 
                mart = mart)#选定的数据库的基因组

# GO: analysis ####

# overlap_GO genes
setwd("E:/goterm")
parallel_male_gene=read.csv('goterm_male.csv',header=T)
temp1=list(parallel_male_gene)

for (i in 1:length(temp1)) {
  temp_3=temp1[[i]]
  name = paste("group", i, sep="")
  #get all go id for DE genes
  group1_xm=merge(as.data.frame(temp_3),combine,by.x="mR",by.y="ga")###这里会失去部分genes
  if(length(is.na(group1_xm))==0){
    next
  }
  ensembl_genes_group1=group1_xm$xm
  group1.gene=getBM(filters = "ensembl_gene_id", 
                    attributes = c("ensembl_gene_id","go_id"), 
                    values = ensembl_genes_group1, 
                    mart = mart)
  #err = try(group1.gene,TRUE)
  
  # Remove blank entries
  Gene_pool <- Gene_pool[Gene_pool$go_id != '',]
  group1.gene <- group1.gene[group1.gene$go_id != '',]
  
  # convert from table format to list format
  geneID2GO <- by(Gene_pool$go_id,
                  Gene_pool$ensembl_gene_id,
                  function(x) as.character(x))
  myInterestingGenes=by(group1.gene$go_id,
                        group1.gene$ensembl_gene_id,
                        function(x) as.character(x))
  # examine result
  head(geneID2GO)
  head(myInterestingGenes)
  
  correction<-"fdr" #多重检验校正
  geneNames = names(geneID2GO)
  myInterestingGenesNames=names(myInterestingGenes)
  geneList = factor(as.integer(geneNames %in% myInterestingGenesNames))#定义背景基因和感兴趣基因
  #有两个levels：0和1，其中1表示感兴趣的基因
  names(geneList) <- geneNames
  head(geneList,3)
  
  ontology=c("MF","BP","CC")
  for (i in 1:length(ontology)) {
    tgData = new("topGOdata", 
                 ontology = ontology[i], 
                 allGenes = geneList, 
                 annot = annFUN.gene2GO, #基因对应的GO注释读取（从geneID2GO中读取）
                 gene2GO = geneID2GO)
    fisherRes = runTest(tgData, algorithm="classic", statistic="fisher")
    fisherResCor = p.adjust(score(fisherRes), method=correction)
    weightRes = runTest(tgData, algorithm="weight01", statistic="fisher")
    weightResCor = p.adjust(score(weightRes), method=correction)
    allRes = GenTable(tgData, 
                      classic=fisherRes, 
                      weight=weightRes, 
                      orderBy="weight", 
                      ranksOf="classic", 
                      topNodes=200)#提取top 200个GO.如何判断显著富集？fisher.COR<0.1
    allRes$fisher.COR = fisherResCor[allRes$GO.ID]
    allRes$weight.COR = weightResCor[allRes$GO.ID]
    write.csv(allRes, paste(name,ontology[i],"csv",sep="."))}
}

#MF
male_MF=read.table("./male.MF.csv",header = T,row.names = 1,sep = ",")
male_MF_gene2GO=merge(male_MF,Gene_pool,by.x="GO.ID",by.y="go_id")
male_MF_gene2GO_final=merge(male_MF_gene2GO,combine,by.x="ensembl_gene_id",by.y="xm")
male_MF_gene2GO_final=male_MF_gene2GO_final[,c("ensembl_gene_id","ga","GO.ID","Term","Annotated","Significant","Expected",
                                                   "classic","weight","fisher.COR","weight.COR" )]

male_MF_gene2GO_final=merge(male_MF_gene2GO_final,as.data.frame(parallel_male_gene),by.x = "ga",by.y = "mR")
male_MF_gene2GO_final=male_MF_gene2GO_final[order(male_MF_gene2GO_final$weight),]
male_MF_gene2GO_final=male_MF_gene2GO_final[which(male_MF_gene2GO_final$weight<0.05),]
length(unique(male_MF_gene2GO_final$ga))#5
length(unique(male_MF_gene2GO_final$GO.ID))#6
write.csv(male_MF_gene2GO_final,"./male_MF_gene2GO_final.csv",row.names = F)

#BP
male_BP=read.table("./male.BP.csv",header = T,row.names = 1,sep = ",")
male_BP_gene2GO=merge(male_BP,Gene_pool,by.x="GO.ID",by.y="go_id")
male_BP_gene2GO_final=merge(male_BP_gene2GO,combine,by.x="ensembl_gene_id",by.y="xm")
male_BP_gene2GO_final=male_BP_gene2GO_final[,c("ensembl_gene_id","ga","GO.ID","Term","Annotated","Significant","Expected",
                                                   "classic","weight","fisher.COR","weight.COR" )]

male_BP_gene2GO_final=merge(male_BP_gene2GO_final,as.data.frame(parallel_male_gene),by.x = "ga",by.y = "mR")
male_BP_gene2GO_final=male_BP_gene2GO_final[order(male_BP_gene2GO_final$weight),]
male_BP_gene2GO_final=male_BP_gene2GO_final[which(male_BP_gene2GO_final$weight<0.05),]
length(unique(male_BP_gene2GO_final$ga))#1
length(unique(male_BP_gene2GO_final$GO.ID))#1
write.csv(male_BP_gene2GO_final,"./male_BP_gene2GO_final.csv",row.names = F)

#CC
male_CC=read.table("./male.CC.csv",header = T,row.names = 1,sep = ",")
male_CC_gene2GO=merge(male_CC,Gene_pool,by.x="GO.ID",by.y="go_id")
male_CC_gene2GO_final=merge(male_CC_gene2GO,combine,by.x="ensembl_gene_id",by.y="xm")
male_CC_gene2GO_final=male_CC_gene2GO_final[,c("ensembl_gene_id","ga","GO.ID","Term","Annotated","Significant","Expected",
                                                   "classic","weight","fisher.COR","weight.COR" )]

male_CC_gene2GO_final=merge(male_CC_gene2GO_final,as.data.frame(parallel_male_gene),by.x = "ga",by.y = "mR")
male_CC_gene2GO_final=male_CC_gene2GO_final[order(male_CC_gene2GO_final$weight),]
male_CC_gene2GO_final=male_CC_gene2GO_final[which(male_CC_gene2GO_final$weight<0.05),]
length(unique(male_CC_gene2GO_final$ga))#1
length(unique(male_CC_gene2GO_final$GO.ID))#1
write.csv(male_CC_gene2GO_final,"./male_CC_gene2GO_final.csv",row.names = F)

# 7 GO analysis: female non-parallel splice gene (regulated by miRNA) #
setwd("E:/exon")
counts=read.csv('filtered_splice_genes.csv',header=T)
rownames(counts) <- counts[, 1]
# import Gambusia annotation, extract filtered genes info
setwd("E:/goterm")
gtf = rtracklayer::import("Gambusia_affinis.ASM309773v1.109.gtf")
gtf_df=as.data.frame(gtf)
gtf_df_gene=gtf_df[gtf$type=="gene",] # only need gene anno
gene_pool=gtf_df_gene[gtf_df_gene$gene_id %in% rownames(counts),]
rownames(gene_pool)=c(1:nrow(gene_pool))
gene_pool_grange=as(gene_pool[,1:3], "GRanges")

# import liftover annotation from platyfish, a.k.a. xm
liftover.gtf=rtracklayer::import("mosquito.ensembl.liftover.gtf")
liftover.gtf_df=as.data.frame(liftover.gtf)
liftover.gtf_df_gene=liftover.gtf_df[liftover.gtf_df$type=="gene",]
rownames(liftover.gtf_df_gene)=c(1:nrow(liftover.gtf_df_gene))
liftover.gtf_df_gene_grange=as(liftover.gtf_df_gene[,1:3], "GRanges")

# find overlaps between filtered genes info and liftover annotation
library(GenomicRanges)
overlap_GO=findOverlaps(gene_pool_grange, liftover.gtf_df_gene_grange, type = c("within"), ignore.strand=T)
overlap_GO_df=data.frame(overlap_GO)

# individually extract ensembl gene names from gene pool and xm
gene_pool_filtered=gene_pool[overlap_GO_df$queryHits,]
liftover.gtf_df_gene_filtered=liftover.gtf_df_gene[overlap_GO_df$subjectHits,]

# this is the key for GO term analysis
combine=data.frame(ga=gene_pool_filtered$gene_id, xm=liftover.gtf_df_gene_filtered$gene_id)

#get GO id for female
library(biomaRt)
mart=useDataset("xmaculatus_gene_ensembl", useMart("ensembl"))

#get all go id for universe
ensembl_genes_universe=combine$xm

#用getBM()函数获取注释（如果报错，可能是网络问题，通过连接vpn解决，或者多试几次）
Gene_pool=getBM(filters = "ensembl_gene_id", #已知的注释类型
                attributes = c("ensembl_gene_id","go_id","name_1006","namespace_1003"), #要获取的注释类型
                values = ensembl_genes_universe,
                mart = mart,
                curl = curl::new_handle(timeout = 10000))#选定的数据库的基因组


setwd("E:/goterm")
cor_female_gene=read.csv('goterm_nonsp_female.csv',header=T)

##female:富集分析
#get all go id for universe
ensembl_genes_universe=combine$xm

#用getBM()函数获取注释（如果报错，可能是网络问题，通过连接vpn解决，或者多试几次）
Gene_pool=getBM(filters = "ensembl_gene_id", #已知的注释类型
                attributes = c("ensembl_gene_id","go_id"), #要获取的注释类型
                values = ensembl_genes_universe,#已知的注释类型的数据，把上面我们通过数据处理得到的ensembl基因序号作为ensembl_gene_id 的值 
                mart = mart)#选定的数据库的基因组

# GO: analysis ####

# overlap_GO genes
setwd("E:/goterm")
cor_female_gene=read.csv('goterm_nonsp_female.csv',header=T)
temp1=list(cor_female_gene)

for (i in 1:length(temp1)) {
  temp_3=temp1[[i]]
  name = paste("group", i, sep="")
  #get all go id for DE genes
  group1_xm=merge(as.data.frame(temp_3),combine,by.x="mR",by.y="ga")###这里会失去部分genes
  if(length(is.na(group1_xm))==0){
    next
  }
  ensembl_genes_group1=group1_xm$xm
  group1.gene=getBM(filters = "ensembl_gene_id", 
                    attributes = c("ensembl_gene_id","go_id"), 
                    values = ensembl_genes_group1, 
                    mart = mart)
  #err = try(group1.gene,TRUE)
  
  # Remove blank entries
  Gene_pool <- Gene_pool[Gene_pool$go_id != '',]
  group1.gene <- group1.gene[group1.gene$go_id != '',]
  
  # convert from table format to list format
  geneID2GO <- by(Gene_pool$go_id,
                  Gene_pool$ensembl_gene_id,
                  function(x) as.character(x))
  myInterestingGenes=by(group1.gene$go_id,
                        group1.gene$ensembl_gene_id,
                        function(x) as.character(x))
  # examine result
  head(geneID2GO)
  head(myInterestingGenes)
  
  correction<-"fdr" #多重检验校正
  geneNames = names(geneID2GO)
  myInterestingGenesNames=names(myInterestingGenes)
  geneList = factor(as.integer(geneNames %in% myInterestingGenesNames))#定义背景基因和感兴趣基因
  #有两个levels：0和1，其中1表示感兴趣的基因
  names(geneList) <- geneNames
  head(geneList,3)
  
  ontology=c("MF","BP","CC")
  for (i in 1:length(ontology)) {
    tgData = new("topGOdata", 
                 ontology = ontology[i], 
                 allGenes = geneList, 
                 annot = annFUN.gene2GO, #基因对应的GO注释读取（从geneID2GO中读取）
                 gene2GO = geneID2GO)
    fisherRes = runTest(tgData, algorithm="classic", statistic="fisher")
    fisherResCor = p.adjust(score(fisherRes), method=correction)
    weightRes = runTest(tgData, algorithm="weight01", statistic="fisher")
    weightResCor = p.adjust(score(weightRes), method=correction)
    allRes = GenTable(tgData, 
                      classic=fisherRes, 
                      weight=weightRes, 
                      orderBy="weight", 
                      ranksOf="classic", 
                      topNodes=200)#提取top 200个GO.如何判断显著富集？fisher.COR<0.1
    allRes$fisher.COR = fisherResCor[allRes$GO.ID]
    allRes$weight.COR = weightResCor[allRes$GO.ID]
    write.csv(allRes, paste(name,ontology[i],"csv",sep="."))}
}

#MF
female_MF=read.table("./group1.MF.csv",header = T,row.names = 1,sep = ",")
female_MF_gene2GO=merge(female_MF,Gene_pool,by.x="GO.ID",by.y="go_id")
female_MF_gene2GO_final=merge(female_MF_gene2GO,combine,by.x="ensembl_gene_id",by.y="xm")
female_MF_gene2GO_final=female_MF_gene2GO_final[,c("ensembl_gene_id","ga","GO.ID","Term","Annotated","Significant","Expected",
                                                   "classic","weight","fisher.COR","weight.COR" )]

female_MF_gene2GO_final=merge(female_MF_gene2GO_final,as.data.frame(cor_female_gene),by.x = "ga",by.y = "mR")
female_MF_gene2GO_final=female_MF_gene2GO_final[order(female_MF_gene2GO_final$weight),]
female_MF_gene2GO_final=female_MF_gene2GO_final[which(female_MF_gene2GO_final$weight<0.05),]
length(unique(female_MF_gene2GO_final$ga))#39
length(unique(female_MF_gene2GO_final$GO.ID))#14
write.csv(female_MF_gene2GO_final,"./nonpa_female_MF_gene2GO_final.csv",row.names = F)

#BP
female_BP=read.table("./group1.BP.csv",header = T,row.names = 1,sep = ",")
female_BP_gene2GO=merge(female_BP,Gene_pool,by.x="GO.ID",by.y="go_id")
female_BP_gene2GO_final=merge(female_BP_gene2GO,combine,by.x="ensembl_gene_id",by.y="xm")
female_BP_gene2GO_final=female_BP_gene2GO_final[,c("ensembl_gene_id","ga","GO.ID","Term","Annotated","Significant","Expected",
                                                   "classic","weight","fisher.COR","weight.COR" )]

female_BP_gene2GO_final=merge(female_BP_gene2GO_final,as.data.frame(cor_female_gene),by.x = "ga",by.y = "mR")
female_BP_gene2GO_final=female_BP_gene2GO_final[order(female_BP_gene2GO_final$weight),]
female_BP_gene2GO_final=female_BP_gene2GO_final[which(female_BP_gene2GO_final$weight<0.05),]
length(unique(female_BP_gene2GO_final$ga))#27
length(unique(female_BP_gene2GO_final$GO.ID))#37
write.csv(female_BP_gene2GO_final,"./nonpa_female_BP_gene2GO_final.csv",row.names = F)

#CC
female_CC=read.table("./group1.CC.csv",header = T,row.names = 1,sep = ",")
female_CC_gene2GO=merge(female_CC,Gene_pool,by.x="GO.ID",by.y="go_id")
female_CC_gene2GO_final=merge(female_CC_gene2GO,combine,by.x="ensembl_gene_id",by.y="xm")
female_CC_gene2GO_final=female_CC_gene2GO_final[,c("ensembl_gene_id","ga","GO.ID","Term","Annotated","Significant","Expected",
                                                   "classic","weight","fisher.COR","weight.COR" )]

female_CC_gene2GO_final=merge(female_CC_gene2GO_final,as.data.frame(cor_female_gene),by.x = "ga",by.y = "mR")
female_CC_gene2GO_final=female_CC_gene2GO_final[order(female_CC_gene2GO_final$weight),]
female_CC_gene2GO_final=female_CC_gene2GO_final[which(female_CC_gene2GO_final$weight<0.05),]
length(unique(female_CC_gene2GO_final$ga))#6
length(unique(female_CC_gene2GO_final$GO.ID))#5
write.csv(female_CC_gene2GO_final,"./nonpa_female_CC_gene2GO_final.csv",row.names = F)


# 8 GO analysis: male non-parallel splice gene (regulated by miRNA) #
setwd("E:/exon")
counts=read.csv('filtered_splice_genes.csv',header=T)
rownames(counts) <- counts[, 1]
# import Gambusia annotation, extract filtered genes info
setwd("E:/goterm")
gtf = rtracklayer::import("Gambusia_affinis.ASM309773v1.109.gtf")
gtf_df=as.data.frame(gtf)
gtf_df_gene=gtf_df[gtf$type=="gene",] # only need gene anno
gene_pool=gtf_df_gene[gtf_df_gene$gene_id %in% rownames(counts),]
rownames(gene_pool)=c(1:nrow(gene_pool))
gene_pool_grange=as(gene_pool[,1:3], "GRanges")

# import liftover annotation from platyfish, a.k.a. xm
liftover.gtf=rtracklayer::import("mosquito.ensembl.liftover.gtf")
liftover.gtf_df=as.data.frame(liftover.gtf)
liftover.gtf_df_gene=liftover.gtf_df[liftover.gtf_df$type=="gene",]
rownames(liftover.gtf_df_gene)=c(1:nrow(liftover.gtf_df_gene))
liftover.gtf_df_gene_grange=as(liftover.gtf_df_gene[,1:3], "GRanges")

# find overlaps between filtered genes info and liftover annotation
library(GenomicRanges)
overlap_GO=findOverlaps(gene_pool_grange, liftover.gtf_df_gene_grange, type = c("within"), ignore.strand=T)
overlap_GO_df=data.frame(overlap_GO)

# individually extract ensembl gene names from gene pool and xm
gene_pool_filtered=gene_pool[overlap_GO_df$queryHits,]
liftover.gtf_df_gene_filtered=liftover.gtf_df_gene[overlap_GO_df$subjectHits,]

# this is the key for GO term analysis
combine=data.frame(ga=gene_pool_filtered$gene_id, xm=liftover.gtf_df_gene_filtered$gene_id)

#get GO id for male
library(biomaRt)
mart=useDataset("xmaculatus_gene_ensembl", useMart("ensembl"))

#get all go id for universe
ensembl_genes_universe=combine$xm

#用getBM()函数获取注释（如果报错，可能是网络问题，通过连接vpn解决，或者多试几次）
Gene_pool=getBM(filters = "ensembl_gene_id", #已知的注释类型
                attributes = c("ensembl_gene_id","go_id","name_1006","namespace_1003"), #要获取的注释类型
                values = ensembl_genes_universe,
                mart = mart,
                curl = curl::new_handle(timeout = 10000))#选定的数据库的基因组


##male:富集分析
#get all go id for universe
ensembl_genes_universe=combine$xm

#用getBM()函数获取注释（如果报错，可能是网络问题，通过连接vpn解决，或者多试几次）
Gene_pool=getBM(filters = "ensembl_gene_id", #已知的注释类型
                attributes = c("ensembl_gene_id","go_id"), #要获取的注释类型
                values = ensembl_genes_universe,#已知的注释类型的数据，把上面我们通过数据处理得到的ensembl基因序号作为ensembl_gene_id 的值 
                mart = mart)#选定的数据库的基因组

# GO: analysis ####

# overlap_GO genes
setwd("E:/goterm")
cor_male_gene=read.csv('goterm_nonsp_male.csv',header=T)

temp1=list(cor_male_gene)

for (i in 1:length(temp1)) {
  temp_3=temp1[[i]]
  name = paste("group", i, sep="")
  #get all go id for DE genes
  group1_xm=merge(as.data.frame(temp_3),combine,by.x="mR",by.y="ga")###这里会失去部分genes
  if(length(is.na(group1_xm))==0){
    next
  }
  ensembl_genes_group1=group1_xm$xm
  group1.gene=getBM(filters = "ensembl_gene_id", 
                    attributes = c("ensembl_gene_id","go_id"), 
                    values = ensembl_genes_group1, 
                    mart = mart)
  #err = try(group1.gene,TRUE)
  
  # Remove blank entries
  Gene_pool <- Gene_pool[Gene_pool$go_id != '',]
  group1.gene <- group1.gene[group1.gene$go_id != '',]
  
  # convert from table format to list format
  geneID2GO <- by(Gene_pool$go_id,
                  Gene_pool$ensembl_gene_id,
                  function(x) as.character(x))
  myInterestingGenes=by(group1.gene$go_id,
                        group1.gene$ensembl_gene_id,
                        function(x) as.character(x))
  # examine result
  head(geneID2GO)
  head(myInterestingGenes)
  
  correction<-"fdr" #多重检验校正
  geneNames = names(geneID2GO)
  myInterestingGenesNames=names(myInterestingGenes)
  geneList = factor(as.integer(geneNames %in% myInterestingGenesNames))#定义背景基因和感兴趣基因
  #有两个levels：0和1，其中1表示感兴趣的基因
  names(geneList) <- geneNames
  head(geneList,3)
  
  ontology=c("MF","BP","CC")
  for (i in 1:length(ontology)) {
    tgData = new("topGOdata", 
                 ontology = ontology[i], 
                 allGenes = geneList, 
                 annot = annFUN.gene2GO, #基因对应的GO注释读取（从geneID2GO中读取）
                 gene2GO = geneID2GO)
    fisherRes = runTest(tgData, algorithm="classic", statistic="fisher")
    fisherResCor = p.adjust(score(fisherRes), method=correction)
    weightRes = runTest(tgData, algorithm="weight01", statistic="fisher")
    weightResCor = p.adjust(score(weightRes), method=correction)
    allRes = GenTable(tgData, 
                      classic=fisherRes, 
                      weight=weightRes, 
                      orderBy="weight", 
                      ranksOf="classic", 
                      topNodes=200)#提取top 200个GO.如何判断显著富集？fisher.COR<0.1
    allRes$fisher.COR = fisherResCor[allRes$GO.ID]
    allRes$weight.COR = weightResCor[allRes$GO.ID]
    write.csv(allRes, paste(name,ontology[i],"csv",sep="."))}
}

#MF
male_MF=read.table("./group1.MF.csv",header = T,row.names = 1,sep = ",")
male_MF_gene2GO=merge(male_MF,Gene_pool,by.x="GO.ID",by.y="go_id")
male_MF_gene2GO_final=merge(male_MF_gene2GO,combine,by.x="ensembl_gene_id",by.y="xm")
male_MF_gene2GO_final=male_MF_gene2GO_final[,c("ensembl_gene_id","ga","GO.ID","Term","Annotated","Significant","Expected",
                                               "classic","weight","fisher.COR","weight.COR" )]

male_MF_gene2GO_final=merge(male_MF_gene2GO_final,as.data.frame(cor_male_gene),by.x = "ga",by.y = "mR")
male_MF_gene2GO_final=male_MF_gene2GO_final[order(male_MF_gene2GO_final$weight),]
male_MF_gene2GO_final=male_MF_gene2GO_final[which(male_MF_gene2GO_final$weight<0.05),]
length(unique(male_MF_gene2GO_final$ga))#72
length(unique(male_MF_gene2GO_final$GO.ID))#37
write.csv(male_MF_gene2GO_final,"./nonpa_male_MF_gene2GO_final.csv",row.names = F)

#BP
male_BP=read.table("./group1.BP.csv",header = T,row.names = 1,sep = ",")
male_BP_gene2GO=merge(male_BP,Gene_pool,by.x="GO.ID",by.y="go_id")
male_BP_gene2GO_final=merge(male_BP_gene2GO,combine,by.x="ensembl_gene_id",by.y="xm")
male_BP_gene2GO_final=male_BP_gene2GO_final[,c("ensembl_gene_id","ga","GO.ID","Term","Annotated","Significant","Expected",
                                               "classic","weight","fisher.COR","weight.COR" )]

male_BP_gene2GO_final=merge(male_BP_gene2GO_final,as.data.frame(cor_male_gene),by.x = "ga",by.y = "mR")
male_BP_gene2GO_final=male_BP_gene2GO_final[order(male_BP_gene2GO_final$weight),]
male_BP_gene2GO_final=male_BP_gene2GO_final[which(male_BP_gene2GO_final$weight<0.05),]
length(unique(male_BP_gene2GO_final$ga))#42
length(unique(male_BP_gene2GO_final$GO.ID))#38
write.csv(male_BP_gene2GO_final,"./nonpa_male_BP_gene2GO_final.csv",row.names = F)

#CC
male_CC=read.table("./group1.CC.csv",header = T,row.names = 1,sep = ",")
male_CC_gene2GO=merge(male_CC,Gene_pool,by.x="GO.ID",by.y="go_id")
male_CC_gene2GO_final=merge(male_CC_gene2GO,combine,by.x="ensembl_gene_id",by.y="xm")
male_CC_gene2GO_final=male_CC_gene2GO_final[,c("ensembl_gene_id","ga","GO.ID","Term","Annotated","Significant","Expected",
                                               "classic","weight","fisher.COR","weight.COR" )]

male_CC_gene2GO_final=merge(male_CC_gene2GO_final,as.data.frame(cor_male_gene),by.x = "ga",by.y = "mR")
male_CC_gene2GO_final=male_CC_gene2GO_final[order(male_CC_gene2GO_final$weight),]
male_CC_gene2GO_final=male_CC_gene2GO_final[which(male_CC_gene2GO_final$weight<0.05),]
length(unique(male_CC_gene2GO_final$ga))#59
length(unique(male_CC_gene2GO_final$GO.ID))#6
write.csv(male_CC_gene2GO_final,"./nonpa_male_CC_gene2GO_final.csv",row.names = F)

### overlap of parallel express gene vs. parallel splice gene
setwd("E:/goterm")
express=read.table("1.all_parallel_express_go_final.csv",header = T,sep = ",")
splice=read.csv("2.all_parallel_splice_go_final.csv",header = T, quote = "")

# 找到重复的GO.ID
common_GO_IDs <- express$GO.ID[express$GO.ID %in% splice$GO.ID]

# 输出重复的GO.ID
print(common_GO_IDs)#"GO:0004463" "GO:0016231" "GO:0008116" "GO:0106256"

### overlap of nonparallel express gene vs. nonparallel splice gene
setwd("E:/goterm")
express=read.table("8.all_nonparallel_express_go_final.csv",header = T,sep = ",")
splice=read.csv("9.all_nonparallel_splice_go_final.csv",header = T, quote = "")

# 找到重复的GO.ID
common_GO_IDs <- express$GO.ID[express$GO.ID %in% splice$GO.ID]#5

# 输出重复的GO.ID
print(common_GO_IDs)#"GO:0005516" "GO:0016525" "GO:0007417" "GO:0051056" "GO:0014069"
