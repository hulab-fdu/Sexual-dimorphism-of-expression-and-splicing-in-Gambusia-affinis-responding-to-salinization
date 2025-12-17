# parallel salinity responsive splice gene
# import ensembl Gambusia annotation, extract filtered genes info
setwd("E:/exon")
exonnames=read.csv("final_filtered_count.csv",header = T)
exonnames <- exonnames[, -1]
rownames(exonnames) <- exonnames[, 1]
exonnames <- exonnames[, -1]
exonnames_genes <- gsub(":.*$", "", rownames(exonnames))
exonnames_genes
exonnames_genes_unique <- unique(exonnames_genes)#15783
write.csv(exonnames_genes_unique, "filtered_splice_genes.csv", row.names = F)

counts=read.csv('filtered_splice_genes.csv',header=T)
rownames(counts) <- counts[, 1]

setwd("E:/liftoff")
gtf = rtracklayer::import("Gambusia_affinis.ASM309773v1.109.gtf")
gtf_df=as.data.frame(gtf)
gtf_df_gene=gtf_df[gtf_df$type=="gene",] # only need gene anno because htseq-count extract counts of mapped reads on genes
gene_pool=gtf_df_gene[gtf_df_gene$gene_id %in% rownames(counts),]
rownames(gene_pool)=c(1:nrow(gene_pool))
gene_pool_grange=as(gene_pool[,1:3], "GRanges")

# import liftover annotation
liftover.gtf=rtracklayer::import("westmosfish.gtf")
liftover.gtf_df=as.data.frame(liftover.gtf)
liftover.gtf_df_gene=liftover.gtf_df[liftover.gtf_df$type=="gene",]
rownames(liftover.gtf_df_gene)=c(1:nrow(liftover.gtf_df_gene))
liftover.gtf_df_gene_grange=as(liftover.gtf_df_gene[,1:3], "GRanges")

# find overlaps between filtered genes info and liftover annotation
library(GenomicRanges)
overlap=findOverlaps(gene_pool_grange, liftover.gtf_df_gene_grange, type = c("any"), ignore.strand=T)
overlap_df=data.frame(overlap)

# individually extract ensembl gene names from gene pool and lift
gene_pool_filtered=gene_pool[overlap_df$queryHits,]
liftover.gtf_df_gene_filtered=liftover.gtf_df_gene[overlap_df$subjectHits,]

# this is the key for GO term analysis (here, ensembl and lift is homologous)
combine=data.frame(ensembl=gene_pool_filtered$gene_id, lift=liftover.gtf_df_gene_filtered$gene_id)

# import ncbi Gambusia annotation
setwd("E:/liftoff")
gtf_ncbi = rtracklayer::import("GCF_019740435.1_SWU_Gaff_1.0_genomic.gtf")
gtf_ncbi_df=as.data.frame(gtf_ncbi)
gtf_ncbi_df_gene=gtf_ncbi_df[gtf_ncbi_df$type=="gene",]
rownames(gtf_ncbi_df_gene)=c(1:nrow(gtf_ncbi_df_gene))
gtf_ncbi_df_gene_grange=as(gtf_ncbi_df_gene[,1:3], "GRanges")

merged_df <- merge(combine, gtf_ncbi_df_gene, by.x = "lift", by.y = "gene_id",
                   all = FALSE, all.x = FALSE, all.y = FALSE,
                   sort = TRUE, suffixes = c(".com", ".ncbi"), no.dups = FALSE,
                   incomparables = NULL)

#Chromosome correspondence
unique_seqnames <- unique(merged_df$seqnames)
print(length(unique_seqnames))
print(unique_seqnames)
chromosome1=read.csv('chromosome.csv',header=T)
merged_df_chro <- merge(chromosome1, merged_df, by.x = "seqnames", by.y = "seqnames", all.x = TRUE)

##1##
#connect parallel (salinity responsive) splice gene with Chromosome
parallel_DEU_genes=read.csv('common_parallel_splicing_genenames.csv',header=T)
merged_DEU_df <- merge(parallel_DEU_genes, merged_df_chro, by.x = "x", by.y = "ensembl", all = FALSE)
merged_DEU_df_unique <- merged_DEU_df[!duplicated(merged_DEU_df$x), ]#782-778

# genes crossing chromosome, these genes were removed in subsequent analyses
rows_not_in_merged <- parallel_DEU_genes[!parallel_DEU_genes$x %in% merged_DEU_df_unique$x, ]
rows_not_in_merged#4 genes, "ENSGAFG00000004741" "ENSGAFG00000015142" "ENSGAFG00000019011" "ENSGAFG00000020450"

write.csv(merged_DEU_df_unique,file="parallel_splice_gene_chromosome.csv")

#Calculate the proportion of chromosomes where parallel DEU genes is located, LG01 is the sex chromosome
LG01_count <- sum(merged_DEU_df_unique$chromosome == "LG01")
LG01_count#43
total_count <- nrow(merged_DEU_df_unique)
total_count#778
LG01_frequency <- LG01_count / total_count
print(LG01_frequency)#0.05526992

##2##
#connect filtered genes with Chromosome
row_names <- rownames(counts)
gene_normal <- data.frame(x = row_names, row.names = NULL)
merged_gene_normal_df <- merge(gene_normal, merged_df_chro, by.x = "x", by.y = "ensembl", all = FALSE)
merged_gene_normal_df_unique <- merged_gene_normal_df[!duplicated(merged_gene_normal_df$x), ]#15783-15477

# filtered genes crossing chromosome, these genes were removed in subsequent analyses
rows_not_in_merged <- gene_normal[!gene_normal$x %in% merged_gene_normal_df_unique$x, ]
rows_not_in_merged#306 genes
write.csv(rows_not_in_merged,file="filtered_splice_genes_crossing_chromosome.csv")

write.csv(merged_gene_normal_df_unique,file="splice_gene_normal_chromosome.csv")

#Calculate the proportion of chromosomes where filtered genes is located, LG01 is the sex chromosome
LG01_count <- sum(merged_gene_normal_df_unique$chromosome == "LG01")
LG01_count#789
total_count <- nrow(merged_gene_normal_df_unique)
total_count#15477
LG01_frequency <- LG01_count / total_count
print(LG01_frequency)#0.05097887

### fisher's exact test ###
# gene normal vs. parallel salinity response splice gene
table <- matrix(c(789, 14688, 43, 735), nrow = 2, byrow = TRUE,
                dimnames = list(c("gene normal", "parallel splice gene"), c("sex", "non-sex")))
table
#                   sex non-sex#
#gene normal          789   14688#
#parallel splice gene  43     735#
result <- fisher.test(table)
result#p-value = 0.560
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:0.6684821 1.2903967
#odds ratio 0.9182213 

#bar chart to visualization
library(ggplot2)
library(tidyr)
library(dplyr)

data <- data.frame(
  Group = c("all filtered genes", "parallelly spliced genes"),
  sexchromosome = c(789/15477, 43/778),
  autosome = c(14688/15477, 735/778)
)

colnames(data) = c("Group", "sex chromosome", "autosome")

data_long <- data %>%
  gather(key = "Type", value = "Value", -Group)

# 计算标签位置并上调autosome标签
data_long <- data_long %>%
  group_by(Group) %>%
  mutate(label_pos = cumsum(Value) - 0.5 * Value) %>%
  # ★ 关键修改：上调autosome标签位置 ★
  mutate(label_pos = ifelse(Type == "autosome", label_pos + 0.44, label_pos)) %>%
  ungroup()

ggplot(data_long, aes(x = Group, y = Value, fill = Type)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(y = label_pos, label = case_when(
    Group == "all filtered genes" & Type == "sex chromosome" ~ "789",
    Group == "all filtered genes" & Type == "autosome" ~ "14688",
    Group == "parallelly spliced genes" & Type == "sex chromosome" ~ "43",
    Group == "parallelly spliced genes" & Type == "autosome" ~ "735"
  )), color = "white", fontface = "bold", size = 5) + 
  scale_fill_manual(values = c("sex chromosome" = "#333333", "autosome" = "#999999")) +
  labs(x = "Alternative splicing", y = "Proportion") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    # ★ 关键修改：调小Group字体 ★
    axis.text.x = element_text(size = 14), # 从18减小到14
    axis.text.y = element_text(size = 18),
    axis.title.x = element_text(size = 20, face = "bold", margin = margin(t = 15) ),
    axis.title.y = element_text(size = 20, face = "bold"),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.border = element_blank(),
    # ★ 可选：调整图例字体 ★
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  ) +
  guides(fill = guide_legend(override.aes = list(y = c(0.15, 0.85))))

ggsave("E:/paper_plot/chromosome.png", device="png",
       #path = "xxx",
       width = 8, height = 8, units="in",#scaling = 0.7,
       dpi=600)


# non-parallelly spliced gene
# import ensembl Gambusia annotation, extract filtered genes info
setwd("E:/exon")
exonnames=read.csv("final_filtered_count.csv",header = T)
exonnames <- exonnames[, -1]
rownames(exonnames) <- exonnames[, 1]
exonnames <- exonnames[, -1]
exonnames_genes <- gsub(":.*$", "", rownames(exonnames))
exonnames_genes
exonnames_genes_unique <- unique(exonnames_genes)#15783

counts=read.csv('filtered_splice_genes.csv',header=T)
rownames(counts) <- counts[, 1]

setwd("E:/liftoff")
gtf = rtracklayer::import("Gambusia_affinis.ASM309773v1.109.gtf")
gtf_df=as.data.frame(gtf)
gtf_df_gene=gtf_df[gtf_df$type=="gene",] # only need gene anno because htseq-count extract counts of mapped reads on genes
gene_pool=gtf_df_gene[gtf_df_gene$gene_id %in% rownames(counts),]
rownames(gene_pool)=c(1:nrow(gene_pool))
gene_pool_grange=as(gene_pool[,1:3], "GRanges")

# import liftover annotation
liftover.gtf=rtracklayer::import("westmosfish.gtf")
liftover.gtf_df=as.data.frame(liftover.gtf)
liftover.gtf_df_gene=liftover.gtf_df[liftover.gtf_df$type=="gene",]
rownames(liftover.gtf_df_gene)=c(1:nrow(liftover.gtf_df_gene))
liftover.gtf_df_gene_grange=as(liftover.gtf_df_gene[,1:3], "GRanges")

# find overlaps between filtered genes info and liftover annotation
library(GenomicRanges)
overlap=findOverlaps(gene_pool_grange, liftover.gtf_df_gene_grange, type = c("any"), ignore.strand=T)
overlap_df=data.frame(overlap)

# individually extract ensembl gene names from gene pool and lift
gene_pool_filtered=gene_pool[overlap_df$queryHits,]
liftover.gtf_df_gene_filtered=liftover.gtf_df_gene[overlap_df$subjectHits,]

# this is the key for GO term analysis (here, ensembl and lift is homologous)
combine=data.frame(ensembl=gene_pool_filtered$gene_id, lift=liftover.gtf_df_gene_filtered$gene_id)

# import ncbi Gambusia annotation
setwd("E:/liftoff")
gtf_ncbi = rtracklayer::import("GCF_019740435.1_SWU_Gaff_1.0_genomic.gtf")
gtf_ncbi_df=as.data.frame(gtf_ncbi)
gtf_ncbi_df_gene=gtf_ncbi_df[gtf_ncbi_df$type=="gene",]
rownames(gtf_ncbi_df_gene)=c(1:nrow(gtf_ncbi_df_gene))
gtf_ncbi_df_gene_grange=as(gtf_ncbi_df_gene[,1:3], "GRanges")

merged_df <- merge(combine, gtf_ncbi_df_gene, by.x = "lift", by.y = "gene_id",
                   all = FALSE, all.x = FALSE, all.y = FALSE,
                   sort = TRUE, suffixes = c(".com", ".ncbi"), no.dups = FALSE,
                   incomparables = NULL)

#Chromosome correspondence
unique_seqnames <- unique(merged_df$seqnames)
print(length(unique_seqnames))
print(unique_seqnames)
chromosome1=read.csv('chromosome.csv',header=T)
merged_df_chro <- merge(chromosome1, merged_df, by.x = "seqnames", by.y = "seqnames", all.x = TRUE)

##1##
#connect parallelly spliced gene with Chromosome
nonparallel_DEU_genes=read.csv('non_pa_splice_genenames.csv',header=T)
merged_DEU_df <- merge(nonparallel_DEU_genes, merged_df_chro, by.x = "x", by.y = "ensembl", all = FALSE)
merged_DEU_df_unique <- merged_DEU_df[!duplicated(merged_DEU_df$x), ]

# genes crossing chromosome, these genes were removed in subsequent analyses
rows_not_in_merged <- nonparallel_DEU_genes[!nonparallel_DEU_genes$x %in% merged_DEU_df_unique$x, ]
rows_not_in_merged#39 genes


#Calculate the proportion of chromosomes where parallelly DEU genes is located, LG01 is the sex chromosome
LG01_count <- sum(merged_DEU_df_unique$chromosome == "LG01")
LG01_count#258
total_count <- nrow(merged_DEU_df_unique)
total_count#4814
LG01_frequency <- LG01_count / total_count
print(LG01_frequency)#0.05359369

##2##
#connect filtered genes with Chromosome
row_names <- rownames(counts)
gene_normal <- data.frame(x = row_names, row.names = NULL)
merged_gene_normal_df <- merge(gene_normal, merged_df_chro, by.x = "x", by.y = "ensembl", all = FALSE)
merged_gene_normal_df_unique <- merged_gene_normal_df[!duplicated(merged_gene_normal_df$x), ]

# filtered genes crossing chromosome, these genes were removed in subsequent analyses
rows_not_in_merged <- gene_normal[!gene_normal$x %in% merged_gene_normal_df_unique$x, ]
rows_not_in_merged#306 genes


#Calculate the proportion of chromosomes where filtered genes is located, LG01 is the sex chromosome
LG01_count <- sum(merged_gene_normal_df_unique$chromosome == "LG01")
LG01_count#789
total_count <- nrow(merged_gene_normal_df_unique)
total_count#15477
LG01_frequency <- LG01_count / total_count
print(LG01_frequency)#0.05097887

### fisher's exact test ###
# gene normal vs. parallel salinity response splice gene
table <- matrix(c(789, 14688, 258, 4556), nrow = 2, byrow = TRUE,
                dimnames = list(c("gene normal", "non-parallel splice gene"), c("sex", "non-sex")))
table
#                         sex non-sex#
#gene normal              789   14688#
#non-parallel splice gene 258    4556#
result <- fisher.test(table)
result#p-value = 0.4785

#bar chart to visualization
library(ggplot2)
library(tidyr)
library(dplyr)

data <- data.frame(
  Group = c("all filtered genes", "non-parallelly spliced genes"),
  sexchromosome = c(789/15477, 258/4814),
  autosome = c(14688/15477, 4556/4814)
)

colnames(data) = c("Group", "sex chromosome", "autosome")

data_long <- data %>%
  gather(key = "Type", value = "Value", -Group)

# 计算标签位置并上调autosome标签
data_long <- data_long %>%
  group_by(Group) %>%
  mutate(label_pos = cumsum(Value) - 0.5 * Value) %>%
  # ★ 关键修改：上调autosome标签位置 ★
  mutate(label_pos = ifelse(Type == "autosome", label_pos + 0.44, label_pos)) %>%
  ungroup()

ggplot(data_long, aes(x = Group, y = Value, fill = Type)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(y = label_pos, label = case_when(
    Group == "all filtered genes" & Type == "sex chromosome" ~ "789",
    Group == "all filtered genes" & Type == "autosome" ~ "14688",
    Group == "non-parallelly spliced genes" & Type == "sex chromosome" ~ "258",
    Group == "non-parallelly spliced genes" & Type == "autosome" ~ "4556"
  )), color = "white", fontface = "bold", size = 5) + 
  scale_fill_manual(values = c("sex chromosome" = "#333333", "autosome" = "#999999")) +
  labs(x = "Alternative splicing", y = "Proportion") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    # ★ 关键修改：调小Group字体 ★
    axis.text.x = element_text(size = 14), # 从18减小到14
    axis.text.y = element_text(size = 18),
    axis.title.x = element_text(size = 20, face = "bold", margin = margin(t = 15) ),
    axis.title.y = element_text(size = 20, face = "bold"),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.border = element_blank(),
    # ★ 可选：调整图例字体 ★
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  ) +
  guides(fill = guide_legend(override.aes = list(y = c(0.15, 0.85))))

ggsave("E:/paper_plot/chro_non_sp.png", device="png",
       #path = "xxx",
       width = 8, height = 8, units="in",#scaling = 0.7,
       dpi=600)
