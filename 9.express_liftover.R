# parallel salinity responsive express gene
# import ensembl Gambusia annotation, extract filtered genes info
setwd("E:/filter")
counts=read.csv('filter_table.txt',header=T,sep="\t")
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
#connect parallel (salinity responsive) express genes with Chromosome
parallel_DEG_genes=read.csv('common_parallel_express_genenames.csv',header=T)
merged_DEG_df <- merge(parallel_DEG_genes, merged_df_chro, by.x = "x", by.y = "ensembl", all = FALSE)
merged_DEG_df_unique <- merged_DEG_df[!duplicated(merged_DEG_df$x), ]#458-455

# genes crossing chromosome, these genes were removed in subsequent analyses
rows_not_in_merged <- parallel_DEG_genes[!parallel_DEG_genes$x %in% merged_DEG_df_unique$x, ]
rows_not_in_merged#3 genes, "ENSGAFG00000010419" "ENSGAFG00000011730" "ENSGAFG00000016716"

write.csv(merged_DEG_df_unique,file="parallel_express_gene_chromosome.csv")

#Calculate the proportion of chromosomes where parallel salinity responsive express gene is located, LG01 is the sex chromosome
LG01_count <- sum(merged_DEG_df_unique$chromosome == "LG01")
LG01_count#28
total_count <- nrow(merged_DEG_df_unique)
total_count#455
LG01_frequency <- LG01_count / total_count
print(LG01_frequency)#0.06153846

##2##
#connect all filtered genes with Chromosome(gene normal)
row_names <- rownames(counts)
gene_normal <- data.frame(x = row_names, row.names = NULL)
merged_gene_normal_df <- merge(gene_normal, merged_df_chro, by.x = "x", by.y = "ensembl", all = FALSE)
merged_gene_normal_df_unique <- merged_gene_normal_df[!duplicated(merged_gene_normal_df$x), ]#16938-16616

# filtered genes crossing chromosome, these genes were removed in subsequent analyses
rows_not_in_merged <- gene_normal[!gene_normal$x %in% merged_gene_normal_df_unique$x, ]
rows_not_in_merged#322 genes
write.csv(rows_not_in_merged,file="filtered_genes_crossing_chromosome.csv")

write.csv(merged_gene_normal_df_unique,file="gene_normal_chromosome.csv")

#Calculate the proportion of chromosomes where filtered genes is located, LG01 is the sex chromosome
LG01_count <- sum(merged_gene_normal_df_unique$chromosome == "LG01")
LG01_count#837
total_count <- nrow(merged_gene_normal_df_unique)
total_count#16616
LG01_frequency <- LG01_count / total_count
print(LG01_frequency)#0.05037313

### fisher's exact test ###
# gene normal vs. parallel express gene
table <- matrix(c(837, 15779, 28, 427), nrow = 2, byrow = TRUE,
                dimnames = list(c("gene normal", "parallel express gene"), c("sex", "non-sex")))
table
#                      sex non-sex#
#gene normal           837   15779#
#parallel express gene  28     427#
result <- fisher.test(table)
result#p-value = 0.278
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:0.5475965 1.2398693
#  odds ratio 0.808951

#bar chart to visualization
library(ggplot2)
library(tidyr)
library(dplyr)

data <- data.frame(
  Group = c("all filtered genes", "parallelly expressed genes"),
  sexchromosome = c(837/16616, 28/455),
  autosome = c(15779/16616, 427/455)
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
    Group == "all filtered genes" & Type == "sex chromosome" ~ "837",
    Group == "all filtered genes" & Type == "autosome" ~ "15779",
    Group == "parallelly expressed genes" & Type == "sex chromosome" ~ "28",
    Group == "parallelly expressed genes" & Type == "autosome" ~ "427"
  )), color = "white", fontface = "bold", size = 5) + 
  scale_fill_manual(values = c("sex chromosome" = "#333333", "autosome" = "#999999")) +
  labs(x = "Gene expression", y = "Proportion") +
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

ggsave("E:/paper_plot/chromosome1.png", device="png",
       #path = "xxx",
       width = 8, height = 8, units="in",#scaling = 0.7,
       dpi=600)


# non-parallelly expressed gene
# import ensembl Gambusia annotation, extract filtered genes info
setwd("E:/filter")
counts=read.csv('filter_table.txt',header=T,sep="\t")
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
#connect non-parallelly expressed genes with Chromosome
non_parallel_DEG_genes=read.csv('non_pa_express_genenames.csv',header=T)
merged_DEG_df <- merge(non_parallel_DEG_genes, merged_df_chro, by.x = "x", by.y = "ensembl", all = FALSE)
merged_DEG_df_unique <- merged_DEG_df[!duplicated(merged_DEG_df$x), ]

# genes crossing chromosome, these genes were removed in subsequent analyses
rows_not_in_merged <- non_parallel_DEG_genes[!non_parallel_DEG_genes$x %in% merged_DEG_df_unique$x, ]
rows_not_in_merged#57 genes


#Calculate the proportion of chromosomes where non-parallelly expressed gene is located, LG01 is the sex chromosome
LG01_count <- sum(merged_DEG_df_unique$chromosome == "LG01")
LG01_count#207
total_count <- nrow(merged_DEG_df_unique)
total_count#4154
LG01_frequency <- LG01_count / total_count
print(LG01_frequency)#0.04983149

##2##
#connect all filtered genes with Chromosome(gene normal)
row_names <- rownames(counts)
gene_normal <- data.frame(x = row_names, row.names = NULL)
merged_gene_normal_df <- merge(gene_normal, merged_df_chro, by.x = "x", by.y = "ensembl", all = FALSE)
merged_gene_normal_df_unique <- merged_gene_normal_df[!duplicated(merged_gene_normal_df$x), ]

# filtered genes crossing chromosome, these genes were removed in subsequent analyses
rows_not_in_merged <- gene_normal[!gene_normal$x %in% merged_gene_normal_df_unique$x, ]
rows_not_in_merged#322 genes

#Calculate the proportion of chromosomes where filtered genes is located, LG01 is the sex chromosome
LG01_count <- sum(merged_gene_normal_df_unique$chromosome == "LG01")
LG01_count#837
total_count <- nrow(merged_gene_normal_df_unique)
total_count#16616
LG01_frequency <- LG01_count / total_count
print(LG01_frequency)#0.05037313
### fisher's exact test ###
# gene normal vs. parallel express gene
table <- matrix(c(837, 15779, 207, 3947), nrow = 2, byrow = TRUE,
                dimnames = list(c("gene normal", "nonparallel express gene"), c("sex", "non-sex")))
table
#                         sex non-sex#
#gene normal              837   15779#
#nonparallel express gene 207    3947#

result <- fisher.test(table)
result#p-value = 0.9053

#bar chart to visualization
library(ggplot2)
library(tidyr)
library(dplyr)

data <- data.frame(
  Group = c("all filtered genes", "non-parallelly expressed genes"),
  sexchromosome = c(837/16616, 207/4154),
  autosome = c(15779/16616, 3947/4154)
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
    Group == "all filtered genes" & Type == "sex chromosome" ~ "837",
    Group == "all filtered genes" & Type == "autosome" ~ "15779",
    Group == "non-parallelly expressed genes" & Type == "sex chromosome" ~ "207",
    Group == "non-parallelly expressed genes" & Type == "autosome" ~ "3947"
  )), color = "white", fontface = "bold", size = 5) + 
  scale_fill_manual(values = c("sex chromosome" = "#333333", "autosome" = "#999999")) +
  labs(x = "Gene expression", y = "Proportion") +
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

ggsave("E:/paper_plot/non_ex_chro.png", device="png",
       #path = "xxx",
       width = 8, height = 8, units="in",#scaling = 0.7,
       dpi=600)
