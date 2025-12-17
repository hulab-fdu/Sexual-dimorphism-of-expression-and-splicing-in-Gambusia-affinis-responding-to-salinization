#Female DESeq2 & WGCNA genes
setwd("E:/filter")
deseq2_genes_female_unique <- read.table("deseq2_genes_female_unique.txt", header = T, sep = "\t", col.names = c("GeneName"))
setwd("E:/filter/wgcna_female")
combined_wgcna_female_file <- read.table("combined_wgcna_female.txt", header = F, sep = "\t", col.names = c("GeneName"))

setwd("E:/filter")
merged_female_genes <- rbind(deseq2_genes_female_unique, combined_wgcna_female_file)
merged_female_genes_unique <- merged_female_genes[!duplicated(merged_female_genes$GeneName), ]
write.table(merged_female_genes_unique, file = "merged_female_genes_unique.txt", sep = "\t", row.names = FALSE, quote = FALSE)
merged_female_genes_unique <- read.table("merged_female_genes_unique.txt", header = T, sep = "\t")


#Male DESeq2 & WGCNA genes
setwd("E:/filter")
deseq2_genes_male_unique <- read.table("deseq2_genes_male_unique.txt", header = T, sep = "\t", col.names = c("GeneName"))
setwd("E:/filter/wgcna_male")
combined_wgcna_male_file <- read.table("combined_wgcna_male.txt", header = F, sep = "\t", col.names = c("GeneName"))

setwd("E:/filter")
merged_male_genes <- rbind(deseq2_genes_male_unique, combined_wgcna_male_file)
merged_male_genes_unique <- merged_male_genes[!duplicated(merged_male_genes$GeneName), ]
write.table(merged_male_genes_unique, file = "merged_male_genes_unique.txt", sep = "\t", row.names = FALSE, quote = FALSE)
merged_male_genes_unique <- read.table("merged_male_genes_unique.txt", header = T, sep = "\t")


#
female_genes <- merged_female_genes_unique$x
male_genes <- merged_male_genes_unique$x


#choose shared genes in female and male
common_genes <- intersect(female_genes, male_genes)
print(common_genes)


write.table(common_genes, file = "shared_genes.txt", sep = "\t", row.names = FALSE, quote = FALSE)
shared_genes <- read.table("shared_genes.txt", header = T, sep = "\t")


#venn for shared genes
library(VennDiagram)
setwd("E:/filter")
list2=list(merged_female_genes_unique$x,merged_male_genes_unique$x)
names(list2) <- c("female","male")

venn = venn.diagram(
  list2,
  filename = "venn_responsive.tiff",
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

#The Hypergeometric Distribution: phyper
phyper(4668, 8861, 8077, 8939, lower.tail = FALSE)#0.5958917


