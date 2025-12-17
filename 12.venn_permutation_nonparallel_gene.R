# permutation test:overlap genes between nonparallel DEG and genes based on nonparallel DEU####
# load files that contain SIGNIFICANT genes in each experiment. each file contains a single column of gene names
library(tidyr)
library(dplyr)
# set minimum overlap threDEUold
min_overlap=2

# observed overlap of genes: Gradually increase group
setwd("E:/liftoff")
nonparallel_DEG_genes=read.csv('non_pa_express_genenames.csv',header=T)
nonparallel_DSG_genes=read.csv('non_pa_splice_genenames.csv',header=T)

temp=c(nonparallel_DEG_genes$x,nonparallel_DSG_genes$x)
temp.tbl=as.data.frame(table(temp))
commongenes=subset(temp.tbl, Freq>=min_overlap);print(commongenes)
N_obs=length(commongenes$temp);print(N_obs)#1549
write.csv(commongenes,"./overlapnonpa_genes.csv",row.names = F)

# perform permutations
setwd("E:/filter")
merged_female_genes_unique <- read.table("merged_female_genes_unique.txt", header = T, sep = "\t")
merged_male_genes_unique <- read.table("merged_male_genes_unique.txt", header = T, sep = "\t")
merged_all_DEG_genes <- rbind(merged_female_genes_unique, merged_male_genes_unique)
merged_all_DEG_genes <- unique(merged_all_DEG_genes)
setwd("E:/exon/female")
merged_female_exons_unique <- read.table("merged_female_exons_unique.txt", header = T, sep = "\t")
setwd("E:/exon/male")
merged_male_exons_unique <- read.table("merged_male_exons_unique.txt", header = T, sep = "\t")
merged_all_DEU_exons <- rbind(merged_female_exons_unique, merged_male_exons_unique)
merged_all_DEU_exons <-unique(merged_all_DEU_exons)
setwd("E:/exon")
shared_exons <- read.table("shared_exons_final.txt", header = T, sep = "\t")
parallel_DEU <- read.csv("common_exonnames.csv", header = T, sep = "\t")
nonpa_DEU <- shared_exons %>%
  filter(!x %in% parallel_DEU$x)   # 只保留 x 值“不在” parallel_DEU$x 里的行

overlappinggenes=c()
iter=1000

for (N in 1:iter){
  # randomly sample as many as candidate genes in population 1
  permset1=sample(merged_all_DEG_genes$x,length(nonparallel_DEG_genes$x)) 
  # randomly sample as many as candidate genes in population 2
  permset2=sample(merged_all_DEU_exons$x,length(nonpa_DEU$x)) 
  permset2=as.data.frame(permset2)
  permset2=separate (permset2,permset2, into= c ("gene","exon"),sep= ":")
  permset2=permset2$gene
  
  
  temp=c(permset1,unique(permset2))#;print(sort(temp))#,permset3,permset4,permset5
  temp.tbl=as.data.frame(table(temp))#;print(temp.tbl)
  randgenes=subset(temp.tbl, Freq>=min_overlap)#;print(randgenes)
  overlappinggenes[N]=length(randgenes$temp) 
  # number of genes overlapping at least min_overlap times in as many populations in this iteration
}

# permutation test
overlappinggenes=as.data.frame(overlappinggenes)
permdf=subset(overlappinggenes,overlappinggenes>=N_obs)
N_perm=length(permdf$overlappinggenes)
p=N_perm/N
print(paste0("Permutation test p-value = ",p))# Permutation test p-value = 1 #
#A small P-value (less than 0.05) indicates that the number of overlapping genes observed is statistically significant and unlikely to have occurred by chance alone

### plot a histogram/density showing the significance, p=1
ggplot(overlappinggenes, aes(x=overlappinggenes)) + geom_histogram(binwidth=1) +
  geom_vline(xintercept = N_obs,col="red",lty=5,linewidth=1) + theme_bw() +
  annotate("text", label="p = 1") +
  theme_bw()+theme(legend.position = "none")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x =element_text(size=12), 
        axis.title.y=element_text(size=12),
        axis.text.x=element_text(vjust=1,size=12),
        axis.text.y=element_text(vjust=1,size=12))+
  xlab("# of Overlapping Genes")+
  ylab("# of Events")+
  scale_y_continuous(expand = c(0,0),limits = c(0,210))

ggsave("E:/liftoff/permutation_non_DEandDS_overlap.png", device="png",
       #path = "xxx",
       width = 4, height = 4, units="in",#scaling = 0.7,
       dpi=600)

#venn plot to look the overlap
library(VennDiagram)
setwd("E:/liftoff")
list2=list(nonparallel_DEG_genes$x,nonparallel_DSG_genes$x)
names(list2) <- c("non-parallelly expressed genes","non-parallelly spliced genes")

setwd("E:/paper_plot")
venn = venn.diagram(
  list2,
  filename = "venn_nonparallel_DE_DS.tiff",
  fontface = "bold",
  fontfamily = "Arial",
  
  # ★★★ 关键修改部分 ★★★
  # 1. 增大数字标签字号
  cex = 2.4,              # 增大数字大小 (原默认1.0)
  
  # 2. 增大分类标签字号
  cat.cex = 2.4,          # 增大分类标签大小 (原0.3)
  cat.fontface = "bold",
  cat.fontfamily = "Arial",
  
  # 3. 可选：调整标签位置避免重叠
  cat.pos = c(0, 0),      # 水平方向 (0=右, 180=左)
  cat.dist = c(0.02, 0.03), # 标签距离圆圈的距离
  
  fill = c("white", "white"),
  ext.text = FALSE,
  
  # 4. 可选：调整图形尺寸以容纳更大文字
  height = 4000,          # 增加高度
  width = 4000,           # 增加宽度
  resolution = 300        # 提高分辨率
)

