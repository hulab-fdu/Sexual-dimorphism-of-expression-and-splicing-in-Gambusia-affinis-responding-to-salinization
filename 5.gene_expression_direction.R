#1h gene expression direction
library(ggplot2)
setwd("E:/filter")
shared_genes <- read.table("shared_genes.txt", header = T, sep = "\t")
Sal1=read.csv("DEgenes_female_1h.csv",header = T)
Sal2=read.csv("DEgenes_male_1h.csv",header = T)


names(Sal1)[names(Sal1) == 'log2FoldChange'] <- 'log2FoldChange.Sal1'
names(Sal2)[names(Sal2) == 'log2FoldChange'] <- 'log2FoldChange.Sal2'
names(Sal1)[names(Sal1) == 'padj'] <- 'padj.Sal1'
names(Sal2)[names(Sal2) == 'padj'] <- 'padj.Sal2'
names(Sal1)[names(Sal1) == 'baseMean'] <- 'baseMean.Sal1'
names(Sal2)[names(Sal2) == 'baseMean'] <- 'baseMean.Sal2'
names(Sal1)[names(Sal1) == 'X'] <- 'genename.Sal1'
names(Sal2)[names(Sal2) == 'X'] <- 'genename.Sal2'
names(Sal1)[names(Sal1) == 'lfcSE'] <- 'lfcSE.Sal1'
names(Sal2)[names(Sal2) == 'lfcSE'] <- 'lfcSE.Sal2'
names(Sal1)[names(Sal1) == 'stat'] <- 'stat.Sal1'
names(Sal2)[names(Sal2) == 'stat'] <- 'stat.Sal2'
names(Sal1)[names(Sal1) == 'pvalue'] <- 'pvalue.Sal1'
names(Sal2)[names(Sal2) == 'pvalue'] <- 'pvalue.Sal2'


Sal1_Sal2_Merged <- cbind(Sal1, Sal2)
Sal1_Sal2_Merged <- Sal1_Sal2_Merged[Sal1_Sal2_Merged$genename.Sal1 %in% shared_genes$x, ]
Sal1_Sal2_Merged <- Sal1_Sal2_Merged[Sal1_Sal2_Merged$genename.Sal2 %in% shared_genes$x, ]


Sal1_Sal2_Merged$direction[which(Sal1_Sal2_Merged$log2FoldChange.Sal1*Sal1_Sal2_Merged$log2FoldChange.Sal2 >0)] <- "same"
Sal1_Sal2_Merged$direction[which(Sal1_Sal2_Merged$log2FoldChange.Sal1*Sal1_Sal2_Merged$log2FoldChange.Sal2 <0)] <- "opposite"


###count the number of gene in different quadrant
Sal1_Sal2_Merged$quarant=ifelse(Sal1_Sal2_Merged$log2FoldChange.Sal1 >0 & Sal1_Sal2_Merged$log2FoldChange.Sal2 >0, "first",
                          ifelse(Sal1_Sal2_Merged$log2FoldChange.Sal1 <0 & Sal1_Sal2_Merged$log2FoldChange.Sal2 >0, "second",
                                 ifelse(Sal1_Sal2_Merged$log2FoldChange.Sal1 <0 & Sal1_Sal2_Merged$log2FoldChange.Sal2 <0, "third", "fourth")))
nrow(Sal1_Sal2_Merged[Sal1_Sal2_Merged$quarant=="first",])#1321
nrow(Sal1_Sal2_Merged[Sal1_Sal2_Merged$quarant=="second",])#1404
nrow(Sal1_Sal2_Merged[Sal1_Sal2_Merged$quarant=="third",])#1131
nrow(Sal1_Sal2_Merged[Sal1_Sal2_Merged$quarant=="fourth",])#813


#1h gene expression direction
p <- ggplot(Sal1_Sal2_Merged, aes(x=log2FoldChange.Sal1, y=log2FoldChange.Sal2, fill=direction)) +
  geom_hline(yintercept=0, col="grey50") +
  geom_vline(xintercept=0, col="grey50") +
  geom_point(colour = "black",shape = 21,size=3, alpha=1)+
  scale_fill_manual(values=c("white","grey"))+
  xlab(expression("log"["2"]~"FC"~"for salinity-responsive expressed genes in female"))+
  ylab(expression("log"["2"]~"FC"~"for salinity-responsive expressed genes in male"))+
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
  geom_text(x=8.6,y=6,label="1321",size=6)+
  geom_text(x=-4.3,y=6,label="1404",size=6)+
  geom_text(x=-4.3,y=-4.3,label="1131",size=6)+
  geom_text(x=8.6,y=-4.3,label="813",size=6)

setwd("E:/paper_plot")
ggsave("plot1h.png", p, width = 10, height = 10, dpi = 300)


#------------------------------------#


#24h gene expression direction
setwd("E:/filter")
shared_genes <- read.table("shared_genes.txt", header = T, sep = "\t")
Sal1=read.csv("DEgenes_female_24h.csv",header = T)
Sal2=read.csv("DEgenes_male_24h.csv",header = T)


names(Sal1)[names(Sal1) == 'log2FoldChange'] <- 'log2FoldChange.Sal1'
names(Sal2)[names(Sal2) == 'log2FoldChange'] <- 'log2FoldChange.Sal2'
names(Sal1)[names(Sal1) == 'padj'] <- 'padj.Sal1'
names(Sal2)[names(Sal2) == 'padj'] <- 'padj.Sal2'
names(Sal1)[names(Sal1) == 'baseMean'] <- 'baseMean.Sal1'
names(Sal2)[names(Sal2) == 'baseMean'] <- 'baseMean.Sal2'
names(Sal1)[names(Sal1) == 'X'] <- 'genename.Sal1'
names(Sal2)[names(Sal2) == 'X'] <- 'genename.Sal2'
names(Sal1)[names(Sal1) == 'lfcSE'] <- 'lfcSE.Sal1'
names(Sal2)[names(Sal2) == 'lfcSE'] <- 'lfcSE.Sal2'
names(Sal1)[names(Sal1) == 'stat'] <- 'stat.Sal1'
names(Sal2)[names(Sal2) == 'stat'] <- 'stat.Sal2'
names(Sal1)[names(Sal1) == 'pvalue'] <- 'pvalue.Sal1'
names(Sal2)[names(Sal2) == 'pvalue'] <- 'pvalue.Sal2'


Sal1_Sal2_Merged <- cbind(Sal1, Sal2)
Sal1_Sal2_Merged <- Sal1_Sal2_Merged[Sal1_Sal2_Merged$genename.Sal1 %in% shared_genes$x, ]
Sal1_Sal2_Merged <- Sal1_Sal2_Merged[Sal1_Sal2_Merged$genename.Sal2 %in% shared_genes$x, ]


Sal1_Sal2_Merged$direction[which(Sal1_Sal2_Merged$log2FoldChange.Sal1*Sal1_Sal2_Merged$log2FoldChange.Sal2 >0)] <- "same"
Sal1_Sal2_Merged$direction[which(Sal1_Sal2_Merged$log2FoldChange.Sal1*Sal1_Sal2_Merged$log2FoldChange.Sal2 <0)] <- "opposite"


###count the number of gene in different quadrant
Sal1_Sal2_Merged$quarant=ifelse(Sal1_Sal2_Merged$log2FoldChange.Sal1 >0 & Sal1_Sal2_Merged$log2FoldChange.Sal2 >0, "first",
                                ifelse(Sal1_Sal2_Merged$log2FoldChange.Sal1 <0 & Sal1_Sal2_Merged$log2FoldChange.Sal2 >0, "second",
                                       ifelse(Sal1_Sal2_Merged$log2FoldChange.Sal1 <0 & Sal1_Sal2_Merged$log2FoldChange.Sal2 <0, "third", "fourth")))
nrow(Sal1_Sal2_Merged[Sal1_Sal2_Merged$quarant=="first",])#1630
nrow(Sal1_Sal2_Merged[Sal1_Sal2_Merged$quarant=="second",])#1545
nrow(Sal1_Sal2_Merged[Sal1_Sal2_Merged$quarant=="third",])#1020
nrow(Sal1_Sal2_Merged[Sal1_Sal2_Merged$quarant=="fourth",])#474


#24h gene expression direction
p <- ggplot(Sal1_Sal2_Merged, aes(x=log2FoldChange.Sal1, y=log2FoldChange.Sal2, fill=direction)) +
  geom_hline(yintercept=0, col="grey50") +
  geom_vline(xintercept=0, col="grey50") +
  geom_point(colour = "black",shape = 21,size=3, alpha=1)+
  scale_fill_manual(values=c("white","grey"))+
  xlab(expression("log"["2"]~"FC"~"for salinity-responsive expressed genes in female"))+
  ylab(expression("log"["2"]~"FC"~"for salinity-responsive expressed genes in male"))+
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
  ggtitle("24 h") +  # 关键添加
  geom_text(x=6.5,y=6.8,label="1630",size=6)+
  geom_text(x=-10.6,y=6.8,label="1545",size=6)+
  geom_text(x=-10.6,y=-10.8,label="1020",size=6)+
  geom_text(x=6.5,y=-10.8,label="474",size=6)

  setwd("E:/paper_plot")
  ggsave("plot24h.png", p, width = 10, height = 10, dpi = 300)


#------------------------------------#


#48h gene expression direction
setwd("E:/filter")
shared_genes <- read.table("shared_genes.txt", header = T, sep = "\t")
Sal1=read.csv("DEgenes_female_48h.csv",header = T)
Sal2=read.csv("DEgenes_male_48h.csv",header = T)


names(Sal1)[names(Sal1) == 'log2FoldChange'] <- 'log2FoldChange.Sal1'
names(Sal2)[names(Sal2) == 'log2FoldChange'] <- 'log2FoldChange.Sal2'
names(Sal1)[names(Sal1) == 'padj'] <- 'padj.Sal1'
names(Sal2)[names(Sal2) == 'padj'] <- 'padj.Sal2'
names(Sal1)[names(Sal1) == 'baseMean'] <- 'baseMean.Sal1'
names(Sal2)[names(Sal2) == 'baseMean'] <- 'baseMean.Sal2'
names(Sal1)[names(Sal1) == 'X'] <- 'genename.Sal1'
names(Sal2)[names(Sal2) == 'X'] <- 'genename.Sal2'
names(Sal1)[names(Sal1) == 'lfcSE'] <- 'lfcSE.Sal1'
names(Sal2)[names(Sal2) == 'lfcSE'] <- 'lfcSE.Sal2'
names(Sal1)[names(Sal1) == 'stat'] <- 'stat.Sal1'
names(Sal2)[names(Sal2) == 'stat'] <- 'stat.Sal2'
names(Sal1)[names(Sal1) == 'pvalue'] <- 'pvalue.Sal1'
names(Sal2)[names(Sal2) == 'pvalue'] <- 'pvalue.Sal2'


Sal1_Sal2_Merged <- cbind(Sal1, Sal2)
Sal1_Sal2_Merged <- Sal1_Sal2_Merged[Sal1_Sal2_Merged$genename.Sal1 %in% shared_genes$x, ]
Sal1_Sal2_Merged <- Sal1_Sal2_Merged[Sal1_Sal2_Merged$genename.Sal2 %in% shared_genes$x, ]


Sal1_Sal2_Merged$direction[which(Sal1_Sal2_Merged$log2FoldChange.Sal1*Sal1_Sal2_Merged$log2FoldChange.Sal2 >0)] <- "same"
Sal1_Sal2_Merged$direction[which(Sal1_Sal2_Merged$log2FoldChange.Sal1*Sal1_Sal2_Merged$log2FoldChange.Sal2 <0)] <- "opposite"


###count the number of gene in different quadrant
Sal1_Sal2_Merged$quarant=ifelse(Sal1_Sal2_Merged$log2FoldChange.Sal1 >0 & Sal1_Sal2_Merged$log2FoldChange.Sal2 >0, "first",
                                ifelse(Sal1_Sal2_Merged$log2FoldChange.Sal1 <0 & Sal1_Sal2_Merged$log2FoldChange.Sal2 >0, "second",
                                       ifelse(Sal1_Sal2_Merged$log2FoldChange.Sal1 <0 & Sal1_Sal2_Merged$log2FoldChange.Sal2 <0, "third", "fourth")))
nrow(Sal1_Sal2_Merged[Sal1_Sal2_Merged$quarant=="first",])#1022
nrow(Sal1_Sal2_Merged[Sal1_Sal2_Merged$quarant=="second",])#880
nrow(Sal1_Sal2_Merged[Sal1_Sal2_Merged$quarant=="third",])#1138
nrow(Sal1_Sal2_Merged[Sal1_Sal2_Merged$quarant=="fourth",])#1629


#48h gene expression direction
p <- ggplot(Sal1_Sal2_Merged, aes(x=log2FoldChange.Sal1, y=log2FoldChange.Sal2, fill=direction)) +
  geom_hline(yintercept=0, col="grey50") +
  geom_vline(xintercept=0, col="grey50") +
  geom_point(colour = "black",shape = 21,size=3, alpha=1)+
  scale_fill_manual(values=c("white","grey"))+
  xlab(expression("log"["2"]~"FC"~"for salinity-responsive expressed genes in female"))+
  ylab(expression("log"["2"]~"FC"~"for salinity-responsive expressed genes in male"))+
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
  ggtitle("48 h") +  # 关键添加
  geom_text(x=4.1,y=5,label="1022",size=6)+
  geom_text(x=-5,y=5,label="880",size=6)+
  geom_text(x=-5,y=-6.2,label="1138",size=6)+
  geom_text(x=4.1,y=-6.2,label="1629",size=6)

setwd("E:/paper_plot")
ggsave("plot48h.png", p, width = 10, height = 10, dpi = 300)


#------------------------------------#


#72h gene expression direction
setwd("E:/filter")
shared_genes <- read.table("shared_genes.txt", header = T, sep = "\t")
Sal1=read.csv("DEgenes_female_72h.csv",header = T)
Sal2=read.csv("DEgenes_male_72h.csv",header = T)


names(Sal1)[names(Sal1) == 'log2FoldChange'] <- 'log2FoldChange.Sal1'
names(Sal2)[names(Sal2) == 'log2FoldChange'] <- 'log2FoldChange.Sal2'
names(Sal1)[names(Sal1) == 'padj'] <- 'padj.Sal1'
names(Sal2)[names(Sal2) == 'padj'] <- 'padj.Sal2'
names(Sal1)[names(Sal1) == 'baseMean'] <- 'baseMean.Sal1'
names(Sal2)[names(Sal2) == 'baseMean'] <- 'baseMean.Sal2'
names(Sal1)[names(Sal1) == 'X'] <- 'genename.Sal1'
names(Sal2)[names(Sal2) == 'X'] <- 'genename.Sal2'
names(Sal1)[names(Sal1) == 'lfcSE'] <- 'lfcSE.Sal1'
names(Sal2)[names(Sal2) == 'lfcSE'] <- 'lfcSE.Sal2'
names(Sal1)[names(Sal1) == 'stat'] <- 'stat.Sal1'
names(Sal2)[names(Sal2) == 'stat'] <- 'stat.Sal2'
names(Sal1)[names(Sal1) == 'pvalue'] <- 'pvalue.Sal1'
names(Sal2)[names(Sal2) == 'pvalue'] <- 'pvalue.Sal2'


Sal1_Sal2_Merged <- cbind(Sal1, Sal2)
Sal1_Sal2_Merged <- Sal1_Sal2_Merged[Sal1_Sal2_Merged$genename.Sal1 %in% shared_genes$x, ]
Sal1_Sal2_Merged <- Sal1_Sal2_Merged[Sal1_Sal2_Merged$genename.Sal2 %in% shared_genes$x, ]


Sal1_Sal2_Merged$direction[which(Sal1_Sal2_Merged$log2FoldChange.Sal1*Sal1_Sal2_Merged$log2FoldChange.Sal2 >0)] <- "same"
Sal1_Sal2_Merged$direction[which(Sal1_Sal2_Merged$log2FoldChange.Sal1*Sal1_Sal2_Merged$log2FoldChange.Sal2 <0)] <- "opposite"


###count the number of gene in different quadrant
Sal1_Sal2_Merged$quarant=ifelse(Sal1_Sal2_Merged$log2FoldChange.Sal1 >0 & Sal1_Sal2_Merged$log2FoldChange.Sal2 >0, "first",
                                ifelse(Sal1_Sal2_Merged$log2FoldChange.Sal1 <0 & Sal1_Sal2_Merged$log2FoldChange.Sal2 >0, "second",
                                       ifelse(Sal1_Sal2_Merged$log2FoldChange.Sal1 <0 & Sal1_Sal2_Merged$log2FoldChange.Sal2 <0, "third", "fourth")))
nrow(Sal1_Sal2_Merged[Sal1_Sal2_Merged$quarant=="first",])#898
nrow(Sal1_Sal2_Merged[Sal1_Sal2_Merged$quarant=="second",])#1488
nrow(Sal1_Sal2_Merged[Sal1_Sal2_Merged$quarant=="third",])#1565
nrow(Sal1_Sal2_Merged[Sal1_Sal2_Merged$quarant=="fourth",])#718


#72h gene expression direction
p <- ggplot(Sal1_Sal2_Merged, aes(x=log2FoldChange.Sal1, y=log2FoldChange.Sal2, fill=direction)) +
  geom_hline(yintercept=0, col="grey50") +
  geom_vline(xintercept=0, col="grey50") +
  geom_point(colour = "black",shape = 21,size=3, alpha=1)+
  scale_fill_manual(values=c("white","grey"))+
  xlab(expression("log"["2"]~"FC"~"for salinity-responsive expressed genes in female"))+
  ylab(expression("log"["2"]~"FC"~"for salinity-responsive expressed genes in male"))+
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
  ggtitle("72 h") +  # 关键添加
  geom_text(x=6.2,y=8.7,label="898",size=6)+
  geom_text(x=-5.2,y=8.7,label="1488",size=6)+
  geom_text(x=-5.2,y=-4.3,label="1565",size=6)+
  geom_text(x=6.2,y=-4.3,label="718",size=6)

setwd("E:/paper_plot")
ggsave("plot72h.png", p, width = 10, height = 10, dpi = 300)

#
#g test
#1h
library(RVAideMemoire)

observed = c(2452, 2217)    # observed frequencies
expected = c(0.5, 0.5)      # expected proportions

G.test(x=observed,
       p=expected)
#G = 11.833, df = 1, p-value = 0.0005819


#24h
observed = c(2650, 2019)    # observed frequencies
expected = c(0.5, 0.5)      # expected proportions

G.test(x=observed,
       p=expected)
#G = 85.539, df = 1, p-value < 2.2e-16


#48h
observed = c(2160, 2509)    # observed frequencies
expected = c(0.5, 0.5)      # expected proportions

G.test(x=observed,
       p=expected)
#G = 26.112, df = 1, p-value = 3.223e-07


#72h
observed = c(2463, 2206)    # observed frequencies
expected = c(0.5, 0.5)      # expected proportions

G.test(x=observed,
       p=expected)
#G = 14.153, df = 1, p-value = 0.0001685
