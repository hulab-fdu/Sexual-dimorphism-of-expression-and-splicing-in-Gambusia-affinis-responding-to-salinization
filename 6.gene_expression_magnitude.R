#gene expression boxplot_1h
library(ggplot2)
library(ggpubr)
library(Rmisc)


Logfc_female_1h <- as.data.frame(cbind(Sal1_Sal2_Merged[["genename.Sal1"]], Sal1_Sal2_Merged[["log2FoldChange.Sal1"]]))
Logfc_male_1h <- as.data.frame(cbind(Sal1_Sal2_Merged[["genename.Sal2"]], Sal1_Sal2_Merged[["log2FoldChange.Sal2"]]))


names(Logfc_female_1h) <- c("genename", "log2FoldChange")
names(Logfc_male_1h) <- c("genename", "log2FoldChange")


#
library(ggplot2)
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
  summarise(mean_value = mean(log2FoldChange, na.rm = TRUE))#平均值 F:0.4438161, M:0.6894090


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
print(wilcox_result)#W = 8301319, p-value < 2.2e-16

ggsave("E:/paper_plot/boxplot_1h.png", device="png",
       #path = "xxx",
       width = 8, height = 4, units="in",#scaling = 0.7,
       dpi=600)
#------------------------------#


#gene expression boxplot_24h
Logfc_female_24h <- as.data.frame(cbind(Sal1_Sal2_Merged[["genename.Sal1"]], Sal1_Sal2_Merged[["log2FoldChange.Sal1"]]))
Logfc_male_24h <- as.data.frame(cbind(Sal1_Sal2_Merged[["genename.Sal2"]], Sal1_Sal2_Merged[["log2FoldChange.Sal2"]]))


names(Logfc_female_24h) <- c("genename", "log2FoldChange")
names(Logfc_male_24h) <- c("genename", "log2FoldChange")


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
  summarise(mean_value = mean(log2FoldChange, na.rm = TRUE))#平均值 F:0.6547953, M:0.8074185


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
print(wilcox_result)#W = 9187287, p-value < 2.2e-16

ggsave("E:/paper_plot/boxplot_24h.png", device="png",
       #path = "xxx",
       width = 8, height = 4, units="in",#scaling = 0.7,
       dpi=600)
#------------------------------#


#gene expression boxplot_48h
Logfc_female_48h <- as.data.frame(cbind(Sal1_Sal2_Merged[["genename.Sal1"]], Sal1_Sal2_Merged[["log2FoldChange.Sal1"]]))
Logfc_male_48h <- as.data.frame(cbind(Sal1_Sal2_Merged[["genename.Sal2"]], Sal1_Sal2_Merged[["log2FoldChange.Sal2"]]))


names(Logfc_female_48h) <- c("genename", "log2FoldChange")
names(Logfc_male_48h) <- c("genename", "log2FoldChange")


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
  summarise(mean_value = mean(log2FoldChange, na.rm = TRUE))#平均值 F:0.4809335, M:0.4024116


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
print(wilcox_result)#W = 11836151, p-value = 6.5e-13

ggsave("E:/paper_plot/boxplot_48h.png", device="png",
       #path = "xxx",
       width = 8, height = 4, units="in",#scaling = 0.7,
       dpi=600)
#------------------------------#


#gene expression boxplot_72h
Logfc_female_72h <- as.data.frame(cbind(Sal1_Sal2_Merged[["genename.Sal1"]], Sal1_Sal2_Merged[["log2FoldChange.Sal1"]]))
Logfc_male_72h <- as.data.frame(cbind(Sal1_Sal2_Merged[["genename.Sal2"]], Sal1_Sal2_Merged[["log2FoldChange.Sal2"]]))


names(Logfc_female_72h) <- c("genename", "log2FoldChange")
names(Logfc_male_72h) <- c("genename", "log2FoldChange")


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
  summarise(mean_value = mean(log2FoldChange, na.rm = TRUE))#平均值 F:0.6521015, M:0.5450619


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
print(wilcox_result)#W = 11408913, p-value = 9.3e-05

ggsave("E:/paper_plot/boxplot_72h.png", device="png",
       #path = "xxx",
       width = 8, height = 4, units="in",#scaling = 0.7,
       dpi=600)
