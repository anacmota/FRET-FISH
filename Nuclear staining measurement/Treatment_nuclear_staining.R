setwd("/Users/anamota/Nuclear staining measurement/Intensities datasets")

require(data.table)
require(ggplot2)
library(tidyverse)
library(ggpubr)

############### Ezh2i H3K27me3
control = read.delim("Ezh2i_control_H3K27me3.csv", header = T, as.is = T, sep =",")
treat = read.delim("Ezh2i_treated_H3K27me3.csv", header = T, as.is = T, sep =",")

pdata <- NULL
pdata$root <- c(rep("control", length(control$Mean)), rep("+Ezh2i", length(treat$Mean[treat$Mean < max(treat$Mean)])))
pdata$mean <- c(control$Mean, treat$Mean[treat$Mean < max(treat$Mean)])
pdata$root <- factor(pdata$root, levels = c("control","+Ezh2i"))
pdata = as.data.frame(pdata)

countX <- NULL
countX$root = c("control","+Ezh2i")
countX$ypos = c(5000, 5000)
countX$count = c(length(which(pdata$root == "control")), length(which(pdata$root == "+Ezh2i")))
countX$root <- factor(countX$root, levels = c("control","+Ezh2i"))
countX = as.data.frame(countX)

my_comparisons <- list( c("control","+Ezh2i"))

ggplot(pdata, aes(x=root, y=mean, fill=root))+
  geom_boxplot(width=0.5, outlier.shape=NA) + #geom_violin(trim=T)+
  scale_fill_manual(values = c("#00AFBB", "#FC4E07")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text=element_text(family="Helvetica", size=7), legend.position="none", axis.text.x=element_text(angle=45, hjust=1)) +
  labs(x="", y="Average nuclear intensity (a.u.)") +
  geom_text(countX, mapping = aes(x = root, y = ypos, label= paste("n = ", count)), size = 2) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", family="Helvetica", size=2) # Add pairwise comparisons p-value

#ggsave(IF, file="IF_Ezh2i_H3K27me3.png", dpi=300, width=3, height=2.5)


############### ATP depletion 1h30
control = read.delim("ATP_control_RNA.csv", header = T, as.is = T, sep =",")
treat = read.delim("ATP_treated_1h30_RNA.csv", header = T, as.is = T, sep =",")

pdata <- NULL
pdata$root <- c(rep("control", length(control$Mean)), rep("ATP dep", length(treat$Mean)))
pdata$mean <- c(control$Mean, treat$Mean)
pdata$root <- factor(pdata$root, levels = c("control","ATP dep"))
pdata = as.data.frame(pdata)

countX <- NULL
countX$root = c("control","ATP dep")
countX$ypos = c(12000, 12000)
countX$count = c(length(which(pdata$root == "control")), length(which(pdata$root == "ATP dep")))
countX$root <- factor(countX$root, levels = c("control","ATP dep"))
countX = as.data.frame(countX)

my_comparisons <- list( c("control","ATP dep"))

ggplot(pdata, aes(x=root, y=mean, fill=root))+
  geom_boxplot(width=0.5, outlier.shape=NA) + #geom_violin(trim=T)+
  scale_fill_manual(values = c("#00AFBB", "#FC4E07")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text=element_text(family="Helvetica", size=7), legend.position="none", axis.text.x=element_text(angle=45, hjust=1)) +
  labs(x="", y="Average nuclear intensity (a.u.)") +
  geom_text(countX, mapping = aes(x = root, y = ypos, label= paste("n = ", count)), size = 2) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", family="Helvetica", size=2) # Add pairwise comparisons p-value


#ggsave(IF, file="ATP_RNA.png", dpi=300, width=3, height=2.5)