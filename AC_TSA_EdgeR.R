############################################################
#####################TSA AC EdgeR DE########################
############################################################
library(EdgeR)
counts <- read.delim("~/GDrive files/foxh1 and hdac1/hdac1 rna seq/ectoderm/ectoderm_raw_counts.txt", header=T)
head(counts)
rownames(counts) <- counts[,1]
counts <- counts[,-1]
colnames(counts) <- c("dms1", "dms2", "tsa1", "tsa2")
group <- factor(substring(colnames(counts), 1, 3))
group <- relevel (group, ref= "dms")
rep <- factor(substring(colnames(counts), 4, 4))
y <- DGEList(counts= counts,group= group)
keep <- filterByExpr(y)
table(keep)
y<- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y$samples
######data explore####
plotMDS(y, col= rep(1:2, each= 2))
########notice the batch effect from PCA##############
##########test log2-fold-change for consistency between reps##########
design <- model.matrix(~rep+rep:group)
logFC <- predFC(y, design, prior.count=1, dispersion=0.05)
cor(logFC[, 2:4])
#########design matrix to an additive linear model##########
design <- model.matrix(~rep+group)
rownames(design) <- colnames(y)
design
#####estimate dispersion##########
y <- estimateDisp(y, design, robust= TRUE)
y$common.dispersion
plotBCV(y)
fit <- glmQLFit(y, design, robust=TRUE)
plotQLDisp(fit)
#########differential expression#########
qlf <- glmQLFTest(fit)
topTags(qlf)
#Here's a closer look at the individual counts-per-million for the top genes. The top genes are
#very consistent across the replicates:
top <- rownames(topTags(qlf))
cpm(y)[top,]
summary(decideTests(qlf))######check how many are less than 0.05 FDR####
FDR <- as.data.frame(p.adjust(qlf$table$PValue, method= "BH"))
colnames(FDR) <- c("FDR")
sum(FDR <0.05)##number should match summary above####
DElist <- as.data.frame(qlf$table)
DElist <- cbind(DElist, FDR)
DElist1 <- DElist[intersect(rownames(DElist[which(abs(DElist$logFC)> 1),]),                                               
                           rownames(DElist[which(DElist$FDR< 0.05),])),]
DEgene <- rownames(DElist1)
downgenes <- DElist1[which(DElist1$logFC<0),]
write.table(downgenes,file="~/GDrive files/foxh1 and hdac1/hdac1 rna seq/ectoderm/ac_tsa_downDEgenes.txt", quote = FALSE, sep = "\t")
upgenes <- DElist1[which(DElist1$logFC>0),]
write.table(upgenes,file="~/GDrive files/foxh1 and hdac1/hdac1 rna seq/ectoderm/ac_tsa_upDEgenes.txt", quote = FALSE, sep = "\t")
#####################################################################
###########################volcano plot##############################
#####################################################################
geneName <- as.data.frame(convert[match(rownames(DElist), convert$V2),4])
geneList <- cbind(geneName, DElist)
b <- c(1,2,5)
geneVPtable <- geneList[,b]
colnames(geneVPtable) <- c("geneName", "logFC", "FDR")
geneOfinterest <- c("foxi4.1", "ventx1.1", "ventx1.2", "ascl1", "bmp4", 
                    "nodal3.1", "sia2", "tbxt.2", "sox17a", "foxa4")
geneOfinterestList <- geneVPtable[match(geneOfinterest, geneVPtable$geneName),]
library(ggplot2)
library(ggrepel)
library(dplyr)
volTable <- geneVPtable %>%
  mutate(threshold = factor(case_when(logFC > 1 & FDR < 0.05 ~ "up-regulated",
                                      logFC < -1 & FDR < 0.05 ~ "down-regulated",
                                      TRUE ~ "n.s")))
vp <- ggplot(data=volTable, aes(x=logFC, y=-log10(FDR))) + 
  geom_point(aes(color=threshold), alpha= 0.5, size= 2)+ 
  geom_vline(xintercept=c(-0.58,0.58), color="green", lty= 2)+ 
  geom_hline(yintercept=1.3, color="green", lty= 2)+ 
  xlab("logFC")+ 
  ylab("-log10(FDR)")+ 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA))+
  theme(aspect.ratio = 1.5)+
  theme(axis.text = element_text(color="black"), text = element_text(size=9))+
  theme(panel.background = element_rect(size= 0.5, fill = "white", colour = "black"))+
  scale_color_manual(name = "Genes",
                     values = c("up-regulated" = "salmon", 
                                "down-regulated" = "cornflowerblue", 
                                "n.s" = "grey"))+
  geom_text_repel(data = geneOfinterestList, 
                  aes(label = geneName), fontface = "italic", size= 3,
                  color = "black", min.segment.length = 0, box.padding = 0.5,
                  segment.linetype = 6, 
                  max.overlaps = 20)
vp
ggsave('C:/Users/latro/Desktop/foxh1 and hdac1/hdac1 rna seq/ectoderm/DE_volcanoPlot.tiff', vp, device = "tiff", dpi = 500)
################################################################
#############DE gene expression in germ layers##################
################################################################
dissect_counts <- read.delim("C:/Users/latro/Desktop/foxh1 and hdac1/hdac1 rna seq/ira_dissect_raw_counts.txt", header=F)
rownames(dissect_counts)= dissect_counts[,1]
dissect_counts <- dissect_counts[,-1]
dissect_counts <- dissect_counts[rowSums(dissect_counts)>0, ]
dc_name <- as.data.frame(convert[match(rownames(dissect_counts), convert$V2),4])
dissect_counts1 <- cbind(dc_name, dissect_counts)
rownames(dissect_counts1)= dissect_counts1[,1]
dissect_counts1 <- dissect_counts1[,-1]
colnames(dissect_counts1) <- c("ac1", "ac2", "dm1", "dm2", "lm1", "lm2",
                               "vm1", "vm2", "vg1", "vg2")
dissect_counts2 <- as.data.frame(matrix(nrow= nrow(dissect_counts1), ncol= 6))
rownames(dissect_counts2) <- rownames(dissect_counts1)
dissect_counts2[,1:2] <- dissect_counts1[,1:2]
dissect_counts2[,5:6] <- dissect_counts1[,9:10]
##########merge 3 mesoderm together################
m1 <- c(3,5,7)
dissect_counts2[,3] <- as.integer(rowMeans(dissect_counts1[,m1]))
m2 <- c(4,6,8)
dissect_counts2[,4] <- as.integer(rowMeans(dissect_counts1[,m2]))
colnames(dissect_counts2) <- c("ac1", "ac2", "me1", "me2", "vg1", "vg2")
###################################################
dissect_y <- DGEList(counts= dissect_counts2)
dissect_y <- calcNormFactors(dissect_y)
dissect_y$samples
logDC <- cpm(dissect_y, log = TRUE, prior.count = 1)
############remove any NA##############
logDC <- logDC[complete.cases(logDC),]
############remove any infinity##########
logDC <- logDC[!is.infinite(rowSums(logDC)),]
############rep1 dissection data#########
rep1 <- c(1,3,5)
logDC_rep1 <- logDC[, rep1]
############remove row std=0#############
logDC_rep1 <- logDC_rep1[!apply(logDC_rep1 , 1 , function(x) sd(x)==0 ), ]
############rep2 dissection data#########
rep2 <- c(2,4,6)
logDC_rep2 <- logDC[, rep2]
############remove row std=0#############
logDC_rep2 <- logDC_rep2[!apply(logDC_rep2 , 1 , function(x) sd(x)==0 ), ]
##########up-regulated genes#############
uprld1 <- logDC_rep1[match(upgeneNames, rownames(logDC_rep1)),]
uprld1 <- uprld1[complete.cases(uprld1),]
uprld2 <- logDC_rep2[match(upgeneNames, rownames(logDC_rep2)),]
uprld2 <- uprld2[complete.cases(uprld2),]
##########heatmap###########
library(pheatmap)
library("RColorBrewer")
ramp<-1:3/3 
cols<-c(rgb(ramp,0,0),rgb(0,ramp,0),rgb(0,0,ramp),rgb(ramp,0,ramp)) 
rep1 <- c(1, 3, 5)
rep2 <- c(2, 4, 6)
rep1av <- c(1, 9)
rep2av <- c(2, 10)
av <- c(1,2,9,10)
hp <- pheatmap(         uprld1, 
                        scale = "row",
                        cutree_rows = 1, 
                        cluster_cols = F, 
                        show_rownames = F,
                        clustering_method = "ward.D2", 
                        col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(255))
ggsave('C:/Users/latro/Desktop/foxh1 and hdac1/hdac1 rna seq/ectoderm/ectoderm_upGene_hp.tiff', hp, device = "tiff", dpi = 500)
##########down-regulated genes##############
downrld1 <- logDC_rep1[match(downgeneNames, rownames(logDC_rep1)),]
downrld1 <- downrld1[complete.cases(downrld1),]
downrld2 <- logDC_rep2[match(downgeneNames, rownames(logDC_rep2)),]
downrld2 <- downrld2[complete.cases(downrld2),]
##########heatmap###########
library(pheatmap)
library("RColorBrewer")
ramp<-1:3/3 
cols<-c(rgb(ramp,0,0),rgb(0,ramp,0),rgb(0,0,ramp),rgb(ramp,0,ramp)) 
rep1 <- c(1, 3, 5)
rep2 <- c(2, 4, 6)
rep1av <- c(1, 9)
rep2av <- c(2, 10)
av <- c(1,2,9,10)
hp <- pheatmap(         downrld1, 
                        scale = "row",
                        cutree_rows = 1, 
                        cluster_cols = F, 
                        show_rownames = F,
                        clustering_method = "ward.D2", 
                        col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(255))
ggsave('C:/Users/latro/Desktop/foxh1 and hdac1/hdac1 rna seq/ectoderm/ectoderm_downGene_hp.tiff', hp, device = "tiff", dpi = 500)
##################################################
###################violin plot####################
##################################################
dissect_tpm <- read.delim("C:/Users/latro/Desktop/foxh1 and hdac1/hdac1 rna seq/TPM_table_ira_dissected_st105_tissues.txt", header=F)
rownames(dissect_tpm)= dissect_tpm[,1]
dissect_tpm <- dissect_tpm[,-1]
dissect_tpm <- dissect_tpm[rowSums(dissect_tpm)>0, ]
dc_name <- as.data.frame(convert[match(rownames(dissect_tpm), convert$V2),4])
dissect_tpm1 <- cbind(dc_name, dissect_tpm)
rownames(dissect_tpm1)= dissect_tpm1[,1]
dissect_tpm1 <- dissect_tpm1[,-1]
colnames(dissect_tpm1) <- c("ac1", "ac2", "dm1", "dm2", "lm1", "lm2",
                               "vm1", "vm2", "vg1", "vg2")
dissect_tpm2 <- as.data.frame(matrix(nrow= nrow(dissect_tpm1), ncol= 6))
rownames(dissect_tpm2) <- rownames(dissect_tpm1)
dissect_tpm2[,1:2] <- dissect_tpm1[,1:2]
dissect_tpm2[,5:6] <- dissect_tpm1[,9:10]
##########merge 3 mesoderm together################
m1 <- c(3,5,7)
dissect_tpm2[,3] <- as.integer(rowMeans(dissect_tpm1[,m1]))
m2 <- c(4,6,8)
dissect_tpm2[,4] <- as.integer(rowMeans(dissect_tpm1[,m2]))
colnames(dissect_tpm2) <- c("ac1", "ac2", "me1", "me2", "vg1", "vg2")
###############log transformation##################
logTPM <- cpm(dissect_tpm2, log = TRUE, prior.count = 1)
############remove any NA##############
logTPM <- logTPM[complete.cases(logTPM),]
############remove any infinity##########
logTPM <- logTPM[!is.infinite(rowSums(logTPM)),]
############rep1 dissection data#########
rep1 <- c(1,3,5)
logTPM_rep1 <- logTPM[, rep1]
############remove row std=0#############
logTPM_rep1 <- logTPM_rep1[!apply(logTPM_rep1 , 1 , function(x) sd(x)==0 ), ]
############rep2 dissection data#########
rep2 <- c(2,4,6)
logTPM_rep2 <- logTPM[, rep2]
############remove row std=0#############
logTPM_rep2 <- logTPM_rep2[!apply(logTPM_rep2 , 1 , function(x) sd(x)==0 ), ]
#########################################
##########up-regulated genes#############
#########################################
uplogtpm1 <- logTPM_rep1[match(upgeneNames, rownames(logTPM_rep1)),]
uplogtpm1 <- uplogtpm1[complete.cases(uplogtpm1),]
uplogtpm2 <- logTPM_rep2[match(upgeneNames, rownames(logTPM_rep2)),]
uplogtpm2 <- uplogtpm2[complete.cases(uplogtpm2),]
############mutate dataset for violin plot##########
uplogtpm1AC <- as.data.frame(uplogtpm1[,1])
colnames(uplogtpm1AC) <- c("logTPM")
uplogtpm1ME <- as.data.frame(uplogtpm1[,2])
colnames(uplogtpm1ME) <- c("logTPM")
uplogtpm1VG <- as.data.frame(uplogtpm1[,3])
colnames(uplogtpm1VG) <- c("logTPM")
library(dplyr)
uplogtpm1AC <- uplogtpm1AC %>% mutate(type= c("AC") )
uplogtpm1ME <- uplogtpm1ME %>% mutate(type= c("MZ") )
uplogtpm1VG <- uplogtpm1VG %>% mutate(type= c("VG") )
uplogtpm1_vlpTable <- rbind(uplogtpm1AC, uplogtpm1ME, uplogtpm1VG)
############violin plot#############
library(ggplot2)
# Basic violin plot
vlp <- ggplot(uplogtpm1_vlpTable, aes(x= type, y= logTPM, fill= type))+ 
  geom_violin(trim = T, scale = "count", color="white")+
  stat_summary(fun=mean, geom="point", size=1, color= "red")+
  scale_fill_manual(values=c("salmon1", "khaki2", "skyblue1"))+
  geom_boxplot(width=0.1, color="grey", alpha=0.5, outlier.shape = NA)+
  theme(aspect.ratio = 1)+
  theme(axis.text = element_text(color="black"), text = element_text(size=9))+
  theme(panel.background = element_rect(size= 0.5, fill = "white", colour = "black"))+
  theme(axis.title.x = element_blank())
vlp
ggsave('C:/Users/latro/Desktop/foxh1 and hdac1/hdac1 rna seq/ectoderm/ectoderm_upGene_violinPlot.tiff', vlp, device = "tiff", dpi = 500)
#############unpaired t test############
x <- uplogtpm1AC$logTPM
y <- uplogtpm1ME$logTPM
z <- uplogtpm1VG$logTPM
t_acvg <- wilcox.test(x, z)
t_acme <- t.test(x, y, alternative = "two.sided", var.equal = FALSE)
t_mevg <- t.test(y, z, alternative = "two.sided", var.equal = FALSE)
t_acvg$p.value
t_acme$p.value
t_mevg$p.value
#########################################
#########down-regulated genes############
#########################################
downlogtpm1 <- logTPM_rep1[match(downgeneNames, rownames(logTPM_rep1)),]
downlogtpm1 <- downlogtpm1[complete.cases(downlogtpm1),]
downlogtpm2 <- logTPM_rep2[match(upgeneNames, rownames(logTPM_rep2)),]
downlogtpm2 <- downlogtpm2[complete.cases(downlogtpm2),]
############mutate dataset for violin plot##########
downlogtpm1AC <- as.data.frame(downlogtpm1[,1])
colnames(downlogtpm1AC) <- c("logTPM")
downlogtpm1ME <- as.data.frame(downlogtpm1[,2])
colnames(downlogtpm1ME) <- c("logTPM")
downlogtpm1VG <- as.data.frame(downlogtpm1[,3])
colnames(downlogtpm1VG) <- c("logTPM")
library(dplyr)
downlogtpm1AC <- downlogtpm1AC %>% mutate(type= c("AC") )
downlogtpm1ME <- downlogtpm1ME %>% mutate(type= c("MZ") )
downlogtpm1VG <- downlogtpm1VG %>% mutate(type= c("VG") )
downlogtpm1_vlpTable <- rbind(downlogtpm1AC, downlogtpm1ME, downlogtpm1VG)
############violin plot#############
library(ggplot2)
# Basic violin plot
vlp <- ggplot(downlogtpm1_vlpTable, aes(x= type, y= logTPM, fill= type))+ 
  geom_violin(trim = T, scale = "count", color="white")+
  stat_summary(fun=mean, geom="point", size=1, color= "red")+
  scale_fill_manual(values=c("salmon1", "khaki2", "skyblue1"))+
  geom_boxplot(width=0.1, color="grey", alpha=0.5, outlier.shape = NA)+
  scale_y_continuous(limits = c(0,12))+
  theme(aspect.ratio = 1)+
  theme(axis.text = element_text(color="black"), text = element_text(size=9))+
  theme(panel.background = element_rect(size= 0.5, fill = "white", colour = "black"))+
  theme(axis.title.x = element_blank())
vlp
ggsave('C:/Users/latro/Desktop/foxh1 and hdac1/hdac1 rna seq/ectoderm/ectoderm_downGene_violinPlot.tiff', vlp, device = "tiff", dpi = 500)
#############t-test--unpaired############
x <- downlogtpm1AC$logTPM
y <- downlogtpm1ME$logTPM
z <- downlogtpm1VG$logTPM
t_acvg <- t.test(x, z, alternative = "two.sided", var.equal = FALSE)
t_acme <- t.test(x, y, alternative = "two.sided", var.equal = FALSE)
t_mevg <- t.test(y, z, alternative = "two.sided", var.equal = FALSE)
t_acvg$p.value
t_acme$p.value
t_mevg$p.value



