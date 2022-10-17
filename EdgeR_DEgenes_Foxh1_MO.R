###############################################################################
######################Differential gene expression#############################
###############################################################################
counts <- read.delim("/Volumes/GoogleDrive/My Drive/foxh1 and hdac1/Foxh1MO ectoderm/raw_counts_ctlMO_vs_foxh1MO.txt")
rownames(counts) <- counts[,1]
counts <- counts[,-1]
colnames(counts) <- c("ctl1", "ctl2", "fox1", "fox2")
library(edgeR)
group <- factor(substring(colnames(counts), 1, 3))
group <- relevel (group, ref= "ctl")
rep <- factor(substring(colnames(counts), 4, 4))
y <- DGEList(counts=counts,group=group)
keep <- filterByExpr(y)
table(keep)
y<- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y$samples
design <- model.matrix(~group)
y <- estimateDisp(y, design)
plotMDS(y)
########notice the batch effect from PCA##############
##########test log2-fold-change for consistency between reps##########
design <- model.matrix(~rep+rep:group)
logFC <- predFC(y,design,prior.count=1,dispersion=0.05)
cor(logFC[,2:4])
#####correlation check is not bad between tsa reps######
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
convert <- read.delim("/Volumes/GoogleDrive/My Drive/foxh1 and hdac1/XENTR_10.0_GeneID_GeneName_converter.txt", header=F)
DEgene1 <- convert[match(DEgene, convert$V2),4]
write.table(as.data.frame(DElist),file="/Volumes/GoogleDrive/My Drive/foxh1 and hdac1/Foxh1MO ectoderm/gene_stats_all.txt", quote = FALSE, sep= "\t")
write.table(as.data.frame(DElist1),file="/Volumes/GoogleDrive/My Drive/foxh1 and hdac1/Foxh1MO ectoderm/DEgene_stats.txt", quote = FALSE, sep= "\t")
write.table(as.data.frame(DEgene1),file="/Volumes/GoogleDrive/My Drive/foxh1 and hdac1/Foxh1MO ectoderm/DEgene_names.txt", quote = FALSE, row.names= FALSE, col.names= FALSE)
downgenes <- DElist1[which(DElist1$logFC<0),]
downgeneNames <- convert[match(rownames(downgenes), convert$V2),4]
write.table(downgeneNames,file="/Volumes/GoogleDrive/My Drive/foxh1 and hdac1/Foxh1MO ectoderm/downDEgene_names.txt", quote = FALSE, row.names= FALSE, col.names= FALSE)
upgenes <- DElist1[which(DElist1$logFC>0),]
upgeneNames <- convert[match(rownames(upgenes), convert$V2),4]
write.table(upgeneNames,file="/Volumes/GoogleDrive/My Drive/foxh1 and hdac1/Foxh1MO ectoderm/upDEgene_names.txt", quote = FALSE, row.names= FALSE, col.names= FALSE)

##################################################################
###################HP#########################
tpm_table <- read.delim("/Volumes/GoogleDrive/My Drive/foxh1 and hdac1/Foxh1MO ectoderm/tpm_ctlMO_vs_foxh1MO.txt")
tpm_degenes <- tpm_table[match(rownames(DElist1), tpm_table$gene_id),]
rownames(tpm_degenes) <- tpm_degenes[,1]
tpm_degenes <- tpm_degenes[, -1]
colnames(tpm_degenes) <- c("ctlMO_1", "ctlMO_2", "foxh1MO_1", "Foxh1MO_2")
library(pheatmap)
library("RColorBrewer")
library(ggplot2)
ramp<-1:3/3 
cols<-c(rgb(ramp,0,0),rgb(0,ramp,0),rgb(0,0,ramp),rgb(ramp,0,ramp)) 
rep1 <- c(1, 3, 5)
rep2 <- c(2, 4, 6)
rep1av <- c(1, 9)
rep2av <- c(2, 10)
av <- c(1,2,9,10)
hp <- pheatmap(         tpm_degenes, 
                        scale = "row",
                        cutree_rows = 2, 
                        cluster_cols = F, 
                        show_rownames = F,
                        treeheight_row = 10,
                        gaps_col = 2,
                        clustering_method = "ward.D2", 
                        col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(255))

ggsave('/Volumes/GoogleDrive/My Drive/foxh1 and hdac1/Foxh1MO ectoderm/tpm_ctlMO_vs_foxh1MO.tiff', hp, device = "tiff", dpi = 200)

##################################################################
###################VP#########################
geneName <- as.data.frame(convert[match(rownames(DElist), convert$V2),4])
geneList <- cbind(geneName, DElist)
b <- c(1,2,5)
geneVPtable <- geneList[,b]
colnames(geneVPtable) <- c("geneName", "logFC", "FDR")
geneOfinterest <- c("ascl1", "chrd.1", "fst", "irx1", "mix1", 
                    "foxh1", "kctd15", "neurog1", "dlx3", "cdx4")
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
  theme(axis.text = element_text(color="black"), text = element_text(size=15))+
  theme(panel.background = element_rect(size= 1, fill = "white", colour = "black"))+
  scale_color_manual(name = "Genes",
                     values = c("up-regulated" = "salmon", 
                                "down-regulated" = "cornflowerblue", 
                                "n.s" = "grey"))+
  geom_text_repel(data = geneOfinterestList, 
                  aes(label = geneName), fontface = "italic", size= 5,
                  color = "black", min.segment.length = 0, box.padding = 0.5,
                  segment.linetype = 6, 
                  max.overlaps = 20)
vp
ggsave('/Volumes/GoogleDrive/My Drive/foxh1 and hdac1/Foxh1MO ectoderm/DE_volcanoPlot.tiff', vp, device = "tiff", dpi = 200)

#################################################################################
###########################Germ layer expression#################################
dissect_tpm <- read.delim("/Volumes/GoogleDrive/My Drive/foxh1 and hdac1/hdac1 rna seq/TPM_table_ira_dissected_st105_tissues.txt", header=F)
rownames(dissect_tpm)= dissect_tpm[,1]
dissect_tpm <- dissect_tpm[,-1]
dissect_tpm <- dissect_tpm[rowSums(dissect_tpm)>0, ]
dc_name <- as.data.frame(convert[match(rownames(dissect_tpm), convert$V2),4])
dissect_tpm1 <- cbind(dc_name, dissect_tpm)
rownames(dissect_tpm1)= dissect_tpm1[,1]
dissect_tpm1 <- dissect_tpm1[,-1]
colnames(dissect_tpm1) <- c("ac1", "ac2", "dm1", "dm2", "lm1", "lm2",
                            "vm1", "vm2", "vg1", "vg2")
dissect_tpm2 <- as.data.frame(matrix(nrow= nrow(dissect_tpm1), ncol= 3))
rownames(dissect_tpm2) <- rownames(dissect_tpm1)
###########average two reps############
dissect_tpm2[,1] <- rowMeans(dissect_tpm1[, 1:2])
##########merge 3 mesoderm together################
dissect_tpm2[,2] <- rowMeans(dissect_tpm1[, 3:8])
dissect_tpm2[,3] <- rowMeans(dissect_tpm1[, 9:10])
###############log transformation##################
logTPM <- cpm(dissect_tpm2, log = TRUE, prior.count = 1)
############remove any NA##############
logTPM <- logTPM[complete.cases(logTPM),]
############remove any infinity##########
logTPM <- logTPM[!is.infinite(rowSums(logTPM)),]
uplogtpm <- data.frame(logTPM[match(upgeneNames, rownames(logTPM)),])
downlogtpm <- data.frame(logTPM[match(downgeneNames, rownames(logTPM)),])
library(dplyr)
######UP######
uplogtpmAC <- as.data.frame(uplogtpm[,1])
colnames(uplogtpmAC) <- c("logTPM")
uplogtpmAC <- uplogtpmAC %>% mutate(type= c("AC") )
uplogtpmME <- as.data.frame(uplogtpm[,2])
colnames(uplogtpmME) <- c("logTPM")
uplogtpmME <- uplogtpmME %>% mutate(type= c("ME") )
uplogtpmVG <- as.data.frame(uplogtpm[,3])
colnames(uplogtpmVG) <- c("logTPM")
uplogtpmVG <- uplogtpmVG %>% mutate(type= c("VG") )
uplogtpm1 <- rbind(uplogtpmAC, uplogtpmME, uplogtpmVG)
uplogtpm1 <- uplogtpm1 %>% mutate(direction= c("up"))
######DOWN######
downlogtpmAC <- as.data.frame(downlogtpm[,1])
colnames(downlogtpmAC) <- c("logTPM")
downlogtpmAC <- downlogtpmAC %>% mutate(type= c("AC") )
downlogtpmME <- as.data.frame(downlogtpm[,2])
colnames(downlogtpmME) <- c("logTPM")
downlogtpmME <- downlogtpmME %>% mutate(type= c("ME") )
downlogtpmVG <- as.data.frame(downlogtpm[,3])
colnames(downlogtpmVG) <- c("logTPM")
downlogtpmVG <- downlogtpmVG %>% mutate(type= c("VG") )
downlogtpm1 <- rbind(downlogtpmAC, downlogtpmME, downlogtpmVG)
downlogtpm1 <- downlogtpm1 %>% mutate(direction= c("down"))

vp_table <- rbind(uplogtpm1, downlogtpm1)
###############plot#####################
library(ggplot2)
library(RColorBrewer)
p <- ggplot(vp_table, aes(x= direction, y= logTPM))+ 
  geom_violin(aes(fill= type), position=position_dodge(.9), trim=F, scale = "area", color= "white")+
  stat_summary(fun=mean, aes(group= type), position=position_dodge(.9), geom="point", size=2, color= "black")+
  scale_fill_manual(values=c("#ECBAC6", "#98E0E2", "#B6E298"))+
  geom_boxplot(width=0.1, aes(fill= type), position=position_dodge(.9), color="grey", alpha=0.5)+
  theme(aspect.ratio = 1/2)+
  theme(axis.text = element_text(color="black"), text = element_text(size=9))+
  theme(panel.background = element_rect(size= 0.5, fill = "white", colour = "black"))+
  theme(axis.title.x = element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p

ggsave('/Volumes/GoogleDrive/My Drive/foxh1 and hdac1/Foxh1MO ectoderm/violinPlot.tiff', p, device = "tiff", dpi = 200)
################statistics################
a <- uplogtpmAC$logTPM
b <- uplogtpmME$logTPM
c <- uplogtpmVG$logTPM

d <- downlogtpmAC$logTPM
e <- downlogtpmME$logTPM
f <- downlogtpmVG$logTPM

t.test(e,f)


################################################################
##################find hdac1 target genes#######################
################################################################
TSA_upGenes <- read.table("/Volumes/GoogleDrive/My Drive/foxh1 and hdac1/hdac1 rna seq/ectoderm/upDEgene_names.txt", quote="\"", comment.char="")
Hdac1_boundGenes <- read.delim("/Volumes/GoogleDrive/My Drive/foxh1 and hdac1/Foxh1MO ectoderm/st105_hdac1_IDR0.05_peaks_associated_genes.txt", header=FALSE)
Hdac1_boundGenes <- data.frame(unique((Hdac1_boundGenes$V14)))
colnames(Hdac1_boundGenes) <- c("V1")
hdac1_supGenes <- Hdac1_boundGenes[match(TSA_upGenes$V1, Hdac1_boundGenes$V1),]
hdac1_supGenes <- data.frame(unique(hdac1_supGenes))
colnames(hdac1_supGenes) <- c("V1")
###########Foxh1 hdac1 overlap up genes#############

foxh1_hdac1_overlap <- hdac1_supGenes[match(upgeneNames, hdac1_supGenes$V1),]
foxh1_hdac1_overlap <- data.frame(unique(foxh1_hdac1_overlap))
colnames(foxh1_hdac1_overlap) <- c("gene")


###########################################################
#####heatmap for spatial expression of overlap genes#######
###########################################################

foxh1_hdac1_gene_TPM <- match(foxh1_hdac1_overlap$gene, rownames(dissect_tpm2))
foxh1_hdac1_gene_TPM <- dissect_tpm2[foxh1_hdac1_gene_TPM, ]
foxh1_hdac1_gene_TPM <- na.omit(foxh1_hdac1_gene_TPM)
colnames(foxh1_hdac1_gene_TPM) <- c("AC", "MZ", "VG")

library(pheatmap)
library("RColorBrewer")
library(ggplot2)
hp <- pheatmap(         foxh1_hdac1_gene_TPM, 
                        scale = "row",
                        cutree_rows = 1, 
                        cluster_cols = F, 
                        show_rownames = T,
                        treeheight_row = 10,
                        clustering_method = "ward.D2", 
                        col=colorRampPalette(brewer.pal(9,"YlOrRd"))(255))

ggsave('/Volumes/GoogleDrive/My Drive/foxh1 and hdac1/Foxh1MO ectoderm/spatial expression of foxh1 hdac1 co-repressed genes.tiff', hp, device = "tiff", dpi = 200)










