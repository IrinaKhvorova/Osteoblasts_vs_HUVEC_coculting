# Transfer and open this file in the directory where you are going to work and where the files with the original raw data and the legend table are located

#### Load the required packages or install if not present ####
if (!require("readxl"))
  install.packages("readxl")
library("readxl")

if (!require("xlsx"))
  install.packages("xlsx")
library("xlsx")

if (!require("RColorBrewer"))
  install.packages("RColorBrewer")
library("RColorBrewer")

if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
  BiocManager::install("limma") 
}
library("limma")

if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
  BiocManager::install("mixOmics")
}
library("mixOmics")

if (!requireNamespace('BiocManager', quietly = TRUE)){
  install.packages('BiocManager')
  BiocManager::install('EnhancedVolcano')
}
library(EnhancedVolcano)

if (!require("ggpubr"))
  install.packages("ggpubr")
library("ggpubr")

if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
  BiocManager::install("rWikiPathways")
}
library("rWikiPathways")

if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
  BiocManager::install("clusterProfiler")
}
library("clusterProfiler")

if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
  BiocManager::install("DESeq2")
}
library("DESeq2")

if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
  BiocManager::install("vsn")
}
library("vsn")

if (!require("stringr"))
  install.packages('stringr')
library(stringr)

if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
  BiocManager::install("apeglm")
}
library(apeglm)

if (!require("ggvenn"))
  install.packages("ggvenn")
library("ggvenn")

BiocManager::install("org.Hs.eg.db")

if (!require("pheatmap"))
  install.packages("pheatmap")
library("pheatmap")

# We start with the endothelial analysis. We will analyze HUVEC proteomic data, HUVEC transcriptomic data, and then move on to the osteoblast analysis

# Transcriptome analysis ####
# read data and legend
tr_legend <- data.frame(read_excel("cocult_trans_legend.xlsx")) 
tr_data <- read.delim("featureCounts_NM+PH.txt")

### Convert ENSEMBL to GENE SYMBOL ####
# remove dots and numbers after them in ID using regular expressions
tr_data$ID <- gsub('\\.\\d*', '', tr_data$ID)
colnames(tr_data)[1] = 'ENSEMBL'
# convert ID to gene symbol
tr_genes <- bitr(tr_data$ENSEMBL, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = 'org.Hs.eg.db')

# not all ids have been mapped and some of them have multiple maps so we solve gaps and multiple matches
tr_genes$ENSEMBL <- make.names(tr_genes$ENSEMBL, unique=TRUE)
# merge two data frames by ID
tr_data <- merge(tr_data, tr_genes, by="ENSEMBL", all.x = TRUE)
# if column with official gene symbols isn't available use ENSEMBL ID
for (i in 1:nrow(tr_data)) {
  if (is.na(tr_data[i, "SYMBOL"])){
    tr_data[i, "SYMBOL"] <- tr_data[i, "ENSEMBL"]
  }
}

# make official gene symbols rownames of data
tr_data$SYMBOL <- make.names(tr_data$SYMBOL, unique=TRUE)
rownames(tr_data) <- tr_data$SYMBOL
str(tr_data)
tr_data <- tr_data[, c(-1, -20)]
# change columns names to tr_legend samples
colnames(tr_data) <- tr_legend[,1]


## HUVEC transcriptome analysis ####
setwd("./HUVEC_transcriptome")
# separate HUVEC's legend and data
fact_huv_tr <- tr_legend[tr_legend$Cell.type == 'HUVEC', ]
huv_tr <- tr_data[, rownames(fact_huv_tr)]
#write.csv(huv_tr, "HUVEC_tr_row_data")

# delete Cell.type and Sample columns
fact_huv_tr <- fact_huv_tr[, c(-1, -2)]
# make factors variables
fact_huv_tr$Condition <- as.factor(fact_huv_tr$Condition)
fact_huv_tr$Donor <- as.factor(fact_huv_tr$Donor)
str(fact_huv_tr)

# Construct DESEQDataSet Object
dds_huv <- DESeqDataSetFromMatrix(countData = huv_tr,
                                  colData = fact_huv_tr,
                                  design= ~ Donor + Condition)
dds_huv

# minimum filtering
keep <- rowSums(counts(dds_huv)) >= 10
dds_huv <- dds_huv[keep,]

# differential expression analysis
dds_huv <- DESeq(dds_huv)
res <- results(dds_huv, alpha=0.05)
resultsNames(dds_huv)

# extracting regularized-logarithm transformed values
rld_huv <- rlog(dds_huv, blind=FALSE)
head(assay(rld_huv), 10)
# dispersion plot
meanSdPlot(assay(rld_huv))
huv_tr <- assay(rld_huv)
#write.csv(huv_tr, 'HUVEC_tr_reg_log_data.csv')


### Сlustering ######
### Sample Principal Component Plot
huv_tr_pca <- pca(t(huv_tr), ncomp = 5, center = TRUE)
#tiff('PCA_HUVEC_tr.tiff', units="in", width=10, height=8, res=300, compression = 'lzw')
plotIndiv(huv_tr_pca, comp = c(1, 2), ind.names = F, group = fact_huv_tr$Condition, legend = TRUE, ellipse = T, title = 'PCA', style = "ggplot2", cex = 6, size.axis = 14, size.xlabel = 18, size.ylabel = 18, size.title = 22, size.legend = 18, size.legend.title = 18)
dev.off()

plotIndiv(huv_tr_pca, comp = c(1, 2), ind.names = T, group = fact_huv_tr$Donor, legend = TRUE, ellipse = T, title = 'PCA')
dev.off()

#PLS-DA 
ordination.optimum.splsda <- splsda(t(huv_tr), fact_huv_tr$Condition, ncomp = 3, keepX = c(15,15,15))
selectVar(ordination.optimum.splsda, comp=1)
selectVar(ordination.optimum.splsda, comp=2)
selectVar(ordination.optimum.splsda, comp=3)


#layout(matrix(c(1, 2), 2, 3))
plotLoadings(ordination.optimum.splsda, comp = 1, size.name = 1, size.title = 1.2, title = "Loadings\n on 1st component", contrib = "max", legend = FALSE, col.ties="black", ndisplay = 15)

plotLoadings(ordination.optimum.splsda, comp = 2, size.name = 1, size.title = 1.2, title = "Loadings\n on 2nd component", contrib = "max",ndisplay = 15,  legend = FALSE, col.ties="black")

#tiff('PLSDA_HUVEC_tr.tiff', units="in", width=10, height=8, res=300, compression = 'lzw')
plotIndiv(ordination.optimum.splsda, ind.names = F, ellipse = T, style = "ggplot2", title = "PLS-DA ordination", legend=TRUE, cex = 6, size.axis = 14, size.xlabel = 18, size.ylabel = 18, size.title = 22, size.legend = 18, size.legend.title = 18)
dev.off()

#tiff('PLSDA_HUVEC_transcr.tiff', units="in", width=14, height=8, res=300, compression = 'lzw')
layout(matrix(c(1, 2, 3, 3, 3, 3), 2, 3))
plotLoadings(ordination.optimum.splsda, comp = 1, size.name = 1.5, size.title = 1.2, 
             title = "Loadings\n on 1st component", contrib = "max", legend = FALSE, 
             col.ties="black", ndisplay = 15)
plotLoadings(ordination.optimum.splsda, comp = 2, size.name = 1.5, size.title = 1.2, 
             title = "Loadings\n on 2nd component", contrib = "max", ndisplay = 15,  
             legend = FALSE, col.ties="black")
plotIndiv(ordination.optimum.splsda, ind.names = F, ellipse = T, style = "graphics",
          title = "PLS-DA ordination",abline = TRUE, cex = 2.5, size.axis = 2.5, 
          size.xlabel = 2.8, size.ylabel = 2.8, size.title = 2.5, legend=TRUE, 
          size.legend = 2, size.legend.title = 2.5)
dev.off()

#tiff('PCA_PLSDA_HUVEC_transcr.tiff', units="in", width=18, height=8, res=300, compression = 'lzw')
layout(matrix(1:2, ncol = 2))
plotIndiv(huv_tr_pca, comp = c(1, 2), ind.names = F, group = fact_huv_tr$Condition, 
          style = "graphics", ellipse = T, title = 'PCA', 
          legend=TRUE, size.title = 2, size.legend = 1.5, size.legend.title = 1.8, 
          abline = TRUE, cex = 2.5, size.axis = 1.5, size.xlabel = 1.8, size.ylabel = 1.8)
plotIndiv(ordination.optimum.splsda, ind.names = F, ellipse = T, style = "graphics",
          title = "PLS-DA ordination", legend=TRUE, size.title = 2,
          size.legend = 1.5, size.legend.title = 1.8, abline = TRUE, cex = 2.5, size.axis = 1.5, size.xlabel = 1.8, size.ylabel = 1.8)
dev.off()


### Differential analysis ####
resultsNames(dds_huv)

#### Condition_direct_vs_сontrol ####
huv_contr_dir_tr <- results(dds_huv, name="Condition_direct_vs_control", alpha=0.05)
#write.csv(huv_contr_dir_tr, file='HUVEC_tr_direct_vs_control.csv')

# shrink log fold changes association with condition:
resLFC_huv_contr_dir_tr <- lfcShrink(dds_huv, coef="Condition_direct_vs_control", 
                               type="apeglm")
# MA-plot
plotMA(huv_contr_dir_tr, ylim=c(-2,2))
plotMA(resLFC_huv_contr_dir_tr, ylim=c(-2,2))
dev.off()

# order the results table by the smallest adjusted p-value:
resLFC_huv_contr_dir_tr
huv_contr_dir_tr <- resLFC_huv_contr_dir_tr[order(resLFC_huv_contr_dir_tr$padj),]
# write the ordered table to a file
#write.csv(huv_contr_dir_tr, file='HUVEC_tr_direct_vs_control_apeglm.csv')

# Vulcano plot
#tiff('HUVEC_tr_direct_vs_control_vulcano.tiff', units="in", width=16, height=12, res=300, compression = 'lzw')
EnhancedVolcano(huv_contr_dir_tr,
                lab = rownames(huv_contr_dir_tr),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1,
                xlim = c(-5.5, 13),
                #xlim = c(-8, 15),
                ylim = c(0, 165),
                title ="Control-direct cocultivation. HUVEC",
                labSize = 7.0,
                pointSize = 3.5,
                axisLabSize = 25,
                titleLabSize = 30,
                subtitleLabSize = 1,
                captionLabSize = 20,
                boxedLabels = F,
                colAlpha = 1)
dev.off()


head(huv_contr_dir_tr)
# counts of reads for a single gene across the Condition
boxplot(as.numeric(huv_tr['POSTN',]) ~ Condition, data = fact_huv_tr, varwidth = TRUE, log = "y", las = 1, ylab = "POSTN")
boxplot(as.numeric(huv_tr['TGFBI',]) ~ Condition, data = fact_huv_tr, varwidth = TRUE, log = "y", las = 1, ylab = "TGFBI")
boxplot(as.numeric(huv_tr['SDC2',]) ~ Condition, data = fact_huv_tr, varwidth = TRUE, log = "y", las = 1, ylab = "SDC2")
boxplot(as.numeric(huv_tr['FOS',]) ~ Condition, data = fact_huv_tr, varwidth = TRUE, log = "y", las = 1, ylab = "FOS")


# control-direct up and down regulated transcripts
# change rows with NA (in p-value)
huv_contr_dir_tr$log2FoldChange[is.na(huv_contr_dir_tr$log2FoldChange)] <- 0
huv_contr_dir_tr$padj[is.na(huv_contr_dir_tr$padj)] <- 1

huv_contr_dir_tr_up <- huv_contr_dir_tr[huv_contr_dir_tr$log2FoldChange >= 1 & huv_contr_dir_tr$padj <= 0.05,]
huv_contr_dir_tr_dn <- huv_contr_dir_tr[huv_contr_dir_tr$log2FoldChange <= -1 & huv_contr_dir_tr$padj <= 0.05,]

#write.csv(huv_contr_dir_tr_up, file='HUVEC_tr_direct_vs_control_up.csv')
#write.csv(huv_contr_dir_tr_dn, file='HUVEC_tr_direct_vs_control_down.csv')


##### Control-Direct enrichment analysis #####

huv_cdir_tr_up <- bitr(rownames(huv_contr_dir_tr_up), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
huv_cdir_tr_dn <- bitr(rownames(huv_contr_dir_tr_dn), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
huv_cdir_tr_bg <- bitr(rownames(huv_contr_dir_tr), fromType = "SYMBOL", toType = "ENTREZID",OrgDb = "org.Hs.eg.db")

###### huv_contr_dir_trans_GO_BP ####
huv_cdir_tr_up_bp <- enrichGO(gene = huv_cdir_tr_up$ENTREZID, universe = huv_cdir_tr_bg$ENTREZID, OrgDb = "org.Hs.eg.db", ont = "BP", pAdjustMethod = "fdr", pvalueCutoff = 0.05, readable = TRUE)
#tiff('HUVEC_tr_direct_vs_control_up_GO_BP.tiff', units="in", width=10, height=8, res=300, compression = 'lzw')
dotplot(huv_cdir_tr_up_bp, showCategory = 20, font.size = 18,  label_format = 50)
dev.off()


huv_cdir_tr_dn_bp <- enrichGO(gene = huv_cdir_tr_dn$ENTREZID, universe = huv_cdir_tr_bg$ENTREZID, OrgDb = "org.Hs.eg.db", ont = "BP", pAdjustMethod = "fdr", pvalueCutoff = 0.05, readable = TRUE)
#tiff('HUVEC_tr_direct_vs_control_down_GO_BP.tiff', units="in", width=11, height=8, res=300, compression = 'lzw')
dotplot(huv_cdir_tr_dn_bp, showCategory = 15, font.size = 18,  label_format = 55)
dev.off()


###### huv_contr_dir_trans_GO_CC ####
huv_cdir_tr_up_cc <- enrichGO(gene = huv_cdir_tr_up$ENTREZID, universe = huv_cdir_tr_bg$ENTREZID, OrgDb = "org.Hs.eg.db", ont = "CC", pAdjustMethod = "fdr", pvalueCutoff = 0.05, readable = TRUE)
#tiff('HUVEC_tr_direct_vs_control_up_GO_CC.tiff', units="in", width=10, height=9, res=300, compression = 'lzw')
dotplot(huv_cdir_tr_up_cc, showCategory = 20, font.size = 18,  label_format = 40)
dev.off()

huv_cdir_tr_dn_cc <- enrichGO(gene = huv_cdir_tr_dn$ENTREZID, universe = huv_cdir_tr_bg$ENTREZID, OrgDb = "org.Hs.eg.db", ont = "CC", pAdjustMethod = "fdr", pvalueCutoff = 0.05, readable = TRUE)
#tiff('HUVEC_tr_direct_vs_control_down_GO_CC.tiff', units="in", width=10, height=9, res=300, compression = 'lzw')
dotplot(huv_cdir_tr_dn_cc, showCategory = 20, font.size = 18,  label_format = 40)
dev.off()

###### huv_contr_dir_trans_KEGG ####
huv_cdir_tr_up_kegg <- enrichKEGG(gene = huv_cdir_tr_up$ENTREZID, universe = huv_cdir_tr_bg$ENTREZID, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05, pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2, use_internal_data = FALSE)
#tiff('HUVEC_tr_direct_vs_control_up_KEGG.tiff', units="in", width=10, height=8, res=300, compression = 'lzw')
dotplot(huv_cdir_tr_up_kegg, showCategory = 20, font.size = 18,  label_format = 50)
dev.off()

huv_cdir_tr_dn_kegg <- enrichKEGG(gene = huv_cdir_tr_dn$ENTREZID, universe = huv_cdir_tr_bg$ENTREZID, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05, pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2, use_internal_data = FALSE)
#tiff('HUVEC_tr_direct_vs_control_down_KEGG.tiff', units="in", width=11, height=4, res=300, compression = 'lzw')
dotplot(huv_cdir_tr_dn_kegg, showCategory = 20, font.size = 18, label_format = 55)
dev.off()


#### Condition_indirect_vs_сontrol ####
huv_contr_indir_tr <- results(dds_huv, name="Condition_indirect_vs_control", alpha=0.05)
#write.csv(huv_contr_indir_tr, file='HUVEC_tr_indirect_vs_control.csv')

# shrink log fold changes association with condition:
resLFC_huv_contr_indir_tr <- lfcShrink(dds_huv, coef="Condition_indirect_vs_control", 
                               type="apeglm")
# MA-plot
plotMA(huv_contr_indir_tr, ylim=c(-2,2))
plotMA(resLFC_huv_contr_indir_tr, ylim=c(-2,2))
dev.off()

# order the results table by the smallest adjusted p-value:
resLFC_huv_contr_indir_tr
huv_contr_indir_tr <- resLFC_huv_contr_indir_tr[order(resLFC_huv_contr_indir_tr$padj),]
# write the ordered table to a file
#write.csv(huv_contr_indir_tr, file='HUVEC_tr_indirect_vs_control_apeglm.csv')


# Vulcano plot
#tiff('HUVEC_tr_indirect_vs_control_vulcano.tiff', units="in", width=16, height=12, res=300, compression = 'lzw')
EnhancedVolcano(huv_contr_indir_tr,
                lab = rownames(huv_contr_indir_tr),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1,
                xlim = c(-10, 7),
                ylim = c(0, 130),
                title ="Control-indirect cocultivation. HUVEC",
                labSize = 7.0,
                pointSize = 3.5,
                axisLabSize = 25,
                titleLabSize = 30,
                subtitleLabSize = 1,
                captionLabSize = 20,
                boxedLabels = F,
                colAlpha = 1)
dev.off()

head(huv_contr_indir_tr)
# counts of reads for a single gene across the Condition
boxplot(as.numeric(huv_tr['FOS',]) ~ Condition, data = fact_huv_tr, varwidth = TRUE, log = "y", las = 1, ylab = "FOS")
boxplot(as.numeric(huv_tr['PCDH12',]) ~ Condition, data = fact_huv_tr, varwidth = TRUE, log = "y", las = 1, ylab = "PCDH12")
boxplot(as.numeric(huv_tr['CXCR4',]) ~ Condition, data = fact_huv_tr, varwidth = TRUE, log = "y", las = 1, ylab = "CXCR4")
boxplot(as.numeric(huv_tr['UNC5B',]) ~ Condition, data = fact_huv_tr, varwidth = TRUE, log = "y", las = 1, ylab = "UNC5B")

# indirect up and down regulated transcripts
# change rows with NA (in p-value)
huv_contr_indir_tr$log2FoldChange[is.na(huv_contr_indir_tr$log2FoldChange)] <- 0
huv_contr_indir_tr$padj[is.na(huv_contr_indir_tr$padj)] <- 1

huv_contr_indir_tr_up <- huv_contr_indir_tr[huv_contr_indir_tr$log2FoldChange >= 1 & huv_contr_indir_tr$padj <= 0.05,]
huv_contr_indir_tr_dn <- huv_contr_indir_tr[huv_contr_indir_tr$log2FoldChange <= -1 & huv_contr_indir_tr$padj <= 0.05,]

#write.csv(huv_contr_indir_tr_up, file='HUVEC_tr_indirect_vs_control_up.csv')
#write.csv(huv_contr_indir_tr_dn, file='HUVEC_tr_indirect_vs_control_down.csv')


##### Control-indirect enrichment analysis #####

huv_cindir_tr_up <- bitr(rownames(huv_contr_indir_tr_up), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
huv_cindir_tr_dn <- bitr(rownames(huv_contr_indir_tr_dn), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
huv_cindir_tr_bg <- bitr(rownames(huv_contr_indir_tr), fromType = "SYMBOL", toType = "ENTREZID",OrgDb = "org.Hs.eg.db")

###### huv_contr_indir_trans_GO_BP ####

huv_cindir_tr_up_bp <- enrichGO(gene = huv_cindir_tr_up$ENTREZID, universe = huv_cindir_tr_bg$ENTREZID, OrgDb = "org.Hs.eg.db", ont = "BP", pAdjustMethod = "fdr", pvalueCutoff = 0.05, readable = TRUE)
#tiff('HUVEC_tr_indirect_vs_control_up_GO_BP.tiff', units="in", width=10, height=8, res=300, compression = 'lzw')
dotplot(huv_cindir_tr_up_bp, showCategory = 20, font.size = 18, label_format = 50)
dev.off()


huv_cindir_tr_dn_bp <- enrichGO(gene = huv_cindir_tr_dn$ENTREZID, universe = huv_cindir_tr_bg$ENTREZID, OrgDb = "org.Hs.eg.db", ont = "BP", pAdjustMethod = "fdr", pvalueCutoff = 0.05, readable = TRUE)
#tiff('HUVEC_tr_indirect_vs_control_down_GO_BP.tiff', units="in", width=10, height=8, res=300, compression = 'lzw')
dotplot(huv_cindir_tr_dn_bp, showCategory = 20, font.size = 18, label_format = 50)
dev.off()


###### huv_contr_indir_trans_GO_CC ####
huv_cindir_tr_up_cc <- enrichGO(gene = huv_cindir_tr_up$ENTREZID, universe = huv_cindir_tr_bg$ENTREZID, OrgDb = "org.Hs.eg.db", ont = "CC", pAdjustMethod = "fdr", pvalueCutoff = 0.05, readable = TRUE)
#tiff('HUVEC_tr_indirect_vs_control_up_GO_CC.tiff', units="in", width=10, height=4, res=300, compression = 'lzw')
dotplot(huv_cindir_tr_up_cc, showCategory = 20, font.size = 18, label_format = 40)
dev.off()

huv_cindir_tr_dn_cc <- enrichGO(gene = huv_cindir_tr_dn$ENTREZID, universe = huv_cindir_tr_bg$ENTREZID, OrgDb = "org.Hs.eg.db", ont = "CC", pAdjustMethod = "fdr", pvalueCutoff = 0.05, readable = TRUE)
#tiff('HUVEC_tr_indirect_vs_control_down_GO_CC.tiff', units="in", width=10, height=3, res=300, compression = 'lzw')
dotplot(huv_cindir_tr_dn_cc, showCategory = 20, font.size = 18, label_format = 40)
dev.off()

###### huv_contr_indir_trans_KEGG ####
huv_cindir_tr_up_kegg <- enrichKEGG(gene = huv_cindir_tr_up$ENTREZID, universe = huv_cindir_tr_bg$ENTREZID, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05, pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2, use_internal_data = FALSE)
#tiff('HUVEC_tr_indirect_vs_control_up_KEGG.tiff', units="in", width=10, height=6, res=300, compression = 'lzw')
dotplot(huv_cindir_tr_up_kegg, showCategory = 20, font.size = 18, label_format = 50)
dev.off()


huv_cindir_tr_dn_kegg <- enrichKEGG(gene = huv_cindir_tr_dn$ENTREZID, universe = huv_cindir_tr_bg$ENTREZID, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05, pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2, use_internal_data = FALSE)
#tiff('HUVEC_tr_indirect_vs_control_down_KEGG.tiff', units="in", width=10, height=4.5, res=300, compression = 'lzw')
dotplot(huv_cindir_tr_dn_kegg, showCategory = 20, font.size = 18, label_format = 55)
dev.off()



#### Condition_indirect_vs_direct ####

# change factors levels to change the baseline for direct-indirect comparison
fact_huv_tr$Condition <- factor(fact_huv_tr$Condition, levels = c("indirect", "direct", "control"))
huv_tr <- tr_data[, rownames(fact_huv_tr)]

# Construct DESEQDataSet Object
dds_huv <- DESeqDataSetFromMatrix(countData = huv_tr,
                                  colData = fact_huv_tr,
                                  design= ~ Donor + Condition)
dds_huv

# minimum filtering
keep <- rowSums(counts(dds_huv)) >= 10
dds_huv <- dds_huv[keep,]

# differential expression analysis
dds_huv <- DESeq(dds_huv)
res <- results(dds_huv, alpha=0.05)
resultsNames(dds_huv)

# extracting regularized-logarithm transformed values
rld_huv <- rlog(dds_huv, blind=FALSE)
head(assay(rld_huv), 30)
huv_tr <- assay(rld_huv)

huv_dir_indir_tr <- results(dds_huv, name="Condition_direct_vs_indirect", alpha=0.05)
#write.csv(huv_dir_indir_tr, file='HUVEC_tr_direct_vs_indirect.csv')

# shrink log fold changes association with condition:
resLFC_huv_dir_indir_tr <- lfcShrink(dds_huv, coef="Condition_direct_vs_indirect", 
                               type="apeglm")
# MA-plot
plotMA(huv_dir_indir_tr, ylim=c(-2,2))
plotMA(resLFC_huv_dir_indir_tr, ylim=c(-2,2))
dev.off()

# order the results table by the smallest adjusted p-value:
resLFC_huv_dir_indir_tr
huv_dir_indir_tr <- resLFC_huv_dir_indir_tr[order(resLFC_huv_dir_indir_tr$padj),]

# write the ordered table to a file
#write.csv(huv_dir_indir_tr, file='HUVEC_tr_direct_vs_indirect_apeglm.csv')


# Vulcano plot
#tiff('HUVEC_tr_direct_vs_indirect_vulcano.tiff', units="in", width=16, height=12, res=300, compression = 'lzw')
EnhancedVolcano(huv_dir_indir_tr,
                lab = rownames(huv_dir_indir_tr),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1,
                xlim = c(-6, 13),
                ylim = c(0, 190),
                title ="Indirect-direct cocultivation. HUVEC",
                labSize = 7.0,
                pointSize = 3.5,
                axisLabSize = 25,
                titleLabSize = 30,
                subtitleLabSize = 1,
                captionLabSize = 20,
                boxedLabels = F,
                colAlpha = 1)
dev.off()

head(huv_dir_indir_tr)
# counts of reads for a single gene across the Condition
boxplot(as.numeric(huv_tr['TGFBI',]) ~ Condition, data = fact_huv_tr, varwidth = TRUE, log = "y", las = 1, ylab = "TGFBI")
boxplot(as.numeric(huv_tr['POSTN',]) ~ Condition, data = fact_huv_tr, varwidth = TRUE, log = "y", las = 1, ylab = "POSTN")
boxplot(as.numeric(huv_tr['IGFBP3',]) ~ Condition, data = fact_huv_tr, varwidth = TRUE, log = "y", las = 1, ylab = "IGFBP3")
boxplot(as.numeric(huv_tr['CHPF',]) ~ Condition, data = fact_huv_tr, varwidth = TRUE, log = "y", las = 1, ylab = "CHPF")


# indirect up and down regulated transcripts
# change rows with NA (in p-value)
huv_dir_indir_tr$log2FoldChange[is.na(huv_dir_indir_tr$log2FoldChange)] <- 0
huv_dir_indir_tr$padj[is.na(huv_dir_indir_tr$padj)] <- 1

huv_dir_indir_tr_up <- huv_dir_indir_tr[huv_dir_indir_tr$log2FoldChange >= 1 & huv_dir_indir_tr$padj <= 0.05,]
huv_dir_indir_tr_dn <- huv_dir_indir_tr[huv_dir_indir_tr$log2FoldChange <= -1 & huv_dir_indir_tr$padj <= 0.05,]

#write.csv(huv_dir_indir_tr_up, file='HUVEC_tr_direct_vs_indirect_up.csv')
#write.csv(huv_dir_indir_tr_dn, file='HUVEC_tr_direct_vs_indirect_down.csv')


##### Indirect-direct enrichment analysis #####

huv_dindir_tr_up <- bitr(rownames(huv_dir_indir_tr_up), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
huv_dindir_tr_dn <- bitr(rownames(huv_dir_indir_tr_dn), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
huv_dindir_tr_bg <- bitr(rownames(huv_dir_indir_tr), fromType = "SYMBOL", toType = "ENTREZID",OrgDb = "org.Hs.eg.db")

###### huv_dir_indir_trans_GO_BP ####
huv_dindir_tr_up_bp <- enrichGO(gene = huv_dindir_tr_up$ENTREZID, universe = huv_dindir_tr_bg$ENTREZID, OrgDb = "org.Hs.eg.db", ont = "BP", pAdjustMethod = "fdr", pvalueCutoff = 0.05, readable = TRUE)
#tiff('HUVEC_tr_direct_vs_indirect_up_GO_BP.tiff', units="in", width=10, height=8, res=300, compression = 'lzw')
dotplot(huv_dindir_tr_up_bp, showCategory = 20, font.size = 18, label_format = 55)
dev.off()


huv_dindir_tr_dn_bp <- enrichGO(gene = huv_dindir_tr_dn$ENTREZID, universe = huv_dindir_tr_bg$ENTREZID, OrgDb = "org.Hs.eg.db", ont = "BP", pAdjustMethod = "fdr", pvalueCutoff = 0.05, readable = TRUE)
#tiff('HUVEC_tr_direct_vs_indirect_down_GO_BP.tiff', units="in", width=10, height=8, res=300, compression = 'lzw')
dotplot(huv_dindir_tr_dn_bp, showCategory = 20, font.size = 18, label_format = 50)
dev.off()


###### huv_dir_indir_trans_GO_CC ####
huv_dindir_tr_up_cc <- enrichGO(gene = huv_dindir_tr_up$ENTREZID, universe = huv_dindir_tr_bg$ENTREZID, OrgDb = "org.Hs.eg.db", ont = "CC", pAdjustMethod = "fdr", pvalueCutoff = 0.05, readable = TRUE)
#tiff('HUVEC_tr_direct_vs_indirect_up_GO_CC.tiff', units="in", width=10, height=9, res=300, compression = 'lzw')
dotplot(huv_dindir_tr_up_cc, showCategory = 20, font.size = 18, label_format = 40)
dev.off()

huv_dindir_tr_dn_cc <- enrichGO(gene = huv_dindir_tr_dn$ENTREZID, universe = huv_dindir_tr_bg$ENTREZID, OrgDb = "org.Hs.eg.db", ont = "CC", pAdjustMethod = "fdr", pvalueCutoff = 0.05, readable = TRUE)
# no result

###### huv_dir_indir_trans_KEGG ####
huv_dindir_tr_up_kegg <- enrichKEGG(gene = huv_dindir_tr_up$ENTREZID, universe = huv_dindir_tr_bg$ENTREZID, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05, pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2, use_internal_data = FALSE)
#tiff('HUVEC_tr_direct_vs_indirect_up_KEGG.tiff', units="in", width=10, height=8, res=300, compression = 'lzw')
dotplot(huv_dindir_tr_up_kegg, showCategory = 20, font.size = 18, label_format = 50)
dev.off()


huv_dindir_tr_dn_kegg <- enrichKEGG(gene = huv_dindir_tr_dn$ENTREZID, universe = huv_dindir_tr_bg$ENTREZID, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05, pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2, use_internal_data = FALSE)
#tiff('HUVEC_tr_direct_vs_indirect_down_KEGG.tiff', units="in", width=10, height=8, res=300, compression = 'lzw')
dotplot(huv_dindir_tr_dn_kegg, showCategory = 20, font.size = 18, label_format = 50)
dev.off()


## HUVEC proteome analysis ####
setwd("../")
# read data
huv_pr <- data.frame(read.table("huv_log_imp_data.csv", sep = ',', header = TRUE))
# data is filtered, NA inserted with log_imp in NAguide
colSums(is.na(huv_pr))

# convert ID to gene symbol
colnames(huv_pr)[1] = 'UNIPROT'
huv_pr_genes <- bitr(huv_pr$UNIPROT, fromType = "UNIPROT", toType = "SYMBOL",
                     OrgDb = 'org.Hs.eg.db')
# not all ids have been mapped and some of them have multiple mapps so we remove gaps and multiple matches
huv_pr_genes$UNIPROT <- make.names(huv_pr_genes$UNIPROT, unique=TRUE)
# merge two data frames by ID
huv_pr <- merge(huv_pr, huv_pr_genes, by="UNIPROT", all.x = TRUE)

# if column with official gene symbols isn't available use UNIPROT ID
for (i in 1:nrow(huv_pr)) {
  if (is.na(huv_pr[i, "SYMBOL"])){
    huv_pr[i, "SYMBOL"] <- huv_pr[i, "UNIPROT"]
  }
}
# make official gene symbols rownames of data
huv_pr$SYMBOL <- make.names(huv_pr$SYMBOL, unique=TRUE)
rownames(huv_pr) <- huv_pr$SYMBOL
str(huv_pr)
huv_pr <- huv_pr[, c(-1, -14)]

# read legend
fact_huv_pr <- data.frame(read_excel("cocult_prot_legend.xlsx"))
# separate HUVEC's legend
fact_huv_pr <- fact_huv_pr[fact_huv_pr$Cell.type == 'HUVEC', ]
# change row names
rownames(fact_huv_pr) <- fact_huv_pr[,1]
# delete Cell.type and Sample columns
fact_huv_pr <- fact_huv_pr[, c(-1, -2)]
# make variables factors
fact_huv_pr$Condition <- as.factor(fact_huv_pr$Condition)
fact_huv_pr$Donor <- as.factor(as.character(fact_huv_pr$Donor))
str(fact_huv_pr)

setwd("./HUVEC_proteome")

#Raw huv_pr
pal <- brewer.pal(n = 3, name = "Set1")
cols <- pal[fact_huv_pr$Condition]
boxplot(huv_pr, outline = FALSE, col = cols, main = "Raw data")
legend("topright", inset=.02, levels(fact_huv_pr$Condition), fill = pal, xpd = T, cex=0.8)
# write row data
#write.csv(huv_pr, 'huv_pr.csv')

#VSN normalization
huv_pr_norm <- normalizeQuantiles(huv_pr)
#write.csv(huv_pr_norm, 'huv_pr_norm.csv')
boxplot(huv_pr_norm, outline = FALSE, col = cols, main = "Normalized data")
legend("topright", levels(fact_huv_pr$Condition), fill = pal, xpd = T, cex=0.8)

#### Сlustering ####
#Sample Principal Component Plot
huv_pr_pca <- pca(t(huv_pr_norm), ncomp = 5, center = TRUE)
#tiff('PCA_HUVEC_prot.tiff', units="in", width=10, height=8, res=300, compression = 'lzw')
plotIndiv(huv_pr_pca, comp = c(1, 2), ind.names = F, group = fact_huv_pr$Condition, 
          legend = TRUE, ellipse = T, title = 'PCA', style = "ggplot2", cex = 6, 
          size.axis = 14, size.xlabel = 18, size.ylabel = 18, size.title = 22, 
          size.legend = 18, size.legend.title = 18)
dev.off()
plotIndiv(huv_pr_pca, comp = c(1, 2), ind.names = T, group = fact_huv_pr$Donor, 
          legend = TRUE, ellipse = T, title = 'PCA')

#PLS-DA 
ordination.optimum.splsda <- splsda(t(huv_pr_norm), fact_huv_pr$Condition, ncomp = 3, keepX = c(15,15,15))
selectVar(ordination.optimum.splsda, comp=1)
selectVar(ordination.optimum.splsda, comp=2)
selectVar(ordination.optimum.splsda, comp=3)

#layout(matrix(c(1, 2), 2, 3))
plotLoadings(ordination.optimum.splsda, comp = 1, size.name = 1, size.title = 1.2, title = "Loadings\n on 1st component", contrib = "max", legend = FALSE, col.ties="black", ndisplay = 15)

plotLoadings(ordination.optimum.splsda, comp = 2, size.name = 1, size.title = 1.2, title = "Loadings\n on 2nd component", contrib = "max",ndisplay = 15,  legend = FALSE, col.ties="black")

#tiff('PLSDA_HUVEC_prot.tiff', units="in", width=10, height=8, res=300, compression = 'lzw')
plotIndiv(ordination.optimum.splsda, ind.names = F, ellipse = T, style = "ggplot2", title = "PLS-DA ordination", legend=TRUE, cex = 6, size.axis = 14, size.xlabel = 18, size.ylabel = 18, size.title = 22, size.legend = 18, size.legend.title = 18)
dev.off()

#tiff('PLSDA_HUVEC_proteom.tiff', units="in", width=16, height=8, res=300, compression = 'lzw')
layout(matrix(c(1, 2, 3, 3, 3, 3), 2, 3))
plotLoadings(ordination.optimum.splsda, comp = 1, size.name = 1, size.title = 1.2, title = "Loadings\n on 1st component", contrib = "max", legend = FALSE, col.ties="black", ndisplay = 15)
plotLoadings(ordination.optimum.splsda, comp = 2, size.name = 1, size.title = 1.2, title = "Loadings\n on 2nd component", contrib = "max",ndisplay = 15,  legend = FALSE, col.ties="black")
plotIndiv(ordination.optimum.splsda, ind.names = F, ellipse = T, style = "graphics", 
          abline = TRUE, cex = 2.5, size.axis = 2.5,           
          size.xlabel = 2.8, size.ylabel = 2.8, size.title = 2.5, legend=TRUE,
          size.legend = 2, size.legend.title = 2.5)
dev.off()

#tiff('PCA_PLSDA_HUVEC_prot.tiff', units="in", width=16, height=8, res=300, compression = 'lzw')
layout(matrix(1:2, ncol = 2))
plotIndiv(huv_pr_pca, comp = c(1, 2), ind.names = F, group = fact_huv_pr$Condition, 
          style = "graphics", ellipse = T, title = 'PCA', legend=TRUE, 
          size.title = 2, size.legend = 1.5, size.legend.title = 1.8, abline = TRUE, 
          cex = 2.5, size.axis = 1.5, size.xlabel = 1.8, size.ylabel = 1.8)
plotIndiv(ordination.optimum.splsda, ind.names = F, ellipse = T, style = "graphics", 
          title = "PLS-DA ordination", legend=TRUE, size.title = 2, size.legend = 1.5, 
          size.legend.title = 1.8, abline = TRUE, cex = 2.5, size.axis = 1.5, 
          size.xlabel = 1.8, size.ylabel = 1.8)
dev.off()


#### Differential expression analysis ####
design <- model.matrix(~ 0 + fact_huv_pr$Condition)
colnames(design) <- c("control","direct","indirect")
fit <- lmFit(huv_pr_norm, design)
contrast.matrix <- makeContrasts(direct-control, indirect-control, direct-indirect,
                                 levels=design)
huv_pr_fit <- contrasts.fit(fit, contrast.matrix)
huv_pr_fit <- eBayes(huv_pr_fit)
results <- decideTests(huv_pr_fit)

# Dif_expr_table
huv_pr_full_list <- topTable(huv_pr_fit, number=nrow(huv_pr_norm))
head(huv_pr_full_list, 10)

# counts of proteins for a single one across the groups
boxplot(as.numeric(huv_pr_norm['HTRA1',]) ~ Condition, data = fact_huv_pr, varwidth = TRUE, log = "y", las = 1, ylab = "HTRA1")
boxplot(as.numeric(huv_pr_norm["TAGLN",]) ~ Condition, data = fact_huv_pr, varwidth = TRUE, log = "y", las = 1, ylab = "TAGLN")
boxplot(as.numeric(huv_pr_norm["COL6A1",]) ~ Condition, data = fact_huv_pr, varwidth = TRUE, log = "y", las = 1, ylab = "COL6A1")
boxplot(as.numeric(huv_pr_norm["MMP1",]) ~ Condition, data = fact_huv_pr, varwidth = TRUE, log = "y", las = 1, ylab = "MMP1")


#### Control-Direct analysis ####
huv_contr_dir_pr <- topTable(huv_pr_fit, coef=1, adjust="BH", number=nrow(huv_pr_norm))
#write.csv(as.data.frame(huv_contr_dir_pr), file="HUVEC_pr_direct_vs_control.csv")

head(huv_contr_dir_pr)
#Vulcano
#tiff('HUVEC_pr_direct_vs_control_volcano.tiff', units="in", width=12, height=12, res=300, compression = 'lzw')
EnhancedVolcano(huv_contr_dir_pr,
                lab = rownames(huv_contr_dir_pr),
                x = 'logFC',
                y = 'adj.P.Val',
                pCutoff = 0.05,
                FCcutoff = 1,
                xlim = c(-3.2, 4.5),
                ylim = c(0, 7.5),
                title ="Control-direct cocultivation. HUVEC",
                labSize = 7.0,
                pointSize = 3.5,
                axisLabSize = 25,
                titleLabSize = 30,
                subtitleLabSize = 1,
                captionLabSize = 20,
                boxedLabels = F,
                colAlpha = 1)
dev.off()

huv_contr_dir_pr_up <- huv_contr_dir_pr[huv_contr_dir_pr$logFC >= 1 & huv_contr_dir_pr$adj.P.Val <= 0.05,]
huv_contr_dir_pr_dn <- huv_contr_dir_pr[huv_contr_dir_pr$logFC <= -1 & huv_contr_dir_pr$adj.P.Val <= 0.05,]

#write.csv(as.data.frame(huv_contr_dir_pr_up), file="HUVEC_pr_direct_vs_control_up.csv")
#write.csv(as.data.frame(huv_contr_dir_pr_dn), file="HUVEC_pr_direct_vs_control_down.csv")


##### Control-Direct enrichment analysis ####
huv_cd_pr_up <- bitr(rownames(huv_contr_dir_pr_up), fromType = "SYMBOL", 
                    toType = "ENTREZID", OrgDb = 'org.Hs.eg.db')
huv_cd_pr_dn <- bitr(rownames(huv_contr_dir_pr_dn), fromType = "SYMBOL", 
                    toType = "ENTREZID", OrgDb = 'org.Hs.eg.db')
huv_cd_pr_bg <- bitr(rownames(huv_contr_dir_pr), fromType = "SYMBOL", 
                    toType = "ENTREZID", OrgDb = 'org.Hs.eg.db')

###### huv_contr_dir_prot_GO_BP ####
huv_contr_dir_pr_up_bp <- enrichGO(gene = huv_cd_pr_up$ENTREZID, universe = huv_cd_pr_bg$ENTREZID, OrgDb = 'org.Hs.eg.db', ont = "BP", pAdjustMethod = "fdr",
                            pvalueCutoff = 0.05, readable = TRUE)
#tiff('HUVEC_pr_direct_vs_control_up_GO_BP.tiff', units="in", width=10, height=8, res=300, compression = 'lzw')
dotplot(huv_contr_dir_pr_up_bp, showCategory = 20, font.size = 18, label_format = 50)
dev.off()


huv_contr_dir_pr_dn_bp <- enrichGO(gene = huv_cd_pr_dn$ENTREZID, universe = huv_cd_pr_bg$ENTREZID,
                            OrgDb = 'org.Hs.eg.db', ont = "BP", pAdjustMethod = "fdr",
                            pvalueCutoff = 0.05, readable = TRUE)
# no result

###### huv_contr_dir_prot_GO_CC ####
huv_contr_dir_pr_up_cc <- enrichGO(gene = huv_cd_pr_up$ENTREZID, 
                            universe = huv_cd_pr_bg$ENTREZID,
                            OrgDb = 'org.Hs.eg.db', ont = "CC", pAdjustMethod = "fdr",
                            pvalueCutoff = 0.05, readable = TRUE)
#tiff('HUVEC_pr_direct_vs_control_up_GO_CC.tiff', units="in", width=10, height=8, res=300, compression = 'lzw')
dotplot(huv_contr_dir_pr_up_cc, showCategory = 20, font.size = 18, label_format = 40)
dev.off()

huv_contr_dir_pr_dn_cc <- enrichGO(gene = huv_cd_pr_dn$ENTREZID, universe = huv_cd_pr_bg$ENTREZID,
                            OrgDb = 'org.Hs.eg.db', ont = "CC", pAdjustMethod = "fdr",
                            pvalueCutoff = 0.05, readable = TRUE)
#tiff('HUVEC_pr_direct_vs_control_down_GO_CC.tiff', units="in", width=10, height=5, res=300, compression = 'lzw')
dotplot(huv_contr_dir_pr_dn_cc, showCategory = 20, font.size = 18, label_format = 40)
dev.off()


###### huv_contr_dir_prot_KEGG ####
huv_contr_dir_pr_up_kegg <- enrichKEGG(gene = huv_cd_pr_up$ENTREZID, universe = huv_cd_pr_bg$ENTREZID, organism = "hsa",
                                keyType = "kegg", pvalueCutoff = 0.05, pAdjustMethod = "BH", 
                                minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2,
                                use_internal_data = FALSE)
#tiff('HUVEC_pr_direct_vs_control_up_KEGG.tiff', units="in", width=10, height=4, res=300, compression = 'lzw')
dotplot(huv_contr_dir_pr_up_kegg, showCategory = 20, font.size = 18, label_format = 50)
dev.off()


huv_contr_dir_pr_dn_kegg <- enrichKEGG(gene = huv_cd_pr_dn$ENTREZID, universe = huv_cd_pr_bg$ENTREZID, organism = "hsa",
                             keyType = "kegg", pvalueCutoff = 0.05, pAdjustMethod = "BH", 
                             minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2,
                             use_internal_data = FALSE)
#tiff('HUVEC_pr_direct_vs_control_down_KEGG.tiff', units="in", width=10, height=3, res=300, compression = 'lzw')
dotplot(huv_contr_dir_pr_dn_kegg, showCategory = 20, font.size = 18, label_format = 50)
dev.off()


#### Control-Indirect analysis ####
huv_contr_indir_pr <- topTable(huv_pr_fit, coef=2, adjust="BH", number=nrow(huv_pr_norm))
#write.csv(as.data.frame(huv_contr_indir_pr), file="HUVEC_pr_indirect_vs_control.csv")

#Vulcano
#tiff('HUVEC_pr_indirect_vs_control_volcano.tiff', units="in", width=12, height=12, res=300, compression = 'lzw')
EnhancedVolcano(huv_contr_indir_pr,
                lab = rownames(huv_contr_indir_pr),
                x = 'logFC',
                y = 'adj.P.Val',
                pCutoff = 0.05,
                FCcutoff = 1,
                xlim = c(-3, 3),
                ylim = c(0, 7.5),
                title ="Control-indirect cocultivation. HUVEC",
                labSize = 7.0,
                pointSize = 3.5,
                axisLabSize = 25,
                titleLabSize = 30,
                subtitleLabSize = 1,
                captionLabSize = 20,                
                boxedLabels = F,                
                colAlpha = 1)
dev.off()

head(huv_contr_indir_pr)

huv_contr_indir_pr_up <- huv_contr_indir_pr[huv_contr_indir_pr$logFC >= 1 & huv_contr_indir_pr$adj.P.Val <= 0.05,]
huv_contr_indir_pr_dn <- huv_contr_indir_pr[huv_contr_indir_pr$logFC <= -1 & huv_contr_indir_pr$adj.P.Val <= 0.05,]

#write.csv(as.data.frame(huv_contr_indir_pr_up), file="HUVEC_pr_indirect_vs_control_up.csv")
#write.csv(as.data.frame(huv_contr_indir_pr_dn), file="HUVEC_pr_indirect_vs_control_down.csv")


##### Control-Indirect enrichment analysis ####
huv_cind_pr_up <- bitr(rownames(huv_contr_indir_pr_up), fromType = "SYMBOL", toType = "ENTREZID",
              OrgDb = 'org.Hs.eg.db')
huv_cind_pr_dn <- bitr(rownames(huv_contr_indir_pr_dn), fromType = "SYMBOL", toType = "ENTREZID",
              OrgDb = 'org.Hs.eg.db')
huv_cind_pr_bg <- bitr(rownames(huv_contr_indir_pr), fromType = "SYMBOL", toType = "ENTREZID",
              OrgDb = 'org.Hs.eg.db')


###### huv_contr_indir_prot_GO_BP ####
huv_contr_indir_pr_up_bp <- enrichGO(gene = huv_cind_pr_up$ENTREZID, universe = huv_cind_pr_bg$ENTREZID,
                            OrgDb = 'org.Hs.eg.db', ont = "BP", pAdjustMethod = "fdr",
                            pvalueCutoff = 0.05, readable = TRUE)
#tiff('HUVEC_pr_indirect_vs_control_up_GO_BP.tiff', units="in", width=10, height=8, res=300, compression = 'lzw')
dotplot(huv_contr_indir_pr_up_bp, showCategory = 20, font.size = 18, label_format = 50)
dev.off()


huv_contr_indir_pr_dn_bp <- enrichGO(gene = huv_cind_pr_dn$ENTREZID, universe = huv_cind_pr_bg$ENTREZID,
                            OrgDb = 'org.Hs.eg.db', ont = "BP", pAdjustMethod = "fdr",
                            pvalueCutoff = 0.05, readable = TRUE)
#tiff('HUVEC_pr_indirect_vs_control_down_GO_BP.tiff', units="in", width=10, height=8, res=300, compression = 'lzw')
dotplot(huv_contr_indir_pr_dn_bp, showCategory = 20, font.size = 18, label_format = 55)
dev.off()


###### huv_contr_indir_prot_GO_CC ####
huv_contr_indir_pr_up_cc <- enrichGO(gene = huv_cind_pr_up$ENTREZID, universe = huv_cind_pr_bg$ENTREZID,
                            OrgDb = 'org.Hs.eg.db', ont = "CC", pAdjustMethod = "fdr",
                            pvalueCutoff = 0.05, readable = TRUE)
#tiff('HUVEC_pr_indirect_vs_control_up_GO_CC.tiff', units="in", width=10, height=8, res=300, compression = 'lzw')
dotplot(huv_contr_indir_pr_up_cc, showCategory = 20, font.size = 18, label_format = 50)
dev.off()

huv_contr_indir_pr_dn_cc <- enrichGO(gene = huv_cind_pr_dn$ENTREZID, universe = huv_cind_pr_bg$ENTREZID,
                            OrgDb = 'org.Hs.eg.db', ont = "CC", pAdjustMethod = "fdr",
                            pvalueCutoff = 0.05, readable = TRUE)
#tiff('HUVEC_pr_indirect_vs_control_down_GO_CC.tiff', units="in", width=10, height=5, res=300, compression = 'lzw')
dotplot(huv_contr_indir_pr_dn_cc, showCategory = 15, font.size = 18, label_format = 50)
dev.off()


###### huv_contr_indir_prot_KEGG ####
huv_contr_indir_pr_up_kegg <- enrichKEGG(gene = huv_cind_pr_up$ENTREZID, universe = huv_cind_pr_bg$ENTREZID, organism = "hsa",
                                keyType = "kegg", pvalueCutoff = 0.05, pAdjustMethod = "BH", 
                                minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2,
                                use_internal_data = FALSE)
#tiff('HUVEC_pr_indirect_vs_control_up_KEGG.tiff', units="in", width=10, height=4, res=300, compression = 'lzw')
dotplot(huv_contr_indir_pr_up_kegg, showCategory = 20, font.size = 18, label_format = 50)
dev.off()


huv_contr_indir_pr_dn_kegg <- enrichKEGG(gene = huv_cind_pr_dn$ENTREZID, universe = huv_cind_pr_bg$ENTREZID, organism = "hsa",
                                keyType = "kegg", pvalueCutoff = 0.05, pAdjustMethod = "BH",
                                minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2,
                                use_internal_data = FALSE)
#tiff('HUVEC_pr_indirect_vs_control_down_KEGG.tiff', units="in", width=10, height=4, res=300, compression = 'lzw')
dotplot(huv_contr_indir_pr_dn_kegg, showCategory = 20, font.size = 18, label_format = 50)
dev.off()


#### Direct-Indirect analysis ####
huv_dir_indir_pr <- topTable(huv_pr_fit, coef=3, adjust="BH", number=nrow(huv_pr_norm))
#write.csv(as.data.frame(huv_dir_indir_pr), file="HUVEC_pr_direct_vs_indirect.csv")

boxplot(as.numeric(huv_pr_norm["TAGLN",]) ~ Condition, data = fact_huv_pr, varwidth = TRUE, log = "y", las = 1, ylab = "TAGLN")

head(huv_dir_indir_pr)
#Vulcano
#tiff('HUVEC_pr_direct_vs_indirect_volcano.tiff', units="in", width=12, height=12, res=300, compression = 'lzw')
EnhancedVolcano(huv_dir_indir_pr,
                lab = rownames(huv_dir_indir_pr),
                x = 'logFC',
                y = 'adj.P.Val',
                pCutoff = 0.05,
                FCcutoff = 1,
                xlim = c(-1.6, 3),
                ylim = c(0, 6.5),
                title ="Direct-indirect cocultivation. HUVEC",
                labSize = 7.0,
                pointSize = 3.5,
                axisLabSize = 25,
                titleLabSize = 30,
                subtitleLabSize = 1,
                captionLabSize = 20,               
                boxedLabels = F,                
                colAlpha = 1)
dev.off()

huv_dir_indir_pr_up <- huv_dir_indir_pr[huv_dir_indir_pr$logFC >= 1 & huv_dir_indir_pr$adj.P.Val <= 0.05,]
huv_dir_indir_pr_dn <- huv_dir_indir_pr[huv_dir_indir_pr$logFC <= -1 & huv_dir_indir_pr$adj.P.Val <= 0.05,]

#write.csv(as.data.frame(huv_dir_indir_pr_up), file="HUVEC_pr_direct_vs_indirect_up.csv")
#write.csv(as.data.frame(huv_dir_indir_pr_dn), file="HUVEC_pr_direct_vs_indirect_down.csv")


##### Direct-Indirect enrichment analysis ####
huv_dind_pr_up <- bitr(rownames(huv_dir_indir_pr_up), fromType = "SYMBOL", toType = "ENTREZID",
              OrgDb = 'org.Hs.eg.db')
huv_dind_pr_dn <- bitr(rownames(huv_dir_indir_pr_dn), fromType = "SYMBOL", toType = "ENTREZID",
              OrgDb = 'org.Hs.eg.db')
huv_dind_pr_bg <- bitr(rownames(huv_dir_indir_pr), fromType = "SYMBOL", toType = "ENTREZID",
              OrgDb = 'org.Hs.eg.db')


###### huv_dir_indir_prot_GO_BP ####
huv_dir_indir_pr_up_bp <- enrichGO(gene = huv_dind_pr_up$ENTREZID, universe = huv_dind_pr_bg$ENTREZID,
                              OrgDb = 'org.Hs.eg.db', ont = "BP", pAdjustMethod = "fdr",
                              pvalueCutoff = 0.05, readable = TRUE)
#tiff('HUVEC_pr_direct_vs_indirect_up_GO_BP.tiff', units="in", width=10, height=8, res=300, compression = 'lzw')
dotplot(huv_dir_indir_pr_up_bp, showCategory = 20, font.size = 18, label_format = 50)
dev.off()


huv_dir_indir_pr_dn_bp <- enrichGO(gene = huv_dind_pr_dn$ENTREZID, universe = huv_dind_pr_bg$ENTREZID,
                              OrgDb = 'org.Hs.eg.db', ont = "BP", pAdjustMethod = "fdr",
                              pvalueCutoff = 0.05, readable = TRUE)
# no result


###### huv_dir_indir_prot_GO_CC ####
huv_dir_indir_pr_up_cc <- enrichGO(gene = huv_dind_pr_up$ENTREZID, universe = huv_dind_pr_bg$ENTREZID,
                              OrgDb = 'org.Hs.eg.db', ont = "CC", pAdjustMethod = "fdr",
                              pvalueCutoff = 0.05, readable = TRUE)
#tiff('HUVEC_pr_direct_vs_indirect_up_GO_CC.tiff', units="in", width=10, height=5, res=300, compression = 'lzw')
dotplot(huv_dir_indir_pr_up_cc, showCategory = 20, font.size = 18, label_format = 50)
dev.off()

huv_dir_indir_pr_dn_cc <- enrichGO(gene = huv_dind_pr_dn$ENTREZID, universe = huv_dind_pr_bg$ENTREZID,
                              OrgDb = 'org.Hs.eg.db', ont = "CC", pAdjustMethod = "fdr",
                              pvalueCutoff = 0.05, readable = TRUE)
# no result


###### huv_dir_indir_prot_KEGG ####
huv_dir_indir_pr_up_kegg <- enrichKEGG(gene = huv_dind_pr_up$ENTREZID, universe = huv_dind_pr_bg$ENTREZID, organism = "hsa",
                                  keyType = "kegg", pvalueCutoff = 0.05, pAdjustMethod = "BH", 
                                  minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2,
                                  use_internal_data = FALSE)
#tiff('HUVEC_pr_direct_vs_indirect_up_KEGG.tiff', units="in", width=10, height=4, res=300, compression = 'lzw')
dotplot(huv_dir_indir_pr_up_kegg, showCategory = 20, font.size = 18, label_format = 50)
dev.off()


huv_dir_indir_pr_dn_kegg <- enrichKEGG(gene = huv_dind_pr_dn$ENTREZID, universe = huv_dind_pr_bg$ENTREZID, organism = "hsa",
                                  keyType = "kegg", pvalueCutoff = 0.05, pAdjustMethod = "BH",
                                  minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2,
                                  use_internal_data = FALSE)
# no result


# HUVEC transcriptome vs proteome analysis ####
setwd("../HUVEC_transcriptome_vs_proteome")

## Venn diagram transcriptome vs proteome all genes ####
venn_col <- c("#0073C2FF", "#CD534CFF")
venn_huv_pr_tr_all <- list(transcriptome = rownames(huv_tr), 
                           proteome = rownames(huv_pr_norm))
#tiff('venn_huv_pr_tr_all.tiff', units="in", width=8, height=5, res=300, compression = 'lzw')
ggvenn(venn_huv_pr_tr_all, fill_color = venn_col, stroke_size = 1, set_name_size = 10, text_size = 10)
dev.off()


## Venn diagram, Correlation and Hetmap ####

### Transcriptome vs proteome direct-control ####

#### Venn diagram up direct-control differentiated genes ####
venn_huv_pr_tr_cdir_up <- list(RNA_up = rownames(huv_contr_dir_tr_up), 
                               Prot.up = rownames(huv_contr_dir_pr_up))
#tiff('venn_huv_pr_tr_cdir_up.tiff', units="in", width=8, height=5, res=300, compression = 'lzw')
ggvenn(venn_huv_pr_tr_cdir_up, fill_color = venn_col, stroke_size = 1, set_name_size = 10, text_size = 10)
dev.off()


#### Venn diagram down direct-control differentiated genes ####
venn_huv_pr_tr_cdir_down <- list(RNA_down = rownames(huv_contr_dir_tr_dn), 
              Prot.down = rownames(huv_contr_dir_pr_dn))
#tiff('venn_huv_pr_tr_cdir_down.tiff', units="in", width=8, height=5, res=300, compression = 'lzw')
ggvenn(venn_huv_pr_tr_cdir_down, fill_color = venn_col, stroke_size = 1, set_name_size = 10, text_size = 10)
dev.off()


#### Correlation of log2 Fold Change between direct-control differentiated genes ####
huv_contr_dir_tr_pr_gene <- intersect(rownames(huv_contr_dir_tr), 
                                      rownames(huv_contr_dir_pr))
huv_contr_dir_tr_FC <- huv_contr_dir_tr[huv_contr_dir_tr_pr_gene, 'log2FoldChange']
huv_contr_dir_pr_FC <- huv_contr_dir_pr[huv_contr_dir_tr_pr_gene, 'logFC']
huv_contr_dir_tr_pr <- data.frame(huv_contr_dir_tr_pr_gene,
                                  huv_contr_dir_tr_FC,
                                  huv_contr_dir_pr_FC)
colnames(huv_contr_dir_tr_pr) <- c('Gene_name', 'Log2FC_tr', 'Log2FC_pr')

#tiff('HUVEC_pr_tr_contr_vs_direct.tiff', units="in", width=12, height=6, res=300, compression = 'lzw')
ggplot(huv_contr_dir_tr_pr, aes(x=Log2FC_tr, y=Log2FC_pr)) +
  geom_point() + 
  geom_vline(xintercept = 1, linetype = 2) + geom_vline(xintercept = -1, linetype = 2) + 
  geom_hline(yintercept = 1, linetype = 2) + geom_hline(yintercept = -1, linetype = 2) + 
  #geom_smooth(method="loess")+
  xlab("RNA log2 Fold Change") +
  ylab("protein log2 Fold Change") +
  ggtitle("correlation of log2 Fold Change between direct and control condition in HUVEC for proteomics and RNA-seq data") + 
  stat_cor(method = "pearson") +
  geom_text(aes(label=ifelse(Log2FC_tr>1 & Log2FC_pr >1,as.character(Gene_name),'')),hjust=0,vjust=-0.5) + 
  geom_text(aes(label=ifelse(Log2FC_pr<(-1) & Log2FC_pr<(-1),as.character(Gene_name),'')),hjust=0,vjust=-0.5) +
  annotate("text", 
           x = 12, # should be optimized manually
           y = 4.1,  # should be optimized manually
           label = paste("Number of genes: ", nrow(huv_contr_dir_tr_pr)))
dev.off()


# Heatmap 
color <- list(sample = c(control = '#FF7370', direct = '#01AFDF', indirect = '#00E0A5'))
tr_hm <- data.frame(sample = fact_huv_tr$Condition)
row.names(tr_hm) <- rownames(fact_huv_tr)

pr_hm <- data.frame(sample = fact_huv_pr$Condition)
row.names(pr_hm) <- rownames(fact_huv_pr)

#### Heatmap direct-control differentiated genes ####

# Find common transcriptome and proteome control-direct DEG
huv_contr_dir_tr_dif <- huv_contr_dir_tr[abs(huv_contr_dir_tr$log2FoldChange) >= 1 & huv_contr_dir_tr$padj <= 0.05,]
huv_contr_dir_pr_dif <- huv_contr_dir_pr[abs(huv_contr_dir_pr$logFC) >= 1 & huv_contr_dir_pr$adj.P.Val <= 0.05,]
huv_contr_dir_tr_pr_dif <- intersect(rownames(huv_contr_dir_tr_dif), rownames(huv_contr_dir_pr_dif))

# Transcriptome heatmap direct-control differentiated genes

# Genes of interest, top common tr-pr DEGs and top DEGs
select <- c('NOTCH3', 'WNT5A', 'WNT5B', 'CRYAB', 'FGF7', 'FGF1', 'FGF16', 'IRS1', 'ADAM12', 'MMP1', huv_contr_dir_tr_pr_dif[1:10], rownames(huv_contr_dir_tr_dif[1:30,]))
select <- unique(select)

#tiff('HUVEC_tr_direct_vs_control_heatmap.tiff', units="in", width=5, height=6, res=300, compression = 'lzw')
pheatmap(huv_tr[select[1:30],], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=tr_hm, 
         annotation_colors = color, show_colnames = FALSE,
         treeheight_col = 20, treeheight_row = 20)
dev.off()


# Proteome heatmap direct-control differentiated genes

# Find common proteome and transcriptome control-direct DEG
huv_contr_dir_tr_pr_dif <- intersect(rownames(huv_contr_dir_pr_dif), rownames(huv_contr_dir_tr_dif))

# Genes of interest, top common tr-pr DEGs and top DEGs
select <- c('POSTN', huv_contr_dir_tr_pr_dif[1:10], rownames(huv_contr_dir_pr_dif[1:30,]))
select <- unique(select)

#tiff('HUVEC_pr_direct_vs_control_heatmap.tiff', units="in", width=6, height=6, res=300, compression = 'lzw')
pheatmap(huv_pr_norm[select[1:30],], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=pr_hm,
         annotation_colors = color, show_colnames = FALSE,
         treeheight_col = 20, treeheight_row = 20)
dev.off()


### Transcriptome vs proteome indirect-control differentiated genes ####

#### Venn diagram up indirect-control differentiated genes ####

venn_huv_pr_tr_cindir_up <- list(RNA_up = rownames(huv_contr_indir_tr_up), 
                               Prot.up = rownames(huv_contr_indir_pr_up))
#tiff('venn_huv_pr_tr_cindir_up.tiff', units="in", width=8, height=5, res=300, compression = 'lzw')
ggvenn(venn_huv_pr_tr_cindir_up, fill_color = venn_col, stroke_size = 1, set_name_size = 10, text_size = 10)
dev.off()


#### Venn diagram down indirect-control differentiated genes ####
venn_huv_pr_tr_cindir_down <- list(RNA_down = rownames(huv_contr_indir_tr_dn), 
                                 Prot.down = rownames(huv_contr_indir_pr_dn))
#tiff('venn_huv_pr_tr_cindir_down.tiff', units="in", width=8, height=5, res=300, compression = 'lzw')
ggvenn(venn_huv_pr_tr_cindir_down, fill_color = venn_col, stroke_size = 1, set_name_size = 10, text_size = 10)
dev.off()


#### Correlation of log2 Fold Change between indirect-control differentiated genes ####
huv_contr_indir_tr_pr_gene <- intersect(rownames(huv_contr_indir_tr), 
                                      rownames(huv_contr_indir_pr))
huv_contr_indir_tr_FC <- huv_contr_indir_tr[huv_contr_indir_tr_pr_gene, 'log2FoldChange']
huv_contr_indir_pr_FC <- huv_contr_indir_pr[huv_contr_indir_tr_pr_gene, 'logFC']
huv_contr_indir_tr_pr <- data.frame(huv_contr_indir_tr_pr_gene,
                                  huv_contr_indir_tr_FC,
                                  huv_contr_indir_pr_FC)
colnames(huv_contr_indir_tr_pr) <- c('Gene_name', 'Log2FC_tr', 'Log2FC_pr')

#tiff('HUVEC_pr_tr_contr_vs_indirect.tiff', units="in", width=12, height=6, res=300, compression = 'lzw')
ggplot(huv_contr_indir_tr_pr, aes(x=Log2FC_tr, y=Log2FC_pr)) +
  geom_point() + 
  geom_vline(xintercept = 1, linetype = 2) + geom_vline(xintercept = -1, linetype = 2) + 
  geom_hline(yintercept = 1, linetype = 2) + geom_hline(yintercept = -1, linetype = 2) + 
  #geom_smooth(method="loess")+
  xlab("RNA log2 Fold Change") +
  ylab("protein log2 Fold Change") +
  ggtitle("correlation of log2 Fold Change between indirect and control condition in HUVEC for proteomics and RNA-seq data") + 
  stat_cor(method = "pearson") +
  geom_text(aes(label=ifelse(Log2FC_tr>=1 & Log2FC_pr>=1, as.character(Gene_name),'')),hjust=0,vjust=-0.5) + 
  geom_text(aes(label=ifelse(Log2FC_pr<=(-1) & Log2FC_pr<=(-1), as.character(Gene_name),'')),hjust=0,vjust=-0.5) +
  annotate("text", 
           x = 2.6, # should be optimized manually
           y = 2.6,  # should be optimized manually
           label = paste("Number of genes: ", nrow(huv_contr_indir_tr_pr)))
dev.off()


#### Heatmap transcriptome and proteome control-indirect differentiated genes ####

# Find common transcriptome and proteome control-indirect DEG
huv_contr_indir_tr_dif <- huv_contr_indir_tr[abs(huv_contr_indir_tr$log2FoldChange) >= 1 & huv_contr_indir_tr$padj <= 0.05,]
huv_contr_indir_pr_dif <- huv_contr_indir_pr[abs(huv_contr_indir_pr$logFC) >= 1 & huv_contr_indir_pr$adj.P.Val <= 0.05,]
huv_contr_indir_tr_pr_dif <- intersect(rownames(huv_contr_indir_tr_dif), rownames(huv_contr_indir_pr_dif))


#### Transcriptome heatmap control-indirect differentiated genes ####

# Genes of interest, top common tr-pr DEGs and top DEGs
select <- c('NOTCH4', 'FGF16', 'IRS2', 'BMPER', 'MYCN', 'TAGLN', 'CRYAB', huv_contr_indir_tr_pr_dif[1:10], rownames(huv_contr_indir_tr_dif[1:30,]))
select <- unique(select)

#tiff('HUVEC_tr_indirect_vs_control_heatmap.tiff', units="in", width=5, height=6, res=300, compression = 'lzw')
pheatmap(huv_tr[select[1:30],], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=tr_hm,
         annotation_colors = color, show_colnames = FALSE,
         treeheight_col = 20, treeheight_row = 20)
dev.off()


#### Proteome heatmap control-indirect differentiated genes ####

# Find common trascriptome and proteome control-indirect DEG
huv_contr_indir_tr_pr_dif <- intersect(rownames(huv_contr_indir_pr_dif), rownames(huv_contr_indir_tr_dif))
huv_contr_indir_pr_dif <- huv_contr_indir_pr_dif[order(huv_contr_indir_pr_dif$logFC),]

# Genes of interest, top common tr-pr DEGs and top DEGs
select <- c('NOTCH3', huv_contr_indir_tr_pr_dif[1:5], rownames(huv_contr_indir_pr_dif[1:30,]))
select <- unique(select)

#tiff('HUVEC_pr_indirect_vs_control_heatmap.tiff', units="in", width=6, height=6, res=300, compression = 'lzw')
pheatmap(huv_pr_norm[select[1:30],], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=pr_hm,
         annotation_colors = color, show_colnames = FALSE,
         treeheight_col = 20, treeheight_row = 20)
dev.off()


### Transcriptome vs proteome direct-indirect differentiated genes ####

#### Venn diagram up direct-indirect differentiated genes ####
venn_huv_pr_tr_dindir_up <- list(RNA_up = rownames(huv_dir_indir_tr_up), 
                                 Prot.up = rownames(huv_dir_indir_pr_up))
#tiff('venn_huv_pr_tr_dindir_up.tiff', units="in", width=8, height=5, res=300, compression = 'lzw')
ggvenn(venn_huv_pr_tr_dindir_up, fill_color = venn_col, stroke_size = 1, set_name_size = 10, text_size = 10)
dev.off()

#### Venn diagram down direct-indirect differentiated genes ####
venn_huv_pr_tr_dindir_down <- list(RNA_down = rownames(huv_dir_indir_tr_dn), 
                                   Prot.down = rownames(huv_dir_indir_pr_dn))
#tiff('venn_huv_pr_tr_dindir_down.tiff', units="in", width=8, height=5, res=300, compression = 'lzw')
ggvenn(venn_huv_pr_tr_dindir_down, fill_color = venn_col, stroke_size = 1, set_name_size = 10, text_size = 10)
dev.off()

#### Correlation of log2 Fold Change between indirect-direct differentiated genes ####
huv_dir_indir_tr_pr_gene <- intersect(rownames(huv_dir_indir_tr), 
                                        rownames(huv_dir_indir_pr))
huv_dir_indir_tr_FC <- huv_dir_indir_tr[huv_dir_indir_tr_pr_gene, 'log2FoldChange']
huv_dir_indir_pr_FC <- huv_dir_indir_pr[huv_dir_indir_tr_pr_gene, 'logFC']
huv_dir_indir_tr_pr <- data.frame(huv_dir_indir_tr_pr_gene,
                                    huv_dir_indir_tr_FC,
                                    huv_dir_indir_pr_FC)
colnames(huv_dir_indir_tr_pr) <- c('Gene_name', 'Log2FC_tr', 'Log2FC_pr')

#tiff('HUVEC_pr_tr_dir_vs_indirect.tiff', units="in", width=12, height=6, res=300, compression = 'lzw')
ggplot(huv_dir_indir_tr_pr, aes(x=Log2FC_tr, y=Log2FC_pr)) +
  geom_point() + 
  geom_vline(xintercept = 1, linetype = 2) + geom_vline(xintercept = -1, linetype = 2) + 
  geom_hline(yintercept = 1, linetype = 2) + geom_hline(yintercept = -1, linetype = 2) + 
  #geom_smooth(method="loess")+
  xlab("RNA log2 Fold Change") +
  ylab("protein log2 Fold Change") +
  ggtitle("correlation of log2 Fold Change between indirect and direct condition in HUVEC for proteomics and RNA-seq data") + 
  stat_cor(method = "pearson") +
  geom_text(aes(label=ifelse(Log2FC_tr>=1 & Log2FC_pr>=1, as.character(Gene_name),'')),hjust=0,vjust=-0.5) + 
  geom_text(aes(label=ifelse(Log2FC_pr<=(-1) & Log2FC_pr<=(-1), as.character(Gene_name),'')),hjust=0,vjust=-0.5) +
  annotate("text", 
           x = 9, # should be optimized manually
           y = 3.5,  # should be optimized manually
           label = paste("Number of genes: ", nrow(huv_dir_indir_tr_pr)))
dev.off()


#### Heatmap direct-indirect differentiated genes ####
# Find common transcriptome and proteome indirect-direct DEG
huv_dir_indir_tr_dif <- huv_dir_indir_tr[abs(huv_dir_indir_tr$log2FoldChange) >= 1 & huv_dir_indir_tr$padj <= 0.05,]
huv_dir_indir_pr_dif <- huv_dir_indir_pr[abs(huv_dir_indir_pr$logFC) >= 1 & huv_dir_indir_pr$adj.P.Val <= 0.05,]
huv_dir_indir_tr_pr_dif <- intersect(rownames(huv_dir_indir_tr_dif), rownames(huv_dir_indir_pr_dif))


# Transcriptome heatmap direct-indirect differentiated genes ####
# Genes of interest, top common tr-pr DEGs and top DEGs
select <- c('NOTCH4', 'FGF7', 'FGF2', 'FGF5', 'ADAM12', 'ALPL', 'ANK2', 'IRS1', 'MYCN', 'TAGLN', 'CRYAB', huv_dir_indir_tr_pr_dif[1:10], rownames(huv_dir_indir_tr_dif[1:30,]))
select <- unique(select)

#tiff('HUVEC_tr_indirect_vs_direct_heatmap.tiff', units="in", width=5, height=6, res=300, compression = 'lzw')
pheatmap(huv_tr[select[1:30],], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=tr_hm,
         annotation_colors = color, show_colnames = FALSE,
         treeheight_col = 20, treeheight_row = 20)
dev.off()


# Proteome heatmap direct-indirect differentiated genes ####
# Find common transcriptome and proteome indirect-direct DEG
huv_dir_indir_tr_pr_dif <- intersect(rownames(huv_dir_indir_pr_dif), rownames(huv_dir_indir_tr_dif))

# Top common tr-pr DEGs and top DEGs
select <- c(huv_dir_indir_tr_pr_dif[1:10], rownames(huv_dir_indir_pr_dif[1:30,]))
select <- unique(select)

#tiff('HUVEC_pr_indirect_vs_indorect_heatmap.tiff', units="in", width=6, height=6, res=300, compression = 'lzw')
pheatmap(huv_pr_norm[select[1:30],], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=pr_hm,
         annotation_colors = color, show_colnames = FALSE,
         treeheight_col = 20, treeheight_row = 20)
dev.off()
