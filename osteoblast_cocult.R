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

# Osteoblast analysis ####
setwd("./Osteoblast_transcriptome")
## Osteoblast transcriptome analysis ####

# separate Osteoblast's legend and data
fact_ost_tr <- tr_legend[tr_legend$Cell.type == 'Osteoblast', ]
ost_tr <- tr_data[, rownames(fact_ost_tr)]
#write.csv(ost_tr, "Osteoblast_tr_row_data")

# delete Cell.type and Sample columns
fact_ost_tr <- fact_ost_tr[, c(-1, -2)]
# make variables factors
fact_ost_tr$Condition <- as.factor(fact_ost_tr$Condition)
fact_ost_tr$Donor <- as.factor(fact_ost_tr$Donor)
str(fact_ost_tr)

# Construct DESEQDataSet Object
dds_ost <- DESeqDataSetFromMatrix(countData = ost_tr,
                                  colData = fact_ost_tr,
                                  design= ~ Donor + Condition)
dds_ost

# minimum filtering
keep <- rowSums(counts(dds_ost)) >= 10
dds_ost <- dds_ost[keep,]

# differential expression analysis
dds_ost <- DESeq(dds_ost)
res <- results(dds_ost, alpha=0.05)
resultsNames(dds_ost)

# extracting regularized-logarithm transformed values
rld_ost <- rlog(dds_ost, blind=FALSE)
head(assay(rld_ost), 30)
# dispersion plot
meanSdPlot(assay(rld_ost))
ost_tr <- assay(rld_ost)
#write.csv(ost_tr, 'Osteoblast_tr_reg_log_data.csv')


### Сlustering ######
### Sample Principal Component Plot
ost_tr_pca <- pca(t(ost_tr), ncomp = 5, center = TRUE)
#tiff('PCA_Osteoblast_tr.tiff', units="in", width=10, height=8, res=300, compression = 'lzw')
plotIndiv(ost_tr_pca, comp = c(1, 2), ind.names = F, group = fact_ost_tr$Condition, legend = TRUE, ellipse = T, title = 'PCA', style = "ggplot2", cex = 6, size.axis = 14, size.xlabel = 18, size.ylabel = 18, size.title = 22, size.legend = 18, size.legend.title = 18)
dev.off()

plotIndiv(ost_tr_pca, comp = c(1, 2), ind.names = T, group = fact_ost_tr$Donor, legend = TRUE, ellipse = T, title = 'PCA')
dev.off()

#PLS-DA 
ordination.optimum.splsda <- splsda(t(ost_tr), fact_ost_tr$Condition, ncomp = 3, keepX = c(15,15,15))
selectVar(ordination.optimum.splsda, comp=1)
selectVar(ordination.optimum.splsda, comp=2)
selectVar(ordination.optimum.splsda, comp=3)


#layout(matrix(c(1, 2), 2, 3))
plotLoadings(ordination.optimum.splsda, comp = 1, size.name = 1, size.title = 1.2, title = "Loadings\n on 1st component", contrib = "max", legend = FALSE, col.ties="black", ndisplay = 15)

plotLoadings(ordination.optimum.splsda, comp = 2, size.name = 1, size.title = 1.2, title = "Loadings\n on 2nd component", contrib = "max",ndisplay = 15,  legend = FALSE, col.ties="black")

#tiff('PLSDA_Osteoblast_tr.tiff', units="in", width=10, height=8, res=300, compression = 'lzw')
plotIndiv(ordination.optimum.splsda, ind.names = F, ellipse = T, style = "ggplot2", title = "PLS-DA ordination", legend=TRUE, cex = 6, size.axis = 14, size.xlabel = 18, size.ylabel = 18, size.title = 22, size.legend = 18, size.legend.title = 18)
dev.off()

#tiff('PLSDA_Osteoblast_transcr.tiff', units="in", width=16, height=8, res=300, compression = 'lzw')
layout(matrix(c(1, 2, 3, 3, 3, 3), 2, 3))
plotLoadings(ordination.optimum.splsda, comp = 1, size.name = 1, size.title = 1.2, 
             title = "Loadings\n on 1st component", contrib = "max", legend = FALSE, 
             col.ties="black", ndisplay = 15)
plotLoadings(ordination.optimum.splsda, comp = 2, size.name = 1, size.title = 1.2, 
             title = "Loadings\n on 2nd component", contrib = "max", ndisplay = 15,  
             legend = FALSE, col.ties="black")
plotIndiv(ordination.optimum.splsda, ind.names = F, ellipse = T, style = "graphics",
          title = "PLS-DA ordination",abline = TRUE, cex = 2.5, size.axis = 2.5, 
          size.xlabel = 2.8, size.ylabel = 2.8, size.title = 2.5, legend=TRUE, 
          size.legend = 2, size.legend.title = 2.5)
dev.off()

#tiff('PCA_PLSDA_Osteoblast_tr.tiff', units="in", width=18, height=8, res=300, compression = 'lzw')
layout(matrix(1:2, ncol = 2))
plotIndiv(ost_tr_pca, comp = c(1, 2), ind.names = F, group = fact_ost_tr$Condition, 
          style = "graphics", ellipse = T, title = 'PCA', 
          legend=TRUE, size.title = 2, size.legend = 1.5, size.legend.title = 1.8, 
          abline = TRUE, cex = 2.5, size.axis = 1.5, size.xlabel = 1.8, size.ylabel = 1.8)
plotIndiv(ordination.optimum.splsda, ind.names = F, ellipse = T, style = "graphics",
          title = "PLS-DA ordination", legend=TRUE, size.title = 2,
          size.legend = 1.5, size.legend.title = 1.8, abline = TRUE, cex = 2.5, size.axis = 1.5, size.xlabel = 1.8, size.ylabel = 1.8)
dev.off()


### Differential analysis ####
resultsNames(dds_ost)

#### Condition_direct_vs_сontrol ####
ost_contr_dir_tr <- results(dds_ost, name="Condition_direct_vs_control", alpha=0.05)
#write.csv(ost_contr_dir_tr, file='Osteoblast_tr_direct_vs_control.csv')

# shrink log fold changes association with condition:
resLFC_ost_contr_dir_tr <- lfcShrink(dds_ost, coef="Condition_direct_vs_control", 
                                     type="apeglm")
# MA-plot
plotMA(ost_contr_dir_tr, ylim=c(-2,2))
plotMA(resLFC_ost_contr_dir_tr, ylim=c(-2,2))
dev.off()

# order the results table by the smallest adjusted p-value:
resLFC_ost_contr_dir_tr
ost_contr_dir_tr <- resLFC_ost_contr_dir_tr[order(resLFC_ost_contr_dir_tr$padj),]
# write the ordered table to a file
#write.csv(ost_contr_dir_tr, "Osteoblast_tr_direct_vs_control_apeglm.csv")


# Vulcano plot
#tiff('Osteoblast_tr_direct_vs_control_vulcano.tiff', units="in", width=16, height=12, res=300, compression = 'lzw')
EnhancedVolcano(ost_contr_dir_tr,
                lab = rownames(ost_contr_dir_tr),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1,
                xlim = c(-8, 12),
                ylim = c(0, 260),
                title ="Control-direct cocultivation. Osteoblast",
                labSize = 7.0,
                pointSize = 3.5,
                axisLabSize = 25,
                titleLabSize = 30,
                subtitleLabSize = 1,
                captionLabSize = 20,
                boxedLabels = F,
                colAlpha = 1)
dev.off()

head(ost_contr_dir_tr)
# counts of reads for a single gene across the Condition
boxplot(as.numeric(ost_tr['COL18A1',]) ~ Condition, data = fact_ost_tr, varwidth = TRUE, log = "y", las = 1, ylab = "COL18A1")
boxplot(as.numeric(ost_tr['HHIP',]) ~ Condition, data = fact_ost_tr, varwidth = TRUE, log = "y", las = 1, ylab = "HHIP")
boxplot(as.numeric(ost_tr['HEYL',]) ~ Condition, data = fact_ost_tr, varwidth = TRUE, log = "y", las = 1, ylab = "HEYL")
boxplot(as.numeric(ost_tr['LAMA5',]) ~ Condition, data = fact_ost_tr, varwidth = TRUE, log = "y", las = 1, ylab = "LAMA5")

# control-direct up and down regulated transcripts
# change rows with NA (in p-value)
ost_contr_dir_tr$log2FoldChange[is.na(ost_contr_dir_tr$log2FoldChange)] <- 0
ost_contr_dir_tr$padj[is.na(ost_contr_dir_tr$padj)] <- 1

ost_contr_dir_tr_up <- ost_contr_dir_tr[ost_contr_dir_tr$log2FoldChange >= 1 & ost_contr_dir_tr$padj <= 0.05,]
ost_contr_dir_tr_dn <- ost_contr_dir_tr[ost_contr_dir_tr$log2FoldChange <= -1 & ost_contr_dir_tr$padj <= 0.05,]

#write.csv(ost_contr_dir_tr_up, file='Osteoblast_tr_direct_vs_control_up.csv')
#write.csv(ost_contr_dir_tr_dn, file='Osteoblast_tr_direct_vs_control_down.csv')


##### Control-Direct enrichment analysis #####

ost_cdir_tr_up <- bitr(rownames(ost_contr_dir_tr_up), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
ost_cdir_tr_dn <- bitr(rownames(ost_contr_dir_tr_dn), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
ost_cdir_tr_bg <- bitr(rownames(ost_contr_dir_tr), fromType = "SYMBOL", toType = "ENTREZID",OrgDb = "org.Hs.eg.db")


###### ost_contr_dir_trans_GO_BP ####

ost_cdir_tr_up_bp <- enrichGO(gene = ost_cdir_tr_up$ENTREZID, universe = ost_cdir_tr_bg$ENTREZID, OrgDb = "org.Hs.eg.db", ont = "BP", pAdjustMethod = "fdr", pvalueCutoff = 0.05, readable = TRUE)
#tiff('Osteoblast_tr_direct_vs_control_up_GO_BP.tiff', units="in", width=10, height=8, res=300, compression = 'lzw')
dotplot(ost_cdir_tr_up_bp, showCategory = 20, font.size = 18, label_format = 50)
dev.off()


ost_cdir_tr_dn_bp <- enrichGO(gene = ost_cdir_tr_dn$ENTREZID, universe = ost_cdir_tr_bg$ENTREZID, OrgDb = "org.Hs.eg.db", ont = "BP", pAdjustMethod = "fdr", pvalueCutoff = 0.05, readable = TRUE)
#tiff('Osteoblast_tr_direct_vs_control_down_GO_BP.tiff', units="in", width=11, height=8, res=300, compression = 'lzw')
dotplot(ost_cdir_tr_dn_bp, showCategory = 18, font.size = 18, label_format = 50)
dev.off()


###### ost_contr_dir_trans_GO_CC ####
ost_cdir_tr_up_cc <- enrichGO(gene = ost_cdir_tr_up$ENTREZID, universe = ost_cdir_tr_bg$ENTREZID, OrgDb = "org.Hs.eg.db", ont = "CC", pAdjustMethod = "fdr", pvalueCutoff = 0.05, readable = TRUE)
#tiff('Osteoblast_tr_direct_vs_control_up_GO_CC.tiff', units="in", width=10, height=8, res=300, compression = 'lzw')
dotplot(ost_cdir_tr_up_cc, showCategory = 20, font.size = 18, label_format = 50)
dev.off()

ost_cdir_tr_dn_cc <- enrichGO(gene = ost_cdir_tr_dn$ENTREZID, universe = ost_cdir_tr_bg$ENTREZID, OrgDb = "org.Hs.eg.db", ont = "CC", pAdjustMethod = "fdr", pvalueCutoff = 0.05, readable = TRUE)
#tiff('Osteoblast_tr_direct_vs_control_down_GO_CC.tiff', units="in", width=10, height=8, res=300, compression = 'lzw')
dotplot(ost_cdir_tr_dn_cc, showCategory = 20, font.size = 18, label_format = 50)
dev.off()

###### ost_contr_dir_trans_KEGG ####
ost_cdir_tr_up_kegg <- enrichKEGG(gene = ost_cdir_tr_up$ENTREZID, universe = ost_cdir_tr_bg$ENTREZID, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05, pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2, use_internal_data = FALSE)
#tiff('Osteoblast_tr_direct_vs_control_up_KEGG.tiff', units="in", width=10, height=8, res=300, compression = 'lzw')
dotplot(ost_cdir_tr_up_kegg, showCategory = 20, font.size = 18, label_format = 45)
dev.off()


ost_cdir_tr_dn_kegg <- enrichKEGG(gene = ost_cdir_tr_dn$ENTREZID, universe = ost_cdir_tr_bg$ENTREZID, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05, pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2, use_internal_data = FALSE)
#tiff('Osteoblast_tr_direct_vs_control_down_KEGG.tiff', units="in", width=11, height=8, res=300, compression = 'lzw')
dotplot(ost_cdir_tr_dn_kegg, showCategory = 20, font.size = 18, label_format = 50)
dev.off()



#### Condition_indirect_vs_сontrol ####
ost_contr_indir_tr <- results(dds_ost, name="Condition_indirect_vs_control", alpha=0.05)
#write.csv(ost_contr_indir_tr, file='Osteoblast_tr_indirect_vs_control.csv')

# shrink log fold changes association with condition:
resLFC_ost_contr_indir_tr <- lfcShrink(dds_ost, coef="Condition_indirect_vs_control", 
                                       type="apeglm")
# MA-plot
plotMA(ost_contr_indir_tr, ylim=c(-2,2))
plotMA(resLFC_ost_contr_indir_tr, ylim=c(-2,2))
dev.off()

# order the results table by the smallest adjusted p-value:
resLFC_ost_contr_indir_tr
ost_contr_indir_tr <- resLFC_ost_contr_indir_tr[order(resLFC_ost_contr_indir_tr$padj),]

# write the ordered table to a file
#write.csv(ost_contr_indir_tr, file='Osteoblast_tr_indirect_vs_control_apeglm.csv')


# Vulcano plot
#tiff('Osteoblast_tr_indirect_vs_control_vulcano.tiff', units="in", width=16, height=12, res=300, compression = 'lzw')
EnhancedVolcano(ost_contr_indir_tr,
                lab = rownames(ost_contr_indir_tr),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1,
                xlim = c(-2.8, 2.5),
                ylim = c(0, 23),
                title ="Control-indirect cocultivation. Osteoblast",
                labSize = 7.0,
                pointSize = 3.5,
                axisLabSize = 25,
                titleLabSize = 30,
                subtitleLabSize = 1,
                captionLabSize = 20,
                boxedLabels = F,
                colAlpha = 1)
dev.off()

head(ost_contr_indir_tr)
# counts of reads for a single gene across the Condition
boxplot(as.numeric(ost_tr['MT1L',]) ~ Condition, data = fact_ost_tr, varwidth = TRUE, log = "y", las = 1, ylab = "MT1L")
boxplot(as.numeric(ost_tr['MT2A',]) ~ Condition, data = fact_ost_tr, varwidth = TRUE, log = "y", las = 1, ylab = "MT2A")
boxplot(as.numeric(ost_tr['HMGA1',]) ~ Condition, data = fact_ost_tr, varwidth = TRUE, log = "y", las = 1, ylab = "HMGA1")
boxplot(as.numeric(ost_tr['TFPI2',]) ~ Condition, data = fact_ost_tr, varwidth = TRUE, log = "y", las = 1, ylab = "TFPI2")

# indirect up and down regulated transcripts
# change rows with NA (in p-value)
ost_contr_indir_tr$log2FoldChange[is.na(ost_contr_indir_tr$log2FoldChange)] <- 0
ost_contr_indir_tr$padj[is.na(ost_contr_indir_tr$padj)] <- 1

ost_contr_indir_tr_up <- ost_contr_indir_tr[ost_contr_indir_tr$log2FoldChange >= 1 & ost_contr_indir_tr$padj <= 0.05,]
ost_contr_indir_tr_dn <- ost_contr_indir_tr[ost_contr_indir_tr$log2FoldChange <= -1 & ost_contr_indir_tr$padj <= 0.05,]

#write.csv(ost_contr_indir_tr_up, file='Osteoblast_tr_indirect_vs_control_up.csv')
#write.csv(ost_contr_indir_tr_dn, file='Osteoblast_tr_indirect_vs_control_down.csv')


##### Control-indirect enrichment analysis #####

ost_cindir_tr_up <- bitr(rownames(ost_contr_indir_tr_up), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
ost_cindir_tr_dn <- bitr(rownames(ost_contr_indir_tr_dn), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
ost_cindir_tr_bg <- bitr(rownames(ost_contr_indir_tr), fromType = "SYMBOL", toType = "ENTREZID",OrgDb = "org.Hs.eg.db")

###### ost_contr_indir_trans_GO_BP ####

ost_cindir_tr_up_bp <- enrichGO(gene = ost_cindir_tr_up$ENTREZID, universe = ost_cindir_tr_bg$ENTREZID, OrgDb = "org.Hs.eg.db", ont = "BP", pAdjustMethod = "fdr", pvalueCutoff = 0.05, readable = TRUE)
#tiff('Osteoblast_tr_indirect_vs_control_up_GO_BP.tiff', units="in", width=11, height=8, res=300, compression = 'lzw')
dotplot(ost_cindir_tr_up_bp, showCategory = 19, font.size = 18, label_format = 50)
dev.off()


ost_cindir_tr_dn_bp <- enrichGO(gene = ost_cindir_tr_dn$ENTREZID, universe = ost_cindir_tr_bg$ENTREZID, OrgDb = "org.Hs.eg.db", ont = "BP", pAdjustMethod = "fdr", pvalueCutoff = 0.05, readable = TRUE)
#tiff('Osteoblast_tr_indirect_vs_control_down_GO_BP.tiff', units="in", width=10, height=8, res=300, compression = 'lzw')
dotplot(ost_cindir_tr_dn_bp, showCategory = 20, font.size = 18, label_format = 50)
dev.off()


###### ost_contr_indir_trans_GO_CC ####
ost_cindir_tr_up_cc <- enrichGO(gene = ost_cindir_tr_up$ENTREZID, universe = ost_cindir_tr_bg$ENTREZID, OrgDb = "org.Hs.eg.db", ont = "CC", pAdjustMethod = "fdr", pvalueCutoff = 0.05, readable = TRUE)
#tiff('Osteoblast_tr_indirect_vs_control_up_GO_CC.tiff', units="in", width=10, height=8, res=300, compression = 'lzw')
dotplot(ost_cindir_tr_up_cc, showCategory = 20, font.size = 18, label_format = 50)
dev.off()

ost_cindir_tr_dn_cc <- enrichGO(gene = ost_cindir_tr_dn$ENTREZID, universe = ost_cindir_tr_bg$ENTREZID, OrgDb = "org.Hs.eg.db", ont = "CC", pAdjustMethod = "fdr", pvalueCutoff = 0.05, readable = TRUE)
#tiff('Osteoblast_tr_indirect_vs_control_down_GO_CC.tiff', units="in", width=10, height=6, res=300, compression = 'lzw')
dotplot(ost_cindir_tr_dn_cc, showCategory = 20, font.size = 18, label_format = 50)
dev.off()

###### ost_contr_indir_trans_KEGG ####
ost_cindir_tr_up_kegg <- enrichKEGG(gene = ost_cindir_tr_up$ENTREZID, universe = ost_cindir_tr_bg$ENTREZID, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05, pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2, use_internal_data = FALSE)
#tiff('Osteoblast_tr_indirect_vs_control_up_KEGG.tiff', units="in", width=10, height=3, res=300, compression = 'lzw')
dotplot(ost_cindir_tr_up_kegg, showCategory = 20, font.size = 18, label_format = 50)
dev.off()


ost_cindir_tr_dn_kegg <- enrichKEGG(gene = ost_cindir_tr_dn$ENTREZID, universe = ost_cindir_tr_bg$ENTREZID, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05, pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2, use_internal_data = FALSE)
# no result


#### Condition_indirect_vs_direct ####

# change factors levels to change the baseline for direct-indirect comparison
fact_ost_tr$Condition <- factor(fact_ost_tr$Condition, levels = c("indirect", "direct", "control"))
ost_tr <- tr_data[, rownames(fact_ost_tr)]

# Construct DESEQDataSet Object
dds_ost <- DESeqDataSetFromMatrix(countData = ost_tr,
                                  colData = fact_ost_tr,
                                  design= ~ Donor + Condition)
dds_ost

# minimum filtering
keep <- rowSums(counts(dds_ost)) >= 10
dds_ost <- dds_ost[keep,]

# differential expression analysis
dds_ost <- DESeq(dds_ost)
res <- results(dds_ost, alpha=0.05)
resultsNames(dds_ost)

# extracting regularized-logarithm transformed values
rld_ost <- rlog(dds_ost, blind=FALSE)
head(assay(rld_ost), 30)
ost_tr <- assay(rld_ost)

ost_dir_indir_tr <- results(dds_ost, name="Condition_direct_vs_indirect", alpha=0.05)
#write.csv(ost_dir_indir_tr, file='Osteoblast_tr_direct_vs_indirect.csv')

# shrink log fold changes association with condition:
resLFC_ost_dir_indir_tr <- lfcShrink(dds_ost, coef="Condition_direct_vs_indirect", 
                                     type="apeglm")
# MA-plot
plotMA(ost_dir_indir_tr, ylim=c(-2,2))
plotMA(resLFC_ost_dir_indir_tr, ylim=c(-2,2))
dev.off()

# order the results table by the smallest adjusted p-value:
resLFC_ost_dir_indir_tr
ost_dir_indir_tr <- resLFC_ost_dir_indir_tr[order(resLFC_ost_dir_indir_tr$padj),]

# write the ordered table to a file
#write.csv(ost_dir_indir_tr, file='Osteoblast_tr_direct_vs_indirect_apeglm.csv')


# Vulcano plot
#tiff('Osteoblast_tr_direct_vs_indirect_vulcano.tiff', units="in", width=16, height=12, res=300, compression = 'lzw')
EnhancedVolcano(ost_dir_indir_tr,
                lab = rownames(ost_dir_indir_tr),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1,
                xlim = c(-7, 16),
                ylim = c(0, 270),
                title ="Indirect-direct cocultivation. Osteoblast",
                labSize = 7.0,
                pointSize = 3.5,
                axisLabSize = 25,
                titleLabSize = 30,
                subtitleLabSize = 1,
                captionLabSize = 20,
                boxedLabels = F,
                colAlpha = 1)
dev.off()

head(ost_dir_indir_tr)
# counts of reads for a single gene across the Condition
boxplot(as.numeric(ost_tr['COL18A1',]) ~ Condition, data = fact_ost_tr, varwidth = TRUE, log = "y", las = 1, ylab = "COL18A1")
boxplot(as.numeric(ost_tr['HHIP',]) ~ Condition, data = fact_ost_tr, varwidth = TRUE, log = "y", las = 1, ylab = "HHIP")
boxplot(as.numeric(ost_tr['HEYL',]) ~ Condition, data = fact_ost_tr, varwidth = TRUE, log = "y", las = 1, ylab = "HEYL")
boxplot(as.numeric(ost_tr['LAMA5',]) ~ Condition, data = fact_ost_tr, varwidth = TRUE, log = "y", las = 1, ylab = "LAMA5")


# indirect up and down regulated transcripts
# change rows with NA (in p-value)
ost_dir_indir_tr$log2FoldChange[is.na(ost_dir_indir_tr$log2FoldChange)] <- 0
ost_dir_indir_tr$padj[is.na(ost_dir_indir_tr$padj)] <- 1

ost_dir_indir_tr_up <- ost_dir_indir_tr[ost_dir_indir_tr$log2FoldChange >= 1 & ost_dir_indir_tr$padj <= 0.05,]
ost_dir_indir_tr_dn <- ost_dir_indir_tr[ost_dir_indir_tr$log2FoldChange <= -1 & ost_dir_indir_tr$padj <= 0.05,]

#write.csv(ost_dir_indir_tr_up, file='Osteoblast_tr_direct_vs_indirect_up.csv')
#write.csv(ost_dir_indir_tr_dn, file='Osteoblast_tr_direct_vs_indirect_down.csv')

##### Indirect-direct enrichment analysis #####

ost_dindir_tr_up <- bitr(rownames(ost_dir_indir_tr_up), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
ost_dindir_tr_dn <- bitr(rownames(ost_dir_indir_tr_dn), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
ost_dindir_tr_bg <- bitr(rownames(ost_dir_indir_tr), fromType = "SYMBOL", toType = "ENTREZID",OrgDb = "org.Hs.eg.db")


###### ost_dir_indir_trans_GO_BP ####

ost_dindir_tr_up_bp <- enrichGO(gene = ost_dindir_tr_up$ENTREZID, universe = ost_dindir_tr_bg$ENTREZID, OrgDb = "org.Hs.eg.db", ont = "BP", pAdjustMethod = "fdr", pvalueCutoff = 0.05, readable = TRUE)
#tiff('Osteoblast_tr_direct_vs_indirect_up_GO_BP.tiff', units="in", width=10, height=8, res=300, compression = 'lzw')
dotplot(ost_dindir_tr_up_bp, showCategory = 20, font.size = 18, label_format = 50)
dev.off()


ost_dindir_tr_dn_bp <- enrichGO(gene = ost_dindir_tr_dn$ENTREZID, universe = ost_dindir_tr_bg$ENTREZID, OrgDb = "org.Hs.eg.db", ont = "BP", pAdjustMethod = "fdr", pvalueCutoff = 0.05, readable = TRUE)
#tiff('Osteoblast_tr_direct_vs_indirect_down_GO_BP.tiff', units="in", width=10, height=8, res=300, compression = 'lzw')
dotplot(ost_dindir_tr_dn_bp, showCategory = 20, font.size = 18, label_format = 50)
dev.off()

###### ost_dir_indir_trans_GO_CC ####
ost_dindir_tr_up_cc <- enrichGO(gene = ost_dindir_tr_up$ENTREZID, universe = ost_dindir_tr_bg$ENTREZID, OrgDb = "org.Hs.eg.db", ont = "CC", pAdjustMethod = "fdr", pvalueCutoff = 0.05, readable = TRUE)
#tiff('Osteoblast_tr_direct_vs_indirect_up_GO_CC.tiff', units="in", width=10, height=8, res=300, compression = 'lzw')
dotplot(ost_dindir_tr_up_cc, showCategory = 20, font.size = 18, label_format = 50)
dev.off()

ost_dindir_tr_dn_cc <- enrichGO(gene = ost_dindir_tr_dn$ENTREZID, universe = ost_dindir_tr_bg$ENTREZID, OrgDb = "org.Hs.eg.db", ont = "CC", pAdjustMethod = "fdr", pvalueCutoff = 0.05, readable = TRUE)
#tiff('Osteoblast_tr_direct_vs_indirect_down_GO_CC.tiff', units="in", width=10, height=9, res=300, compression = 'lzw')
dotplot(ost_dindir_tr_dn_cc, showCategory = 20, font.size = 18, label_format = 50)
dev.off()

###### ost_dir_indir_trans_KEGG ####
ost_dindir_tr_up_kegg <- enrichKEGG(gene = ost_dindir_tr_up$ENTREZID, universe = ost_dindir_tr_bg$ENTREZID, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05, pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2, use_internal_data = FALSE)
#tiff('Osteoblast_tr_direct_vs_indirect_up_KEGG.tiff', units="in", width=10, height=8, res=300, compression = 'lzw')
dotplot(ost_dindir_tr_up_kegg, showCategory = 20, font.size = 18, label_format = 50)
dev.off()

ost_dindir_tr_dn_kegg <- enrichKEGG(gene = ost_dindir_tr_dn$ENTREZID, universe = ost_dindir_tr_bg$ENTREZID, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05, pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2, use_internal_data = FALSE)
#tiff('Osteoblast_tr_direct_vs_indirect_down_KEGG.tiff', units="in", width=10, height=4, res=300, compression = 'lzw')
dotplot(ost_dindir_tr_dn_kegg, showCategory = 20, font.size = 18, label_format = 50)
dev.off()


## Osteoblast proteome analysis ####
setwd("../")
# read data
ost_pr <- data.frame(read.table("ost_log_imp_data.csv", sep = ',', header = TRUE))
# data is filtered, NA inserted with log_imp in NAguide
colSums(is.na(ost_pr))

# convert ID to gene symbol
colnames(ost_pr)[1] = 'UNIPROT'
ost_pr_genes <- bitr(ost_pr$UNIPROT, fromType = "UNIPROT", toType = "SYMBOL",
                     OrgDb = 'org.Hs.eg.db')
# not all ids have been mapped and some of them have multiple mapps so we remove gaps and multiple matches
ost_pr_genes$UNIPROT <- make.names(ost_pr_genes$UNIPROT, unique=TRUE)
# merge two data frames by ID
ost_pr <- merge(ost_pr, ost_pr_genes, by="UNIPROT", all.x = TRUE)

# if column with official gene symbols isn't available use UNIPROT ID
for (i in 1:nrow(ost_pr)) {
  if (is.na(ost_pr[i, "SYMBOL"])){
    ost_pr[i, "SYMBOL"] <- ost_pr[i, "UNIPROT"]
  }
}
# make official gene symbols rownames of data
ost_pr$SYMBOL <- make.names(ost_pr$SYMBOL, unique=TRUE)
rownames(ost_pr) <- ost_pr$SYMBOL
str(ost_pr)
ost_pr <- ost_pr[, c(-1, -14)]

# read legend
fact_ost_pr <- data.frame(read_excel("cocult_prot_legend.xlsx"))
# separate Osteoblast's legend
fact_ost_pr <- fact_ost_pr[fact_ost_pr$Cell.type == 'Osteoblast', ]
# change row names
rownames(fact_ost_pr) <- fact_ost_pr[,1]
# delete Cell.type and Sample columns
fact_ost_pr <- fact_ost_pr[, c(-1, -2)]
# make variables factors
fact_ost_pr$Condition <- as.factor(fact_ost_pr$Condition)
fact_ost_pr$Donor <- as.factor(as.character(fact_ost_pr$Donor))
str(fact_ost_pr)


setwd("./Osteoblast_proteome")


#Raw ost_pr
pal <- brewer.pal(n = 3, name = "Set1")
cols <- pal[fact_ost_pr$Condition]
boxplot(ost_pr, outline = FALSE, col = cols, main = "Raw data")
legend("topright", inset=.02, levels(fact_ost_pr$Condition), fill = pal, xpd = T, cex=0.8)
# write row data
#write.csv(ost_pr, 'ost_pr.csv')


#VSN normalization
ost_pr_norm <- normalizeQuantiles(ost_pr)
#write.csv(ost_pr_norm, 'ost_pr_norm.csv')
boxplot(ost_pr_norm, outline = FALSE, col = cols, main = "Normalized data")
legend("topright", levels(fact_ost_pr$Condition), fill = pal, xpd = T, cex=0.8)

#### Сlustering ####
#Sample Principal Component Plot
ost_pr_pca <- pca(t(ost_pr_norm), ncomp = 5, center = TRUE)

#tiff('PCA_Osteoblast_pr.tiff', units="in", width=10, height=8, res=300, compression = 'lzw')
plotIndiv(ost_pr_pca, comp = c(1, 2), ind.names = F, group = fact_ost_pr$Condition, 
          legend = TRUE, ellipse = T, title = 'PCA', style = "ggplot2", cex = 6, size.axis = 14, size.xlabel = 18, size.ylabel = 18, size.title = 22, size.legend = 18, size.legend.title = 18)
dev.off()
plotIndiv(ost_pr_pca, comp = c(1, 2), ind.names = T, group = fact_ost_pr$Donor, 
          legend = TRUE, ellipse = T, title = 'PCA')

#PLS-DA 
ordination.optimum.splsda <- splsda(t(ost_pr_norm), fact_ost_pr$Condition, ncomp = 3, keepX = c(15,15,15))
selectVar(ordination.optimum.splsda, comp=1)
selectVar(ordination.optimum.splsda, comp=2)
selectVar(ordination.optimum.splsda, comp=3)

#layout(matrix(c(1, 2), 2, 3))
plotLoadings(ordination.optimum.splsda, comp = 1, size.name = 1, size.title = 1.2, title = "Loadings\n on 1st component", contrib = "max", legend = FALSE, col.ties="black", ndisplay = 15)

plotLoadings(ordination.optimum.splsda, comp = 2, size.name = 1, size.title = 1.2, title = "Loadings\n on 2nd component", contrib = "max",ndisplay = 15,  legend = FALSE, col.ties="black")

#tiff('PLSDA_Osteoblast_pr.tiff', units="in", width=10, height=8, res=300, compression = 'lzw')
plotIndiv(ordination.optimum.splsda, ind.names = F, ellipse = T, style = "graphics", abline = TRUE, cex = 2, size.axis = 1.0, size.xlabel = 1.2, size.ylabel = 1.2, title = "PLS-DA ordination", size.title = 1.5, legend=TRUE)
dev.off()

#tiff('PLSDA_Osteoblast_prot.tiff', units="in", width=16, height=8, res=300, compression = 'lzw')
layout(matrix(c(1, 2, 3, 3, 3, 3), 2, 3))
plotLoadings(ordination.optimum.splsda, comp = 1, size.name = 1, size.title = 1.2, title = "Loadings\n on 1st component", contrib = "max", legend = FALSE, col.ties="black", ndisplay = 15)
plotLoadings(ordination.optimum.splsda, comp = 2, size.name = 1, size.title = 1.2, title = "Loadings\n on 2nd component", contrib = "max",ndisplay = 15,  legend = FALSE, col.ties="black")
plotIndiv(ordination.optimum.splsda, ind.names = F, ellipse = T, style = "graphics", abline = TRUE, cex = 1.5, size.axis = 1.2, size.xlabel = 1.5, size.ylabel = 1.5, title = "PLS-DA ordination", size.title = 1.5, legend=TRUE)
dev.off()

#tiff('PCA_PLSDA_Osteoblast_pr.tiff', units="in", width=16, height=8, res=300, compression = 'lzw')
layout(matrix(1:2, ncol = 2))
plotIndiv(ost_pr_pca, comp = c(1, 2), ind.names = F, group = fact_ost_pr$Condition, 
          style = "graphics", size.title = 1.4, ellipse = T, title = 'PCA', legend=TRUE)
plotIndiv(ordination.optimum.splsda, ind.names = F, ellipse = T, style = "graphics", size.title = 1.4,
          title = "PLS-DA ordination", legend=TRUE)
dev.off()


#### Differential expression analysis ####
design <- model.matrix(~ 0 + fact_ost_pr$Condition)
colnames(design) <- c("control","direct","indirect")
fit <- lmFit(ost_pr_norm, design)
contrast.matrix <- makeContrasts(direct-control, indirect-control, direct-indirect,
                                 levels=design)
ost_pr_fit <- contrasts.fit(fit, contrast.matrix)
ost_pr_fit <- eBayes(ost_pr_fit)
results <- decideTests(ost_pr_fit)

# Dif_expr_table
ost_pr_full_list <- topTable(ost_pr_fit, number=nrow(ost_pr_norm))
head(ost_pr_full_list, 10)

# counts of proteins for a single one across the groups
boxplot(as.numeric(ost_pr_norm['HTRA1',]) ~ Condition, data = fact_ost_pr, varwidth = TRUE, log = "y", las = 1, ylab = "HTRA1")
boxplot(as.numeric(ost_pr_norm["VASN",]) ~ Condition, data = fact_ost_pr, varwidth = TRUE, log = "y", las = 1, ylab = "VASN")
boxplot(as.numeric(ost_pr_norm["CCN1",]) ~ Condition, data = fact_ost_pr, varwidth = TRUE, log = "y", las = 1, ylab = "CCN1")
boxplot(as.numeric(ost_pr_norm["FBLN1",]) ~ Condition, data = fact_ost_pr, varwidth = TRUE, log = "y", las = 1, ylab = "FBLN1")


#### Control-Direct analysis ####
ost_contr_dir_pr <- topTable(ost_pr_fit, coef=1, adjust="BH", number=nrow(ost_pr_norm))
#write.csv(as.data.frame(ost_contr_dir_pr), file="Osteoblast_pr_control_vs_direct.csv")

head(ost_contr_dir_pr)
#Vulcano
#tiff('Osteoblast_pr_control_vs_direct_volcano.tiff', units="in", width=12, height=12, res=300, compression = 'lzw')
EnhancedVolcano(ost_contr_dir_pr,
                lab = rownames(ost_contr_dir_pr),
                x = 'logFC',
                y = 'adj.P.Val',
                pCutoff = 0.05,
                FCcutoff = 1,
                xlim = c(-3.2, 1.7),
                ylim = c(0, 5),
                title ="Control with direct cocultivation. Osteoblast",
                labSize = 7.0,
                pointSize = 3.5,
                axisLabSize = 25,
                titleLabSize = 30,
                subtitleLabSize = 1,
                captionLabSize = 20,
                boxedLabels = F,
                colAlpha = 1)
dev.off()

boxplot(as.numeric(ost_pr_norm["HTRA1",]) ~ Condition, data = fact_ost_pr, varwidth = TRUE, log = "y", las = 1, ylab = "HTRA1")

ost_contr_dir_pr_up <- ost_contr_dir_pr[ost_contr_dir_pr$logFC >= 1 & ost_contr_dir_pr$adj.P.Val <= 0.05,]
ost_contr_dir_pr_dn <- ost_contr_dir_pr[ost_contr_dir_pr$logFC <= -1 & ost_contr_dir_pr$adj.P.Val <= 0.05,]

#write.csv(as.data.frame(ost_contr_dir_pr_up), file="Osteoblast_pr_control_vs_direct_up.csv")
#write.csv(as.data.frame(ost_contr_dir_pr_dn), file="Osteoblast_pr_control_vs_direct_down.csv")


##### Control-Direct enrichment analysis ####
ost_cd_pr_up <- bitr(rownames(ost_contr_dir_pr_up), fromType = "SYMBOL", 
                     toType = "ENTREZID", OrgDb = 'org.Hs.eg.db')
ost_cd_pr_dn <- bitr(rownames(ost_contr_dir_pr_dn), fromType = "SYMBOL", 
                     toType = "ENTREZID", OrgDb = 'org.Hs.eg.db')
ost_cd_pr_bg <- bitr(rownames(ost_contr_dir_pr), fromType = "SYMBOL", 
                     toType = "ENTREZID", OrgDb = 'org.Hs.eg.db')

###### ost_contr_dir_prot_GO_BP ####
ost_contr_dir_pr_up_bp <- enrichGO(gene = ost_cd_pr_up$ENTREZID, universe = ost_cd_pr_bg$ENTREZID, OrgDb = 'org.Hs.eg.db', ont = "BP", pAdjustMethod = "fdr", pvalueCutoff = 0.05, readable = TRUE)
# no result


ost_contr_dir_pr_dn_bp <- enrichGO(gene = ost_cd_pr_dn$ENTREZID, universe = ost_cd_pr_bg$ENTREZID,
                                   OrgDb = 'org.Hs.eg.db', ont = "BP", pAdjustMethod = "fdr",
                                   pvalueCutoff = 0.05, readable = TRUE)
#tiff('Osteoblast_pr_control_vs_direct_down_GO_BP.tiff', units="in", width=10, height=7, res=300, compression = 'lzw')
dotplot(ost_contr_dir_pr_dn_bp, showCategory = 20, font.size = 18, label_format = 50)
dev.off()


###### ost_contr_dir_prot_GO_CC ####
ost_contr_dir_pr_up_cc <- enrichGO(gene = ost_cd_pr_up$ENTREZID, universe = ost_cd_pr_bg$ENTREZID, OrgDb = 'org.Hs.eg.db', ont = "CC", pAdjustMethod = "fdr", pvalueCutoff = 0.05, readable = TRUE)
# no result

ost_contr_dir_pr_dn_cc <- enrichGO(gene = ost_cd_pr_dn$ENTREZID, universe = ost_cd_pr_bg$ENTREZID, OrgDb = 'org.Hs.eg.db', ont = "CC", pAdjustMethod = "fdr", pvalueCutoff = 0.05, readable = TRUE)
#tiff('Osteoblast_pr_control_vs_direct_down_GO_CC.tiff', units="in", width=10, height=4, res=300, compression = 'lzw')
dotplot(ost_contr_dir_pr_dn_cc, showCategory = 20, font.size = 18, label_format = 50)
dev.off()


###### ost_contr_dir_prot_KEGG ####
ost_contr_dir_pr_up_kegg <- enrichKEGG(gene = ost_cd_pr_up$ENTREZID, universe = ost_cd_pr_bg$ENTREZID, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05, pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2, use_internal_data = FALSE)
# no result

ost_contr_dir_pr_dn_kegg <- enrichKEGG(gene = ost_cd_pr_dn$ENTREZID, universe = ost_cd_pr_bg$ENTREZID, organism = "hsa",
                                       keyType = "kegg", pvalueCutoff = 0.05, pAdjustMethod = "BH", 
                                       minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2,
                                       use_internal_data = FALSE)
#tiff('Osteoblast_pr_control_vs_direct_down_KEGG.tiff', units="in", width=10, height=3, res=300, compression = 'lzw')
dotplot(ost_contr_dir_pr_dn_kegg, showCategory = 20, font.size = 18, label_format = 50)
dev.off()


#### Control-Indirect analysis ####
ost_contr_indir_pr <- topTable(ost_pr_fit, coef=2, adjust="BH", number=nrow(ost_pr_norm))
#write.csv(as.data.frame(ost_contr_indir_pr), file="Osteoblast_pr_control_vs_indirect.csv")

#Vulcano
#tiff('Osteoblast_pr_control_vs_indirect_volcano.tiff', units="in", width=12, height=12, res=300, compression = 'lzw')
EnhancedVolcano(ost_contr_indir_pr,
                lab = rownames(ost_contr_indir_pr),
                x = 'logFC',
                y = 'P.Value',
                pCutoff = 0.05,
                FCcutoff = 1,
                xlim = c(-3, 3),
                ylim = c(0, 7.5),
                title ="Control with indirect cocultivation. Osteoblast",
                labSize = 7.0,
                pointSize = 3.5,
                axisLabSize = 25,
                titleLabSize = 30,
                subtitleLabSize = 1,
                captionLabSize = 20,                
                boxedLabels = F,                
                colAlpha = 1)
dev.off()

head(ost_contr_indir_pr)
# no DEG


#### Direct-Indirect analysis ####
ost_dir_indir_pr <- topTable(ost_pr_fit, coef=3, adjust="BH", number=nrow(ost_pr_norm))
#write.csv(as.data.frame(ost_dir_indir_pr), file="Osteoblast_pr_direct_vs_indirect.csv")

#Vulcano
#tiff('Osteoblast_pr_direct_vs_indirect_volcano.tiff', units="in", width=12, height=12, res=300, compression = 'lzw')
EnhancedVolcano(ost_dir_indir_pr,
                lab = rownames(ost_dir_indir_pr),
                x = 'logFC',
                y = 'adj.P.Val',
                pCutoff = 0.05,
                FCcutoff = 1,
                xlim = c(-3, 2),
                ylim = c(0, 4.5),
                title ="Direct with indirect cocultivation. Osteoblast",
                labSize = 7.0,
                pointSize = 3.5,
                axisLabSize = 25,
                titleLabSize = 30,
                subtitleLabSize = 1,
                captionLabSize = 20,                
                boxedLabels = F,                
                colAlpha = 1)
dev.off()

head(ost_dir_indir_pr)

boxplot(as.numeric(ost_pr_norm["HTRA1",]) ~ Condition, data = fact_ost_pr, varwidth = TRUE, log = "y", las = 1, ylab = "HTRA1")


ost_dir_indir_pr_up <- ost_dir_indir_pr[ost_dir_indir_pr$logFC >= 1 & ost_dir_indir_pr$adj.P.Val <= 0.05,]
ost_dir_indir_pr_dn <- ost_dir_indir_pr[ost_dir_indir_pr$logFC <= -1 & ost_dir_indir_pr$adj.P.Val <= 0.05,]

#write.csv(as.data.frame(ost_dir_indir_pr_up), file="Osteoblast_pr_direct_vs_indirect_up.csv")
#write.csv(as.data.frame(ost_dir_indir_pr_dn), file="Osteoblast_pr_direct_vs_indirect_down.csv")


##### Direct-Indirect enrichment analysis ####
ost_dind_pr_up <- bitr(rownames(ost_dir_indir_pr_up), fromType = "SYMBOL", toType = "ENTREZID",
                       OrgDb = 'org.Hs.eg.db')
ost_dind_pr_dn <- bitr(rownames(ost_dir_indir_pr_dn), fromType = "SYMBOL", toType = "ENTREZID",
                       OrgDb = 'org.Hs.eg.db')
ost_dind_pr_bg <- bitr(rownames(ost_dir_indir_pr), fromType = "SYMBOL", toType = "ENTREZID",
                       OrgDb = 'org.Hs.eg.db')


###### ost_dir_indir_prot_GO_BP ####
ost_dir_indir_pr_up_bp <- enrichGO(gene = ost_dind_pr_up$ENTREZID, universe = ost_dind_pr_bg$ENTREZID, OrgDb = 'org.Hs.eg.db', ont = "BP", pAdjustMethod = "fdr",
                                   pvalueCutoff = 0.05, readable = TRUE)
# no result

ost_dir_indir_pr_dn_bp <- enrichGO(gene = ost_dind_pr_dn$ENTREZID, universe = ost_dind_pr_bg$ENTREZID, OrgDb = 'org.Hs.eg.db', ont = "BP", pAdjustMethod = "fdr", pvalueCutoff = 0.05, readable = TRUE)
#tiff('Osteoblast_pr_direct_vs_indirect_down_GO_BP.tiff', units="in", width=11, height=9, res=300, compression = 'lzw')
dotplot(ost_dir_indir_pr_dn_bp, showCategory = 20, font.size = 18, label_format = 50)
dev.off()


###### ost_dir_indir_prot_GO_CC ####
ost_dir_indir_pr_up_cc <- enrichGO(gene = ost_dind_pr_up$ENTREZID, universe = ost_dind_pr_bg$ENTREZID, OrgDb = 'org.Hs.eg.db', ont = "CC", pAdjustMethod = "fdr",
                                   pvalueCutoff = 0.05, readable = TRUE)
# no result

ost_dir_indir_pr_dn_cc <- enrichGO(gene = ost_dind_pr_dn$ENTREZID, universe = ost_dind_pr_bg$ENTREZID,
                                   OrgDb = 'org.Hs.eg.db', ont = "CC", pAdjustMethod = "fdr",
                                   pvalueCutoff = 0.05, readable = TRUE)
#tiff('Osteoblast_pr_direct_vs_indirect_down_GO_CC.tiff', units="in", width=10, height=6, res=300, compression = 'lzw')
dotplot(ost_dir_indir_pr_dn_cc, showCategory = 15)
dev.off()


###### ost_dir_indir_prot_KEGG ####
ost_dir_indir_pr_up_kegg <- enrichKEGG(gene = ost_dind_pr_up$ENTREZID, universe = ost_dind_pr_bg$ENTREZID, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05, pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2, use_internal_data = FALSE)
# no result

ost_dir_indir_pr_dn_kegg <- enrichKEGG(gene = ost_dind_pr_dn$ENTREZID, universe = ost_dind_pr_bg$ENTREZID, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05, pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2,
                                       use_internal_data = FALSE)
# no result


# Osteoblast transcriptome vs proteome analysis ####
setwd("../Osteoblast_transcriptome_vs_proteome")


## Venn diagram transcriptome vs proteome all genes ####
venn_col <- c("#0073C2FF", "#CD534CFF")
venn_ost_pr_tr_all <- list(transcriptome = rownames(ost_tr), 
                           proteome = rownames(ost_pr_norm))
#tiff('venn_ost_pr_tr_all.tiff', units="in", width=8, height=5, res=300, compression = 'lzw')
ggvenn(venn_ost_pr_tr_all, fill_color = venn_col, stroke_size = 1, set_name_size = 10, text_size = 10)
dev.off()


## Venn diagram, Correlation and Hetmap ####

### Transcriptome vs proteome direct-control ####

#### Venn diagram up direct-control differentiated genes ####
venn_ost_pr_tr_cdir_up <- list(RNA_up = rownames(ost_contr_dir_tr_up), 
                               Prot.up = rownames(ost_contr_dir_pr_up))
#tiff('venn_ost_pr_tr_cdir_up.tiff', units="in", width=8, height=5, res=300, compression = 'lzw')
ggvenn(venn_ost_pr_tr_cdir_up, fill_color = venn_col, stroke_size = 1, set_name_size = 10, text_size = 10)
dev.off()


#### Venn diagram down direct-control differentiated genes ####
venn_ost_pr_tr_cdir_down <- list(RNA_down = rownames(ost_contr_dir_tr_dn), 
                                 Prot.down = rownames(ost_contr_dir_pr_dn))
#tiff('venn_ost_pr_tr_cdir_down.tiff', units="in", width=8, height=5, res=300, compression = 'lzw')
ggvenn(venn_ost_pr_tr_cdir_down, fill_color = venn_col, stroke_size = 1, set_name_size = 10, text_size = 10)
dev.off()


#### Correlation of log2 Fold Change between direct-control differentiated genes ####
ost_contr_dir_tr_pr_gene <- intersect(rownames(ost_contr_dir_tr), 
                                      rownames(ost_contr_dir_pr))
ost_contr_dir_tr_FC <- ost_contr_dir_tr[ost_contr_dir_tr_pr_gene, 'log2FoldChange']
ost_contr_dir_pr_FC <- ost_contr_dir_pr[ost_contr_dir_tr_pr_gene, 'logFC']
ost_contr_dir_tr_pr <- data.frame(ost_contr_dir_tr_pr_gene,
                                  ost_contr_dir_tr_FC,
                                  ost_contr_dir_pr_FC)
colnames(ost_contr_dir_tr_pr) <- c('Gene_name', 'Log2FC_tr', 'Log2FC_pr')

#tiff('Osteoblasts_pr_tr_contr_vs_direct.tiff', units="in", width=12, height=6, res=300, compression = 'lzw')
ggplot(ost_contr_dir_tr_pr, aes(x=Log2FC_tr, y=Log2FC_pr)) +
  geom_point() + 
  geom_vline(xintercept = 1, linetype = 2) + geom_vline(xintercept = -1, linetype = 2) + 
  geom_hline(yintercept = 1, linetype = 2) + geom_hline(yintercept = -1, linetype = 2) + 
  #geom_smooth(method="loess")+
  xlab("RNA log2 Fold Change") +
  ylab("protein log2 Fold Change") +
  ggtitle("correlation of log2 Fold Change between direct and control condition in Osteoblasts for proteomics and RNA-seq data") + 
  stat_cor(method = "pearson") +
  geom_text(aes(label=ifelse(Log2FC_tr>1 & Log2FC_pr >1,as.character(Gene_name),'')),hjust=0,vjust=-0.5) + 
  geom_text(aes(label=ifelse(Log2FC_pr<(-1) & Log2FC_pr<(-1),as.character(Gene_name),'')),hjust=0,vjust=-0.5) +
  annotate("text", 
           x = 8.9, # should be optimized manually
           y = 1.6,  # should be optimized manually
           label = paste("Number of genes: ", nrow(ost_contr_dir_tr_pr)))
dev.off()


# Heatmap
color <- list(sample = c(control = '#FF7370', 
                         direct = '#01AFDF', 
                         indirect = '#00E0A5'))
tr_hm <- data.frame(sample = fact_ost_tr$Condition)
row.names(tr_hm) <- rownames(fact_ost_tr)

pr_hm <- data.frame(sample = fact_ost_pr$Condition)
row.names(pr_hm) <- rownames(fact_ost_pr)


#### Heatmap direct-control differentiated genes ####

# Find common transcriptome and proteome control_direct DEG
ost_contr_dir_tr_dif <- ost_contr_dir_tr[abs(ost_contr_dir_tr$log2FoldChange) >= 1 & ost_contr_dir_tr$padj <= 0.05,]
ost_contr_dir_pr_dif <- ost_contr_dir_pr[abs(ost_contr_dir_pr$logFC) >= 1 & ost_contr_dir_pr$adj.P.Val <= 0.05,]
ost_contr_dir_tr_pr_dif <- intersect(rownames(ost_contr_dir_tr_dif), rownames(ost_contr_dir_pr_dif))


# Transcriptome heatmap direct-control differentiated genes

# Genes of interest, top common tr-pr DEGs and top DEGs
ost_contr_indir_tr_dif <- ost_contr_indir_tr[abs(ost_contr_indir_tr$log2FoldChange) >= 1 & ost_contr_indir_tr$padj <= 0.05,]
ost_contr_dir_trrr_dif <- intersect(rownames(ost_contr_dir_tr_dif), rownames(ost_contr_indir_tr_dif))

select <- c('ADAM19', 'DLL4', 'TGFB1', 'LIF', 'ADAMTS14', 'JAG1', 'JAG2', 'HEYL', 'VASN', 'ALCAM', 'NOTUM', ost_contr_dir_trrr_dif[1:8],
            rownames(ost_contr_dir_tr_dif[1:30,]))
select <- unique(select)


#tiff('Osteoblast_tr_direct_vs_control_heatmap.tiff', units="in", width=5, height=6, res=300, compression = 'lzw')
pheatmap(ost_tr[select[1:30],], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=tr_hm, 
         annotation_colors = color, show_colnames = FALSE,
         treeheight_col = 20, treeheight_row = 20)
dev.off()


# Proteome heatmap direct-control differentiated genes

# Usw only proteome control-direct DEG
select <- rownames(ost_contr_dir_pr_dif)

#tiff('Osteoblast_pr_direct_vs_control_heatmap.tiff', units="in", width=6, height=4.5, res=300, compression = 'lzw')
pheatmap(ost_pr_norm[select,], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=pr_hm,
         annotation_colors = color, show_colnames = FALSE,
         treeheight_col = 20, treeheight_row = 20)
dev.off()


### Transcriptome vs proteome indirect-control differentiated genes ####

#we don't have DEG for proteome, use only transcriptome

venn_ost_tr_cindir_up_down <- list(RNA_up = rownames(ost_contr_indir_tr_up), 
                                   RNA_down = rownames(ost_contr_indir_tr_dn))
#tiff('venn_ost_tr_cindir_up_down.tiff', units="in", width=8, height=5, res=300, compression = 'lzw')
ggvenn(venn_ost_tr_cindir_up_down, fill_color = venn_col, stroke_size = 1, set_name_size = 10, text_size = 10)
dev.off()


#### Correlation of log2 Fold Change between indirect-control differentiated genes ####
ost_contr_indir_tr_pr_gene <- intersect(rownames(ost_contr_indir_tr), 
                                      rownames(ost_contr_indir_pr))
ost_contr_indir_tr_FC <- ost_contr_indir_tr[ost_contr_indir_tr_pr_gene, 'log2FoldChange']
ost_contr_indir_pr_FC <- ost_contr_indir_pr[ost_contr_indir_tr_pr_gene, 'logFC']
ost_contr_indir_tr_pr <- data.frame(ost_contr_indir_tr_pr_gene,
                                  ost_contr_indir_tr_FC,
                                  ost_contr_indir_pr_FC)
colnames(ost_contr_indir_tr_pr) <- c('Gene_name', 'Log2FC_tr', 'Log2FC_pr')

#tiff('Osteoblasts_pr_tr_contr_vs_indirect.tiff', units="in", width=12, height=6, res=300, compression = 'lzw')
ggplot(ost_contr_indir_tr_pr, aes(x=Log2FC_tr, y=Log2FC_pr)) +
  geom_point() + 
  geom_vline(xintercept = 1, linetype = 2) + geom_vline(xintercept = -1, linetype = 2) + 
  geom_hline(yintercept = 1, linetype = 2) + geom_hline(yintercept = -1, linetype = 2) + 
  #geom_smooth(method="loess")+
  xlab("RNA log2 Fold Change") +
  ylab("protein log2 Fold Change") +
  ggtitle("correlation of log2 Fold Change between indirect and control condition in Osteoblasts for proteomics and RNA-seq data") + 
  stat_cor(method = "pearson") +
  geom_text(aes(label=ifelse(Log2FC_tr>1 & Log2FC_pr >1,as.character(Gene_name),'')),hjust=0,vjust=-0.5) + 
  geom_text(aes(label=ifelse(Log2FC_pr<(-1) & Log2FC_pr<(-1),as.character(Gene_name),'')),hjust=0,vjust=-0.5) +
  annotate("text", 
           x = 1.0, # should be optimized manually
           y = 1.05,  # should be optimized manually
           label = paste("Number of genes: ", nrow(ost_contr_indir_tr_pr)))
dev.off()


#### Heatmap transcriptome and proteome control-indirect differentiated genes ####

# Find common transcriptome and proteome control-indirect DEG
ost_contr_indir_tr_dif <- ost_contr_indir_tr[abs(ost_contr_indir_tr$log2FoldChange) >= 1 & ost_contr_indir_tr$padj <= 0.05,]


#### Transcriptome heatmap control-indirect differentiated genes ####
# Genes top DEGs
select <- rownames(ost_contr_indir_tr_dif)


#tiff('Osteoblast_tr_indirect_vs_control_heatmap.tiff', units="in", width=5, height=6, res=300, compression = 'lzw')
pheatmap(ost_tr[select[1:30],], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=tr_hm,
         annotation_colors = color, show_colnames = FALSE,
         treeheight_col = 20, treeheight_row = 20)
dev.off()


### Transcriptome vs proteome direct-indirect differentiated genes ####

#### Venn diagram up direct-indirect differentiated genes ####
venn_ost_pr_tr_dindir_up <- list(RNA_up = rownames(ost_dir_indir_tr_up), 
                                 Prot.up = rownames(ost_dir_indir_pr_up))
#tiff('venn_ost_pr_tr_dindir_up.tiff', units="in", width=8, height=5, res=300, compression = 'lzw')
ggvenn(venn_ost_pr_tr_dindir_up, fill_color = venn_col, stroke_size = 1, set_name_size = 10, text_size = 10)
dev.off()

venn_ost_pr_tr_dindir_dn <- list(RNA_down = rownames(ost_dir_indir_tr_dn), 
                                 Prot.down = rownames(ost_dir_indir_pr_dn))
#tiff('venn_ost_pr_tr_dindir_down.tiff', units="in", width=8, height=5, res=300, compression = 'lzw')
ggvenn(venn_ost_pr_tr_dindir_dn, fill_color = venn_col, stroke_size = 1, set_name_size = 10, text_size = 10)
dev.off()

# Correlation of log2 Fold Change between indirect and direct
ost_dir_indir_tr_pr_gene <- intersect(rownames(ost_dir_indir_tr), 
                                        rownames(ost_dir_indir_pr))
ost_dir_indir_tr_FC <- ost_dir_indir_tr[ost_dir_indir_tr_pr_gene, 'log2FoldChange']
ost_dir_indir_pr_FC <- ost_dir_indir_pr[ost_dir_indir_tr_pr_gene, 'logFC']
ost_dir_indir_tr_pr <- data.frame(ost_dir_indir_tr_pr_gene,
                                    ost_dir_indir_tr_FC,
                                    ost_dir_indir_pr_FC)
colnames(ost_dir_indir_tr_pr) <- c('Gene_name', 'Log2FC_tr', 'Log2FC_pr')

tiff('Osteoblasts_pr_tr_dir_vs_indirect.tiff', units="in", width=12, height=6, res=300, compression = 'lzw')
ggplot(ost_dir_indir_tr_pr, aes(x=Log2FC_tr, y=Log2FC_pr)) +
  geom_point() + 
  geom_vline(xintercept = 1, linetype = 2) + geom_vline(xintercept = -1, linetype = 2) + 
  geom_hline(yintercept = 1, linetype = 2) + geom_hline(yintercept = -1, linetype = 2) + 
  #geom_smooth(method="loess")+
  xlab("RNA log2 Fold Change") +
  ylab("protein log2 Fold Change") +
  ggtitle("correlation of log2 Fold Change between indirect and direct condition in Osteoblasts for proteomics and RNA-seq data") + 
  stat_cor(method = "pearson") +
  geom_text(aes(label=ifelse(Log2FC_tr>1 & Log2FC_pr >1,as.character(Gene_name),'')),hjust=0,vjust=-0.5) + 
  geom_text(aes(label=ifelse(Log2FC_pr<(-1) & Log2FC_pr<(-1),as.character(Gene_name),'')),hjust=0,vjust=-0.5) +
  annotate("text", 
           x = 9.0, # should be optimized manually
           y = 1.7,  # should be optimized manually
           label = paste("Number of genes: ", nrow(ost_dir_indir_tr_pr)))
dev.off()


#### Heatmap direct-indirect differentiated genes ####
# Find common transcriptome and proteome indirect-direct DEG
ost_dir_indir_tr_dif <- ost_dir_indir_tr[abs(ost_dir_indir_tr$log2FoldChange) >= 1 & ost_dir_indir_tr$padj <= 0.05,]
ost_dir_indir_pr_dif <- ost_dir_indir_pr[abs(ost_dir_indir_pr$logFC) >= 1 & ost_dir_indir_pr$adj.P.Val <= 0.05,]
ost_dir_indir_tr_pr_dif <- intersect(rownames(ost_dir_indir_tr_dif), rownames(ost_dir_indir_pr_dif))

ost_dir_indir_trrr_dif <- intersect(rownames(ost_dir_indir_tr_dif), rownames(ost_contr_indir_tr_dif))

# Transcriptome heatmap direct-indirect differentiated genes ####
# Genes of interest, top common tr-pr DEGs and top DEGs
select <- c(ost_dir_indir_tr_pr_dif, ost_dir_indir_trrr_dif[1:5], 
            rownames(ost_dir_indir_tr_dif[1:30,]))
select <- unique(select)

#tiff('Osteoblast_tr_direct_vs_indirect_heatmap.tiff', units="in", width=5, height=6, res=300, compression = 'lzw')
pheatmap(ost_tr[select[1:30],], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=tr_hm,
         annotation_colors = color, show_colnames = FALSE,
         treeheight_col = 20, treeheight_row = 20)
dev.off()


# Proteome heatmap direct-indirect differentiated genes ####
select <- rownames(ost_dir_indir_pr_dif)

#tiff('Osteoblast_pr_indirect_vs_control_heatmap.tiff', units="in", width=6, height=6, res=300, compression = 'lzw')
pheatmap(ost_pr_norm[select,], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=pr_hm,
         annotation_colors = color, show_colnames = FALSE,
         treeheight_col = 20, treeheight_row = 20)
dev.off()
