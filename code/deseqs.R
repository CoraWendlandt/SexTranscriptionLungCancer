CPTAC3_n_cts <- read.csv('data/CPTAC3_NORMAL_deseqs_counts.csv',row.names=1)
CPTAC3_n_meta <- read.csv('data/CPTAC3_NORMAL_meta.csv',row.names = 1)
install.packages("BiocManager")
library(BiocManager)
BiocManager::install("DESeq2", force= TRUE)
library("DESeq2")
CPTAC3_n_dds <- DESeq2::DESeqDataSetFromMatrix(countData = CPTAC3_n_cts,
                              colData = CPTAC3_n_meta,
                              design = ~ gender)
CPTAC3_n_dds <- DESeq(CPTAC3_n_dds)
CPTAC3_n_res <- results(CPTAC3_n_dds)
CPTAC3_n_res
CPTAC3_n_res <- CPTAC3_n_res[order(CPTAC3_n_res$log2FoldChange,decreasing=TRUE),]
BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)
CPTAC3_n_res_noNA <- na.omit(CPTAC3_n_res)
C_n_selectLab <- CPTAC3_n_res_noNA[abs(CPTAC3_n_res_noNA$log2FoldChange) > 2 & CPTAC3_n_res_noNA$padj < 0.05,]@rownames
C_n_selectLab

CPTAC3_t_cts <- read.csv('data/CPTAC3_TUMOR_deseqs_counts.csv',row.names=1)
CPTAC3_t_meta <- read.csv('data/CPTAC3_TUMOR_meta.csv',row.names = 1)
CPTAC3_t_dds <- DESeq2::DESeqDataSetFromMatrix(countData = CPTAC3_t_cts,
                                               colData = CPTAC3_t_meta,
                                               design = ~ gender)
CPTAC3_t_dds <- DESeq(CPTAC3_t_dds)
CPTAC3_t_res <- results(CPTAC3_t_dds)
CPTAC3_t_res
CPTAC3_t_res <- CPTAC3_t_res[order(CPTAC3_t_res$log2FoldChange,decreasing=TRUE),]
CPTAC3_t_res_noNA <- na.omit(CPTAC3_t_res)
C_t_selectLab <- CPTAC3_t_res_noNA[abs(CPTAC3_t_res_noNA$log2FoldChange) > 2 & CPTAC3_t_res_noNA$padj < 0.05,]@rownames
C_t_selectLab

C_t_only = C_t_selectLab[(C_t_selectLab) %in% C_n_selectLab == FALSE]



LUAD_n_cts <- read.csv('data/LUAD_NORMAL_deseqs_counts.csv',row.names=1)
LUAD_n_meta <- read.csv('data/LUAD_NORMAL_meta.csv',row.names = 1)
LUAD_n_dds <- DESeq2::DESeqDataSetFromMatrix(countData = LUAD_n_cts,
                                               colData = LUAD_n_meta,
                                               design = ~ gender)
LUAD_n_dds <- DESeq(LUAD_n_dds)
LUAD_n_res <- results(LUAD_n_dds)
LUAD_n_res
LUAD_n_res <- LUAD_n_res[order(LUAD_n_res$log2FoldChange,decreasing=TRUE),]

LUAD_n_res_noNA <- na.omit(LUAD_n_res)
LUAD_n_selectLab <- LUAD_n_res_noNA[abs(LUAD_n_res_noNA$log2FoldChange) > 2 & LUAD_n_res_noNA$padj < 0.05,]@rownames
LUAD_n_selectLab

LUAD_t_cts <- read.csv('data/LUAD_TUMOR_deseqs_counts.csv',row.names=1)
LUAD_t_meta <- read.csv('data/LUAD_TUMOR_meta.csv',row.names = 1)
LUAD_t_dds <- DESeq2::DESeqDataSetFromMatrix(countData = LUAD_t_cts,
                                               colData = LUAD_t_meta,
                                               design = ~ gender)
LUAD_t_dds <- DESeq(LUAD_t_dds)
LUAD_t_res <- results(LUAD_t_dds)
LUAD_t_res
LUAD_t_res <- LUAD_t_res[order(LUAD_t_res$log2FoldChange,decreasing=TRUE),]
LUAD_t_res_noNA <- na.omit(LUAD_t_res)
LUAD_t_selectLab <- LUAD_t_res_noNA[abs(LUAD_t_res_noNA$log2FoldChange) > 2 & LUAD_t_res_noNA$padj < 0.05,]@rownames
LUAD_t_selectLab

LUAD_t_only = LUAD_t_selectLab[(LUAD_t_selectLab) %in% LUAD_n_selectLab == FALSE]




LUSC_n_cts <- read.csv('data/LUSC_NORMAL_deseqs_counts.csv',row.names=1)
LUSC_n_meta <- read.csv('data/LUSC_NORMAL_meta.csv',row.names = 1)
LUSC_n_dds <- DESeq2::DESeqDataSetFromMatrix(countData = LUSC_n_cts,
                                             colData = LUSC_n_meta,
                                             design = ~ gender)
LUSC_n_dds <- DESeq(LUSC_n_dds)
LUSC_n_res <- results(LUSC_n_dds)
LUSC_n_res
LUSC_n_res <- LUSC_n_res[order(LUSC_n_res$log2FoldChange,decreasing=TRUE),]

LUSC_n_res_noNA <- na.omit(LUSC_n_res)
LUSC_n_selectLab <- LUSC_n_res_noNA[abs(LUSC_n_res_noNA$log2FoldChange) > 2 & LUSC_n_res_noNA$padj < 0.05,]@rownames
LUSC_n_selectLab

LUSC_t_cts <- read.csv('data/LUSC_TUMOR_deseqs_counts.csv',row.names=1)
LUSC_t_meta <- read.csv('data/LUSC_TUMOR_meta.csv',row.names = 1)
LUSC_t_dds <- DESeq2::DESeqDataSetFromMatrix(countData = LUSC_t_cts,
                                             colData = LUSC_t_meta,
                                             design = ~ gender)
LUSC_t_dds <- DESeq(LUSC_t_dds)
LUSC_t_res <- results(LUSC_t_dds)
LUSC_t_res
LUSC_t_res <- LUSC_t_res[order(LUSC_t_res$log2FoldChange,decreasing=TRUE),]
LUSC_t_res_noNA <- na.omit(LUSC_t_res)
LUSC_t_selectLab <- LUSC_t_res_noNA[abs(LUSC_t_res_noNA$log2FoldChange) > 2 & LUSC_t_res_noNA$padj < 0.05,]@rownames
LUSC_t_selectLab

LUSC_t_only = LUSC_t_selectLab[(LUSC_t_selectLab) %in% LUSC_n_selectLab == FALSE]


t_only = intersect(intersect(C_t_only,LUAD_t_only),LUSC_t_only)

EnhancedVolcano(LUSC_t_res_noNA,
                lab = rownames(LUSC_t_res_noNA),
                x = 'log2FoldChange',
                y = 'padj',
                selectLab = t_only,
                xlab = bquote(~Log[2]~ 'fold change'),
                title = 'LUSC Male vs Female Tumors',
                pCutoff = 0.05,
                FCcutoff = 2.0,
                pointSize = 4.0,
                labSize = 6.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                colAlpha = 4/5,
                legendPosition = 'right',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black')

EnhancedVolcano(LUAD_t_res_noNA,
                lab = rownames(LUAD_t_res_noNA),
                x = 'log2FoldChange',
                y = 'padj',
                selectLab = t_only,
                xlab = bquote(~Log[2]~ 'fold change'),
                title = 'LUAD Male vs Female Tumors',
                pCutoff = 0.05,
                FCcutoff = 2.0,
                pointSize = 4.0,
                labSize = 6.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                colAlpha = 4/5,
                legendPosition = 'right',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black')

EnhancedVolcano(CPTAC3_t_res_noNA,
                lab = rownames(CPTAC3_t_res_noNA),
                x = 'log2FoldChange',
                y = 'padj',
                selectLab = t_only,
                xlab = bquote(~Log[2]~ 'fold change'),
                title = 'CPTAC3 Male vs Female Tumors',
                pCutoff = 0.05,
                FCcutoff = 2.0,
                pointSize = 4.0,
                labSize = 6.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                colAlpha = 4/5,
                legendPosition = 'right',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black')
