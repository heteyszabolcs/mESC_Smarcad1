
data <- read.table("../results/rna_seq_deseq2/TEtranscripts/ERVK_IAPez_extended/TEtranscripts_SMARCAD1_WTKO_IAPEz_extended.cntTable",header=T,row.names=1)
groups <- factor(c(rep("TGroup",3),rep("CGroup",3)))
min_read <- 1
data <- data[apply(data,1,function(x){max(x)}) > min_read,]
sampleInfo <- data.frame(groups,row.names=colnames(data))
suppressPackageStartupMessages(library(DESeq2))
dds <- DESeqDataSetFromMatrix(countData = data, colData = sampleInfo, design = ~ groups)
dds$groups = relevel(dds$groups,ref="CGroup")
dds <- DESeq(dds)

# nornalized values
norm_counts = counts(dds, normalized = TRUE)
write.table(norm_counts, file = "../results/rna_seq_deseq2/TEtranscripts/ERVK_IAPez_extended/TEtranscripts_SMARCAD1_KOWT_IAPEz_normalized_expr.tsv", quote = FALSE, sep = "\t")

res <- results(dds)
write.table(res, file="TEtranscripts_SMARCAD1_WTKO_IAPEz_extended_gene_TE_analysis.txt", sep="\t",quote=F)
resSig <- res[(!is.na(res$padj) & (res$padj < 0.050000) &         (abs(res$log2FoldChange)> 0.000000)), ]
write.table(resSig, file="TEtranscripts_SMARCAD1_WTKO_IAPEz_extended_sigdiff_gene_TE.txt",sep="\t", quote=F)
