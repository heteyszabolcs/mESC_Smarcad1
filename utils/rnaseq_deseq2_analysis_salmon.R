# packages
suppressPackageStartupMessages({
  library("DESeq2")
  library("data.table")
  library("tidyverse")
  library("EnhancedVolcano")
  library("glue")
  library("enrichR")
  library("scales")
  library("patchwork")
})

# folders
result_folder = "../results/rna_seq_deseq2/"

# prepare tables
counts = fread("../data/RNA-Seq/nextflow_output/salmon/salmon.merged.gene_counts_scaled.tsv")
counts = counts %>% group_by(gene_name) %>% summarise_all(mean) %>% dplyr::select(-gene_id)
counts = counts %>%
  mutate_if(is.numeric, round)
counts = as.data.frame(counts)
rownames(counts) = counts$gene_name
counts = counts[-1]

coldata = fread("../data/coldata.txt", sep = "\t", header = TRUE, stringsAsFactors = TRUE)
coldata = as.data.frame(coldata)
rownames(coldata) = coldata$sample
all(rownames(coldata) == colnames(counts))

## DESeq2 protocol 
dds = DESeqDataSetFromMatrix(countData = counts,
                             colData = coldata,
                             design = ~ treatment)
# keep = rowSums(counts(dds)) >= 10
# dds = dds[keep, ]
dds$condition = factor(dds$treatment, levels = c("WT", "treatment"))
dds$condition = relevel(dds$treatment, ref = "WT")
dds = DESeq(dds)

# normalized read counts (median of ratios)
norm_counts = counts(dds, normalized = TRUE)
wts = norm_counts[,c("Ph001", "Ph006", "Ph011")]
highly_expressed_genes = tibble(gene = rownames(wts)[order(rowMeans(wts), decreasing = TRUE)][1:200])
write_tsv(highly_expressed_genes, "../data/top200_highly_expr_gene-wt.txt", col_names = FALSE)

# result table
# SMARCAD1 KO vs. WT
ko_res = results(dds, contrast = c("treatment", "SMARCAD1_KO", "WT"))
ko_res = as.data.frame(ko_res)
ko_res["gene_name"] = rownames(ko_res)
write_tsv(ko_res, glue("{result_folder}DESeq2_res-ko_vs_wt_salmon.tsv"))
# SMARCAD1 KO WT vs. WT
kowt_res = results(dds, contrast = c("treatment", "SMARCAD1_KO_WT", "WT"))
kowt_res = as.data.frame(kowt_res)
kowt_res["gene_name"] = rownames(kowt_res)
write_tsv(kowt_res, glue("{result_folder}DESeq2_res-kowt_vs_wt_salmon.tsv"))
# SMARCAD1 KO d425 vs. WT
kod425_res = results(dds, contrast = c("treatment", "SMARCAD1_KO_d425", "WT"))
kod425_res = as.data.frame(kod425_res)
kod425_res["gene_name"] = rownames(kod425_res)
write_tsv(kod425_res, glue("{result_folder}DESeq2_res-kod425_vs_wt_salmon.tsv"))
# SMARCAD1 KO K523R vs. WT
kok523r_res = results(dds, contrast = c("treatment", "SMARCAD1_KO_K523R", "WT"))
kok523r_res = as.data.frame(kok523r_res)
kok523r_res["gene_name"] = rownames(kok523r_res)
write_tsv(kok523r_res, glue("{result_folder}DESeq2_res-kok523r_vs_wt_salmon.tsv"))

ko_volc = EnhancedVolcano(
  ko_res,
  lab = ko_res$gene_name,
  labSize = 3.5,
  title = "SMARCAD1 KO - WT",
  pCutoff = 0.05,
  FCcutoff = 1,
  colAlpha = 0.5,
  legendLabels = c(" ", " ", " ", " "),
  legendLabSize = 14,
  legendIconSize = 0,
  legendDropLevels = TRUE,
  legendPosition = "right",
  x = 'log2FoldChange',
  y = 'pvalue',
  subtitle = "pCutoff: 0.05, FC: 1",
  col = c('grey', 'grey', 'grey', '#2ca25f')
)
ko_volc

kowt_volc = EnhancedVolcano(
  kowt_res,
  lab = kowt_res$gene_name,
  labSize = 3.5,
  title = "SMARCAD1 KO WT - WT",
  pCutoff = 0.05,
  FCcutoff = 1,
  colAlpha = 0.5,
  legendLabels = c(" ", " ", " ", " "),
  legendLabSize = 14,
  legendIconSize = 0,
  legendDropLevels = TRUE,
  legendPosition = "right",
  x = 'log2FoldChange',
  y = 'pvalue',
  subtitle = "pCutoff: 0.05, FC: 1",
  col = c('grey', 'grey', 'grey', '#2ca25f')
)
kowt_volc

kok523r_volc = EnhancedVolcano(
  kok523r_res,
  lab = kok523r_res$gene_name,
  labSize = 3.5,
  title = "SMARCAD1 KO K523R - WT",
  pCutoff = 0.05,
  FCcutoff = 1,
  colAlpha = 0.5,
  legendLabels = c(" ", " ", " ", " "),
  legendLabSize = 14,
  legendIconSize = 0,
  legendDropLevels = TRUE,
  legendPosition = "right",
  x = 'log2FoldChange',
  y = 'pvalue',
  subtitle = "pCutoff: 0.05, FC: 1",
  col = c('grey', 'grey', 'grey', '#2ca25f')
)
kok523r_volc

kod425_volc = EnhancedVolcano(
  kod425_res,
  lab = kod425_res$gene_name,
  labSize = 3.5,
  title = "SMARCAD1 KO d425 - WT",
  pCutoff = 0.05,
  FCcutoff = 1,
  colAlpha = 0.5,
  legendLabels = c(" ", " ", " ", " "),
  legendLabSize = 14,
  legendIconSize = 0,
  legendDropLevels = TRUE,
  legendPosition = "right",
  x = 'log2FoldChange',
  y = 'pvalue',
  subtitle = "pCutoff: 0.05, FC: 1",
  col = c('grey', 'grey', 'grey', '#2ca25f')
)
kod425_volc

# volcano plot
png(
  file = glue("{result_folder}RNA-Seq_volcanos_salmon.png"),
  width = 10,
  height = 10,
  units = 'in',
  res = 500
)
ko_volc + kowt_volc + kok523r_volc + kod425_volc
dev.off()

# volcano plot
pdf(
  file = glue("{result_folder}RNA-Seq_volcanos_salmon.pdf"),
  width = 10,
  height = 10
)
ko_volc + kowt_volc + kok523r_volc + kod425_volc
dev.off()



