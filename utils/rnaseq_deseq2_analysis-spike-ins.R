suppressPackageStartupMessages({
  library("tidyverse")
  library("data.table")
  library("DESeq2")
  library("pheatmap")
  library("EnhancedVolcano")
  library("glue")
  library("patchwork")
})

# folders
result_folder = "../results/rna_seq_deseq2/"

# prepare tables
counts = fread("../data/RNA-Seq/nextflow_output/star_rsem_spikein/rsem.merged.transcript_counts_spikein.tsv")
coldata = data.frame(fread("../data/coldata.txt", header = TRUE, sep = "\t"), row.names = 1)

# spike-in count table (aligned on dm6)
spikein_counts = counts %>% select(-"transcript_id")
spikein_counts = as.data.frame(spikein_counts)
# rownames(spikein_counts) = spikein_counts$gene_id
spikein_counts = spikein_counts %>% select(-"gene_id")
spikein_counts = as.matrix(spikein_counts)

# check order
all(rownames(coldata) == colnames(spikein_counts))

spikeins = DESeqDataSetFromMatrix(
  countData = round(spikein_counts),
  colData = coldata,
  design = ~ treatment
)
spikeins

# estimate size factor based on the Drosophila whole-cell spike-ins
spikeins = estimateSizeFactors(spikeins)
sizeFactors(spikeins)

spikein_sf = tibble(normalization = "Dros. spike-in", 
                    sample = names(sizeFactors(spikeins)), 
                    size_factor = sizeFactors(spikeins))
write_tsv(spikein_sf, glue("{result_folder}spikein_size_factors.tsv"))

# prepare tables
counts = fread("../data/RNA-Seq/nextflow_output/star_rsem/rsem.merged.transcript_counts.tsv")
counts = counts %>% group_by(gene_id) %>% summarise_all(mean) %>% dplyr::select(-transcript_id)
counts = counts %>%
  mutate_if(is.numeric, round)
counts = as.data.frame(counts)
rownames(counts) = counts$gene_id
counts = counts[-1]

coldata = fread("../data/coldata.txt", sep = "\t", header = TRUE, stringsAsFactors = TRUE)
coldata = as.data.frame(coldata)
rownames(coldata) = coldata$sample
all(rownames(coldata) == colnames(counts))

## DESeq2 protocol 
dds = DESeqDataSetFromMatrix(countData = counts,
                             colData = coldata,
                             design = ~ treatment)

# add spike-in size factors for normalization
sizeFactors(dds) = sizeFactors(spikeins)
sizeFactors(dds)

# keep = rowSums(counts(dds)) >= 10
# dds = dds[keep, ]
dds$condition = factor(dds$treatment, levels = c("WT", "treatment"))
dds$condition = relevel(dds$treatment, ref = "WT")
dds = DESeq(dds)

# result table
# SMARCAD1 KO vs. WT
ko_res = results(dds, contrast = c("treatment", "SMARCAD1_KO", "WT"))
ko_res = as.data.frame(ko_res)
ko_res["gene_name"] = rownames(ko_res)
write_tsv(ko_res, glue("{result_folder}DESeq2_res-ko_vs_wt_spikein.tsv"))
# SMARCAD1 KO WT vs. WT
kowt_res = results(dds, contrast = c("treatment", "SMARCAD1_KO_WT", "WT"))
kowt_res = as.data.frame(kowt_res)
kowt_res["gene_name"] = rownames(kowt_res)
write_tsv(kowt_res, glue("{result_folder}DESeq2_res-kowt_vs_wt_spikein.tsv"))
# SMARCAD1 KO d425 vs. WT
kod425_res = results(dds, contrast = c("treatment", "SMARCAD1_KO_d425", "WT"))
kod425_res = as.data.frame(kod425_res)
kod425_res["gene_name"] = rownames(kod425_res)
write_tsv(kod425_res, glue("{result_folder}DESeq2_res-kod425_vs_wt_spikein.tsv"))
# SMARCAD1 KO K523R vs. WT
kok523r_res = results(dds, contrast = c("treatment", "SMARCAD1_KO_K523R", "WT"))
kok523r_res = as.data.frame(kok523r_res)
kok523r_res["gene_name"] = rownames(kok523r_res)
write_tsv(kok523r_res, glue("{result_folder}DESeq2_res-kok523r_vs_wt_spikein.tsv"))

# result table
# SMARCAD1 KO vs. SMARCAD1 KO WT
ko_vs_kowt_res = results(dds, contrast = c("treatment", "SMARCAD1_KO", "SMARCAD1_KO_WT"))
ko_vs_kowt_res = as.data.frame(ko_vs_kowt_res)
ko_vs_kowt_res["gene_name"] = rownames(ko_vs_kowt_res)

ko_vs_kowt_deg = ko_vs_kowt_res %>% 
  filter(abs(log2FoldChange) > 1) %>% 
  filter(padj < 0.05)

write_tsv(ko_vs_kowt_deg, glue("{result_folder}DESeq2_res_ko-vs-kowt_padj0.05_fc1.tsv"))

kod425_vs_kowt_res = results(dds, contrast = c("treatment", "SMARCAD1_KO_d425", "SMARCAD1_KO_WT"))
kod425_vs_kowt_res = as.data.frame(kod425_vs_kowt_res)
kod425_vs_kowt_res["gene_name"] = rownames(kod425_vs_kowt_res)

koK523R_vs_kowt_res = results(dds, contrast = c("treatment", "SMARCAD1_KO_K523R", "SMARCAD1_KO_WT"))
koK523R_vs_kowt_res = as.data.frame(koK523R_vs_kowt_res)
koK523R_vs_kowt_res["gene_name"] = rownames(koK523R_vs_kowt_res)

wt_vs_kowt_res = results(dds, contrast = c("treatment", "WT", "SMARCAD1_KO_WT"))
wt_vs_kowt_res = as.data.frame(wt_vs_kowt_res)
wt_vs_kowt_res["gene_name"] = rownames(wt_vs_kowt_res)

kowt_foldchanges = tibble(gene_name = wt_vs_kowt_res$gene_name, 
                          "KO_vs_KOWT" = ko_vs_kowt_res$log2FoldChange,
                          "KO-d425_vs_KOWT" = kod425_vs_kowt_res$log2FoldChange,
                          "KO-K523R_vs_KOWT" = koK523R_vs_kowt_res$log2FoldChange,
                          "WT_vs_KOWT" = wt_vs_kowt_res$log2FoldChange)

write_tsv(kowt_foldchanges, glue("{result_folder}DESeq2_res-kowt_foldchanges.tsv"))

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
  col = c('grey', 'grey', 'grey', '#f03b20')
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
  col = c('grey', 'grey', 'grey', '#f03b20')
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
  col = c('grey', 'grey', 'grey', '#f03b20')
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
  col = c('grey', 'grey', 'grey', '#f03b20')
)
kod425_volc

# volcano plot
png(
  file = glue("{result_folder}RNA-Seq_spikein-norm_volcanos.png"),
  width = 10,
  height = 10,
  units = 'in',
  res = 500
)
ko_volc + kowt_volc + kok523r_volc + kod425_volc
dev.off()

# normalized read counts (median of ratios)
norm_counts = counts(dds, normalized = TRUE)
genes = rownames(norm_counts)
norm_counts = as_tibble(norm_counts)
norm_counts = norm_counts %>% mutate(gene_symbol = genes)
write_tsv(norm_counts, glue("{result_folder}median_of_ratios_counts_spikein.tsv"))

rlog_counts = rlog(dds, blind = FALSE)

plotPCA(rlog_counts, intgroup = "treatment", ntop = 1000, returnData = FALSE) +
  geom_point(colour = "black", size = 7, alpha = 0.4) +
  geom_point(
    aes(colour = group),
    size = 6,
    alpha = 0.4
  ) +
  labs(title = "PCA on rlog norm. expression", color = " ") +
  theme_classic() +
  theme(
    text = element_text(size = 18),
    plot.title = element_text(size = 15),
    axis.text.x = element_text(size = 13, color = "black")
  ) +
  scale_color_brewer(palette = "Set2")
ggsave(
  glue("{result_folder}RNA-Seq_PCA_plot_spikein.png"),
  width = 10,
  height = 10,
  dpi = 500,
)

rlog_counts = assay(rlog_counts)
gene_names = rownames(rlog_counts)
rlog_counts = as_tibble(rlog_counts)
rlog_counts = rlog_counts %>% mutate(gene_symbol = gene_names)
write_tsv(rlog_counts, glue("{result_folder}reg_log_norm_counts_spikein.tsv"))