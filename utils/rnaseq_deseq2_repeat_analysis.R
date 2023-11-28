suppressPackageStartupMessages({
  library("DESeq2")
  library("EnhancedVolcano")
  library("tidyverse")
  library("wigglescout")
  library("data.table")
  library("ggplot2")
  library("glue")
  library("patchwork")
})

# result / peak folder
result_folder = "../results/rna_seq_deseq2/"
bigwigs = "../data/RNA-Seq/bigwig_rpgc/"
bigwigs = list.files(bigwigs, full.names = TRUE)

repeats = "../data/bed/UCSC_RepeatMasker_ERVK_mm10.bed"

coldata = fread("../data/coldata_repeat.txt", sep = "\t", header = TRUE, stringsAsFactors = TRUE)
coldata = as.data.frame(coldata)

repeat_cov = bw_loci(bigwigs, loci = repeats)
repeat_cov = as_tibble(repeat_cov)
repeat_cov = repeat_cov %>% drop_na()
counts = repeat_cov %>% select_if(is.numeric) %>% select(starts_with("Ph")) %>% as.matrix %>% round
rownames(counts) = repeat_cov$name

## DESeq2 protocol 
dds = DESeqDataSetFromMatrix(countData = counts,
                             colData = coldata,
                             design = ~ treatment)
dds$condition = factor(dds$treatment, levels = c("WT", "treatment"))
dds$condition = relevel(dds$treatment, ref = "WT")
dds = DESeq(dds)

# result table
# SMARCAD1 KO vs. WT
ko_res = results(dds, contrast = c("treatment", "SMARCAD1_KO", "WT"))
ko_res = as.data.frame(ko_res)
ko_res["gene_name"] = rownames(ko_res)
repeat_cov["gene_name_with_id"] = ko_res["gene_name"]
write_tsv(ko_res, glue("{result_folder}DESeq2_repeats-res-ko_vs_wt.tsv"))

# repetitive element with increased expr. upon SMARCAD1 KO
top = ko_res %>% drop_na() %>% dplyr::filter(padj < 0.05) %>% 
  dplyr::filter(log2FoldChange > 0) %>% 
  top_n(wt = log2FoldChange, 50)
regions = top %>% inner_join(repeat_cov, by = c("gene_name" = "gene_name_with_id")) %>% 
  dplyr::select(seqnames, start, end)

write_tsv(regions, glue("{result_folder}DESeq2_repeats-res-ko_vs_wt-top_fc_regions.bed"), 
          col_names = FALSE)
write_tsv(regions, "../data/bed/DESeq2_repeats-res-ko_vs_wt-top_fc_regions.bed", 
          col_names = FALSE)

# SMARCAD1 KO WT vs. WT
kowt_res = results(dds, contrast = c("treatment", "SMARCAD1_KO_WT", "WT"))
kowt_res = as.data.frame(kowt_res)
kowt_res["gene_name"] = rownames(kowt_res)
write_tsv(kowt_res, glue("{result_folder}DESeq2_repeats-res-kowt_vs_wt.tsv"))
# SMARCAD1 KO d425 vs. WT
kod425_res = results(dds, contrast = c("treatment", "SMARCAD1_KO_d425", "WT"))
kod425_res = as.data.frame(kod425_res)
kod425_res["gene_name"] = rownames(kod425_res)
write_tsv(kod425_res, glue("{result_folder}DESeq2_repeats-kod425_vs_wt.tsv"))
# SMARCAD1 KO K523R vs. WT
kok523r_res = results(dds, contrast = c("treatment", "SMARCAD1_KO_K523R", "WT"))
kok523r_res = as.data.frame(kok523r_res)
kok523r_res["gene_name"] = rownames(kok523r_res)
write_tsv(kok523r_res, glue("{result_folder}DESeq2_repeats-kok523r_vs_wt.tsv"))

# volcano plots
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

ko_res_aggr = ko_res %>%  drop_na() %>% mutate(direction = case_when(log2FoldChange > 1 & pvalue < 0.05 ~ "upregulated",
                                                                     log2FoldChange < -1 & pvalue < 0.05 ~ "downregulated")) %>% 
  group_by(direction) %>% summarise(count = n()) %>% drop_na() %>% mutate(contrast = "KO vs. WT") 
  
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

kowt_res_aggr = kowt_res %>%  drop_na() %>% mutate(direction = case_when(log2FoldChange > 1 & pvalue < 0.05 ~ "upregulated",
                                                                     log2FoldChange < -1 & pvalue < 0.05 ~ "downregulated")) %>% 
  group_by(direction) %>% summarise(count = n()) %>% drop_na() %>% mutate(contrast = "KO-WT vs. WT") 

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

kok523r_res_aggr = kok523r_res %>%  drop_na() %>% mutate(direction = case_when(log2FoldChange > 1 & pvalue < 0.05 ~ "upregulated",
                                                                       log2FoldChange < -1 & pvalue < 0.05 ~ "downregulated")) %>% 
  group_by(direction) %>% summarise(count = n()) %>% drop_na() %>% mutate(contrast = "KO K523R vs. WT") 

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

kod425_res_aggr = kod425_res %>%  drop_na() %>% mutate(direction = case_when(log2FoldChange > 1 & pvalue < 0.05 ~ "upregulated",
                                                                               log2FoldChange < -1 & pvalue < 0.05 ~ "downregulated")) %>% 
  group_by(direction) %>% summarise(count = n()) %>% drop_na() %>% mutate(contrast = "KO d425 vs. WT") 

# volcano plot
png(
  file = glue("{result_folder}RNA-Seq_repeats_volcanos.png"),
  width = 10,
  height = 10,
  units = 'in',
  res = 500
)
ko_volc + kowt_volc + kok523r_volc + kod425_volc
dev.off()

pdf(
  file = glue("{result_folder}RNA-Seq_repeats_volcanos.pdf"),
  width = 10,
  height = 10
)
ko_volc + kowt_volc + kok523r_volc + kod425_volc
dev.off()

aggrs = rbind(ko_res_aggr, kowt_res_aggr,  kok523r_res_aggr, kod425_res_aggr)

# order = aggrs %>% filter(direction == "upregulated") %>% arrange(desc(count)) %>% pull(contrast)
order = factor(aggrs$contrast, levels = c("KO vs. WT", "KO-WT vs. WT", "KO K523R vs. WT", "KO d425 vs. WT"))

ggplot(data = aggrs, aes(x = order, y = count, fill = direction)) +
  geom_bar(stat="identity", position=position_dodge(), color = "black") +
  ylim(0, 300) +
  labs(
    title = "repmasker_clusters_lt2kb",
    x = " ",
    y = "# of diff. expressed repeats",
    fill = " "
  ) +
  scale_fill_manual(values = c("#bdbdbd", "#fc9272")) +
  theme_classic() +
  theme(
    text = element_text(size = 20),
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 17, color = "black"),
    axis.text.y = element_text(size = 20, color = "black")
  )

ggsave(
  glue("{result_folder}RNA-Seq_repeats_bars.pdf"),
  plot = last_plot(),
  width = 12,
  height = 10,
  device = "pdf",
)

ggsave(
  glue("{result_folder}RNA-Seq_repeats_bars.png"),
  plot = last_plot(),
  width = 12,
  height = 10,
  dpi = 300
)

png(
  file = glue("{result_folder}RNA-Seq_repeats_volcanos.png"),
  width = 10,
  height = 10,
  units = 'in',
  res = 500
)
ko_volc + kowt_volc + kok523r_volc + kod425_volc
dev.off()

pdf(
  file = glue("{result_folder}RNA-Seq_repeats_volcanos.pdf"),
  width = 10,
  height = 10
)
ko_volc + kowt_volc + kok523r_volc + kod425_volc
dev.off()




