suppressPackageStartupMessages({
  library("glue")
  library("tidyverse")
  library("data.table")
  library("ggrepel")
  library("ggpubr")
  library("RColorBrewer")
  library("biomaRt")
})

result_folder = "../results/rna_seq_deseq2/"

# IAPEz-int
te = fread("../results/rna_seq_deseq2/TEtranscripts/ERVK_IAPez_extended/TEtranscripts_SMARCAD1_KO_UCSC_Repeatmasker_sigdiff_gene_TE.txt")
te_all = fread("../results/rna_seq_deseq2/TEtranscripts/ERVK_IAPez_extended/TEtranscripts_SMARCAD1_KO_UCSC_Repeatmasker_gene_TE_analysis.txt")
# ERVKs
te_all_iapez = fread("../results/rna_seq_deseq2/TEtranscripts/ERVK_IAPez_extended/TEtranscripts_SMARCAD1_KO_IAPEz-int_extended_gene_TE_analysis.txt")

ensembl = useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")

attr = listAttributes(ensembl)
refseq_annotations = tibble(refseq = unlist(lapply(te$V1, function(x) {strsplit(x, ".", fixed = TRUE)[[1]][1]})))
xms = refseq_annotations[grep(refseq_annotations$refseq, pattern = "XM"),]
res = getBM(values = refseq_annotations$refseq, attributes=c("refseq_mrna", "mgi_symbol"), filters= "refseq_mrna", mart = ensembl)
xm_res =  getBM(values = xms$refseq, attributes=c("refseq_mrna_predicted", "mgi_symbol"), filters= "refseq_mrna_predicted", mart = ensembl)


te_annot = te %>% 
  mutate(refseq = refseq_annotations$refseq) %>% 
  left_join(., res, by = c("refseq" = "refseq_mrna")) %>% 
  dplyr::select(refseq_id = V1, refseq, gene_symbol = mgi_symbol, everything())

for(x in seq(te_annot$refseq)) {
  if(te_annot$refseq[x] %in% xm_res$refseq_mrna_predicted) {
    gene = xm_res %>% dplyr::filter(refseq_mrna_predicted == te_annot$refseq[x]) %>% 
      pull(mgi_symbol) %>% unique %>% as.character
    te_annot$gene_symbol[x] = gene
    print(te_annot$refseq[x])
  }
}

leftover = te_annot[is.na(te_annot$gene_symbol)]

ncbi_all = fread("C:/Szabolcs/Karolinska/Data/reference_data/NCBI_RefSeq-all_mm10.txt")
ncbi_all = ncbi_all %>% mutate(refseq = unlist(lapply(ncbi_all$name, function(x) {strsplit(x, ".", fixed = TRUE)[[1]][1]})))
leftover = leftover %>% left_join(., ncbi_all, by = c("refseq" = "refseq")) %>% dplyr::select(refseq, gene_symbol = name2)

for(x in seq(te_annot$refseq)) {
  if(te_annot$refseq[x] %in% leftover$refseq) {
    gene = leftover %>% dplyr::filter(refseq == te_annot$refseq[x]) %>% 
      pull(gene_symbol) %>% unique %>% as.character
    te_annot$gene_symbol[x] = gene
  }
}

te_annot = drop_na(te_annot)

# volcano of ERVKs, IAPez
reps = te_all[!grep(te_all$V1, pattern = "XR|XM|NR|NM")]
iapez = reps[grep(reps$V1, pattern = "IAPEz")]
ervks = reps[grep(reps$V1, pattern = "ERV")]
iapez_ervks = rbind(iapez, ervks)
iapez_ervks = iapez_ervks %>% separate(V1, sep = ":", into = "te_name")
iapez_ervks = drop_na(iapez_ervks)

volc_input = iapez_ervks %>% 
  # group_by(gene_symbol) %>%
  # dplyr::filter(avg_log2FC == max(abs(avg_log2FC), na.rm=TRUE)) %>% 
  mutate(group = case_when(
    log2FoldChange > 0.1 & padj < 0.05 ~ "up",
    log2FoldChange < -0.1 & padj < 0.05 ~ "down",
    log2FoldChange >= -0.1 & log2FoldChange <= 0.1 ~ "unaltered",
    TRUE ~ "unaltered"
  ))
volc_input = volc_input %>% mutate(sign_label = case_when(
  log2FoldChange > 0.25 & padj < 0.05 ~ te_name,
  log2FoldChange < -0.25 & padj < 0.05 ~ te_name,
  log2FoldChange >= -0.25 & log2FoldChange <= 0.25 ~ ""
))

labels = volc_input %>% pull(sign_label)

# add color, size, alpha indications
cols = c("up" = "#fc9272", "down" = "#a1d99b", "unaltered" = "grey")
sizes = c("up" = 4, "down" = 4, "unaltered" = 2)
alphas = c("up" = 1, "down" = 1, "unaltered" = 0.5)

ggplot_volc = volc_input %>%
  ggplot(aes(x = log2FoldChange,
             y = -log10(padj),
             size = group,
             fill = group,
             alpha = group)) +
  geom_point(shape = 21, # Specify shape and colour as fixed local parameters    
             colour = "black") +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") +
  geom_vline(xintercept = c(-0.1, 0.1),
             linetype = "dashed") +
  scale_fill_manual(values = cols) +
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  scale_x_continuous(breaks = c(seq(-2, 2, 0.5)),  	 
                     limits = c(-2, 2)) +
  labs(
    title = "TEtranscripts (DESeq2) output - Smarcad1 KO vs. WT",
    subtitle = "UCSC RepeatMasker: IAPez, ERVK repeat classes",
    x = "log2FoldChange",
    y = "-log10 adj.p-value",
    fill = " "
  ) +
  ylim(0, 100) +
  guides(alpha = FALSE, size = FALSE, fill = guide_legend(override.aes = list(size = 5))) +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 20, color = "black")
  ) +
  geom_text_repel(label = labels, size = 6, max.overlaps = Inf) # add labels
ggplot_volc

ggsave(
  glue("{result_folder}TEtranscripts-IAPEz_ERVKs_volc.pdf"),
  plot = last_plot(),
  width = 12,
  height = 8,
  dpi = 300,
)

# volcano of only, individual IAPez
reps = te_all_iapez[!grep(te_all$V1, pattern = "XR|XM|NR|NM")]
iapez = reps[grep(reps$V1, pattern = "IAPEz")]
iapez = iapez %>% separate(V1, sep = ":", into = "te_name")
iapez = drop_na(iapez)

volc_input = iapez %>% 
  # group_by(gene_symbol) %>%
  # dplyr::filter(avg_log2FC == max(abs(avg_log2FC), na.rm=TRUE)) %>% 
  mutate(group = case_when(
    log2FoldChange > 0.1 & padj < 0.05 ~ "up",
    log2FoldChange < -0.1 & padj < 0.05 ~ "down",
    log2FoldChange >= -0.1 & log2FoldChange <= 0.1 ~ "unaltered",
    TRUE ~ "unaltered"
  ))
volc_input = volc_input %>% mutate(sign_label = case_when(
  log2FoldChange > 3 & padj < 0.001 ~ te_name,
  log2FoldChange < -3 & padj < 0.001 ~ te_name,
  padj < 1e-10 ~ te_name,
  log2FoldChange >= -3 & log2FoldChange <= 3 & padj > 1e-10 ~ ""
))

labels = volc_input %>% pull(sign_label)

# add color, size, alpha indications
cols = c("up" = "#fc9272", "down" = "#a1d99b", "unaltered" = "grey")
sizes = c("up" = 4, "down" = 4, "unaltered" = 2)
alphas = c("up" = 1, "down" = 1, "unaltered" = 0.5)

ggplot_volc_iapez = volc_input %>%
  ggplot(aes(x = log2FoldChange,
             y = -log10(padj),
             size = group,
             fill = group,
             alpha = group)) +
  geom_point(shape = 21, # Specify shape and colour as fixed local parameters    
             colour = "black") +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") +
  geom_vline(xintercept = c(-0.1, 0.1),
             linetype = "dashed") +
  scale_fill_manual(values = cols) +
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  scale_x_continuous(breaks = c(seq(-5, 5, 1)),  	 
                     limits = c(-5, 5)) +
  labs(
    title = "TEtranscripts (DESeq2) output - Smarcad1 KO vs. WT",
    subtitle =  "IAPEz-int elements", 
    x = "log2FoldChange",
    y = "-log10 adj.p-value",
    fill = " "
  ) +
  ylim(0, 50) +
  guides(alpha = FALSE, size = FALSE, fill = guide_legend(override.aes = list(size = 5))) +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 20, color = "black")
  ) +
  geom_text_repel(label = labels, size = 6, max.overlaps = Inf) # add labels
ggplot_volc_iapez

ggsave(
  glue("{result_folder}TEtranscripts-IAPEz_volc.pdf"),
  plot = last_plot(),
  width = 12,
  height = 8,
  dpi = 300,
)

iapez_gtf = fread("C:/Szabolcs/Karolinska/Data/reference_data/repeatmasker_mm10-IAPEz-int-gene_name_mod.gtf")
sign_pos_iapez = volc_input %>% dplyr::filter(log2FoldChange > 2 & padj < 0.05) %>% pull(te_name)
sign_pos_iapez = iapez_gtf[unname(unlist(lapply(sign_pos_iapez, function(x) {grep(iapez_gtf$V9, pattern = x)}))),]
sign_pos_iapez = sign_pos_iapez %>% dplyr::select(V1, V4, V5) %>% mutate(V1 = paste0("chr", V1)) 
write_tsv(sign_pos_iapez, "../data/bed/RNA-Seq-sign_upreg_IAPez-int_regions_uponSmarcad1_KO.bed", col_names = FALSE)


