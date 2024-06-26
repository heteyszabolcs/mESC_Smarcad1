if (!require("pacman"))
  install.packages("pacman")
pacman::p_load("data.table",
               "glue",
               "tidyverse",
               "data.table",
               "ggrepel",
               "ggpubr",
               "RColorBrewer",
               "biomaRt",
               "ggrastr"
)

# export
result_folder = "../results/rna_seq_deseq2/"

## TEtranscripts output
# RepeatMasker UCSC 
te = fread(
  "../results/rna_seq_deseq2/TEtranscripts/ERVK_IAPez_extended/TEtranscripts_SMARCAD1_KO_UCSC_Repeatmasker_sigdiff_gene_TE.txt"
)
te_all = fread(
  "../results/rna_seq_deseq2/TEtranscripts/ERVK_IAPez_extended/TEtranscripts_SMARCAD1_KO_UCSC_Repeatmasker_gene_TE_analysis.txt"
)

# ERVKs
te_all_ervk = fread(
  "../results/rna_seq_deseq2/TEtranscripts/ERVK_extended/TEtranscripts_SMARCAD1_KO_ERVK_extended_gene_TE_analysis.txt"
)

# IAPEz-int
te_all_iapez = fread(
  "../results/rna_seq_deseq2/TEtranscripts/ERVK_IAPez_extended/TEtranscripts_SMARCAD1_KO_IAPEz-int_extended_gene_TE_analysis.txt"
)
te_all_iapez = te_all_iapez %>% separate(V1, sep = "\\:", into = c("TE", "rest"), 
                                         remove = FALSE) %>% 
  dplyr::select(-rest) %>% 
  mutate(TE = ifelse(str_detect(TE, "IAP"), TE, "non_TE"))

te_iapez_norm = fread("../results/rna_seq_deseq2/TEtranscripts/ERVK_IAPez_extended/TEtranscripts_SMARCAD1_KO_IAPEz_normalized_expr.tsv")

# diff. expressed IAPez-int elements after Smarcad1 KO
te_sign_iapez = fread(
  "../results/rna_seq_deseq2/TEtranscripts/ERVK_IAPez_extended/TEtranscripts_SMARCAD1_KO_IAPEz-int_extended_sigdiff_gene_TE.txt"
)
te_sign_iapez = te_sign_iapez[grep("IAPEz-int", x = te_sign_iapez$V1), ]
unname(sapply(te_sign_iapez %>% arrange(desc(log2FoldChange)) %>% top_n(5, wt = log2FoldChange) %>% pull(V1), function(x) {
  strsplit(x, ":")[[1]][1]
}))

iapez_gtf = fread("C:/Szabolcs/Karolinska/Data/reference_data/repeatmasker_mm10-IAPEz-int.gtf")
te_sign_iapez_mod = te_sign_iapez %>% separate(V1, sep = "\\:", into = c("IAPEz_name", "rest")) %>% 
  dplyr::select(-rest)

te_sign_iapez_regions = iapez_gtf %>% separate(V9, sep = "transcript_id ", into = c("V9", "rest")) %>% 
  separate(rest, sep = "\\; family_id", into = c("rest2", "rest3")) %>% 
  dplyr::select(V1, V4, V5, rest2) %>% mutate(rest2 = str_replace_all(rest2, "\"", "")) %>% 
  mutate(V1 = paste0("chr", V1)) %>% 
  rename(seq = V1, start = V4, end = V5, name = rest2) %>% 
  inner_join(., te_sign_iapez_mod, by = c("name" = "IAPEz_name")) %>% 
  dplyr::select("baseMean":"padj", "name", "seq", "start", "end")

te_all_iapez = te_all_iapez %>% mutate(diff_expressed = ifelse(TE %in% te_sign_iapez_regions$name,
                                       TRUE, FALSE))

write_tsv(te_sign_iapez_regions,
          glue("{result_folder}TEtranscripts_SMARCAD1_KO_IAPEz-int_extended_sigdiff_TE-regions.txt")
)
write_tsv(te_all_iapez,
          glue("{result_folder}TEtranscripts_SMARCAD1_KO_IAPEz-int_extended_gene_TE-labeled.txt")
)


# diff. expressed IAPez-int elements after rescue
te_sign_iapez_kowt = fread(
  "../results/rna_seq_deseq2/TEtranscripts/ERVK_IAPez_extended/TEtranscripts_SMARCAD1_WTKO_IAPEz_extended_sigdiff_gene_TE.txt"
)

te_all_iapez_kowt = fread(
  "../results/rna_seq_deseq2/TEtranscripts/ERVK_IAPez_extended/TEtranscripts_SMARCAD1_WTKO_IAPEz_extended_gene_TE_analysis.txt"
)

te_iapez_norm_kowt = fread("../results/rna_seq_deseq2/TEtranscripts/ERVK_IAPez_extended/TEtranscripts_SMARCAD1_KOWT_IAPEz_normalized_expr.tsv")

# annotation
# ensembl = useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")
# 
# attr = listAttributes(ensembl)
# refseq_annotations = tibble(refseq = unlist(lapply(te$V1, function(x) {
#   strsplit(x, ".", fixed = TRUE)[[1]][1]
# })))
# xms = refseq_annotations[grep(refseq_annotations$refseq, pattern = "XM"), ]
# res = getBM(
#   values = refseq_annotations$refseq,
#   attributes = c("refseq_mrna", "mgi_symbol"),
#   filters = "refseq_mrna",
#   mart = ensembl
# )
# xm_res =  getBM(
#   values = xms$refseq,
#   attributes = c("refseq_mrna_predicted", "mgi_symbol"),
#   filters = "refseq_mrna_predicted",
#   mart = ensembl
# )
# 
# 
# te_annot = te %>%
#   mutate(refseq = refseq_annotations$refseq) %>%
#   left_join(., res, by = c("refseq" = "refseq_mrna")) %>%
#   dplyr::select(refseq_id = V1,
#                 refseq,
#                 gene_symbol = mgi_symbol,
#                 everything())
# 
# for (te_iapez_norm in seq(te_annot$refseq)) {
#   if (te_annot$refseq[te_iapez_norm] %in% xm_res$refseq_mrna_predicted) {
#     gene = xm_res %>% dplyr::filter(refseq_mrna_predicted == te_annot$refseq[te_iapez_norm]) %>%
#       pull(mgi_symbol) %>% unique %>% as.character
#     te_annot$gene_symbol[te_iapez_norm] = gene
#     print(te_annot$refseq[te_iapez_norm])
#   }
# }
# 
# leftover = te_annot[is.na(te_annot$gene_symbol)]
# 
# # annotation by whole RefSeq table (NCBI)
# ncbi_all = fread("C:/Szabolcs/Karolinska/Data/reference_data/NCBI_RefSeq-all_mm10.txt")
# ncbi_all = ncbi_all %>% mutate(refseq = unlist(lapply(ncbi_all$name, function(x) {
#   strsplit(x, ".", fixed = TRUE)[[1]][1]
# })))
# leftover = leftover %>% left_join(., ncbi_all, by = c("refseq" = "refseq")) %>% dplyr::select(refseq, gene_symbol = name2)
# 
# for (te_iapez_norm in seq(te_annot$refseq)) {
#   if (te_annot$refseq[te_iapez_norm] %in% leftover$refseq) {
#     gene = leftover %>% dplyr::filter(refseq == te_annot$refseq[te_iapez_norm]) %>%
#       pull(gene_symbol) %>% unique %>% as.character
#     te_annot$gene_symbol[te_iapez_norm] = gene
#   }
# }
# 
# te_annot = drop_na(te_annot)

## plots
# volcano of ERVKs, IAPez
reps = te_all[!grep(te_all$V1, pattern = "XR|XM|NR|NM")]
iapez = reps[grep(reps$V1, pattern = "IAPEz")]
ervks = reps[grep(reps$V1, pattern = "ERV")]
iapez_ervks = rbind(iapez, ervks)
iapez_ervks = iapez_ervks %>% separate(V1, sep = ":", into = "te_name")
iapez_ervks = drop_na(iapez_ervks)

# ERVK volcano
volc_input_classes = iapez_ervks  %>%
  mutate(
    group = case_when(
      log2FoldChange > 1 & padj < 0.05 ~ "up",
      log2FoldChange < -1 & padj < 0.05 ~ "down",
      log2FoldChange >= -1 & log2FoldChange <= 1 ~ "unaltered",
      TRUE ~ "unaltered"
    )
  )
volc_input_classes = volc_input_classes %>% mutate(
  sign_label = case_when(
    log2FoldChange > 1 & padj < 0.001 ~ te_name,
    log2FoldChange < -1 & padj < 0.001 ~ te_name,
    log2FoldChange >= -1 & log2FoldChange <= 1 | padj > 0.05 ~ ""
  )
)

labels = volc_input_classes %>% pull(sign_label)

# add color, size, alpha indications
cols = c(
  "up" = "#fc9272",
  "down" = "#a1d99b",
  "unaltered" = "grey"
)
sizes = c("up" = 4,
          "down" = 4,
          "unaltered" = 2)
alphas = c("up" = 1,
           "down" = 1,
           "unaltered" = 0.5)

ggplot_volc_classes = volc_input_classes %>%
  ggplot(aes(
    x = log2FoldChange,
    y = -log10(padj),
    size = group,
    fill = group,
    alpha = group
  )) +
  geom_point(shape = 21, # Specify shape and colour as fixed local parameters
             colour = "black") +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1),
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
  guides(alpha = FALSE,
         size = FALSE,
         fill = guide_legend(override.aes = list(size = 5))) +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 10),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 20, color = "black")
  ) +
  geom_text_repel(label = labels, size = 6) # add labels
ggplot_volc_classes

ggsave(
  glue("{result_folder}TEtranscripts-IAPEz_ERVKs_volc.pdf"),
  plot = last_plot(),
  width = 12,
  height = 8,
  dpi = 300,
)

# ERVKs (elementwise)
ervks_ind = te_all_ervk[grep(te_all_ervk$V1, pattern = "MMERV")]
ervks_ind = ervks_ind %>% separate(V1, sep = ":", into = "te_name")
ervks_ind = drop_na(ervks_ind)

volc_input_ervks = ervks_ind %>%
  mutate(
    group = case_when(
      log2FoldChange > 1 & padj < 0.05 ~ "up",
      log2FoldChange < -1 & padj < 0.05 ~ "down",
      log2FoldChange >= -1 & log2FoldChange <= 1 ~ "unaltered",
      TRUE ~ "unaltered"
    )
  )
volc_input_ervks = volc_input_ervks %>% mutate(
  sign_label = case_when(
    log2FoldChange > 1 & padj < 0.001 ~ te_name,
    log2FoldChange < -1 & padj < 0.001 ~ te_name,
    log2FoldChange >= -1 & log2FoldChange <= 1 | padj > 0.05 ~ ""
  )
)

labels = volc_input_ervks %>% pull(sign_label)

# add color, size, alpha indications
cols = c(
  "up" = "#fc9272",
  "down" = "#a1d99b",
  "unaltered" = "grey"
)
sizes = c("up" = 4,
          "down" = 4,
          "unaltered" = 2)
alphas = c("up" = 1,
           "down" = 1,
           "unaltered" = 0.5)

ggplot_volc_ervks = volc_input_ervks %>%
  ggplot(aes(
    x = log2FoldChange,
    y = -log10(padj),
    size = group,
    fill = group,
    alpha = group
  )) +
  geom_point(shape = 21, # Specify shape and colour as fixed local parameters
             colour = "black") +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1),
             linetype = "dashed") +
  scale_fill_manual(values = cols) +
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  scale_x_continuous(breaks = c(seq(-2, 2, 0.5)),
                     limits = c(-2, 2)) +
  labs(
    title = "TEtranscripts (DESeq2) output - Smarcad1 KO vs. WT",
    subtitle = "UCSC RepeatMasker: ERVK elements",
    x = "log2FoldChange",
    y = "-log10 adj.p-value",
    fill = " "
  ) +
  ylim(0, 25) +
  guides(alpha = FALSE,
         size = FALSE,
         fill = guide_legend(override.aes = list(size = 5))) +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 10),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 20, color = "black")
  ) +
  geom_text_repel(label = labels, size = 6) # add labels
ggplot_volc_ervks

ggsave(
  glue("{result_folder}TEtranscripts-ERVKs_volc.pdf"),
  plot = last_plot(),
  width = 12,
  height = 8,
  dpi = 300,
)


# IAPez-int volcano
# volcano of individual IAPez elements
reps = te_all_iapez[!grep(te_all_iapez$V1, pattern = "XR|XM|NR|NM")]
iapez = reps[grep(reps$V1, pattern = "IAPEz")]
iapez = iapez %>% separate(V1, sep = ":", into = "te_name")
iapez = drop_na(iapez)

volc_input_iapez = iapez %>%
  mutate(
    group = case_when(
      log2FoldChange > 1 & padj < 0.05 ~ "up",
      log2FoldChange < -1 & padj < 0.05 ~ "down",
      log2FoldChange >= -1 & log2FoldChange <= 1 ~ "unaltered",
      TRUE ~ "unaltered"
    )
  )
volc_input_iapez = volc_input_iapez %>% mutate(
  sign_label = case_when(
    log2FoldChange > 1 & padj < 10e-12 ~ te_name,
    log2FoldChange < -1 & padj < 10e-12 ~ te_name,
    log2FoldChange >= -1 & log2FoldChange <= 1 | padj > 0.05 ~ ""
  )
)

labels = volc_input_iapez %>% pull(sign_label)

# add color, size, alpha indications
cols = c(
  "up" = "#fc9272",
  "down" = "#a1d99b",
  "unaltered" = "grey"
)
sizes = c("up" = 4,
          "down" = 4,
          "unaltered" = 2)
alphas = c("up" = 1,
           "down" = 1,
           "unaltered" = 0.5)

ggplot_volc_iapez = volc_input_iapez %>%
  ggplot(aes(
    x = log2FoldChange,
    y = -log10(padj),
    size = group,
    fill = group,
    alpha = group
  )) +
  geom_point(shape = 21, # Specify shape and colour as fixed local parameters
             colour = "black") +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1),
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
  guides(alpha = FALSE,
         size = FALSE,
         fill = guide_legend(override.aes = list(size = 5))) +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 10),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 20, color = "black")
  ) +
  geom_text_repel(label = labels,
                  size = 6,
                  max.overlaps = Inf) # add labels
ggplot_volc_iapez

ggsave(
  glue("{result_folder}TEtranscripts-IAPEz_volc.pdf"),
  plot = last_plot(),
  width = 12,
  height = 8,
  dpi = 300,
)

# Smarcad1 KO WT rescue
reps = te_all_iapez_kowt[!grep(te_all_iapez_kowt$V1, pattern = "XR|XM|NR|NM")]
iapez = reps[grep(reps$V1, pattern = "IAPEz")]
iapez = iapez %>% separate(V1, sep = ":", into = "te_name")
iapez = drop_na(iapez)

volc_input_iapez_kowt = iapez %>%
  mutate(
    group = case_when(
      log2FoldChange > 1 & padj < 0.05 ~ "up",
      log2FoldChange < -1 & padj < 0.05 ~ "down",
      log2FoldChange >= -1 & log2FoldChange <= 1 ~ "unaltered",
      TRUE ~ "unaltered"
    )
  )
volc_input_iapez_kowt = volc_input_iapez_kowt %>% mutate(
  sign_label = case_when(
    log2FoldChange > 1 & padj < 10e-12 ~ te_name,
    log2FoldChange < -1 & padj < 10e-12 ~ te_name,
    log2FoldChange < -3 & padj < 0.05 ~ te_name,
    log2FoldChange > 3 & padj < 0.05 ~ te_name,
    log2FoldChange >= -1 & log2FoldChange <= 1 | padj > 0.05 ~ ""
  )
)

labels = volc_input_iapez_kowt %>% pull(sign_label)

# add color, size, alpha indications
cols = c(
  "up" = "#fc9272",
  "down" = "#a1d99b",
  "unaltered" = "grey"
)
sizes = c("up" = 4,
          "down" = 4,
          "unaltered" = 2)
alphas = c("up" = 1,
           "down" = 1,
           "unaltered" = 0.5)

ggplot_volc_iapez_kowt = volc_input_iapez_kowt %>%
  ggplot(aes(
    x = log2FoldChange,
    y = -log10(padj),
    size = group,
    fill = group,
    alpha = group
  )) +
  geom_point(shape = 21, # Specify shape and colour as fixed local parameters
             colour = "black") +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1),
             linetype = "dashed") +
  scale_fill_manual(values = cols) +
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  scale_x_continuous(breaks = c(seq(-5, 5, 1)),
                     limits = c(-5, 5)) +
  labs(
    title = "TEtranscripts (DESeq2) output - Smarcad1 KO vs. KO-WT",
    subtitle =  "IAPEz-int elements",
    x = "log2FoldChange",
    y = "-log10 adj.p-value",
    fill = " "
  ) +
  ylim(0, 50) +
  guides(alpha = FALSE,
         size = FALSE,
         fill = guide_legend(override.aes = list(size = 5))) +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 10),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 20, color = "black")
  ) +
  geom_text_repel(label = labels,
                  size = 6,
                  max.overlaps = Inf) # add labels
ggplot_volc_iapez_kowt

ggsave(
  glue("{result_folder}TEtranscripts-IAPEz_volc_KOWT.pdf"),
  plot = last_plot(),
  width = 12,
  height = 8,
  dpi = 300,
)

# arrange
ggarrange(ggplot_volc_iapez, ggplot_volc_iapez_kowt)

ggsave(
  glue("{result_folder}TEtranscripts-IAPEz_Smarcad1KO_volcanos.pdf"),
  plot = last_plot(),
  width = 14,
  height = 13,
  dpi = 300,
)

# arrange
ggarrange(ggplot_volc_classes, ggplot_volc_ervks, ggplot_volc_iapez)

# retrieve regions of upregulated IAPez-int elements
iapez_gtf = fread(
  "C:/Szabolcs/Karolinska/Data/reference_data/repeatmasker_mm10-IAPEz-int-gene_name_mod.gtf"
) # repeatmasker table
sign_pos_iapez = volc_input_iapez %>% dplyr::filter(log2FoldChange > 2 &
                                                padj < 0.05) %>% 
  arrange(desc(log2FoldChange)) %>% pull(te_name) 
sign_pos_iapez = iapez_gtf[unname(unlist(lapply(sign_pos_iapez, function(x) {
  grep(iapez_gtf$V9, pattern = x)
}))), ]
sign_pos_iapez = sign_pos_iapez %>% dplyr::select(V1, V4, V5) %>% mutate(V1 = paste0("chr", V1))
write_tsv(
  sign_pos_iapez,
  "../data/bed/RNA-Seq-sign_upreg_IAPez-int_regions_uponSmarcad1_KO.bed",
  col_names = FALSE
)

ervk_gtf = fread(
  "C:/Szabolcs/Karolinska/Data/reference_data/repeatmasker_mm10-ERVK-int-gene_name_mod.gtf"
)
sign_pos_ervks = volc_input_ervks %>% dplyr::filter(log2FoldChange > 2 &
                                                padj < 0.05) %>% 
  arrange(desc(log2FoldChange)) %>% pull(te_name) 
sign_pos_ervks = ervk_gtf[unname(unlist(lapply(sign_pos_ervks, function(x) {
  grep(ervk_gtf$V9, pattern = x)
}))), ]
sign_pos_ervks = sign_pos_ervks %>% dplyr::select(V1, V4, V5) %>% mutate(V1 = paste0("chr", V1))
write_tsv(
  sign_pos_ervks,
  "../data/bed/RNA-Seq-sign_upreg_ERVK_regions_uponSmarcad1_KO.bed",
  col_names = FALSE
)

volc_input_ervks_iapez = rbind(volc_input_ervks, volc_input_iapez)
iapez_ervk_gtf = rbind(iapez_gtf, ervk_gtf)
sign_pos_ervks_iapez = volc_input_ervks_iapez %>% dplyr::filter(log2FoldChange > 2 &
                                                      padj < 0.05) %>% 
  arrange(desc(log2FoldChange)) %>% pull(te_name) 

sign_pos_ervks_iapez = iapez_ervk_gtf[unname(unlist(lapply(sign_pos_ervks_iapez, function(x) {
  grep(iapez_ervk_gtf$V9, pattern = x)
}))), ]
sign_pos_ervks_iapez = sign_pos_ervks_iapez %>% dplyr::select(V1, V4, V5) %>% mutate(V1 = paste0("chr", V1))
write_tsv(
  sign_pos_ervks_iapez,
  "../data/bed/RNA-Seq-sign_upreg_ERVK_IAPEz_regions_uponSmarcad1_KO.bed",
  col_names = FALSE
)

### scatters
# norm expression of KOWT rescue
te_iapez_norm_kowt = te_iapez_norm_kowt %>% mutate(., KOWT_mean = round(rowMeans(dplyr::select(., starts_with("SMARCAD1_KOWT")), na.rm = TRUE)),
                                         Smarcad1_KO_mean = round(rowMeans(dplyr::select(., starts_with("SMARCAD1_KO_rep")), na.rm = TRUE))) %>% 
  dplyr::select(gene = V1, Smarcad1_KO_mean, KOWT_mean) %>% 
  separate(gene, sep = ":", into = "gene")

# norm expression of KO cell lines
te_iapez_norm = te_iapez_norm %>% mutate(., WT_mean = round(rowMeans(dplyr::select(., starts_with("WT")), na.rm = TRUE)),
                             Smarcad1_KO_mean = round(rowMeans(dplyr::select(., starts_with("SMARC")), na.rm = TRUE))) %>% 
  dplyr::select(gene = V1, starts_with("Smarcad1_KO_m"), starts_with("WT_mean")) %>% 
  separate(gene, sep = ":", into = "gene")

sign_alt_iapez = volc_input_iapez %>% dplyr::filter(abs(log2FoldChange) > 2 &
                                                      padj < 0.05) %>% 
  arrange(desc(log2FoldChange)) %>% pull(te_name) 

unalt_iapez = volc_input_iapez %>% dplyr::filter(log2FoldChange < 2 | log2FoldChange > -2 |
                                                      padj > 0.05) %>% 
  arrange(desc(log2FoldChange)) %>% pull(te_name) 

sign_alt_iapez_rescue = volc_input_iapez_kowt %>% dplyr::filter(abs(log2FoldChange) > 2 &
                                                                  padj < 0.05) %>%
  arrange(desc(log2FoldChange)) %>% pull(te_name)

te_iapez_norm = te_iapez_norm %>% mutate(IAPez = case_when(gene %in% sign_alt_iapez ~ "alt. IAPEz-int",
                                   gene %in% unalt_iapez ~ "unaltered IAPEz-int", 
                                   .default = "non-IAPEz"))

scatter = ggplot(te_iapez_norm, aes(x = Smarcad1_KO_mean, y = WT_mean, color = IAPez)) +
  ggrastr::geom_point_rast() +
  xlim(0, 1000) +
  ylim(0, 1000) +
  scale_color_manual(values = c("#f0f0f0", "#fdae6b", "#de2d26")) +
  theme_minimal() +
  theme(
    text = element_text(size = 18),
    legend.text = element_text(size = 18),
    plot.title = element_text(size = 18, face = "bold"),
    axis.text.x = element_text(size = 16, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.title.x = element_text(size = 16, color = "black"),
    axis.title.y = element_text(size = 16, color = "black")
  ) +
  guides(fill = "none", color = guide_legend(title = "type")) +
  # geom_abline(slope = 1,
  #             intercept = 0,
  #             color = "#de2d26",
  #             linetype = "dashed") +
  labs(title = "DESeq2 normalized expr.",
       x = "Smarcad1 KO",
       y = "Smarcad1 WT")

te_all_iapez_mod = te_all_iapez %>% separate(V1, sep = ":", into = "gene") 
te_all_iapez_mod = drop_na(te_all_iapez_mod)
te_all_iapez_mod = te_all_iapez_mod %>% mutate(IAPez = case_when(gene %in% sign_alt_iapez ~ "alt. IAPEz-int",
                                                                                gene %in% unalt_iapez ~ "unaltered IAPEz-int", 
                                                                                .default = "non-IAPEz"))
                                               

# MA plot
ma = ggplot(te_all_iapez_mod, aes(x = log2(baseMean), y = log2FoldChange, color = IAPez)) +
  ggrastr::geom_point_rast() +
  xlim(0, 15) +
  ylim(-10, 10) +
  scale_color_manual(values = c("#de2d26", "#f0f0f0", "orange")) +
  theme_minimal() +
  theme(
    text = element_text(size = 18),
    legend.text = element_text(size = 18),
    plot.title = element_text(size = 18, face = "bold"),
    axis.text.x = element_text(size = 16, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.title.x = element_text(size = 16, color = "black"),
    axis.title.y = element_text(size = 16, color = "black")
  ) +
  guides(fill = "none") +
  # geom_abline(slope = 1,
  #             intercept = 0,
  #             color = "#de2d26",
  #             linetype = "dashed") +
  labs(title = "MA-plot (Smarcad1 KO vs. WT)",
       x = "log2(baseMean)",
       y = "log2FoldChange)")

# arrange
ggarrange(ma, scatter)

ggsave(
  glue("{result_folder}TEtranscripts-IAPEz_MAplot.pdf"),
  plot = last_plot(),
  width = 14,
  height = 7,
  dpi = 300,
)

te_all_iapez_kowt_mod = te_all_iapez_kowt %>% separate(V1, sep = ":", into = "gene") 
te_all_iapez_kowt_mod = te_all_iapez_kowt_mod %>% mutate(IAPez = ifelse(gene %in% sign_pos_iapez, "upreg. IAPEz-int", ""))
te_all_iapez_kowt_mod = te_all_iapez_kowt_mod %>% mutate(IAPez = case_when(gene %in% sign_alt_iapez ~ "alt. IAPEz-int",
                                                                                          gene %in% unalt_iapez ~ "unaltered IAPEz-int", 
                                                                                          .default = "non-IAPEz"))

# MA plot, rescue
ma_kowt = ggplot(te_all_iapez_kowt_mod, aes(x = log2(baseMean), y = log2FoldChange, color = IAPez)) +
  ggrastr::geom_point_rast() +
  xlim(0, 15) +
  ylim(-10, 10) +
  scale_color_manual(values = c("#de2d26", "#f0f0f0", "orange")) +
  theme_minimal() +
  theme(
    text = element_text(size = 18),
    legend.text = element_text(size = 18),
    plot.title = element_text(size = 18, face = "bold"),
    axis.text.x = element_text(size = 16, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.title.x = element_text(size = 16, color = "black"),
    axis.title.y = element_text(size = 16, color = "black")
  ) +
  guides(fill = "none") +
  # geom_abline(slope = 1,
  #             intercept = 0,
  #             color = "#de2d26",
  #             linetype = "dashed") +
  labs(title = "MA-plot (Smarcad1 KO vs. Smarcad1 KOWT)",
       x = "log2(baseMean)",
       y = "log2FoldChange)")

# arrange
ggarrange(ma, ma_kowt)

ggsave(
  glue("{result_folder}TEtranscripts-IAPEz_MAplot_KOWT.pdf"),
  plot = last_plot(),
  width = 14,
  height = 7,
  dpi = 300,
)

# MA plot, protein coding genes
te_all_iapez_mod3 = te_all_iapez %>% separate(V1, sep = ":", into = "gene") 
te_all_iapez_mod3 = drop_na(te_all_iapez_mod3)

sign_alt_genes = te_all_iapez_mod3 %>% dplyr::filter(abs(log2FoldChange) > 2 &
                                                       padj < 0.05) %>% 
  dplyr::filter(!str_detect(gene, "IAPEz")) %>% pull(gene)

te_all_iapez_mod3 = te_all_iapez_mod3 %>% mutate(alteration = case_when(gene %in% sign_alt_genes ~ "alt. non-rep. gene",
                                                                           .default = "")) %>% arrange(alteration)

ma3 = ggplot(te_all_iapez_mod3, aes(x = log2(baseMean), y = log2FoldChange, color = alteration)) +
  ggrastr::geom_point_rast() + 
  xlim(0, 15) +
  ylim(-10, 10) +
  scale_color_manual(values = c("#f0f0f0", "#2ca25f")) +
  theme_minimal() +
  theme(
    text = element_text(size = 18),
    legend.text = element_text(size = 18),
    plot.title = element_text(size = 18, face = "bold"),
    axis.text.x = element_text(size = 16, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.title.x = element_text(size = 16, color = "black"),
    axis.title.y = element_text(size = 16, color = "black")
  ) +
  guides(fill = "none", color = guide_legend(title = "")) +
  # geom_abline(slope = 1,
  #             intercept = 0,
  #             color = "#de2d26",
  #             linetype = "dashed") +
  labs(title = "MA-plot (Smarcad1 KO vs. Smarcad1 WT)",
       x = "log2(baseMean)",
       y = "log2FoldChange)")

te_all_iapez_kowt_mod3 = te_all_iapez_kowt %>% separate(V1, sep = ":", into = "gene") 
te_all_iapez_kowt_mod3 = drop_na(te_all_iapez_kowt_mod3)

te_all_iapez_kowt_mod3 = te_all_iapez_kowt_mod3 %>% mutate(alteration = case_when(gene %in% sign_alt_genes ~ "alt. non-rep. gene",
                                                                   .default = "")) %>% arrange(alteration)



ma_kowt3 = ggplot(te_all_iapez_kowt_mod3, aes(x = log2(baseMean), y = log2FoldChange, color = alteration)) +
  ggrastr::geom_point_rast() +
  xlim(0, 15) +
  ylim(-10, 10) +
  scale_color_manual(values = c("#f0f0f0", "#2ca25f")) +
  theme_minimal() +
  theme(
    text = element_text(size = 18),
    legend.text = element_text(size = 18),
    plot.title = element_text(size = 18, face = "bold"),
    axis.text.x = element_text(size = 16, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.title.x = element_text(size = 16, color = "black"),
    axis.title.y = element_text(size = 16, color = "black")
  ) +
  guides(fill = "none", color = guide_legend(title = "")) +
  # geom_abline(slope = 1,
  #             intercept = 0,
  #             color = "#de2d26",
  #             linetype = "dashed") +
  labs(title = "MA-plot (Smarcad1 KO vs. Smarcad1 KOWT)",
       x = "log2(baseMean)",
       y = "log2FoldChange)")

# arrange
ggarrange(ma3, ma_kowt3)

ggsave(
  glue("{result_folder}TEtranscripts-IAPEz_MAplot_KOWT-genes.pdf"),
  plot = last_plot(),
  width = 14,
  height = 7,
  dpi = 300,
)


# heatmaps
# DESeq2 normalized expression levels
sign_pos_iapez = volc_input_iapez %>% dplyr::filter(log2FoldChange > 2 &
                                                      padj < 0.05) %>% 
  arrange(desc(log2FoldChange)) %>% pull(te_name) 

sign_pos_iapez_norm = te_iapez_norm %>% dplyr::filter(gene %in% sign_pos_iapez)
sign_pos_iapez_norm_kowt = te_iapez_norm_kowt %>% dplyr::filter(gene %in% sign_pos_iapez) %>% dplyr::select(gene, KOWT_mean)
sign_pos_norm = sign_pos_iapez_norm %>% inner_join(., sign_pos_iapez_norm_kowt, by = "gene") %>% dplyr::select(gene, WT_mean, Smarcad1_KO_mean, KOWT_mean)
colnames(sign_pos_norm) = c("gene", "WT", "Smarcad1 KO", "Smarcad1 KOWT")
sign_pos_norm = sign_pos_norm %>% pivot_longer(cols = "WT":"Smarcad1 KOWT", names_to = "condition", values_to = "norm_expr") 

y_order = sign_pos_norm %>% dplyr::filter(condition == "WT" | condition == "Smarcad1 KO") %>% 
  pivot_wider(names_from = condition, values_from = "norm_expr") %>% mutate(diff = `Smarcad1 KO` - WT) %>% 
  arrange(diff) %>% pull(gene)
y_order = factor(sign_pos_norm$gene, levels = unique(y_order))
x_order = factor(sign_pos_norm$condition, levels = c("WT", "Smarcad1 KO", "Smarcad1 KOWT"))

hm = ggplot(sign_pos_norm, aes(x = x_order, y = y_order, fill = log2(norm_expr))) +
  geom_tile(color = "white",
            lwd = 0.5,
            linetype = 1) +
  scale_fill_gradient2(
    low = "#075AFF",
    mid = "#FFFFCC",
    high = "#FF0000",
    midpoint = 4,
    limits = c(0, 8)
  ) +
  xlab(label = "") +
  ylab(label = "upregulated IAPEz-int upon Smarcad1 KO") +
  labs(fill = "log(norm. expr.)", title = "") +
  theme(
    axis.text.x = element_text(
      color = "black",
      size = 10,
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = element_text(color = "black", size = 10),
    axis.title = element_text(size = 14)
  ) +
  coord_fixed()
print(hm)

ggsave(
  glue("{result_folder}upreg_IAPez_upon_Smarcad1KO-norm_expr_hm.pdf"),
  plot = last_plot(),
  width = 5,
  height = 7,
  dpi = 300,
)
