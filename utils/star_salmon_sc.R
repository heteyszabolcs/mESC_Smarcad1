# packages
suppressPackageStartupMessages({
  library("DESeq2")
  library("data.table")
  library("tidyverse")
  library("glue")
  library("ggplot2")
  library("ggpointdensity")
  library("ggpubr")
})

# folders
result_folder = "../results/rna_seq_deseq2/"

fc_1 = fread(glue("{result_folder}DESeq2_res-kowt_vs_wt_salmon.tsv"))
salmon_sign_genes = fc_1 %>% dplyr::filter(abs(log2FoldChange > 1) & padj < 0.05) %>% pull(gene_name)
fc_2 = fread(glue("{result_folder}DESeq2_res-kowt_vs_wt.tsv"))
starrsem_sign_genes = fc_2 %>% dplyr::filter(abs(log2FoldChange > 1) & padj < 0.05) %>% pull(gene_name)

both = intersect(starrsem_sign_genes, salmon_sign_genes)

fc_joined = fc_1 %>% inner_join(., fc_2, by = c("gene_name" = "gene_name")) %>%
  dplyr::select(starts_with("log2FoldChange"), gene_name) %>% drop_na() %>%
  rename("KOWT_vs_WT_Salmon" = log2FoldChange.x,
         "KOWT_vs_WT_STAR_RSEM" = log2FoldChange.y)

sc1 = ggplot(
  data = fc_joined,
  mapping = aes(x = KOWT_vs_WT_Salmon, y = KOWT_vs_WT_STAR_RSEM)
) +
  geom_pointdensity(show.legend = FALSE) +
  scale_colour_gradient(low = "#f03b20",
                        high = "yellow") +
  labs(colour = " ", title = "Salmon vs. STAR-RSEM fold changes") +
  xlab("KO-WT vs WT - Salmon") +
  ylab("KO-WT vs WT - STAR-RSEM") +
  theme_minimal() +
  # xlim(-10, 3) +
  # ylim(-10, 3) +
  theme(
    axis.line = element_line(size = 0.5, colour = "black"),
    text = element_text(size = 20),
    plot.title = element_text(size = 15),
    axis.text.x = element_text(
      size = 13,
      color = "black",
      angle = 0,
      hjust = 1,
      colour = "black"
    ),
    axis.text.y = element_text(size = 13, color = "black")
  ) +
  stat_cor(
    method = "pearson",
    label.x = -20.0,
    label.y = 20,
    size = 6
  )
sc1

fc_joined = fc_joined %>% 
  mutate(sign = ifelse(gene_name %in% salmon_sign_genes, "sign. in Salmon", "non-sign")) %>%
  mutate(sign = ifelse(gene_name %in% starrsem_sign_genes, "sign. in STAR-RSEM", sign)) %>% 
  mutate(sign = ifelse(gene_name %in% both, "sign. in both", sign))
  
sc2 = ggplot(
  data = fc_joined, aes(x = KOWT_vs_WT_Salmon, y = KOWT_vs_WT_STAR_RSEM, color = sign)
) +
  geom_point() +
  scale_color_manual(values = c("#bdbdbd", "#31a354", "#f03b20", "#2b8cbe")) +
  labs(color = " ", title = "Salmon vs. STAR-RSEM fold changes") +
  xlab("KO-WT vs WT - Salmon") +
  ylab("KO-WT vs WT - STAR-RSEM") +
  theme_minimal() +
  theme(
    axis.line = element_line(size = 0.5, colour = "black"),
    text = element_text(size = 20),
    plot.title = element_text(size = 15),
    axis.text.x = element_text(
      size = 13,
      color = "black",
      angle = 0,
      hjust = 1,
      colour = "black"
    ),
    axis.text.y = element_text(size = 13, color = "black")
  ) 
sc2  

scs = list(sc1, sc2)
ggarrange(plotlist = scs)

ggsave(
  plot = last_plot(),
  glue("{result_folder}STAR-RSEM_vs._Salmon_sc.pdf"),
  width = 12,
  height = 6
)

ggsave(
  plot = last_plot(),
  glue("{result_folder}STAR-RSEM_vs._Salmon_sc.png"),
  width = 12,
  height = 6,
  dpi = 300
)






