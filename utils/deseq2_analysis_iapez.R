suppressPackageStartupMessages({
  library("tidyverse")
  library("data.table")
  library("glue")
  library("elsasserlib")
  library("apeglm")
  library("GenomeInfoDb")
  library("ggpubr")
  library("wigglescout")
  library("ggpointdensity")
})

# folders
result_folder = "../results/bw_diff/"
data_folder = "../data/"
bigwig_folder = "../data/ATAC-Seq_SMARCAD1_KO/bigwig/"
rnaseq_folder = "../results/rna_seq_deseq2/"

# for subsetting
bed = "../data/bed/repmasker_IAPEz-int.mm10.bed"

wt_bws = c(
  glue("{bigwig_folder}SMARCAD1_WT_R1.mLb.clN.bigWig"),
  glue("{bigwig_folder}SMARCAD1_WT_R2.mLb.clN.bigWig")
)

ko_wt_bws = c(
  glue("{bigwig_folder}SMARCAD1_KO_WT_R1.mLb.clN.bigWig"),
  glue("{bigwig_folder}SMARCAD1_KO_WT_R2.mLb.clN.bigWig")
)

ko_bws = c(
  glue("{bigwig_folder}SMARCAD1_KO_R1.mLb.clN.bigWig"),
  glue("{bigwig_folder}SMARCAD1_KO_R2.mLb.clN.bigWig")
)

# DESeq2 wrapper
ko_wt_iapez_fc = bw_bed_diff_analysis(
  bwfiles_c1 = ko_bws,
  bwfiles_c2 = wt_bws,
  bed = bed,
  label_c1 = "SMARCAD1_KO",
  label_c2 = "SMARCAD1_WT"
)

ko_wt_iapez_fc = as.data.frame(ko_wt_iapez_fc)
ko_wt_iapez_fc = tibble(ko_wt_iapez_fc)
ko_wt_iapez_fc = ko_wt_iapez_fc %>% mutate(contrast = "SMARCAD1 KO vs. SMARCAD1")

kowt_wt_iapez_fc = bw_bed_diff_analysis(
  bwfiles_c1 = ko_wt_bws,
  bwfiles_c2 = wt_bws,
  bed = bed,
  label_c1 = "SMARCAD1_KO_WT",
  label_c2 = "SMARCAD1_WT"
)

kowt_wt_iapez_fc = as.data.frame(kowt_wt_iapez_fc)
kowt_wt_iapez_fc = tibble(kowt_wt_iapez_fc)
kowt_wt_iapez_fc = kowt_wt_iapez_fc %>% mutate(contrast = "SMARCAD1 KO - rescued vs. SMARCAD1")

iapez_fc = rbind(kowt_wt_iapez_fc, ko_wt_iapez_fc)

boxplot = iapez_fc %>% ggplot(., aes(
  x = contrast,
  y = log2FoldChange
)) +
  geom_boxplot(fill = "#fc9272") +
  theme_minimal() +
  geom_density(alpha = 1.0) +
  labs(title = "ATAC-Seq - Diff. accessible regions over IAPez ",
       x = "",
       y = "DESeq2 fold change",
       fill = "") +
  ylim(-1.5, 1.5) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 13),
    axis.text.x = element_text(size = 8, color = "black", angle = 45, vjust = 1, hjust = 1),
    axis.text.y = element_text(size = 10)
  ) +
  stat_compare_means(label.y = 1, label.x = 1.3, size = 3)
boxplot

ggsave(
  glue("{result_folder}ATAC-Seq_SMARCAD1-diff_acc_over_IAPez.png"),
  plot = boxplot,
  width = 5,
  height = 5,
  dpi = 300,
)

rep_masker = fread(bed)
rep_masker = GRanges(
  seqnames = rep_masker$V1,
  ranges = IRanges(start = rep_masker$V2,
                   end = rep_masker$V3))

# RNA-Seq RPGC between conditions
# wild-type Vs. Smarcad1 KO-WT
kowt_vs_wt = plot_bw_bins_scatter(
  "../data/RNA-Seq/bigwig_rpgc/Ph001_RPGC.bigWig",
  "../data/RNA-Seq/bigwig_rpgc/Ph008_RPGC.bigwig",
  bin_size = 2000,
  verbose = FALSE,
  remove_top = 0.01,
  genome = "mm10",
  selection = rep_masker
) +
  geom_pointdensity(show.legend = FALSE) +
  scale_colour_gradient(low = "#f03b20",
                        high = "yellow") +
  labs(colour = " ", title = "RNA-Seq - RPGC") +
  xlab("Smarcad1 WT") +
  ylab("Smarcad1 KO-WT") +
  theme_minimal() +
  xlim(0, 6) +
  ylim(0, 6) +
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
    label.x = 3,
    label.y = 5,
    size = 3
  )
  
# wild-type Vs. Smarcad1 KO-WT
wt_vs_wt = plot_bw_bins_scatter(
  "../data/RNA-Seq/bigwig_rpgc/Ph001_RPGC.bigWig",
  "../data/RNA-Seq/bigwig_rpgc/Ph006_RPGC.bigwig",
  bin_size = 2000,
  verbose = FALSE,
  remove_top = 0.01,
  genome = "mm10",
  selection = rep_masker
) +
  geom_pointdensity(show.legend = FALSE) +
  scale_colour_gradient(low = "#f03b20",
                        high = "yellow") +
  labs(colour = " ", title = "RNA-Seq - RPGC") +
  xlab("Smarcad1 WT, rep 1") +
  ylab("Smarcad1 WT, rep 2") +
  theme_minimal() +
  xlim(0, 6) +
  ylim(0, 6) +
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
    label.x = 3,
    label.y = 5,
    size = 3
  )
# wild-type Vs. Smarcad1 KO-WT
ko_vs_wt = plot_bw_bins_scatter(
  "../data/RNA-Seq/bigwig_rpgc/Ph001_RPGC.bigWig",
  "../data/RNA-Seq/bigwig_rpgc/Ph007_RPGC.bigwig",
  bin_size = 2000,
  verbose = FALSE,
  remove_top = 0.01,
  genome = "mm10",
  selection = rep_masker
) +
  geom_pointdensity(show.legend = FALSE) +
  scale_colour_gradient(low = "#f03b20",
                        high = "yellow") +
  labs(colour = " ", title = "RNA-Seq - RPGC") +
  xlab("Smarcad1 WT") +
  ylab("Smarcad1 KO") +
  theme_minimal() +
  xlim(0, 6) +
  ylim(0, 6) +
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
    label.x = 3,
    label.y = 5,
    size = 3
  )

kod425_vs_wt = plot_bw_bins_scatter(
  "../data/RNA-Seq/bigwig_rpgc/Ph001_RPGC.bigWig",
  "../data/RNA-Seq/bigwig_rpgc/Ph005_RPGC.bigwig",
  bin_size = 2000,
  verbose = FALSE,
  remove_top = 0.01,
  genome = "mm10",
  selection = rep_masker
) +
  geom_pointdensity(show.legend = FALSE) +
  scale_colour_gradient(low = "#f03b20",
                        high = "yellow") +
  labs(colour = " ", title = "RNA-Seq - RPGC") +
  xlab("Smarcad1 WT") +
  ylab("Smarcad1 KO d425") +
  theme_minimal() +
  xlim(0, 6) +
  ylim(0, 6) +
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
    label.x = 3,
    label.y = 5,
    size = 3
  )

scs = list(wt_vs_wt, kowt_vs_wt, ko_vs_wt, kod425_vs_wt)
ggarrange(plotlist = scs)

ggsave(
  plot = last_plot(),
  glue("{rnaseq_folder}bigwig_rpgc_scatters.pdf"),
  width = 10,
  height = 8
)

ggsave(
  plot = last_plot(),
  glue("{rnaseq_folder}bigwig_rpgc_scatters.png"),
  width = 10,
  height = 8,
  dpi = 300
)
