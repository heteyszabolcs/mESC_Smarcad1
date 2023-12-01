suppressPackageStartupMessages({
  library("glue")
  library("tidyverse")
  library("data.table")
  library("wigglescout")
  library("GenomicRanges")
  library("TxDb.Mmusculus.UCSC.mm10.knownGene")
  library("ggrepel")
  library("ggpubr")
  library("RColorBrewer")
})

# annotation script
source("C:/Szabolcs/Karolinska/Data/scripts/annotation.R")

# result folder
result_folder = "../results/wigglescout/"

# folders
bw = list.files("../data/MINUTE/bw/", full.names = TRUE)
dip_bws = bw[grep(x = bw, pattern = "DIP")]
dip_bws = dip_bws[grep(x = dip_bws, pattern = ".mm10.scaled.bw")]
k9me3_bws = bw[grep(x = bw, pattern = "H3K9me3")]
k9me3_bws = k9me3_bws[grep(x = k9me3_bws, pattern = ".mm10.scaled.bw")]
dip_mapq_bws = bw[grep(x = bw, pattern = "DIP")]
dip_mapq_bws = dip_mapq_bws[grep(x = dip_mapq_bws, pattern = "*.mm10.scaled.mapq.bw")]
k9me3_mapq_bws = bw[grep(x = bw, pattern = "H3K9me3")]
k9me3_mapq_bws = k9me3_mapq_bws[grep(x = k9me3_mapq_bws, pattern = "*.mm10.scaled.mapq.bw")]
rnaseq_bws = list.files("../data/RNA-Seq/STAR_output/", full.names = TRUE, pattern = "*.bigwig")
atac_bws = list.files("../data/ATAC-Seq_SMARCAD1_KO/bigwig/", full.names = TRUE)
atac_bws = atac_bws[grep(x = atac_bws, pattern = "SMARCAD1_KO_R|SMARCAD1_WT")]

# regions
ervk = fread("../data/bed/UCSC_RepeatMasker_ERVK_mm10.bed")
ervk$V4 = "ERVK"
ervk = GRanges(
  seqnames = ervk$V1,
  ranges = IRanges(
    start = ervk$V2,
    end = ervk$V3,
    names = ervk$V4
  )
)

regions_top_fc_upon_Smarcad1KO = fread("../data/bed/DESeq2_repeats-res-ko_vs_wt-top_fc_regions.bed")
regions_top_fc_upon_Smarcad1KO$V4 = "tops"
regions_top_fc_upon_Smarcad1KO = GRanges(
  seqnames = regions_top_fc_upon_Smarcad1KO$V1,
  ranges = IRanges(
    start = regions_top_fc_upon_Smarcad1KO$V2,
    end = regions_top_fc_upon_Smarcad1KO$V3,
    names = regions_top_fc_upon_Smarcad1KO$V4
  )
)

rltr4 = fread("../data/bed/UCSC_RepeatMasker_RLTR4_mm10.bed")
rltr4$V4 = "RLTR4"
rltr4 = GRanges(
  seqnames = rltr4$V1,
  ranges = IRanges(
    start = rltr4$V2,
    end = rltr4$V3,
    names = rltr4$V4
  )
)

rltrs = fread("../data/bed/UCSC_RepeatMasker_RLTRs_mm10.bed")
rltrs$V4 = "RLTRs"
rltrs = GRanges(
  seqnames = rltrs$V1,
  ranges = IRanges(
    start = rltrs$V2,
    end = rltrs$V3,
    names = rltrs$V4
  )
)

imprinteds = fread("../data/bed/mm10_imprinted_genes.bed")
imprinteds$V4 = "imprint"
imprinteds = GRanges(
  seqnames = imprinteds$V1,
  ranges = IRanges(
    start = imprinteds$V2,
    end = imprinteds$V3,
    names = imprinteds$V4
  )
)

high_expr = fread("../data/bed/top200_highly_expr_gene-wt.bed")
high_expr$V4 = "high_expr"
high_expr = GRanges(
  seqnames = high_expr$V1,
  ranges = IRanges(
    start = high_expr$V2,
    end = high_expr$V3,
    names = high_expr$V4
  )
)

upreg_iapez_upon_smarcad1ko = fread("../data/bed/RNA-Seq-sign_upreg_IAPez-int_regions_uponSmarcad1_KO.bed")
upreg_iapez_upon_smarcad1ko$V7 = "IAPEz"
upreg_iapez_upon_smarcad1ko = GRanges(
  seqnames = upreg_iapez_upon_smarcad1ko$V1,
  ranges = IRanges(
    start = upreg_iapez_upon_smarcad1ko$V2,
    end = upreg_iapez_upon_smarcad1ko$V3,
    names = upreg_iapez_upon_smarcad1ko$V7
  )
)

## functions
# heatmap function
plot_heatmap = function(bigwig, axis_label, loci, title, zmax = 200, palette = "Greens") {
  plot_bw_heatmap(
    mode = "stretch",
    zmax = zmax,
    bigwig,
    cmap = palette,
    loci = loci,
    verbose = FALSE,
    upstream = 1000,
    downstream = 1000,
  ) +
    labs(title = title,
         x = "",
         y = axis_label) +
    theme(
      text = element_text(size = 8),
      plot.title = element_text(size = 7),
      axis.text.x = element_text(size = 8, color = "black")
    )
}

# function for calculate signal differences between bigwigs
calc_bigwig_fc = function(bigwig,
                          bg_bigwig,
                          subset,
                          output_file) {
  # wigglescout's bw_loci
  fc = bw_loci(bigwig, bg_bwfiles = bg_bigwig, loci = subset)
  fc = as.data.frame(fc)
  fc = fc[!is.infinite(fc[, 6]),]
  fc = fc[!is.nan(fc[, 6]), ]
  
  feature = colnames(fc)[6]
  
  # annotation
  fc = mm10_annotation(
    regions = fc,
    start_col = "start",
    end_col = "end",
    seqname_col = "seqnames",
    feature_1 = feature,
    feature_2 = feature
  )
  
  fc = fc %>%
    dplyr::select(
      seqnames,
      start,
      end,
      fold_change = feature_1,
      gene_symbol = SYMBOL,
      distanceToTSS
    ) %>%
    dplyr::filter(fold_change >= 2)
  
  #write_tsv(fc, glue("{result_folder}{output_file}"))
  
  return(fc)
}

# function wigglescout profile plot
plot_profile_plot = function(bigwig,
                        gr,
                        labels,
                        loci,
                        title,
                        gr_title,
                        ymax = 10,
                        colors) {
  plot = plot_bw_profile(
    upstream = 3000,
    downstream = 3000,
    bin_size = 50,
    bigwig,
    loci = gr,
    labels = c(labels),
    verbose = FALSE, default_na = 0
  ) +
    ylim(0, ymax) +
    theme(
      text = element_text(size = 16),
      legend.text = element_text(size = 14),
      plot.title = element_text(size = 15, face = "bold"),
      axis.text.x = element_text(size = 12, color = "black"),
      axis.text.y = element_text(size = 12, color = "black"),
      axis.title.x = element_text(size = 14, color = "black"),
      axis.title.y = element_text(size = 14, color = "black")
    ) +
    scale_color_manual(values = colors) +
    guides(fill = "none") +
    labs(title = title,
         x = gr_title,
         y = "RPGC")
  print(plot)
  return(plot)
}

# function for boxplot with statistics
plot_boxplot = function(bigwig,
                        gr,
                        list_for_comparisons = list(
                          c(
                            "DIPh_KO_P6_rep1.mm10.scaled",
                            "DIPh_WT_P6_rep2.mm10.scaled"
                          ),
                          c(
                            "DIPh_KO_P6_rep1.mm10.scaled",
                            "DIPh_WT_P6_rep1.mm10.scaled"
                          ),
                          c(
                            "DIPh_KO_P6_rep1.mm10.scaled",
                            "DIPh_KO_P6_rep2.mm10.scaled"
                          )
                        ),
                        plot_title,
                        colors,
                        ymax = 15, label.y = c(5, 7, 9)) {
  read_dens = bw_loci(bigwig, loci = gr)
  read_dens = as.data.frame(read_dens)
  col1 = colnames(read_dens)[6]
  coln = colnames(read_dens)[ncol(read_dens)]
  read_dens = read_dens %>% pivot_longer(.,
                                         all_of(col1):all_of(coln),
                                         names_to = "minute",
                                         values_to = "read_dens") %>%
    dplyr::select(minute, read_dens)
  
  plot = ggplot(read_dens, aes(x = minute, y = read_dens, fill = minute)) +
    geom_boxplot(color = "black") +
    scale_fill_manual(values = colors) +
    guides(fill = "none") +
    ylim(0, ymax) +
    labs(
      title = plot_title,
      x = "",
      y = "RPGC",
      fill = ""
    ) +
    theme_classic() +
    theme(
      text = element_text(size = 12),
      plot.title = element_text(size = 17),
      axis.text.x = element_text(
        size = 17,
        color = "black",
        angle = 90,
        vjust = 0.5,
        hjust = 1
      ),
      axis.text.y = element_text(size = 17, color = "black")
    ) + stat_compare_means(
      comparisons = list_for_comparisons,
      label = "p.signif",
      label.y = label.y,
      tip.length = 0.005
    )
  print(plot)
  return(plot)
}

### repetitive regions with higher expression upon Smarcad1 KO
# methylation status after Smarcad1 KO
p1 = plot_heatmap(bigwig = dip_bws[1], loci = regions_top_fc_upon_Smarcad1KO, 
             axis_label = "ERVs with highest fold change upon Smarcad1 KO",
             title = "MeDIP Smarcad1 KO, rep 1")
p2 = plot_heatmap(bigwig = dip_bws[2], loci = regions_top_fc_upon_Smarcad1KO, 
             axis_label = "ERVs with highest fold change upon Smarcad1 KO",
             title = "MeDIP Smarcad1 KO, rep 2")
p3 = plot_heatmap(bigwig = dip_bws[4], loci = regions_top_fc_upon_Smarcad1KO, 
             axis_label = "ERVs with highest fold change upon Smarcad1 KO",
             title = "MeDIP Smarcad1 WT, rep 1")
p4 = plot_heatmap(bigwig = dip_bws[5], loci = regions_top_fc_upon_Smarcad1KO, 
             axis_label = "ERVs with highest fold change upon Smarcad1 KO",
             title = "MeDIP Smarcad1 WT, rep 2")

ggsave(
  plot = ggarrange(plotlist = list(p1, p2, p3, p4), nrow = 2, ncol = 2),
  glue("{result_folder}top_expressing_ERVs-MeDIP_signal_Smarcad1_KO.png"),
  width = 10,
  height = 10,
  dpi = 500,
)
ggsave(
  plot = ggarrange(plotlist = list(p1, p2, p3, p4), nrow = 2, ncol = 2),
  glue("{result_folder}top_expressing_ERVs-MeDIP_signal_Smarcad1_KO.pdf"),
  width = 10,
  height = 10
)


# H3K9me3 status after Smarcad1 KO
p1 = plot_heatmap(bigwig = k9me3_bws[1], loci = regions_top_fc_upon_Smarcad1KO, 
             axis_label = "ERVs with highest fold change upon Smarcad1 KO",
             title = "H3K9me3 Smarcad1 KO, rep 1", zmax = 100)
p2 = plot_heatmap(bigwig = k9me3_bws[2], loci = regions_top_fc_upon_Smarcad1KO, 
                  axis_label = "ERVs with highest fold change upon Smarcad1 KO",
                  title = "H3K9me3 Smarcad1 KO, rep 2", zmax = 100)
p3 = plot_heatmap(bigwig = k9me3_bws[5], loci = regions_top_fc_upon_Smarcad1KO, 
             axis_label = "ERVs with highest fold change upon Smarcad1 KO",
             title = "H3K9me3 Smarcad1 WT, rep 1", zmax = 100)
p4 = plot_heatmap(bigwig = k9me3_bws[6], loci = regions_top_fc_upon_Smarcad1KO, 
                  axis_label = "ERVs with highest fold change upon Smarcad1 KO",
                  title = "H3K9me3 Smarcad1 WT, rep 2", zmax = 100)

ggsave(
  plot = ggarrange(plotlist = list(p1, p2, p3, p4), nrow = 2, ncol = 2),
  glue("{result_folder}top_expressing_ERVs-H3K9me3_signal_Smarcad1_KO.png"),
  width = 10,
  height = 10,
  dpi = 500,
)
ggsave(
  plot = ggarrange(plotlist = list(p1, p2, p3, p4), nrow = 2, ncol = 2),
  glue("{result_folder}top_expressing_ERVs-H3K9me3_signal_Smarcad1_KO.pdf"),
  width = 10,
  height = 10
)

### ERVKs with elevated demethylation, after Smarcad1 KO
high_demethly_ervk = calc_bigwig_fc(bigwig = "../data/MINUTE/bw/DIPh_WT_P6_rep1.mm10.scaled.bw",
                            bg_bigwig = "../data/MINUTE/bw/DIPh_KO_P6_rep1.mm10.scaled.bw",
                            subset = ervk)

highest_demethly_ervks = high_demethly_ervk %>% drop_na() %>% 
  top_n(wt = fold_change, 50) %>% dplyr::select(seqnames, start, end)
highest_demethly_ervks$V4 = "ERVK"
highest_demethly_ervks = GRanges(
  seqnames = highest_demethly_ervks$seqnames,
  ranges = IRanges(
    start = highest_demethly_ervks$start,
    end = highest_demethly_ervks$end,
    names = highest_demethly_ervks$V4
  )
)

# methylation heatmaps
p1 = plot_heatmap(bigwig = dip_bws[1], loci = highest_demethly_ervks, zmax = 100,
             axis_label = "ERVs with highest demethylation upon Smarcad1 KO",
             title = "MeDIP Smarcad1 KO, rep 1")
p2 = plot_heatmap(bigwig = dip_bws[2], loci = highest_demethly_ervks, zmax = 100,
             axis_label = "ERVs with highest demethylation upon Smarcad1 KO",
             title = "MeDIP Smarcad1 KO, rep 2")
p3 = plot_heatmap(bigwig = dip_bws[4], loci = highest_demethly_ervks, zmax = 100,
             axis_label = "ERVs with highest demethylation upon Smarcad1 KO",
             title = "MeDIP Smarcad1 WT, rep 1")
p4 = plot_heatmap(bigwig = dip_bws[5], loci = highest_demethly_ervks, zmax = 100,
             axis_label = "ERVs with highest demethylation upon Smarcad1 KO",
             title = "MeDIP Smarcad1 WT, rep 2")

ggsave(
  plot = ggarrange(plotlist = list(p1, p2, p3, p4), nrow = 2, ncol = 2),
  glue("{result_folder}highest_demethy_ERVs-MeDIP_signal_Smarcad1_KO.png"),
  width = 10,
  height = 10,
  dpi = 500,
)
ggsave(
  plot = ggarrange(plotlist = list(p1, p2, p3, p4), nrow = 2, ncol = 2),
  glue("{result_folder}highest_demethy_ERVs-MeDIP_signal_Smarcad1_KO.pdf"),
  width = 10,
  height = 10
)

# H3K9me3 heatmaps
p1 = plot_heatmap(bigwig = k9me3_bws[1], loci = highest_demethly_ervks, zmax = 100,
             axis_label = "ERVs with highest demethylation upon Smarcad1 KO",
             title = "H3K9me3 Smarcad1 KO, rep 1")
p2 = plot_heatmap(bigwig = k9me3_bws[2], loci = highest_demethly_ervks, zmax = 100,
             axis_label = "ERVs with highest demethylation upon Smarcad1 KO",
             title = "H3K9me3 Smarcad1 KO, rep 2")
p3 = plot_heatmap(bigwig = k9me3_bws[5], loci = highest_demethly_ervks, zmax = 100,
             axis_label = "ERVs with highest demethylation upon Smarcad1 KO",
             title = "H3K9me3 Smarcad1 WT, rep 1")
p4 = plot_heatmap(bigwig = k9me3_bws[6], loci = highest_demethly_ervks, zmax = 100,
             axis_label = "ERVs with highest demethylation upon Smarcad1 KO",
             title = "H3K9me3 Smarcad1 WT, rep 2")

ggsave(
  plot = ggarrange(plotlist = list(p1, p2, p3, p4), nrow = 2, ncol = 2),
  glue("{result_folder}highest_demethy_ERVs-H3K9me3_signal_Smarcad1_KO.png"),
  width = 10,
  height = 10,
  dpi = 500,
)
ggsave(
  plot = ggarrange(plotlist = list(p1, p2, p3, p4), nrow = 2, ncol = 2),
  glue("{result_folder}highest_demethy_ERVs-H3K9me3_signal_Smarcad1_KO.pdf"),
  width = 10,
  height = 10
)


### RLTRs
# RLTR4
# methylation heatmaps
p1 = plot_heatmap(bigwig = dip_bws[1], loci = rltr4, zmax = 100,
                  axis_label = "RLTR4 elements",
                  title = "MeDIP Smarcad1 KO, rep 1", palette = "Blues")
p2 = plot_heatmap(bigwig = dip_bws[2], loci = rltr4, zmax = 100,
                  axis_label = "RLTR4 elements",
                  title = "MeDIP Smarcad1 KO, rep 2", palette = "Blues")
p3 = plot_heatmap(bigwig = dip_bws[4], loci = rltr4, zmax = 100,
                  axis_label = "RLTR4 elements",
                  title = "MeDIP Smarcad1 WT, rep 1", palette = "Blues")
p4 = plot_heatmap(bigwig = dip_bws[5], loci = rltr4, zmax = 100,
                  axis_label = "RLTR4 elements",
                  title = "MeDIP Smarcad1 WT, rep 2", palette = "Blues")

ggsave(
  plot = ggarrange(plotlist = list(p1, p2, p3, p4), nrow = 2, ncol = 2),
  glue("{result_folder}RLTR4-MeDIP_signal_Smarcad1_KO.png"),
  width = 10,
  height = 10,
  dpi = 500,
)
ggsave(
  plot = ggarrange(plotlist = list(p1, p2, p3, p4), nrow = 2, ncol = 2),
  glue("{result_folder}RLTR4-MeDIP_signal_Smarcad1_KO.pdf"),
  width = 10,
  height = 10
)

# H3K9me3 heatmaps
p1 = plot_heatmap(bigwig = k9me3_bws[1], loci = rltr4, zmax = 100,
                  axis_label = "RLTR4 elements",
                  title = "H3K9me3 Smarcad1 KO, rep 1", palette = "Blues")
p2 = plot_heatmap(bigwig = k9me3_bws[2], loci = rltr4, zmax = 100,
                  axis_label = "RLTR4 elements",
                  title = "H3K9me3 Smarcad1 KO, rep 2", palette = "Blues")
p3 = plot_heatmap(bigwig = k9me3_bws[5], loci = rltr4, zmax = 100,
                  axis_label = "RLTR4 elements",
                  title = "H3K9me3 Smarcad1 WT, rep 1", palette = "Blues")
p4 = plot_heatmap(bigwig = k9me3_bws[6], loci = rltr4, zmax = 100,
                  axis_label = "RLTR4 elements",
                  title = "H3K9me3 Smarcad1 WT, rep 2", palette = "Blues")

ggsave(
  plot = ggarrange(plotlist = list(p1, p2, p3, p4), nrow = 2, ncol = 2),
  glue("{result_folder}RLTR4-H3K9me3_signal_Smarcad1_KO.png"),
  width = 10,
  height = 10,
  dpi = 500,
)
ggsave(
  plot = ggarrange(plotlist = list(p1, p2, p3, p4), nrow = 2, ncol = 2),
  glue("{result_folder}RLTR4-H3K9me3_signal_Smarcad1_KO.pdf"),
  width = 10,
  height = 10
)

# all RLTR elements
# methylation heatmaps
p1 = plot_heatmap(bigwig = dip_bws[1], loci = rltrs, zmax = 50,
                  axis_label = "all RLTR elements",
                  title = "MeDIP Smarcad1 KO, rep 1", palette = "Blues")
p2 = plot_heatmap(bigwig = dip_bws[2], loci = rltrs, zmax = 50,
                  axis_label = "all RLTR elements",
                  title = "MeDIP Smarcad1 KO, rep 2", palette = "Blues")
p3 = plot_heatmap(bigwig = dip_bws[4], loci = rltrs, zmax = 50,
                  axis_label = "all RLTR elements",
                  title = "MeDIP Smarcad1 WT, rep 1", palette = "Blues")
p4 = plot_heatmap(bigwig = dip_bws[5], loci = rltrs, zmax = 50,
                  axis_label = "all RLTR elements",
                  title = "MeDIP Smarcad1 WT, rep 2", palette = "Blues")

ggsave(
  plot = ggarrange(plotlist = list(p1, p2, p3, p4), nrow = 2, ncol = 2),
  glue("{result_folder}all_RLTRs-MeDIP_signal_Smarcad1_KO.png"),
  width = 10,
  height = 10,
  dpi = 500,
)
ggsave(
  plot = ggarrange(plotlist = list(p1, p2, p3, p4), nrow = 2, ncol = 2),
  glue("{result_folder}all_RLTRs-MeDIP_signal_Smarcad1_KO.pdf"),
  width = 10,
  height = 10
)

# H3K9me3 heatmaps
p1 = plot_heatmap(bigwig = k9me3_bws[1], loci = rltrs, zmax = 100,
                  axis_label = "all RLTR elements",
                  title = "H3K9me3 Smarcad1 KO, rep 1", palette = "Blues")
p2 = plot_heatmap(bigwig = k9me3_bws[2], loci = rltrs, zmax = 100,
                  axis_label = "all RLTR elements",
                  title = "H3K9me3 Smarcad1 KO, rep 2", palette = "Blues")
p3 = plot_heatmap(bigwig = k9me3_bws[5], loci = rltrs, zmax = 100,
                  axis_label = "all RLTR elements",
                  title = "H3K9me3 Smarcad1 WT, rep 1", palette = "Blues")
p4 = plot_heatmap(bigwig = k9me3_bws[6], loci = rltrs, zmax = 100,
                  axis_label = "all RLTR elements",
                  title = "H3K9me3 Smarcad1 WT, rep 2", palette = "Blues")

ggsave(
  plot = ggarrange(plotlist = list(p1, p2, p3, p4), nrow = 2, ncol = 2),
  glue("{result_folder}all_RLTRs-H3K9me3_signal_Smarcad1_KO.png"),
  width = 10,
  height = 10,
  dpi = 500,
)
ggsave(
  plot = ggarrange(plotlist = list(p1, p2, p3, p4), nrow = 2, ncol = 2),
  glue("{result_folder}all_RLTRs-H3K9me3_signal_Smarcad1_KO.pdf"),
  width = 10,
  height = 10
)

# ERVKs with elevated K9me3, less methylation after Smarcad1 KO
high_k9me3_ervk = calc_bigwig_fc(bigwig = "../data/MINUTE/bw/H3K9me3_KO_P6_rep1.mm9.scaled.bw",
              bg_bigwig = "../data/MINUTE/bw/DIPh_KO_P6_rep1.mm9.scaled.bw",
              subset = ervk)

### mm10 imprinted regions
colors = c("#fc9272", "#fc9272", "#a1d99b", "#a1d99b", "#9ecae1", "#9ecae1")

plot_profile_plot(
  bigwig = k9me3_bws,
  gr = imprinteds,
  gr_title = "mm10 imprinted genes",
  title = "H3K9me3",
  labels = c(
    "KO (rep1)",
    "KO (rep2)",
    "KO, res (rep1)",
    "KO, res (rep2)",
    "WT (rep1)",
    "WT (rep2)"
  ),
  ymax = 10,
  colors = colors)

ggsave(
  glue("{result_folder}mm10_imprinted_genes-H3K9me3-profiles.png"),
  width = 10,
  height = 7,
  dpi = 500,
)
ggsave(
  glue("{result_folder}mm10_imprinted_genes-H3K9me3-profiles.pdf"),
  width = 10,
  height = 7
)

plot_profile_plot(
  bigwig = dip_bws,
  gr = imprinteds,
  gr_title = "mm10 imprinted genes",
  title = "hDIP",
  labels = c(
    "KO (rep1)",
    "KO (rep2)",
    "KO, res (rep1)",
    "WT (rep1)",
    "WT (rep2)"
  ),
  ymax = 10,
  colors = colors)

ggsave(
  glue("{result_folder}mm10_imprinted_genes-MeDIP-profiles.png"),
  width = 10,
  height = 7,
  dpi = 500,
)
ggsave(
  glue("{result_folder}mm10_imprinted_genes-MeDIP-profiles.pdf"),
  width = 10,
  height = 7
)

plot_boxplot(
  bigwig = dip_bws,
  gr = imprinteds,
  plot_title = "hDIP at imprinted genes (mm10)",
  list_for_comparisons = list(c(
    "DIPh_KO_P6_rep1.mm10.scaled",
    "DIPh_WT_P6_rep2.mm10.scaled"
  ),
  c(
    "DIPh_KO_P6_rep1.mm10.scaled",
    "DIPh_WT_P6_rep1.mm10.scaled"
  ),
  c(
    "DIPh_KO_P6_rep1.mm10.scaled",
    "DIPh_KO_P6_rep2.mm10.scaled"
  )
  ), colors = colors, ymax = 15, label.y = c(5, 7, 9))

ggsave(
  glue("{result_folder}mm10_imprinted_genes-MeDIP-bp.png"),
  width = 10,
  height = 10,
  dpi = 500,
)
ggsave(
  glue("{result_folder}mm10_imprinted_genes-MeDIP-bp.pdf"),
  width = 10,
  height = 10
)

plot_boxplot(
  bigwig = k9me3_bws,
  gr = imprinteds,
  plot_title = "H3K9me3 at imprinted genes (mm10)",
  list_for_comparisons = list(c(
    "H3K9me3_KO_P6_rep1.mm10.scaled",
    "H3K9me3_WT_P6_rep2.mm10.scaled"
  ),
  c(
    "H3K9me3_KO_P6_rep1.mm10.scaled",
    "H3K9me3_WT_P6_rep1.mm10.scaled"
  ),
  c(
    "H3K9me3_KO_P6_rep1.mm10.scaled",
    "H3K9me3_KO_P6_rep2.mm10.scaled"
  )
), colors = colors, ymax = 15, label.y = c(5, 7, 9))

ggsave(
  glue("{result_folder}mm10_imprinted_genes-H3K9me3-bp.png"),
  width = 10,
  height = 10,
  dpi = 500,
)
ggsave(
  glue("{result_folder}mm10_imprinted_genes-H3K9me3-bp.pdf"),
  width = 10,
  height = 10
)

## highly expressed genes (based on Salmon and DESeq2 normalized values) 
plot_profile_plot(
  bigwig = k9me3_bws,
  gr = high_expr,
  gr_title = "highly expressed genes (n = 200)",
  title = "H3K9me3",
  labels = c(
    "KO (rep1)",
    "KO (rep2)",
    "KO, res (rep1)",
    "KO, res (rep2)",
    "WT (rep1)",
    "WT (rep2)"
  ),
  ymax = 10,
  colors = colors)

ggsave(
  glue("{result_folder}highly_expr_genes-H3K9me3-profiles.png"),
  width = 10,
  height = 7,
  dpi = 500,
)
ggsave(
  glue("{result_folder}highly_expr_genes-H3K9me3-profiles.pdf"),
  width = 10,
  height = 7
)

plot_profile_plot(
  bigwig = dip_bws,
  gr = high_expr,
  gr_title = "highly expressed genes (n = 200)",
  title = "hDIP",
  labels = c(
    "KO (rep1)",
    "KO (rep2)",
    "KO, res (rep1)",
    "WT (rep1)",
    "WT (rep2)"
  ),
  ymax = 10,
  colors = colors)

ggsave(
  glue("{result_folder}highly_expr_genes-MeDIP-profiles.png"),
  width = 10,
  height = 7,
  dpi = 500,
)
ggsave(
  glue("{result_folder}highly_expr_genes-MeDIP-profiles.pdf"),
  width = 10,
  height = 7
)

plot_boxplot(
  bigwig = dip_bws,
  gr = high_expr,
  plot_title = "hDIP at highly expressed genes",
  list_for_comparisons = list(
    c(
      "DIPh_KO_P6_rep1.mm10.scaled",
      "DIPh_WT_P6_rep2.mm10.scaled"
    ),
    c(
      "DIPh_KO_P6_rep1.mm10.scaled",
      "DIPh_WT_P6_rep1.mm10.scaled"
    ),
    c(
      "DIPh_KO_P6_rep1.mm10.scaled",
      "DIPh_KO_P6_rep2.mm10.scaled"
    )
  ),
  colors = colors,
  ymax = 15,
  label.y = c(5, 7, 9)
)

ggsave(
  glue("{result_folder}highly_expr_genes-MeDIP-bp.png"),
  width = 10,
  height = 10,
  dpi = 500,
)
ggsave(
  glue("{result_folder}highly_expr_genes-MeDIP-bp.pdf"),
  width = 10,
  height = 10
)

### highly expressed IAPez (based on TEtranscripts analysis) 
plot_profile_plot(
  bigwig = k9me3_bws,
  gr = upreg_iapez_upon_smarcad1ko,
  gr_title = "upregulated IAPEz-int upon Smarcad1 KO",
  title = "H3K9me3",
  labels = c(
    "KO (rep1)",
    "KO (rep2)",
    "KO, res (rep1)",
    "KO, res (rep2)",
    "WT (rep1)",
    "WT (rep2)"
  ),
  ymax = 50,
  colors = colors)

ggsave(
  glue("{result_folder}upreg_IAPez_upon_Smarcad1KO-H3K9me3-profile.png"),
  width = 10,
  height = 10,
  dpi = 500,
)
ggsave(
  glue("{result_folder}upreg_IAPez_upon_Smarcad1KO-H3K9me3-profile.pdf"),
  width = 10,
  height = 10
)

plot_profile_plot(
  bigwig = dip_bws,
  gr = upreg_iapez_upon_smarcad1ko,
  gr_title = "upregulated IAPEz-int upon Smarcad1 KO",
  title = "hDIP",
  labels = c(
    "KO (rep1)",
    "KO (rep2)",
    "KO, res (rep1)",
    "WT (rep1)",
    "WT (rep2)"
  ),
  ymax = 50,
  colors = colors)

ggsave(
  glue("{result_folder}upreg_IAPez_upon_Smarcad1KO-MeDIP-profile.png"),
  width = 10,
  height = 10,
  dpi = 500,
)
ggsave(
  glue("{result_folder}upreg_IAPez_upon_Smarcad1KO-MeDIP-profile.pdf"),
  width = 10,
  height = 10
)

## mapq bigwigs. highly expressed IAPez (based on TEtranscripts analysis) 
plot_profile_plot(
  bigwig = k9me3_mapq_bws,
  gr = upreg_iapez_upon_smarcad1ko,
  gr_title = "upregulated IAPEz-int upon Smarcad1 KO",
  title = "H3K9me3 (MAPQ > 10)",
  labels = c(
    "KO (rep1)",
    "KO (rep2)",
    "KO, res (rep1)",
    "KO, res (rep2)",
    "WT (rep1)",
    "WT (rep2)"
  ),
  ymax = 50,
  colors = colors)

ggsave(
  glue("{result_folder}upreg_IAPez_upon_Smarcad1KO-mapq_H3K9me3-profile.png"),
  width = 10,
  height = 10,
  dpi = 500,
)
ggsave(
  glue("{result_folder}upreg_IAPez_upon_Smarcad1KO-mapq_H3K9me3-profile.pdf"),
  width = 10,
  height = 10
)

plot_profile_plot(
  bigwig = dip_mapq_bws,
  gr = upreg_iapez_upon_smarcad1ko,
  gr_title = "upregulated IAPEz-int upon Smarcad1 KO",
  title = "hDIP (MAPQ > 10)",
  labels = c(
    "KO (rep1)",
    "KO (rep2)",
    "KO, res (rep1)",
    "KO, res (rep2)",
    "WT (rep1)",
    "WT (rep2)"
  ),
  ymax = 50,
  colors = colors)

ggsave(
  glue("{result_folder}upreg_IAPez_upon_Smarcad1KO-mapq_MeDIP-profile.png"),
  width = 10,
  height = 10,
  dpi = 500,
)
ggsave(
  glue("{result_folder}upreg_IAPez_upon_Smarcad1KO-mapq_MeDIP-profile.pdf"),
  width = 10,
  height = 10
)



