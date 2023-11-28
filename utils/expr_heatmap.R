suppressPackageStartupMessages({
  library("tidyverse")
  library("data.table")
  library("ggplot2")
  library("glue")
  library("GenomicRanges")
  library("wigglescout")
  library("ComplexHeatmap")
  library("circlize")
  library("cowplot")
})

# folders
deseq2_folder = "../results/rna_seq_deseq2/"

# RNA-Seq data
rlog_counts = fread(glue("{deseq2_folder}reg_log_norm_counts_spikein.tsv"))
ko_fcs = fread(glue("{deseq2_folder}DESeq2_res-kowt_foldchanges.tsv"))

# fold change table
fc = fread(glue("{deseq2_folder}DESeq2_res_ko-vs-kowt_padj0.05_fc1.tsv"))

# heatmap of upregulated genes 
sign_up_genes = fc %>% filter(padj < 0.01) %>% 
  filter(log2FoldChange > 2) %>%
  arrange(desc(log2FoldChange))

mat = rlog_counts %>% filter(gene_symbol %in% sign_up_genes$gene_name) 
genes = mat$gene_symbol
mat = mat %>% select("Ph001":"Ph015") %>% as.matrix
rownames(mat) = genes
mat = mat[sign_up_genes$gene_name,]

col_fun = colorRamp2(c(2, 6, 10), c("#9ecae1", "white", "#fc9272"))
ha = HeatmapAnnotation(
  condition = c(colnames(mat)),
  col = list(condition = c("Ph001" = "#a6bddb", "Ph006" = "#a6bddb", "Ph011" = "#a6bddb",
                           "Ph002" = "#e6550d", "Ph007" = "#e6550d", "Ph012" = "#e6550d",
                           "Ph003" = "#fdae6b", "Ph008" = "#fdae6b", "Ph013" = "#fdae6b",
                           "Ph004" = "#31a354", "Ph009" = "#31a354", "Ph014" = "#31a354",
                           "Ph005" = "#a1d99b", "Ph010" = "#a1d99b", "Ph015" = "#a1d99b")
  ),
  gp = gpar(col = "black", fontsize = 4),
  show_legend = FALSE
)

png(
  file = glue("{deseq2_folder}rlog_heatmap-sign_up_genes_KO-vs-KOWT.png"),
  width = 10,
  height = 10,
  units = 'in',
  res = 500
)
hm = Heatmap(
  mat,
  name = "exp",
  clustering_distance_rows = "pearson",
  col = col_fun,
  rect_gp = gpar(col = "black", lwd = 1),
  top_annotation = ha,
  cluster_columns = TRUE,
  row_title = "KO vs. KO-WT DEG - adj.p < 0.01, log2FC > 2",
  column_title = "rlog norm. expression",
  cluster_rows = FALSE,
  width = unit(14, "cm"),
  height = unit(14, "cm"),
  row_names_gp = gpar(fontsize = 8)
) 

lgd = Legend(at = c("WT", "KO", "KO-WT", "KO-d425", "KO-K523R"), 
             title = "SMARCAD1", legend_gp = gpar(fill = c("#a6bddb", "#e6550d", "#fdae6b", "#31a354", "#a1d99b")))
hm = draw(hm, heatmap_legend_list = lgd)
hm

dev.off()

# heatmap of downregulated genes 
sign_down_genes = fc %>% filter(padj < 0.01) %>% 
  filter(log2FoldChange < -2) %>%
  arrange(desc(log2FoldChange))

mat = rlog_counts %>% filter(gene_symbol %in% sign_down_genes$gene_name) 
genes = mat$gene_symbol
mat = mat %>% select("Ph001":"Ph015") %>% as.matrix
rownames(mat) = genes
mat = mat[sign_down_genes$gene_name,]

col_fun = colorRamp2(c(2, 6, 10), c("#9ecae1", "white", "#fc9272"))
ha = HeatmapAnnotation(
  condition = c(colnames(mat)),
  col = list(condition = c("Ph001" = "#a6bddb", "Ph006" = "#a6bddb", "Ph011" = "#a6bddb",
                           "Ph002" = "#e6550d", "Ph007" = "#e6550d", "Ph012" = "#e6550d",
                           "Ph003" = "#fdae6b", "Ph008" = "#fdae6b", "Ph013" = "#fdae6b",
                           "Ph004" = "#31a354", "Ph009" = "#31a354", "Ph014" = "#31a354",
                           "Ph005" = "#a1d99b", "Ph010" = "#a1d99b", "Ph015" = "#a1d99b")
  ),
  gp = gpar(col = "black", fontsize = 4),
  show_legend = FALSE
)

png(
  file = glue("{deseq2_folder}rlog_heatmap-sign_down_genes_KO-vs-KOWT.png"),
  width = 10,
  height = 10,
  units = 'in',
  res = 500
)
hm = Heatmap(
  mat,
  name = "exp",
  clustering_distance_rows = "pearson",
  col = col_fun,
  rect_gp = gpar(col = "black", lwd = 1),
  top_annotation = ha,
  cluster_columns = TRUE,
  row_title = "KO vs. KO-WT DEG - adj.p < 0.01, log2FC < -2",
  column_title = "rlog norm. expression",
  cluster_rows = FALSE,
  width = unit(14, "cm"),
  height = unit(14, "cm"),
  row_names_gp = gpar(fontsize = 8)
) 

lgd = Legend(at = c("WT", "KO", "KO-WT", "KO-d425", "KO-K523R"), 
             title = "SMARCAD1", legend_gp = gpar(fill = c("#a6bddb", "#e6550d", "#fdae6b", "#31a354", "#a1d99b")))
hm = draw(hm, heatmap_legend_list = lgd)
hm

dev.off()

# fold change heatmap of upregulated genes 
ko_up_fcs = ko_fcs %>% filter(gene_name %in% sign_up_genes$gene_name)
genes = ko_up_fcs$gene_name
vars = ko_up_fcs %>% select("KO_vs_KOWT":"WT_vs_KOWT") %>% as.matrix %>% rowVars
ko_up_fcs = ko_up_fcs %>% mutate(variance = vars) %>% arrange(desc(variance)) %>% select(-variance) 
rows = ko_up_fcs$gene_name
ko_up_fcs = ko_up_fcs %>% select(-gene_name)
ko_up_fcs = as.matrix(ko_up_fcs)
rownames(ko_up_fcs) = rows

png(
  file = glue("{deseq2_folder}upregulated_fc_heatmap-KOWT_contrasts.png"),
  width = 7,
  height = 10,
  units = 'in',
  res = 500
)

col_fun = colorRamp2(c(2, 5, 8), c("#9ecae1", "white", "#fc9272"))
hm = Heatmap(
  ko_up_fcs,
  name = "log2 fc",
  clustering_distance_rows = "pearson",
  col = col_fun,
  rect_gp = gpar(col = "black", lwd = 1),
  cluster_columns = TRUE,
  row_title = "KO vs. KO-WT DEG - adj.p < 0.01, log2FC > 2",
  column_title = "upregulated genes of KO vs. WT",
  cluster_rows = FALSE,
  width = unit(1, "cm"),
  height = unit(10, "cm"),
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 6)
) 
hm

dev.off()

# fold change heatmap of downregulated genes 
ko_down_fcs = ko_fcs %>% filter(gene_name %in% sign_down_genes$gene_name)
genes = ko_down_fcs$gene_name
vars = ko_down_fcs %>% select("KO_vs_KOWT":"WT_vs_KOWT") %>% as.matrix %>% rowVars
ko_down_fcs = ko_down_fcs %>% mutate(variance = vars) %>% arrange(desc(variance)) %>% select(-variance) 
rows = ko_down_fcs$gene_name
ko_down_fcs = ko_down_fcs %>% select(-gene_name)
ko_down_fcs = as.matrix(ko_down_fcs)
rownames(ko_down_fcs) = rows

png(
  file = glue("{deseq2_folder}downregulated_fc_heatmap-KOWT_contrasts.png"),
  width = 7,
  height = 10,
  units = 'in',
  res = 500
)

col_fun = colorRamp2(c(-8, -5, -2), c("#fc9272", "white", "#9ecae1"))
hm = Heatmap(
  ko_down_fcs,
  name = "log2 fc",
  clustering_distance_rows = "pearson",
  col = col_fun,
  rect_gp = gpar(col = "black", lwd = 1),
  cluster_columns = TRUE,
  row_title = "KO vs. KO-WT DEG - adj.p < 0.01, log2FC < -2",
  column_title = "downregulated genes of KO vs. WT",
  cluster_rows = FALSE,
  width = unit(1, "cm"),
  height = unit(10, "cm"),
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 6)
) 
hm

dev.off()

