suppressPackageStartupMessages({
  library("tidyverse")
  library("data.table")
  library("glue")
  library("dplyr")
  library("circlize")
  library("Hmisc")
})

result_folder = "../results/rna_seq_deseq2/"

# source: https://medium.com/@putri.a.purwono/visualizing-correlation-analysis-results-through-a-chord-diagram-using-the-circlize-package-on-r-6a2d76f65a6d

te_ko = fread(
  "../results/rna_seq_deseq2/TEtranscripts/ERVK_IAPez_extended/TEtranscripts_SMARCAD1_KO_IAPEz_normalized_expr.tsv"
)
te_kowt = fread(
  "../results/rna_seq_deseq2/TEtranscripts/ERVK_IAPez_extended/TEtranscripts_SMARCAD1_KOWT_IAPEz_normalized_expr.tsv"
)
te = te_ko %>% inner_join(., te_kowt, by = "V1") %>%
  filter(str_detect(V1, "IAPEz")) %>%
  select(
    "WT_1" = "WT_rep1Aligned.sortedByCoord.out.bam.C",
    "WT_2" = "WT_rep2Aligned.sortedByCoord.out.bam.C",
    "WT_3" = "WT_rep3Aligned.sortedByCoord.out.bam.C",
    "SMARCAD1_KO_1" =   "SMARCAD1_KO_rep1Aligned.sortedByCoord.out.bam.T.x",
    "SMARCAD1_KO_2" =   "SMARCAD1_KO_rep2Aligned.sortedByCoord.out.bam.T.x",
    "SMARCAD1_KO_3" =   "SMARCAD1_KO_rep3Aligned.sortedByCoord.out.bam.T.x",
    "SMARCAD1_KOWT_1" = "SMARCAD1_KOWT_rep1Aligned.sortedByCoord.out.bam.C",
    "SMARCAD1_KOWT_2" = "SMARCAD1_KOWT_rep2Aligned.sortedByCoord.out.bam.C",
    "SMARCAD1_KOWT_3" = "SMARCAD1_KOWT_rep3Aligned.sortedByCoord.out.bam.C"
  )

result = (rcorr(as.matrix(te), type = c("spearman")))
result = result$r[-c(1:3),-c(4:9)]
df = data.frame(
  from = rep(rownames(result), times = ncol(result)),
  to = rep(colnames(result), each = nrow(result)),
  value = as.vector(result),
  stringsAsFactors = FALSE
)

grid.col = c(rep("#9ecae1", 3), rep("#fc9272", 3), rep("#addd8e", 3))
col_fun = colorRamp2(c(0.50, 0.75, 1),
                     c("#a6bddb", "#fc9272", "#de2d26"),
                     transparency = 0.5)

pdf(
  file = glue(
    "{result_folder}Smarcad1_KO_RNA-Seq-IAPEz_chord_diagram.pdf"
  ),
  width = 7,
  height = 7
)

chordDiagram(
  df,
  col = col_fun,
  grid.col = grid.col,
  annotationTrack = "grid",
  preAllocateTracks = list(track.height = max(strwidth(unlist(
    dimnames(df)
  ))))
)



circos.track(
  track.index = 1,
  panel.fun = function(x, y) {
    circos.text(
      CELL_META$xcenter,
      CELL_META$ylim[1],
      CELL_META$sector.index,
      facing = "clockwise",
      niceFacing = TRUE,
      adj = c(0, 0.5)
    )
  },
  bg.border = NA
)

title(main = "Correlations of normalized IAPEz expressions")
dev.off()
