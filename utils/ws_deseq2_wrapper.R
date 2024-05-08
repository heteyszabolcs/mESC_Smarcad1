suppressPackageStartupMessages({
  library("tidyverse")
  library("data.table")
  library("glue")
  library("elsasserlib")
  library("apeglm")
  library("wigglescout")
  library("GenomeInfoDb")
  library("plotly")
})

# folders
result_folder = "../results/wigglescout/"
data_folder = "../data/"
bigwig_folder = "../data/MINUTE/bw/"
rnaseq_folder = "../data/RNA-Seq/STAR_output/"

# for subsetting
bed = "../data/bed/RNA-Seq-sign_upreg_ERVK_IAPEz_regions_uponSmarcad1_KO.bed"
bed_wo_chr = "../data/bed/RNA-Seq-sign_upreg_ERVK_IAPEz_regions_uponSmarcad1_KO_wo_chr.bed"

# bws
k9me3_wt = list.files(bigwig_folder, pattern = "H3K9me3_WT", full.names = TRUE)
k9me3_wt = k9me3_wt[grep(k9me3_wt, pattern = ".mm10.scaled.mapq.bw")]
k9me3_smarcad1ko = list.files(bigwig_folder, pattern = "DIPh_KO_P", full.names = TRUE)
k9me3_smarcad1ko = k9me3_smarcad1ko[grep(k9me3_smarcad1ko, pattern = ".mm10.scaled.mapq.bw")]
diph_wt = list.files(bigwig_folder, pattern = "DIPh_WT", full.names = TRUE)
diph_wt = diph_wt[grep(diph_wt, pattern = ".mm10.scaled.mapq.bw")]
diph_smarcad1ko = list.files(bigwig_folder, pattern = "DIPh_KO_P", full.names = TRUE)
diph_smarcad1ko = diph_smarcad1ko[grep(diph_smarcad1ko, pattern = ".mm10.scaled.mapq.bw")]

# RNA-Seq bigwigs
rna_wt = list.files(rnaseq_folder, pattern = "^WT_rep", full.names = TRUE)
rna_smarcad1ko = list.files(rnaseq_folder, pattern = "^SMARCAD1_KO_rep", full.names = TRUE)

# DESeq2 wrapper
k9me3_fc = bw_bed_diff_analysis(
   bwfiles_c1 = k9me3_wt,
   bwfiles_c2 = k9me3_smarcad1ko,
   bed = bed,
   label_c1 = "KO",
   label_c2 = "WT", p_cutoff = 0.1, shrink = TRUE
)

# export
k9me3_fc = as.data.frame(k9me3_fc)
k9me3_fc = k9me3_fc %>% mutate(id = paste(as.character(seqnames), as.character(start), as.character(end), sep = "_"))

diph_fc = bw_bed_diff_analysis(
  bwfiles_c1 = diph_wt,
  bwfiles_c2 = diph_smarcad1ko,
  bed = bed,
  label_c1 = "KO",
  label_c2 = "WT", p_cutoff = 0.1, shrink = TRUE
)

# export
diph_fc = as.data.frame(diph_fc)
diph_fc = diph_fc %>% mutate(id = paste(as.character(seqnames), as.character(start), as.character(end), sep = "_"))

rna_fc = bw_bed_diff_analysis(
  bwfiles_c1 = rna_wt,
  bwfiles_c2 = rna_smarcad1ko,
  bed = bed_wo_chr,
  label_c1 = "KO",
  label_c2 = "WT", p_cutoff = 0.1, shrink = TRUE
)
rna_fc = as.data.frame(rna_fc)
rna_fc = rna_fc %>% mutate(seqnames = paste0("chr", as.character(seqnames)))
rna_fc = rna_fc %>% mutate(id = paste(as.character(seqnames), as.character(start), as.character(end), sep = "_"))

fc = rna_fc %>% inner_join(., k9me3_fc, by = "id") %>% dplyr::select(RNASeq_log2FC = log2FoldChange.x, K9me3_log2FC = log2FoldChange.y, id) %>% 
  inner_join(., diph_fc, by = "id") %>% dplyr::select(DIPh_log2FC = log2FoldChange, RNASeq_log2FC, K9me3_log2FC) %>% drop_na()

fc_long = fc %>% pivot_longer("RNASeq_log2FC":"K9me3_log2FC", names_to = "comparison", values_to = "log2FC")
write_csv(fc, glue("{result_folder}wigglescout_DESeq2_output.csv"))

p = plot_ly(data = fc, 
  x = ~RNASeq_log2FC,
  y = ~K9me3_log2FC,
  z = ~DIPh_log2FC,
  type = "scatter3d",
  mode = "markers",
  color = fc$RNASeq_log2FC
)  %>%
  layout(title = "Smarcad1 KO vs. WT", plot_bgcolor = "white", xaxis = list(title = "RNA-Seq log2FC"), 
  yaxis = list(title = "H3K9me3 log2FC"), xaxis = list(title = "hDIP log2FC"))
p





