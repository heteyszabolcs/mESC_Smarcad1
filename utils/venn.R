suppressPackageStartupMessages({
  library("tidyverse")
  library("data.table")
  library("Vennerable")
})

# folders
result_folder = "../results/rna_seq_deseq2/"

# Venn diagram function (fc > 2, adj p < 0.05)
create_venn = function(x, 
                       y,
                       filename) {
  x = fread(glue("{result_folder}{x}"))
  y = fread(glue("{result_folder}{y}"))
  x = x %>% filter(padj < 0.05) %>% filter(log2FoldChange > 2) %>% pull(gene_name)
  y = y %>% filter(padj < 0.05) %>% filter(log2FoldChange > 2) %>% pull(gene_name)
  z = list(x, y)
  venn = Venn(z)
  colnames(venn@IndicatorWeight) = c("spike-in", "DESeq2", ".Weight")
  
  png(
    file = glue("{result_folder}{filename}"),
    # The directory you want to save the file in
    width = 8,
    height = 8,
    units = "in",
    res = 500
  )
  plot(
    venn,
    doWeights = FALSE,
    show = list(
      SetLabels = TRUE,
      Faces = FALSE,
      DarkMatter = FALSE
    )
  )
  dev.off() 
}

## overlap between DESeq2 and DESeq2 spike-in gene sets
# create Venns for all fold change table
create_venn(x = "DESeq2_res-kod425_vs_wt_spikein.tsv", 
            y = "DESeq2_res-kod425_vs_wt.tsv",
            filename = "KO-d425_DESEq2-spikein_Venn_padj0.05_fc2.png")
create_venn(x = "DESeq2_res-kok523r_vs_wt_spikein.tsv", 
            y = "DESeq2_res-kok523r_vs_wt.tsv",
            filename = "KO-K523R_DESEq2-spikein_Venn_padj0.05_fc2.png")
create_venn(x = "DESeq2_res-ko_vs_wt_spikein.tsv", 
            y = "DESeq2_res-ko_vs_wt.tsv",
            filename = "KO_DESEq2-spikein_Venn_padj0.05_fc2.png")
create_venn(x = "DESeq2_res-kowt_vs_wt_spikein.tsv", 
            y = "DESeq2_res-kowt_vs_wt.tsv",
            filename = "KOWT_DESEq2-spikein_Venn_padj0.05_fc2.png")


