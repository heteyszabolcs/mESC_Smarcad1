if (!require("pacman")) install.packages("pacman")
pacman::p_load("tidyverse",
               "data.table",
               "wigglescout",
               "ggpubr",
               "glue",
               "RColorBrewer",
               "GenomicRanges",
               "magick",
               "ggrastr") 

# output folder
result_folder = "../results/wigglescout/"

# Minute-ChIP scaled bigwigs
minutes = list.files("../data/MINUTE/bw/",
                   pattern = "*bw",
                   full.names = TRUE)

# regions
# Sachs et al peaks
smarcad1_peaks = fread("../data/bed/Sachs_et_al_WT_SMARCAD1_R1_peaks.narrowPeak")
smarcad1_peaks$V11 = "Sachs et al, Smarcad1"
smarcad1_peaks = GRanges(
  seqnames = smarcad1_peaks$V1,
  ranges = IRanges(
    start = smarcad1_peaks$V2,
    end = smarcad1_peaks$V3,
    names = smarcad1_peaks$V11
  )
)

iapez = fread("../data/bed/Navarro_2020_IAPEz_consensus.bed")
iapez$V7 = "IAPEz"
iapez = GRanges(
  seqnames = iapez$V1,
  ranges = IRanges(
    start = iapez$V2,
    end = iapez$V3,
    names = iapez$V7
  )
)



create_scatter = function(bigwig1 = minutes[1], bigwig2 = minutes[13]) {
  name1 = strsplit(strsplit(bigwig1, "../data/MINUTE/bw/")[[1]][2],
                   ".mm9.scaled.bw")[[1]][1]
  name2 = strsplit(strsplit(bigwig2, "../data/MINUTE/bw/")[[1]][2],
                   ".mm9.scaled.bw")[[1]][1]
  
  done = list.files(glue("{result_folder}"), full.names = FALSE)
  if (glue("{name1}-{name2}_sc.pdf") %in% done) {
    return(print(glue("{name1}-{name2} is done!")))
  } else if (glue("{name2}-{name1}_sc.pdf") %in% done) {
    return(print(glue("{name2}-{name1} is done!")))
  } else {
    print(glue("Working on: {name1}-{name2}"))
  }
  
  p = plot_bw_loci_scatter(
    bigwig1,
    bigwig2,
    loci = iapez,
    remove_top = 0.0001,
    verbose = FALSE
  ) + 
    xlim(0, 80) +
    ylim(0, 80) +
    theme(
      text = element_text(size = 16),
      legend.text = element_text(size = 14),
      plot.title = element_text(size = 7, face = "bold"),
      axis.text.x = element_text(size = 12, color = "black"),
      axis.text.y = element_text(size = 12, color = "black"),
      axis.title.x = element_text(size = 14, color = "black"),
      axis.title.y = element_text(size = 14, color = "black")
    ) +
    guides(color = "none", fill = "none") +
    geom_abline(slope = 1,
                intercept = 0,
                color = "#de2d26",
                linetype = "dashed") +
    labs(title = glue("{name1} - {name2}"),
         x = name1,
         y = name2) +
    stat_cor(method = "pearson", label.x = 45, label.y = 75, size = 2)
  # rasterize!
  p = rasterize(p, layers='Point', dpi = 150)
  # save
  ggsave(
    glue("{result_folder}{name1}-{name2}_sc.pdf"),
    plot = p,
    width = 3,
    height = 3
  )
  return(p)
}

lapply(minutes, function(x) {
    lapply(minutes, function(y) {
      create_scatter(x, y)
    })
  })
