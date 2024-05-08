library("data.table")
library("tidyverse")

gtf = fread(
  "C:/Szabolcs/Karolinska/Data/reference_data/elsasser_repeat_atlas_v1.0.mm39.gtf"
)
mod = gtf %>% separate(
  V9,
  sep = ";",
  into = c(
    "gene_id",
    "transcript_id",
    "family_id",
    "class_id",
    "gene_name"
  )
)
mod = mod %>% mutate(gene_id = transcript_id)
mod = mod %>% mutate(gene_id = str_replace(gene_id, " transcript_id", "gene_id"))
mod = mod %>% mutate(V9 = paste(gene_id, transcript_id, family_id, class_id, gene_name, sep = ";")) %>%
  dplyr::select(-c(gene_id, transcript_id, family_id, class_id, gene_name))


write.table(
  mod,
  "C:/Szabolcs/Karolinska/Data/reference_data/elsasser_repeat_atlas_v1.0.mm39_mod.gtf",
  col.names = FALSE,
  quote = FALSE,
  row.names = FALSE,
  sep = "\t"
)
