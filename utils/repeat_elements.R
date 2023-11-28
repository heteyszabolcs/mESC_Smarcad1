repeats = fread("../data/bed/UCSC_RepeatMasker_mm9.tsv")
families = unique(repeats$repFamily)
families = families[grep("\\?", families, invert = TRUE)]
families = families[grep("DNA", families, invert = TRUE)]
families = families[grep("RNA", families, invert = TRUE)]
families = families[grep("centr", families, invert = TRUE)]
families = families[grep("Other", families, invert = TRUE)]
families = families[grep("Unknown", families, invert = TRUE)]
