#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("glue"))

# create parser object
parser <- ArgumentParser()

# add arguments and specify our desired options 
parser$add_argument("-bw", "--bigwig", type = "character",
                    help = "bigwig to be visualized")
parser$add_argument("-b", "--bed", type = "character",
                    help = "bed to subset")
parser$add_argument("-l", "--label", type = "character",
                    help = "y axis label")
parser$add_argument("-o", "--output", type = "character",
                    help = "output filename")
parser$add_argument("-ymax", "--ymax", type = "integer",
                    help = "ymax of profile plot", default = 150)
parser$add_argument("-yzmax", "--yzmax", type = "integer",
                    help = "z max of heatmap", default = 150)

args <- parser$parse_args()

bigwig = args$bigwig
bed = args$bed
label = args$label
output = args$output
ymax = args$ymax
zmax = args$yzmax

# deeptools
system(glue("computeMatrix scale-regions -S {bigwig} \\
-R {bed} \\
--beforeRegionStartLength 3000 \\
--regionBodyLength 5000 \\
--afterRegionStartLength 3000 \\
--skipZeros -o ../results/deeptools/matrix.mat.gz"))

system(glue("plotProfile -m ../results/deeptools/matrix.mat.gz \\
-out {output}_profile \\
--numPlotsPerRow 2 \\
--plotTitle \"UCSC repmasker mm10\""))

# system(glue("computeMatrix reference-point -S {bigwig} \\
# -R  {bed} \\
# -b 3000 \\
# -a 3000 \\
# --samplesLabel {label} \\
# --skipZeros -o ../results/deeptools/matrix.mat.gz"))
# 
# system(glue("plotHeatmap -m ../results/deeptools/matrix.mat.gz \\
# -out {output} \\
# --refPointLabel \"repeat\" \\
# --heatmapHeight 14 \\
# --colorMap \"Greens\" \\
# --yMin 0 \\
# --yMax {ymax} \\
# --zMax {zmax} \\
# --yAxisLabel \"\" \\
# --xAxisLabel \"\" \\
# --regionsLabel \"UCSC repmasker mm10\" \\
# --legendLocation \"none\""))