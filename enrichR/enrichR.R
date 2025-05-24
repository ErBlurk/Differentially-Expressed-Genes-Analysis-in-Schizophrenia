library(enrichR)
library(ggplot2)
library(forcats)
library(stringr)

setwd("C:\\Users\\lnrdm\\Desktop\\Differentially Expressed Genes Analysis in Schizophrenia")

# load helper functions
source("enrichR/getEnrichment.R")
source("enrichR/getEnrichmentPlot.R")

################################################
# 0) set up paths & parameters
################################################
dataset   <- "GEO_Schizophrenia"
dirRes    <- "Results/"
dirDataset<- file.path(dirRes, dataset, "/")
dirEnrich <- file.path(dirRes, "Functional_Enrichment", "/")
for (d in c(dirRes, dirDataset, dirEnrich)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

top_term <- 10
thr_pval <- 0.15

################################################
# 1) load DEG list
################################################
file_input_list <- file.path(dirDataset, "DEG.txt")
deg        <- read.table(file_input_list, sep = "\t", header = TRUE, check.names = FALSE, quote = "")
up_genes   <- deg$genes[deg$direction == "UP"]
down_genes <- deg$genes[deg$direction == "DOWN"]

################################################
# 2) choose databases & run enrichr
################################################
dbs <- c("DisGeNET", "GO_Molecular_Function_2021", "GO_Biological_Process_2021", "KEGG_2021_Human")

enrich_up   <- enrichr(up_genes,   dbs)
enrich_down <- enrichr(down_genes, dbs)

################################################
# 3) plot results
################################################
if (length(enrich_up) > 0) {
  getEnrichment(enrich_up, "UP")
} else {
  message("No enrichment results for UP genes.")
}

if (length(enrich_down) > 0) {
  getEnrichment(enrich_down, "DOWN")
} else {
  message("No enrichment results for DOWN genes.")
}