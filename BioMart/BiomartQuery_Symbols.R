# source("https://bioconductor.org/biocLite.R")
# biocLite("biomaRt")
# https://bioconductor.org/packages/release/bioc/html/biomaRt.html

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("biomaRt")

library(biomaRt)

setwd("C:\\Users\\lnrdm\\Desktop\\Differentially Expressed Genes Analysis in Schizophrenia")

# Results folder
dirBMart    <- "BioMart/"
if (!dir.exists(dirBMart)) {
  dir.create(dirBMart)
} else {
  message("the directory ", dirBMart, " already exists")
}

options(stringsAsFactors=FALSE)
###########################
# input
# use DEG.txt from the schizophrenia case study
###########################
input_list <- read.table(
  file.path("Results", "GEO_Schizophrenia", "DEG.txt"),
  header = TRUE, sep = "\t", stringsAsFactors = FALSE
)
colnames(input_list) <- "GeneSymbol"

# download DB
ensembl <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
attributes <- listAttributes(ensembl)  # list of possibile output elements
filters    <- listFilters(ensembl)     # list of possibile input elements

# download attributes
df_entrez <- getBM(
  attributes = c("hgnc_symbol", "entrezgene_id"),
  filters    = "hgnc_symbol",
  values     = input_list$GeneSymbol,
  mart       = ensembl
)

df_Ensembl <- getBM(
  attributes = c("hgnc_symbol", "ensembl_gene_id", "ensembl_transcript_id"),
  filters    = "hgnc_symbol",
  values     = input_list$GeneSymbol,
  mart       = ensembl
)

df_refseq <- getBM(
  attributes = c("hgnc_symbol", "refseq_mrna"),
  filters    = "hgnc_symbol",
  values     = input_list$GeneSymbol,
  mart       = ensembl
)
ind <- which(df_refseq$refseq_mrna == "")
df_refseq <- df_refseq[-ind, ]

df_gene_type <- getBM(
  attributes = c("hgnc_symbol", "gene_biotype"),
  filters    = "hgnc_symbol",
  values     = input_list$GeneSymbol,
  mart       = ensembl
)

rm(attributes, filters, ensembl)

# file output (to Results/GEO_Schizophrenia folder)
write.table(
  df_entrez,
  file      = file.path("Results", "BioMart", "biomart_entrez.txt"),
  row.names = FALSE,
  col.names = c("GeneSymbol", "EntrezGene_ID"),
  sep       = "\t",
  quote     = FALSE
)

write.table(
  df_gene_type,
  file      = file.path("Results", "BioMart", "biomart_geneBiotype.txt"),
  row.names = FALSE,
  col.names = c("GeneSymbol", "GeneBiotype"),
  sep       = "\t",
  quote     = FALSE
)

write.table(
  df_Ensembl,
  file      = file.path("Results", "BioMart", "biomart_Ensembl.txt"),
  row.names = FALSE,
  col.names = c("GeneSymbol", "Ensembl_gene_ID", "Ensembl_transcript_ID"),
  sep       = "\t",
  quote     = FALSE
)

write.table(
  df_refseq,
  file      = file.path("Results", "BioMart", "biomart_refseq.txt"),
  row.names = FALSE,
  col.names = c("GeneSymbol", "RefSeq"),
  sep       = "\t",
  quote     = FALSE
)
