rm(list = ls())
options(stringsAsFactors = FALSE)

library(biomaRt)

setwd("C:\\Users\\lnrdm\\Desktop\\Differentially Expressed Genes Analysis in Schizophrenia")

################################################
# 0) Define paths
################################################
dirRes       <- "Results/"
if (!dir.exists(dirRes)) dir.create(dirRes)

degFolder  <- "GEO_Schizophrenia"
dirDataset <- file.path(dirRes, degFolder, "/")
if (!dir.exists(dirDataset)) stop("Dataset folder not found: ", dirDataset)

bioMartFolder <- "BioMart"
dirBMart      <- file.path(dirRes, bioMartFolder, "/")
if (!dir.exists(dirBMart)) dir.create(dirBMart)

################################################
# 1) Read DEG gene list
################################################
degFile <- file.path(dirDataset, "DEG.txt") 
if (!file.exists(degFile)) stop("DEG file not found: ", degFile)

deg      <- read.table(
  degFile,
  header    = TRUE,
  sep       = "\t",
  quote     = "",
  stringsAsFactors = FALSE
)
geneList <- deg$genes  # first column

################################################
# 2) Connect to Ensembl BioMart
################################################
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

################################################
# 3) Query HGNC symbol + GO term names
################################################
bmResults <- getBM(
  attributes = c("hgnc_symbol", "name_1006"),
  filters    = "hgnc_symbol",
  values     = geneList,
  mart       = ensembl
)
colnames(bmResults) <- c("GeneSymbol", "GO_term")

################################################
# 4) Find transcription factors
################################################
ind_TF   <- grep(
  "transcription factor activity",
  bmResults$GO_term,
  ignore.case = TRUE
)
TF_genes <- unique(bmResults$GeneSymbol[ind_TF])

################################################
# 5) Save outputs
################################################
annotFile <- file.path(dirBMart, "DEG_with_GO.txt")
write.table(
  bmResults,
  file      = annotFile,
  sep       = "\t",
  row.names = FALSE,
  quote     = FALSE
)

tfFile <- file.path(dirBMart, "TF.txt")
write.table(
  TF_genes,
  file      = tfFile,
  row.names = FALSE,
  col.names = FALSE,
  quote     = FALSE
)

message("Annotated DEG written to: ", annotFile)
message("Transcription factor list written to: ", tfFile)
