rm(list = ls())
options(stringsAsFactors = FALSE)

library(R.utils)

setwd("C:\\Users\\lnrdm\\Desktop\\Differentially Expressed Genes Analysis in Schizophrenia")


################################################
# 0) Define paths
################################################
dirRes      <- "Results/"
if (!dir.exists(dirRes)) dir.create(dirRes)

degFolder   <- "GEO_Schizophrenia"
dirDataset  <- file.path(dirRes, degFolder, "/")
if (!dir.exists(dirDataset)) stop("Dataset folder not found: ", dirDataset)

# Input DEG file
degFile     <- file.path(dirDataset, "DEG.txt")
if (!file.exists(degFile)) stop("Cannot find DEG file at: ", degFile)

# Output NCBI info folder
ncbiFolder  <- file.path(dirRes, "NCBI_gene_info/")
if (!dir.exists(ncbiFolder)) dir.create(ncbiFolder)

################################################
# 1) Read DEG.txt and extract geneList
################################################
degFile  <- file.path(dirDataset, "DEG.txt")
deg       <- read.table(
  degFile,
  header       = TRUE,
  sep          = "\t",
  quote        = "",
  stringsAsFactors = FALSE,
  check.names  = FALSE
)
# Confirm tehre really is a column named "genes"
if (!"genes" %in% colnames(deg)) {
  stop("DEG.txt has columns: ", paste(colnames(deg), collapse = ", "),
       "; but expected a 'genes' column.")
}
geneList  <- deg$genes
message("→ Read ", length(geneList), " genes from DEG.txt")

################################################
# 2) Download and unzip NCBI Homo sapiens gene_info
################################################
fileURL   <- "ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz"
destGz    <- file.path(ncbiFolder, "Homo_sapiens_gene_info.gz")
destTxt   <- file.path(ncbiFolder, "Homo_sapiens_gene_info.txt")

download.file(fileURL, destGz, mode = "wb")
gunzip(destGz, destname = destTxt, overwrite = TRUE)

################################################
# 3) Read NCBI gene_info and merge
################################################

# 3a) Grab the header line (it starts with '#') and parse column names
hdr <- readLines(destTxt, n = 1)
# remove the leading '#' and split on tabs
col_names <- strsplit(sub("^#", "", hdr), "\t")[[1]]

# 3b) Read the data skipping the first line, using those column names
NCBI <- read.table(
  destTxt,
  sep           = "\t",
  header        = FALSE,
  quote         = "",
  comment.char  = "",      # already handled header, so don’t skip #
  skip          = 1,       # skip the '#...' line
  fill          = TRUE,    # pad short rows
  stringsAsFactors = FALSE,
  check.names   = FALSE,
  col.names     = col_names
)
message("→ Loaded NCBI with ", nrow(NCBI), " rows and ", ncol(NCBI), " columns")

# 3c) Verify 'Symbol' is now a real column
if (!"Symbol" %in% names(NCBI)) {
  stop("Still no 'Symbol' column after setting col.names; found: ",
       paste(names(NCBI), collapse = ", "))
}

# 3d) Merge with gene list
geneDF <- data.frame(GeneSymbol = geneList, stringsAsFactors = FALSE)
merged <- merge(
  geneDF,
  NCBI,
  by.x    = "GeneSymbol",
  by.y    = "Symbol",
  all.x   = TRUE
)

################################################
# 4) Write output and clean up
################################################
outFile <- file.path(ncbiFolder, paste0("NCBI_gene_info_DEG.txt"))
write.table(merged, outFile,
            sep       = "\t",
            quote     = FALSE,
            row.names = FALSE)

# Remove the unzipped txt
file.remove(destTxt)

message("NCBI gene info merged and saved to: ", outFile)
