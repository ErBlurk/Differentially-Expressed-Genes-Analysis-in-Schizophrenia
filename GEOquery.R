# to install the package
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("GEOquery")
# a

readr::local_edition(1)
library("GEOquery")

##### Set Results folder##### 
dirRes <- "Results/"
if (!dir.exists(dirRes)) {
  dir.create(dirRes)
} else {
  message("The directory ", dirRes, " already exists")
}

# Project folder
dataset <- "GEO_Schizophrenia"
dirDataset <- file.path(dirRes, dataset)
if (!dir.exists(dirDataset)) {
  dir.create(dirDataset)
} else {
  message("The directory ", dirDataset, " already exists")
}

###########################################
# STEP 1: Downloading data
###########################################
# schizophrenia whole blood expression dataset - already log scaled and quantile normalized
series <- "GSE38484"
tmp <- getGEO(GEO = series)

# If there are multiple platforms, pick the first one:
set <- tmp[[1]]

pData <- phenoData(set)       # information about the samples
metadata <- pData@data
aData <- assayData(set)
exprs_mat <- data.frame(aData$exprs)
# exprs_mat <- normalizeBetweenArrays(exprs_mat, method="quantile")             # data already normalized
rm(pData, aData, tmp)

###########################################
# STEP 2: Preparing data
###########################################
annotation <- fData(set)
geneSymbol <- annotation$Symbol

# reorder rows so they match annotation IDs
keep <- match(annotation$ID, rownames(exprs_mat))
exprs_mat <- exprs_mat[keep, , drop = FALSE]
geneSymbol <- geneSymbol[keep]

# collapse probes to genes by mean
exprs_mat <- aggregate(exprs_mat, by = list(geneSymbol), FUN = mean)

# remove rows without a gene symbol
empty_genes <- which(exprs_mat$Group.1 == "")
exprs_mat <- exprs_mat[-empty_genes, ]

# set gene symbols as row names
rownames(exprs_mat) <- exprs_mat$Group.1
exprs_mat <- exprs_mat[, -1]
rm(set, annotation, empty_genes, geneSymbol, keep)

############################################
# STEP 3: Extracting case/control samples
###########################################
# keep only the relevant metadata columns
metadata <- metadata[, c("geo_accession", "status:ch1", "age:ch1", "gender:ch1")]

# check the actual levels
print(unique(metadata$`status:ch1`))
# [1] "schizophrenia (SCZ)" "CONTROL"

# split into case vs control by the exact 'status:ch1' values
sampleList <- split(metadata, metadata$`status:ch1`)
control    <- sampleList$`CONTROL`
case       <- sampleList$`schizophrenia (SCZ)`

# extract the GEO IDs
control_ids <- control$geo_accession
case_ids    <- case$geo_accession

# subset the expression matrix
data_all     <- exprs_mat[, c(control_ids, case_ids), drop = FALSE]
data_control <- exprs_mat[, control_ids,    drop = FALSE]
data_case    <- exprs_mat[, case_ids,       drop = FALSE]

# reorder metadata to match the columns of data_all
metadata <- metadata[match(colnames(data_all), metadata$geo_accession), ]
rownames(metadata) <- metadata$geo_accession

# now write out
write.table(control_ids, file = file.path(dirDataset, "control.txt"),
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(case_ids,    file = file.path(dirDataset, "case.txt"),
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

rm(sampleList)

############################################
# STEP 4: Export data
############################################
# now data_all is ready for differential expression analysis
write.table(
  data_all,
  file      = file.path(dirDataset, "matrix.txt"),
  sep       = "\t",
  col.names = NA,
  row.names = TRUE,
  quote     = FALSE
)
write.table(
  control_ids,
  file      = file.path(dirDataset, "control.txt"),
  sep       = "\t",
  col.names = FALSE,
  row.names = FALSE,
  quote     = FALSE
)
write.table(
  case_ids,
  file      = file.path(dirDataset, "case.txt"),
  sep       = "\t",
  col.names = FALSE,
  row.names = FALSE,
  quote     = FALSE
)
write.table(
  metadata,
  file      = file.path(dirDataset, "metadata.txt"),
  sep       = "\t",
  col.names = TRUE,
  row.names = FALSE,
  quote     = FALSE
)

# Remove downloaded GEO data (temp files)
unlink(series, recursive = TRUE)

