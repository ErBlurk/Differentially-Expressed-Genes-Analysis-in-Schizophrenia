rm(list = ls())
options(stringsAsFactors = FALSE)

library(stringr)
library(pheatmap)    # heatmap

# Results folder
dirRes    <- "Results/"
if (!dir.exists(dirRes)) {
  dir.create(dirRes)
} else {
  message("the directory ", dirRes, " already exists")
}

# Load downloaded data (from GEOquery.R script)
dataset    <- "GEO_Schizophrenia"
dirDataset <- file.path(dirRes, dataset, "/")
if (!dir.exists(dirDataset)) {
  dir.create(dirDataset)
} else {
  message("the directory ", dirDataset, " already exists")
}

######################
# DEG parameters
prc_IQR  <- 0.1
thr_fc   <- 1.22
thr_pval <- 0.05
paired   <- FALSE

######################
# Output files after DEG analysis
filename_DEG         <- file.path(dirDataset, "DEG.txt")  
filename_matrix_DEG  <- file.path(dirDataset, "matrix_DEG.txt")
filename_heatmap     <- file.path(dirDataset, "heatmap.png")

#################
# Step 1: Data loading
#################

# Import data matrix
data        <- read.table(
  file.path(dirDataset, "matrix.txt"),
  header    = TRUE,
  sep       = "\t",
  quote     = "",
  check.names = FALSE,
  row.names = 1
)

# Import control cases
control_ids <- read.table(
  file.path(dirDataset, "control.txt"),
  header    = FALSE,
  sep       = "\t",
  quote     = "",
  stringsAsFactors = FALSE
)[,1]

# Import case cases
case_ids <- read.table(
  file.path(dirDataset, "case.txt"),
  header    = FALSE,
  sep       = "\t",
  quote     = "",
  stringsAsFactors = FALSE
)[,1]

# Import the extracted metadata
metadata <- read.table(
  file.path(dirDataset, "metadata.txt"),
  header    = TRUE,
  sep       = "\t",
  quote     = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Store the data in a useful format
dataN    <- data[, control_ids]
dataC    <- data[, case_ids]
data     <- data[, c(control_ids, case_ids)]

#################
# step 2: ANALYSIS
#################

gene_counts <- list()
gene_counts$Init <- nrow(data)

# Preprocessing 1.1: overall mean - filtering out 0 mean values
genes <- rownames(data)
overall_mean <- apply(data, 1, mean) # we delete those with overall 0
ind <- which(overall_mean == 0)
if (length(ind) > 0) {
  dataN <- dataN[-ind, , drop = FALSE]
  dataC <- dataC[-ind, , drop = FALSE]
  data  <- data[-ind, , drop = FALSE]
  genes <- genes[-ind]
}
rm(ind)
gene_counts$ZeroMean <- nrow(data)



# Preprocessing 1.2: log transformation, sipped due to already transformed data
# step 2.1 : logarithmic transformation
# dataN <- log2(dataN + 1)
# dataC <- log2(dataC + 1)
# data  <- log2(data  + 1)



# Preprocessing 1.3: IQR value for each gene and deletion of out of IQR range genes
# step 2.2 : pre-processing
variation <- apply(data, 1, IQR) # IQR helps measure variation around the median
# IQR filtering
thr_prc <- quantile(variation, prc_IQR)
ind <- which(variation <= thr_prc)
dataN <- dataN[-ind, , drop = FALSE]
dataC <- dataC[-ind, , drop = FALSE]
data  <- data[-ind, , drop = FALSE]
genes <- genes[-ind]
rm(ind)
# IQR distribution plot
hist(
  variation,
  main   = "IQR Frequency distribution",
  breaks = 100,
  xlab   = "IQR value",
  ylab   = "Frequency",
  col    = "seagreen"
)
abline(v = thr_prc, lty = 2, lwd = 4, col = "darkgoldenrod1")
gene_counts$IQR <- nrow(data)
# The result histogram is skewed towards the origin, this should indicate small IQR
# and small variation in genes expression -> successfull filtering



# Filtering 2.0: setup
# data is already log scaled, log transform also the threeshold fold change
thr_logFC <- log2(thr_fc)



# Filtering 2.4: volcano plot over whole data
# Save originals for whole data volcano plot
data_all   <- data
dataN_all  <- dataN
dataC_all  <- dataC
# Compute logFC & p-values on all genes
N_all      <- ncol(dataN_all)
M_all      <- ncol(dataC_all)
logFC_all  <- rowMeans(dataC_all) - rowMeans(dataN_all)
pval_all   <- apply(data_all, 1, function(x)
  t.test(x[1:N_all], x[(N_all+1):(N_all+M_all)], paired = paired)$p.value
)
pval_adj_all <- p.adjust(pval_all, method = "fdr")
logP_all     <- -log10(pval_adj_all)
# Color mapping
colori_all <- ifelse(
  logFC_all >  thr_logFC & pval_adj_all < thr_pval, "red",
  ifelse(logFC_all < -thr_logFC & pval_adj_all < thr_pval, "yellow", "black")
)
# Full volcano (all points)
plot(
  logFC_all, logP_all,
  col   = colori_all,
  pch   = 20,
  main  = "Volcano Plot",
  xlab  = "log2 Fold Change",
  ylab  = "-log10 P-value",
  xlim  = c(min(logFC_all) - 0.2, max(logFC_all) + 0.2),
  ylim  = c(0,              max(logP_all) + 2)
)
abline(h = -log10(thr_pval),      col = "blue4", lty = 2)
abline(v = c(-thr_logFC, thr_logFC), col = "blue4", lty = 2)
# In black, points to be filtered out



# Filtering 2.1: compute logFC on log data + data filtering
# step 2.3 : Log FC
logFC <- rowMeans(dataC) - rowMeans(dataN)   # remember data is already log scaled
print(log2(thr_fc))
hist(
  logFC,
  main   = "FC (logarithmic) frequency distribution",
  breaks = 100,
  xlab   = "logFC",
  ylab   = "frequency",
  col    = "orange"
)
abline(v = c(-log2(thr_fc), log2(thr_fc)), lty = 2, lwd = 4, col = "red")
# gene probes up/down regulation indicator - well distributed
# far threshold lines -> more stringent condition for filtering and analysis



# Filtering 2.2: p-value computation
N <- ncol(dataN)
M <- ncol(dataC)
pval <- apply(data, 1, function(x) {
  t.test(x[1:N], x[(N+1):(M+N)], paired = paired)$p.value
})
# Welch's unpaired test - can't do a student test due to unpaired cases-control
# Applies a two-sample t-test to each gene to compare control (dataN) vs. case (dataC) groups



# Filtering 2.3: p-value adjustment
# adjustment of p-values
pval_adj <- p.adjust(pval, method = "fdr")
# Controls the false discovery rate using the Benjamini-Hochberg method



# Filtering 2.5: Remove genes with low fold change (|log2FC| < threshold)
# Keeps only genes with strong up/down regulation
# Filter data
ind <- which(abs(logFC) < log2(thr_fc))
if (length(ind) > 0) {
  dataN    <- dataN[-ind, , drop = FALSE]
  dataC    <- dataC[-ind, , drop = FALSE]
  data     <- data[-ind, , drop = FALSE]
  genes    <- genes[-ind]
  logFC    <- logFC[-ind]
  pval     <- pval[-ind]
  pval_adj <- pval_adj[-ind]
}
rm(ind)
gene_counts$FC <- nrow(data)



# Filtering 2.6: Remove genes with non-significant adjusted p-values (FDR > threshold)
# Keeps only genes passing the statistical significance threshold
# p-value filtering
ind <- which(pval_adj > thr_pval)
if (length(ind) > 0) {
  dataN    <- dataN[-ind, , drop = FALSE]
  dataC    <- dataC[-ind, , drop = FALSE]
  data     <- data[-ind, , drop = FALSE]
  genes    <- genes[-ind]
  logFC    <- logFC[-ind]
  pval     <- pval[-ind]
  pval_adj <- pval_adj[-ind]
}
rm(ind)
gene_counts$Pval <- nrow(data)



# Plot the number of genes at each step
steps <- names(gene_counts)
counts <- unlist(gene_counts)

# Plot the changes in the number of genes
bar_heights <- barplot(
  counts, 
  names.arg = steps, 
  col = "skyblue", 
  main = "Changes in Gene Counts After Filtering", 
  xlab = "Filtering Step", 
  ylab = "Number of Genes", 
  las = 1  # Rotate x-axis labels for better readability
)
# Add the gene counts inside the bars
text(
  x = bar_heights,
  y = counts / 2, 
  labels = counts,
  col = "black",  
  cex = 0.8,      
  font = 2        
)

##############
# STEP 3 : exporting results
##############

# Create a data frame with the final DEG list
results   <- data.frame(genes = genes, logFC = logFC, pval = pval, pval_adj = pval_adj)
direction <- ifelse(logFC > 0, "UP", "DOWN")
results   <- cbind(results, direction)

# Export to disk
write.table(results,   file = filename_DEG,        row.names = FALSE, sep = "\t", quote = FALSE)
write.table(data,      file = filename_matrix_DEG, row.names = TRUE,  col.names = NA, sep = "\t", quote = FALSE)

# After exporting the full DEG table, also export separate up/down gene lists
write.table(
  results$genes[results$direction == "UP"],
  file = file.path(dirDataset, "up.txt"),
  row.names = FALSE, col.names = FALSE, quote = FALSE
)
write.table(
  results$genes[results$direction == "DOWN"],
  file = file.path(dirDataset, "down.txt"),
  row.names = FALSE, col.names = FALSE, quote = FALSE
)

###########################
# STEP 4 : PLOTS 
###########################
# 4.0 Describe cohort: age distribution
hist(
  metadata$`age:ch1`,
  main   = "Age Distribution",
  xlab   = "Age (years)",
  ylab   = "Count",
  breaks = 10,
  col    = "steelblue"
)

# 4.1 Describe cohort: gender breakdown
# Count the number of males and females
gender_counts <- table(metadata$`gender:ch1`)
# bar plot
bar_heights <- barplot(
  gender_counts,
  main = "Gender Distribution",
  ylab = "Count",
  col = c("#BB0044", "#5555AA")
)
# Add the gender counts inside the bars
text(
  x = bar_heights,
  y = gender_counts / 2,
  labels = gender_counts,
  col = "black",
  cex = 0.8,
  font = 2
)

# Subset data to include only 'case' samples
case_data <- subset(metadata, `status:ch1` == "schizophrenia (SCZ)")
gender_counts <- table(case_data$`gender:ch1`)
# xreate a bar plot for gender distribution among cases
bar_heights <- barplot(
  gender_counts,
  main = "Gender Distribution (Cases)",
  ylab = "Count",
  col = c("#BB0044", "#5555AA")
)
# Add the gender counts inside the bars
text(
  x = bar_heights,
  y = gender_counts / 2,
  labels = gender_counts,
  col = "black",
  cex = 0.8,
  font = 2
)


# 4.2 Pie chart of disease distribution
status_counts <- table(metadata$`status:ch1`)
pie(
  status_counts,
  labels = paste0( names(status_counts), " ", round(100 * status_counts / sum(status_counts), 2), "%" ),
  col = c("grey", "red2"),
  main = "Case vs Control"
)

# 4.3 Final volcano plot (adjusted p-values used)
volcano_cols <- ifelse(results$direction == "UP", "red", "blue")
plot(
  results$logFC, -log10(results$pval_adj),
  main  = "Volcano Plot",
  xlab  = "log2 fold change",
  ylab  = "-log10 adjusted p-value",
  ylim = c(0, max(-log10(results$pval_adj))),
  col   = volcano_cols,
  pch   = 16,
  cex   = 0.7
)
abline(h = -log10(thr_pval), col = "grey", lty = 2, lwd = 2)
abline(v = c(-thr_logFC, thr_logFC), col = "grey", lty = 2, lwd = 2)

# 4.4 Boxplot for top up- and down-regulated genes
if (nrow(dataN) > 0 && nrow(dataC) > 0) {
  ind <- which.max(results$logFC)
  gene_id <- results$genes[ind]
  boxplot(
    as.numeric(dataN[ind, ]), as.numeric(dataC[ind, ]),
    main   = paste(gene_id, ", adj. p-value =", format(results$pval_adj[ind], digits = 2)),
    notch  = FALSE,
    ylab   = "Gene expression value",
    xlab   = "Most up-regulated gene",
    col    = c("yellow", "red"),
    names  = c("Normal", "Case"),
    pars   = list(boxwex = 0.3, staplewex = 0.6)
  )
  stripchart(
    list(as.numeric(dataN[ind, ]), as.numeric(dataC[ind, ])),
    method = "jitter",
    pch    = 1,
    col    = "blue4",
    vertical = TRUE,
    add     = TRUE
  )
  
  
  ind <- which.min(results$logFC)
  gene_id <- results$genes[ind]
  boxplot(
    as.numeric(dataN[ind, ]), as.numeric(dataC[ind, ]),
    main   = paste(gene_id, ", adj. p-value =", format(results$pval_adj[ind], digits = 2)),
    notch  = FALSE,
    ylab   = "Gene expression value",
    xlab   = "Most down-regulated gene",
    col    = c("yellow", "red"),
    names  = c("Normal", "Case"),
    pars   = list(boxwex = 0.3, staplewex = 0.6)
  )
  stripchart(
    list(as.numeric(dataN[ind, ]), as.numeric(dataC[ind, ])),
    method = "jitter",
    pch    = 1,
    col    = "blue4",
    vertical = TRUE,
    add     = TRUE
  )
  rm(ind)
}

# 4.5 Pie chart of UP vs DOWN
count <- table(results$direction)
pie(
  count,
  labels = paste0(names(count), " ", round(100 * count / sum(count), 2), "%"),
  col    = c("blue", "gold")
)

# 4.6 Heatmap of DEGs
annotation <- data.frame(status = metadata$'status:ch1')
rownames(annotation) <- metadata$geo_accession
vect_color <- c("CONTROL" = "purple", "schizophrenia (SCZ)" = "red")
annotation_colors <- list(status = vect_color)

pheatmap(
  data,
  scale                   = "row",
  border_color            = NA,
  cluster_cols            = TRUE,
  cluster_rows            = TRUE,
  clustering_distance_rows= "correlation",
  clustering_distance_cols= "correlation",
  clustering_method       = "average",
  annotation_col          = annotation,
  annotation_colors       = annotation_colors,
  color                   = colorRampPalette(c("blue","blue3","black","yellow3","yellow"))(100),
  show_rownames           = FALSE,
  show_colnames           = FALSE,
  cutree_cols             = 2,
  cutree_rows             = 2,
  width                   = 10,
  height                  = 10,
  filename                = filename_heatmap
)
