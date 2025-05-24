rm(list = ls())
options(stringsAsFactors = FALSE)

library(ggbiplot)
library(qcc)
library(ggpubr)
library(factoextra)
library(corrplot)
library(FactoMineR)
library(RColorBrewer)
library(survminer)

################################################
# 0. Setup: results folder and dataset naming
################################################
dirRes       <- "Results/"
if (!dir.exists(dirRes)) dir.create(dirRes)

dataset      <- "GEO_Schizophrenia"
dirDataset   <- file.path(dirRes, dataset, "/")
if (!dir.exists(dirDataset)) dir.create(dirDataset)

plots      <- "Plots" 
dirPlots   <- file.path(dirRes, plots, "/")
if (!dir.exists(dirPlots)) dir.create(dirPlots)

# Output file names
file_score_plot          <- file.path(dirPlots, "score_plot.pdf")
file_pareto_scree_plot   <- file.path(dirPlots, "pareto_scree_plot.pdf")
file_loading_plot        <- file.path(dirPlots, "loading_plot.pdf")
file_contribution_plot   <- file.path(dirPlots, "PC_contribution_plot.pdf")

################################################
# 1. Importing data
################################################
# DEG matrix produced in your previous step (genes × samples)
data_raw <- read.table(
  file.path(dirDataset, "matrix_DEG.txt"),
  header     = TRUE, sep = "\t", quote = "",
  check.names = FALSE, row.names = 1
)

# Which samples are controls / cases?
control_ids <- scan( file.path(dirDataset, "control.txt"), what = character(), sep = "\n" )
case_ids    <- scan( file.path(dirDataset, "case.txt"), what = character(), sep = "\n" )

# Transpose so rows = samples, columns = genes and put cases first, then controls
data <- t(data_raw[, c(case_ids, control_ids)])

# Create a grouping factor for plotting
groups <- c( rep("schizophrenia (SCZ)", length(case_ids)), rep("CONTROL", length(control_ids)) )

################################################
# 2. Apply PCA
# Rows = observations, columns = variables
################################################
pca <- prcomp(data,
              center = TRUE,
              scale. = TRUE,
              retx   = TRUE)

################################################
# 3. Score plot (PC1 vs PC2)
################################################
pdf(file_score_plot, width = 5, height = 5)
g <- ggbiplot(pca,
              obs.scale = 1,
              var.axes  = FALSE,
              ellipse   = TRUE,
              groups    = groups) +
  theme_minimal() +
  labs(title = "PCA score plot: cases vs controls")
print(g)
dev.off()

# Alternative with factoextra
fviz_pca_ind(pca,
             col.ind     = groups,
             addEllipses = TRUE,
             repel       = TRUE)

################################################
# 4. Pareto / Scree plot
################################################
eigenvalue  <- pca$sdev^2
varS        <- round(eigenvalue / sum(eigenvalue) * 100, 2)
names(varS) <- paste0("PC", seq_along(varS))

pdf(file_pareto_scree_plot, width = 5, height = 5)
pareto.chart(varS[1:10], main = "Pareto chart: % variance explained (first 10 PCs)")
fviz_eig(pca, addlabels = TRUE) +
  ggtitle("Scree plot")
dev.off()

################################################
# 5. Loadings plot
################################################
pdf(file_loading_plot, width = 5, height = 5)
fviz_pca_var(pca,
             col.var     = "cos2", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel       = TRUE) +
  ggtitle("Variable loadings (cos2) on PC1–PC2")
dev.off()

################################################
# 6. Contributions of variables to PCs
################################################
contrib_var <- get_pca_var(pca)$contrib
colnames(contrib_var) <- paste0("PC", seq_len(ncol(contrib_var)))
contrib_var <- contrib_var[order(contrib_var[,"PC1"], decreasing = TRUE), ]

pdf(file_contribution_plot, width = 10, height = 10)
corrplot(contrib_var[1:10, 1:10],
         is.corr     = FALSE,
         tl.col      = "black",
         method      = "color",
         col         = brewer.pal(n = 9, name = "BuPu"),
         addCoef.col = "black",
         tl.cex      = 0.8)

# Bar‐plots of top contributors
fviz_contrib(pca, choice = "var", axes = 1, top = 10) +
  ggtitle("Top 10 contributors to PC1")
fviz_contrib(pca, choice = "var", axes = 2, top = 10) +
  ggtitle("Top 10 contributors to PC2")
fviz_contrib(pca, choice = "var", axes = 3, top = 10) +
  ggtitle("Top 10 contributors to PC3")
fviz_contrib(pca, choice = "var", axes = 1:3, top = 15) +
  ggtitle("Top 15 contributors to PC1–PC3")
dev.off()
