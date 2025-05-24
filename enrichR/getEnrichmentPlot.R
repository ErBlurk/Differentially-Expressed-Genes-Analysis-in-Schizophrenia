getEnrichmentPlot <- function(annotation, type, top_term, thr_pval, dirEnrich) {
  
  # Filter by p-value
  annotation <- annotation[annotation$Adjusted.P.value < thr_pval, ]
  annotation <- annotation[order(annotation$Adjusted.P.value, decreasing = FALSE), ]
  
  # exit if no meaningful values are found
  if (nrow(annotation) == 0) {
    message(paste0("No enrichment terms passed the threshold for ", type))
    return(NULL)
  }
  
  # number of genes
  annotation$Gene_count <- sapply(annotation$Genes, function(x) {
    tmp <- unlist(strsplit(x, split = ";"))
    length(tmp)
  })
  
  # Gene_ratio
  annotation$Gene_ratio <- unlist(lapply(annotation$Overlap, function(x) {
    parts <- as.numeric(unlist(strsplit(x, "/")))
    parts[1] / parts[2]
  }))
  
  # top ranking terms
  if (length(top_term) != 0 && top_term <= nrow(annotation)) {
    annotation_top <- annotation[1:top_term, ]
  } else {
    annotation_top <- annotation
  }
  
  # bar plot
  g1 <- ggplot(annotation_top, 
               aes(x = Gene_count, 
                   y = fct_reorder(Term, Gene_count), 
                   fill = Adjusted.P.value)) +
    geom_bar(stat = "identity") +
    scale_fill_continuous(low = "red", high = "blue", name = "Adjusted.P.value",
                          guide = guide_colorbar(reverse = TRUE)) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
    theme_bw(base_size = 10) +
    ylab(NULL)
  
  pdf(file.path(dirEnrich, paste0(type, "_barplot.pdf")))
  print(g1)
  dev.off()
  
  # dot plot
  g2 <- ggplot(annotation_top, 
               aes(x = Gene_count, 
                   y = fct_reorder(Term, Gene_count))) +
    geom_point(aes(size = Gene_ratio, color = Adjusted.P.value)) +
    scale_colour_gradient(low = "red", high = "blue") +
    theme_bw(base_size = 10) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
    ylab(NULL)
  
  pdf(file.path(dirEnrich, paste0(type, "_dotplot.pdf")))
  print(g2)
  dev.off()
  
  # Export results
  write.table(annotation[, c("Term", "Overlap", "P.value", "Adjusted.P.value", 
                             "Gene_count", "Gene_ratio", "Genes")],
              file.path(dirEnrich, paste0(type, "_adj_pval_", thr_pval, ".txt")),
              sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
}
