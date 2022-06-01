





## Utils
MatVar <- function(x) {
  "
  Function to calculate the individual varaince for samples in each row

  Input:
    x = count matrix

  Output:
    Returns matrix of z-scores by row
  "
  a <- rowSds(x)
  b <- rowMeans(as.matrix(x))
  v <- (x - b) / a
  return(v)
}


geneset_score <- function(counts_matrix, geneset, stat = "Mean") {
  "
  Function to calculate the mean z-score of a gene set for each cell

  Inputs:
    counts_matrix = rows of genes by columns of cells with gene names as rownames and cells as colnames
    geneset = list of genes that you want to calculate the mean score for
    stat = Return the mean or the sum of the gene z-scores

  Output:
    out = A list of the geneset expression for each cell
  "

  ol <- which(tolower(rownames(counts_matrix)) %in% tolower(geneset))
  if (length(ol) > 0) {
    counts <- as.matrix(counts_matrix[ol, ])
  }else{
    stop("No genes from this gene set are found in your counts matrix")
  }

  matvar <- MatVar(counts)

  n <- which(is.na(matvar[,1]))
  if (length(n) > 0) {
    matvar <- matvar[-n,]
  }

  if (tolower(stat) == "mean") {
    out <- colMeans(matvar)
  }else if (tolower(stat) == "sum"){
    out <- colSums(matvar)
  }else{
    stop("stat must be one of the following; c('Mean', 'Sum')")
  }

  return(out)
}

compare_geneset_signatures_seurat <- function(seurat_obj, geneset_1, geneset_2, assay = "RNA",stat = "mean") {
  "
  Function to compare two geneset signatures in a Seurat scRNA seq data set

  Inputs:
    seurat_obj = scRNA seq object with dimensional reduction completed
    geneset_1 = list of gene names for a given gene set signature
    geneset_2 = list of gene names for a given gene set signature
    assay = With assay from seurat object to pull @data from
    stat = Passing argument to geneset_score() for what statistic to calculate

  Outputs:
    Returns seurat_obj with signatures and rank in meta.data
  "

  if (tolower(stat) %in% c("mean", "stat")) {
    if (assay %in% names(seurat_obj@assays)) {
      t1_sig <- geneset_score(seurat_obj@assays[[assay]]@data, geneset_1, stat = stat)
      t2_sig <- geneset_score(seurat_obj@assays[[assay]]@data, geneset_2, stat = stat)
    }else{
      stop("assay must be in names(seurat_object@assays)")
    }
  }else{
    stop("stat must be one of the following; c('mean', 'stat')")
  }

  rdf <- data.frame(T1_sig = t1_sig,
                    T2_sig = t2_sig,
                    Z_diff = t1_sig - t2_sig)

  rdf <- rdf[order(rdf$Z_diff, decreasing = F), ]
  rdf$Rank <- c(1:nrow(rdf))

  rep_c <- c("T1_sig", "T2_sig", "Z_diff", "Rank")
  for (k in rep_c) {
    if (k %in% colnames(seurat_obj@meta.data)) {
      seurat_obj@meta.data <- seurat_obj@meta.data[, -which(colnames(seurat_obj@meta.data) %in% rep_c)]
    }
  }

  rdf <- rdf[match(rownames(seurat_obj@meta.data), rownames(rdf)),]
  seurat_obj@meta.data <- cbind(seurat_obj@meta.data, rdf)

  return(seurat_obj)
}



plot_z_diff_seurat <- function(seurat_obj, add_rug = "", color_values = c("")) {
  "
  Function to plot individual signatures and Z-diff

  Input:
    seurat_obj =  scRNA seq object after running compare_geneset_signatures_seurat()
    add_rug = Meta.data value to plot on geom_rug bottom
    color_values = color values associated to the levels of 'add_rug'

  Returns:
    Ggplot formatted plot of the individual signatures and Z-diff
  "

  b <- ggplot(seurat_obj@meta.data, aes(x = Rank, y = Z_diff)) +
    geom_point() +
    theme_classic() +
    ggtitle("T1 - T2 Signature") +
    geom_hline(yintercept = 0, linetype = "dotted") +
    theme(text = element_text(face = "bold", size = 15), plot.title = element_text(hjust = 0.5), line = element_line(size = 1))

  if (add_rug != "") {
    if (add_rug %in% colnames(seurat_obj@meta.data)) {
      b <- b + geom_rug(sides = "b", aes_string(color = add_rug))
    } else{
      stop("add_rug must be in colnames(seurat_obj@meta.data")
    }
  }

  if ((add_rug != "") & (length(color_values) > 1)) {
    if (length(color_values) == length(unique(seurat_obj@meta.data[[add_rug]]))) {
      b <- b + scale_color_manual(values = color_values)
    }else{
      stop("Length of color_values must be the same as number of levels being plotted from add_rug")
    }
  }

  a <- ggplot(seurat_obj@meta.data, aes(x = Rank)) +
    geom_point(aes(y = T1_sig, color = "T1_sig")) +
    geom_point(aes(y = T2_sig, color = "T2_sig")) +
    scale_color_manual(values = c("orange", "dark green")) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    theme_classic() +
    theme(text = element_text(face = "bold", size = 15), plot.title = element_text(hjust = 0.5), line = element_line(size = 1))

  plot <- ggarrange(a,b, ncol = 1)
  return(plot)
}



plot_dimred_seurat <- function(seurat_obj, dimred = "umap", color = "", color_values = c("")) {
  "
  Function to plot dimensional reduction colored by meta.data column and rank

  Inputs:
    seurat_obj = scRNA seq object after running compare_geneset_signatures_seurat()
    dimred = Which dimensional reduction to plot; c('pca', 'tsne', 'umap')
    color = Meta.data value to color the points
    color_values = color values associated to the levels of 'color'

  Outputs:
    Returns ggplot formatted plot of dimensional reduction colored by meta.data column and Z_diff rank
  "

  if (dimred %in% c('pca', 'tsne', 'umap')) {
    if (dimred %in% names(seurat_obj@reductions)) {
      red_coords <- as.data.frame(seurat_obj@reductions[[dimred]]@cell.embeddings)
      red_coords <- cbind(red_coords, seurat_obj@meta.data[match(rownames(red_coords), rownames(seurat_obj@meta.data)),])
    }else{
      stop("Run dimensional reduction on Seurat object before trying to plot these coordinates")
    }
  }else{
    stop("dimred must be one of the following; c('pca', 'tsne', 'umap')")
  }


  if (color == "") {
    c <- ggplot(red_coords, aes_string("UMAP_1", "UMAP_2")) +
      geom_point() +
      theme_classic() +
      theme(text = element_text(face = "bold", size = 15), plot.title = element_text(hjust = 0.5), line = element_line(size = 1))
  }else if (color %in% colnames(red_coords)) {
    c <- ggplot(red_coords, aes_string("UMAP_1", "UMAP_2", color = color)) +
      geom_point() +
      theme_classic() +
      theme(text = element_text(face = "bold", size = 15), plot.title = element_text(hjust = 0.5), line = element_line(size = 1))
  }else{
    stop("color must be in colnames(seurat_obj@meta.data) in order to plot")
  }


  if (length(color_values) != 1) {
    if (length(color_values) == length(unique(red_coords[[color]]))) {
      c <- c + scale_color_manual(values = color_values)
    }else{
      stop("Length of color_values must be the same as number of levels being plotted from color")
    }
  }

  d <- ggplot(red_coords, aes(UMAP_1, UMAP_2, color = Rank)) +
    geom_point() +
    theme_classic() +
    scale_color_viridis_c(option = "inferno") +
    theme(text = element_text(face = "bold", size = 15), plot.title = element_text(hjust = 0.5), line = element_line(size = 1))

  plot <- ggarrange(c,d, ncol = 1)
  return(plot)
}




plot_enrichment_seurat <- function (seurat_obj, metadata_col, figure_header = "", color_values = c(""), abs = F, pval = 0.05)
{
  "
  Function to calculate enrichment of a metadata level for signature 1 vs signature 2

  Inputs:
    seurat_obj = scRNA seq object with dimensional reduction completed
    metadata_col = name of meta.data column to be used for geneset enrichment analysis
    figure_header = character string to add figure_header to output figure

  Output:
    prints statistical results of fgsea and will plot enrichment curves if statistically significant enrichment is found.
  "

  if (("Z_diff" %in% colnames(seurat_obj@meta.data)) == F) {
    stop("Run compare_geneset_signatures_seurat() to calculate Z_diff before running enrichment analysis")
  }

  RNK_order <- order(abs(seurat_obj@meta.data$Z_diff), decreasing = T)
  RNK <- seurat_obj@meta.data$Z_diff[RNK_order]
  names(RNK) <- rownames(seurat_obj@meta.data)[RNK_order]
  if ((metadata_col %in% colnames(seurat_obj@meta.data)) ==
      F) {
    stop("metadata_col must be in colnames(seurat_obj@meta.data)")
  }
  genesets <- list()
  for (i in unique(seurat_obj@meta.data[[metadata_col]])) {
    genesets[[i]] <- rownames(seurat_obj@meta.data)[which(seurat_obj@meta.data[[metadata_col]] ==
                                                            i)]
  }
  sets <- cust_sets(genesets)
  res <- fgsea(genesets, RNK, nperm = 10000)
  print(res)
  plot <- enrichment_analysis_dev(geneset_list = genesets, fgsea_RNK = RNK,
                                  msigdb_sets = sets, figure_header = figure_header, color_values = color_values, abs = abs, pval = pval)
  return(plot)
}


























































