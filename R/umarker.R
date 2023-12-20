#' Extract expression matrix from Seurat object or matrix
#'
#' This function extracts the expression matrix from either a Seurat object or a matrix object. If a Seurat object is provided, the function can optionally specify the assay and slot to be used for extracting the expression matrix. If a matrix is provided, the function will create a Seurat object with the provided matrix and default metadata fields. The function can also specify the minimum number of cells and features required to keep genes in the expression matrix.
#'
#' @name assay2mtx
#' @author Qiong Zhang
#' @param object A Seurat object or matrix containing expression data
#' @param assay A character string specifying the assay name to extract expression data from the Seurat object
#' @param slot A character string specifying the slot name to extract expression data from the Seurat object
#' @param min.cells An integer specifying the minimum number of cells required to keep genes in the expression matrix
#' @param min.features An integer specifying the minimum number of features required to keep genes in the expression matrix
#' @return A matrix containing the expression data, with gene names as row names and cell barcodes as column names
#' @examples
#' expmtx <- assay2mtx(object = seurat_object, assay = "RNA", slot = "counts", min.cells = 10, min.features = 5)
#' expmtx <- assay2mtx(object = expression_matrix, min.cells = 10, min.features = 5)
#'
assay2mtx <- function(object, assay = NULL, slot = "counts", min.cells = 3, min.features = 1, update = F, excludeRM = T) {
  if (any(methods::is(object) %in% "Seurat")) {
    assay <- if (is.null(assay)) Seurat::DefaultAssay(object) else assay
    object <- if (isTRUE(update)) Seurat::UpdateSeuratObject(object) else object
  } else {
    object <- Seurat::CreateSeuratObject(counts = object, project = "uMarker", assay = assay, min.cells = min.cells, min.features = min.features)
  }
  .expMat <- Seurat::GetAssayData(object, assay = assay, slot = slot)
  if (isTRUE(excludeRM)){
    .genes.all <- rownames(.expMat)
    .qual_qc_terms.h <- tibble::tibble(word_reg = c("^RP[SL]", "^MT-", "^HB[^(P)]"), term_name = c("percent_ribo", "percent_mito", "percent_hb"))
    .qual_qc_terms.m <- tibble::tibble(word_reg = c("^Rp[sl]", "^Mt-", "^Hb[^(p)]"), term_name = c("percent_ribo", "percent_mito", "percent_hb"))
    .genes.filter.preset <- grep(stringr::str_c(c(.qual_qc_terms.h$word_reg, .qual_qc_terms.m$word_reg), collapse = "|"), .genes.all, value = T)
    .expMat <- .expMat[!(.genes.all %in% .genes.filter.preset), ]
  }
  return(.expMat)
}


#' Transform expression matrix to binary type based on give mode or customize value
#'
#' @param mat vector
#' @param pattern transform pattern
#' @name expMat2Bin
#' @return binary matrix
#' @author Qiong Zhang
#'
expMat2Bin <- function(mat, pattern = "min", cutoff = NULL) {
  .stabe2b <- function(mat, pattern, cutoff = NULL) {
    if (pattern == "rowMean") {
      .x.m <- Matrix::rowMeans(mat)
      return(list(mat = mat > .x.m[rownames(mat)], cutoff = cutoff))
    }
    .x <- mat@x
    if (is.null(cutoff)) {
      .cutoff <- switch(pattern,
                        "min" = min(.x),
                        "max" = max(.x),
                        "mean" = mean(.x),
                        "median" = median(.x),
                        stop("Error: pattern need to be one of 'min', 'max', 'median', 'mean' and 'adpmed'. ")
      )
    } else {
      .cutoff <- as.numeric(cutoff)
    }
    .index <- which(.x > .cutoff)
    .dp <- diff(mat@p)
    .col.index <- (rep(seq_along(.dp), .dp))[.index]
    .row.index <- (mat@i + 1)[.index]
    .mat.new <- Matrix::sparseMatrix(.row.index[1], .col.index[1], dims = mat@Dim)
    .mat.new[cbind(.row.index, .col.index)] <- 1
    .mat.new <- as(.mat.new, "dgCMatrix")
    rownames(.mat.new) <- rownames(mat)
    return(list(mat = .mat.new, cutoff = .cutoff))
  }
  .adpmed <- function(mat) {
    .m.multi <- Mclust(Matrix::rowSums(mat), 1:30, modelNames = c("V"), verbose = FALSE)
    if (.m.multi$G == 1) {
      return(.stabe2b(mat, type = "median"))
    }
    .mat.new <- list()
    .cutoff <- rep(0, .m.multi$G)
    for (i in seq_along(1:.m.multi$G)) {
      .x.clu <- mat[which(.m.multi$classification == i), ]
      .cutoff[i] <- median(.x.clu@x)
      .index <- which(.x.clu@x >= .cutoff[i])
      .dp <- diff(.x.clu@p)
      .col.index <- (rep(seq_along(.dp), .dp))[.index]
      .row.index <- (.x.clu@i + 1)[.index]
      .mat.new[[i]] <- Matrix::sparseMatrix(.row.index[1], .col.index[1], dims = .x.clu@Dim)
      .mat.new[[i]][cbind(.row.index, .col.index)] <- 1
      .mat.new[[i]] <- as(.mat.new[[i]], "dgCMatrix")
      rownames(.mat.new[[i]]) <- rownames(.x.clu)
    }
    .mat.new <- do.call(rbind, .mat.new)
    return(list(mat = .mat.new, cutoff = .cutoff))
  }
  if (pattern == "adpmed") {
    .exp2bin <- .adpmed(mat)
  } else {
    .exp2bin <- .stabe2b(mat = mat, pattern = pattern, cutoff = cutoff)
  }
  return(.exp2bin)
}

#' rescale a vector between 0 (lowest) and 1 (highest)
#'
#' @param x vector
#' @name scaleZO
#' @return rescale vector with value between 0 and 1
#' @author Qiong Zhang
#'
scaleZO <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

#' calculate prior probabilities for clusters
#'
#' @param x vector
#' @name calProb4C
#' @return probabilities of clusters
#' @author Qiong Zhang
#'
calProb4C <- function(x) {
  table(x) / length(x)
}

#' calculate prior probabilities for genes
#'
#' @param mat expression matrix
#' @name calProb4G
#' @return probabilities for genes
#' @author Qiong Zhang
#'
calProb4G <- function(mat) {
  .prob4G <- as.vector((Matrix::rowSums(mat)) / ncol(mat))
  names(.prob4G) <- rownames(mat)
  return(.prob4G)
}

#' calculate probabilities genes in clusters
#'
#' @param mat expression matrix
#' @param cluProb cluster
#' @param cluID
#' @param cpu
#' @name calProbGIC
#' @return probabilities of genes in clusters
#' @author Qiong Zhang
#'
# get gene probability conditional on cluster
calProbGIC <- function(mat, cluProb, cluID) {
  .clusterID <- names(cluProb)
  .GIC <- do.call(cbind, lapply(1:length(.clusterID), function(i) Matrix::rowMeans(mat[, which(cluID == .clusterID[i])])))
  .GIC <- as(.GIC, "dgCMatrix")
  colnames(.GIC) <- .clusterID
  return(.GIC)
}

#' calculate probabilities of genes in given condition
#'
#' @param mat expression matrix
#' @param cluProb
#' @param cluID
#' @param cpu
#' @name calProbCIG
#' @return probabilities of genes in given condition
#' @author Qiong Zhang
#'
calProbCIG <- function(binmat, probG, ProbC) {
  .prob <- t(apply(log(binmat), 1, function(x) x + log(ProbC)))
  .prob <- exp(apply(.prob, 2, function(x) x - log(probG)))
  colnames(.prob) <- names(ProbC)
  return(as(.prob, "dgCMatrix"))
}

#' calculate the feature score of genes
#'
#' @param ginc gene expression probability
#' @param binmat geneXcell binary matrix
#' @param scale scale the scores between groups
#' @name calGeneFea
#' @return feature score of genes
#' @author Qiong Zhang
#'
calGeneFea <- function(ginc, binmat, scale = F) {
  if (isTRUE(scale)) {
    .genefeature <- apply(ginc * binmat, 2, scaleZO)
  } else {
    .genefeature <- exp(log(ginc) + log(binmat))
  }
  colnames(.genefeature) <- colnames(binmat)
  return(.genefeature)
}


#' QC for clusterIDs
#'
#' @param mat  gene expression matrix
#' @param clusterIDs a vector contains clusterID for cells
#' @name clusidQC
#' @return validated clusterIDs
#' @author Qiong Zhang
#'
clusidQC <- function(mat, clusterIDs) {
  if (any(table(clusterIDs) < 2)) stop("Several clusters comprised of only one cell. Please check clusterIDs")
  if (length(levels(clusterIDs)) < 2) stop("Only one cluster exists in clusterIDs!")
  .cluIDs.len <- length(clusterIDs)
  if (ncol(mat) != .cluIDs.len) stop("Inequal cell number and clusterIDs.")
  return(droplevels(clusterIDs))
}


#' Calculate marker specificity
#'
#' @param mat expression matrix
#' @param clusterIDs cluster IDs for each cell.
#' @param binmethod Method used to binarize the matrix
#' @param expR The proportion of feature expressed cells in a cluster.
#' @param cpu How many cores will be used for detection of markers
#' @param marker The number of markers returned.
#' @name LLRMarker
#' @author Qiong Zhang
#' @return data.frame with marker features of each cluster.
#'
LLRMarker <- function(mat, clusterIDs, binmethod = "min", expR = 0.1, model = 1, marker = 20, cutoff = 3, scType = 'scRNA') {
  clusterIDs <- clusidQC(mat, clusterIDs)
  .clusterIDs.lev <- levels(clusterIDs)
  .mat <- as(mat, "dgCMatrix")
  .cutoff <- cutoff
  if (scType == 'scATAC') .cutoff <- 0
  .mat.bin <- expMat2Bin(.mat, pattern = binmethod, cutoff = .cutoff)
  .clu.prob <- calProb4C(as.integer(clusterIDs))
  names(.clu.prob) <- .clusterIDs.lev
  .binmat <- calProbGIC(.mat.bin$mat, .clu.prob, clusterIDs)
  .binmat.f.l <- apply(.binmat, 1, function(x) any(x > expR))
  .binmat.f <- .binmat[.binmat.f.l, ]
  .gene.prob <- calProb4G(.mat.bin$mat[.binmat.f.l, ])
  .ginc <- calProbCIG(.binmat.f, .gene.prob, .clu.prob)
  .gfeat <- calGeneFea(.ginc, .binmat.f)
  .w <- replicate(ncol(.gfeat), Matrix::rowSums(.gfeat))
  .gfeat.n <- .gfeat^3 / .w^2
  .gfeat.s <- sqrt(.gfeat / (.w - .gfeat))
  .marker.n <- markerComb(matbin = .mat.bin$mat, gfeat = switch(model, .gfeat.n, .gfeat.s), clusterIDs = clusterIDs, marker = marker)
  return(list(marker = .marker.n, pct = .binmat))
}


#' calulate the coverage of cell for marker gene combination in given cluster
#'
#' @author Qiong Zhang
#' @name markerComb
#' @param matbin binary matrix
#' @param gfeat feature scores
#' @param clusterIDs Identity of cells
#' @param marker marker number
#' @return
#'
markerComb <- function(matbin, gfeat, clusterIDs, marker) {
  .clusters <- levels(clusterIDs)
  .cell.in.cluster.index.lst <- split(1:length(clusterIDs), clusterIDs)
  .marker.info <- lapply(1:length(.clusters), function(x, marker) {
    .cand.marker <- names(sort(gfeat[, x], decreasing = T))[1:100]
    .cell.index <- .cell.in.cluster.index.lst[[x]]
    .cell.n.cluster <- length(.cell.index)
    .gene.index <- match(.cand.marker, rownames(matbin))
    .mat.bin.cand <- matbin[.gene.index, .cell.index]
    .uncov.marker <- sapply(Reduce(function(x,y) x + y, Matrix::t(.mat.bin.cand), accumulate = T), function(s,l) sum(!s)/l, l = .cell.n.cluster)
    if (identical(marker, 'Auto')){
      .p <- which(.uncov.marker < 0.05)[1]
      .p <- ifelse(.p < 10, 10, .p)
      marker <- .p
    }
    return(data.frame(markers = .cand.marker[1:marker], loss = .uncov.marker[1:marker]))
  }, marker = marker)
  names(.marker.info) <- .clusters
  return(.marker.info)
}


#' Recommend marker genes for expression matrix or Seurat object
#'
#' @param obj Seurat object or expression matrix
#' @param clusterIDs The cluster ID for each cell or the column name which contains cluster ID information in Seurat object.
#' @param binmethod Method used to binarize the matrix
#' @param expR The proportion of feature expressed cells in a cluster.
#' @param cpu How many cores will be used for detection of markers
#' @param marker The number of markers returned.
#' @name markerREC
#' @export
#' @author Qiong Zhang
#' @return data.frame of marker features in each cluster.
#' @examples marker.mr <- markerREC(pbmc3k.final)
#'
umarker <- function(obj, scType = 'scRNA', assay = "RNA", slot = "counts", binmethod = "min", expR = 0.3, model = 1, marker = 20, clusterIDs = NULL, cutoff = NULL, excludeRM = T) {
  .mat <- assay2mtx(object = obj, assay = assay, slot = slot, excludeRM = excludeRM)
  if (is.null(clusterIDs)) {
    .clusterIDs <- Seurat::Idents(obj)
  } else {
    if (length(clusterIDs) != ncol(.mat)) .clusterIDs <- obj@meta.data[[clusterIDs]]
  }
  .markers <- LLRMarker(mat = .mat, clusterIDs = .clusterIDs, binmethod = binmethod, expR = expR, model = model, marker = marker, cutoff = cutoff, scType = scType)
  return(.markers)
}



