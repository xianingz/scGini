#' Neighborhood Smoothing of the scGini
#'
#' This functions implements a neighborhood difussion map to improve the scGini prediction
#' by utilizing information from neighboring cells.
#'
#' @param dat Single-cell counts matrix with genes as rows and cells as columns.
#' @param scgini Corrected Gini index, results from `scGini`
#' @param nGenes Number of genes used for neighborhood smoothing
#' @return A vector with smoothed scGini index is returned.
#' @export
scGiniSmooth <- function(mat, scgini, nGenes=1000){
  expmat <- normalize(mat)
  #Filter out cells not expressing any of the 1000 most variable genes
  hvgmat <- hvg(expmat, nGenes)
  rm1 <- colSums(hvgmat) == 0
  hvgmat <- hvgmat[, !rm1]
  scgini <- scgini[!rm1]
  D <- similarity_matrix_cleaned(HiClimR::fastCor(hvgmat))
  cgini_smoothed <- diffused(D, scgini)
  return(cgini_smoothed)
}

