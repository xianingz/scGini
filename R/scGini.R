#' Calculate corrected Gin index
#'
#' This function calculates corrected Gini index from single cell RNA-seq counts data
#' without any preprocessings.
#'
#' @param mat Single-cell counts matrix with genes as rows and cells as columns.
#' @param ncore Number of cores for parallel computing, default is 1.
#' @return A data.frame with estimated parameters and directly calculated Gini index (Gini)
#' and corrected Gini index (cGini).
#' @export
scGini <- function(mat, ncore=1){
  tl <- dim(mat)[2]
  lr=1e-2
  res <- pbmcapply::pbmclapply(c(1:tl), function(i){
    exp = data.frame(mat[,i])
    colnames(exp) <- colnames(mat)[i]
    gnf.nb(exp, lr)
  }, mc.cores = ncore, mc.preschedule = TRUE)
  df <- c()
  for(i in res){
    df <- rbind(df, i)
  }
  return(df)
}

