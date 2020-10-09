#' Plot the goodness of fitting
#'
#' This function fit the Negative Binomial model on one single cell and plot the
#' frequecy of oberseved counts and the theoretical probability from the NB model.
#' This function helps evaluate the goodness of fitting of the model to the dataset
#' before applying.
#'
#' @param dat Expression counts vector for a single cell.
#' @param seed Random seed for the simulated curve from theoretical models.
#' @return A plot with observed and theoretical distribution of the expression counts.
#' @export
NBFitPlot <- function(dat,seed=1){
  dat.fre <- as.data.frame(table(dat))
  dat.fre$prob <- dat.fre$Freq/sum(dat.fre$Freq)
  dat.est <- est.nb(dat)
  dat.fre$Exp <- as.numeric(as.character(dat.fre$dat))
  est.prob <- unlist(lapply(as.numeric(as.character(dat.fre$dat)), function(x) gampoi(x, dat.est[2], dat.est[1])))
  dat.fre$est.prob <- est.prob
  set.seed(seed)
  sim = rnbinom(length(dat), dat.est[2], 1-dat.est[1])
  sim.fre <- as.data.frame(table(sim))
  sim.fre$est.prob <- sim.fre$Freq/sum(sim.fre$Freq)
  sim.fre$Exp <- as.numeric(as.character(sim.fre$sim))
  p <- ggplot2::ggplot(dat.fre) +
    ggplot2::theme(panel.background = element_rect(fill = "white", colour = "black")) +
    ggplot2::geom_point(aes(log(Exp+1), prob), col="blue",alpha=0.5) +
    ggplot2::geom_point(data=sim.fre, aes(log(Exp+1),est.prob),col="red", alpha=0.5)+
    ggplot2::stat_smooth(method=lm, formula = y ~ splines::bs(x, 3), aes(log(Exp+1),est.prob),col="red", alpha=0.3, size=0.3) +
    ggplot2::scale_y_log10(limits=c(1e-6,1)) + ggplot2::ylab("Probability")
  return(p)
}
