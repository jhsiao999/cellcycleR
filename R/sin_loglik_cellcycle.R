#' @title Loglikelihood of sinusoidal model
#'
#' @param cycle_data a N x G matrix, where N is number of cells, G number of genes
#' @param cell_times a N x 1 vector of cell times
#' @param amp the amplitude vector (G x 1) over the genes
#' @param phi the G x 1 vector of phase values over genes
#' @param sigma: the G x 1 vector of gene variation
#' @param freq: frequency of the signal of the gene. Defaults to 1.
#'
#' @description Computes the loglikelihood of all the cells in the cycle under sinusoidal
#'              gene expression patterns.
#'
#'  @author  Kushal K Dey
#'
#'  @examples
#'  G <- 500;
#'  num_cells <- 400;
#'  amp_genes <- rep(10, G);
#'  phi_genes <- runif(G, 0, 2*pi)
#'  sigma_genes <- rchisq(G, 4);
#'  cell_times_sim <- sample(seq(0,2*pi, 2*pi/(num_cells-1)), num_cells, replace=FALSE);
#'  cycle_data <- sim_sinusoidal_cycle(G, amp_genes, phi_genes, sigma_genes, cell_times_sim);
#'
#'  sing_loglik_cellcycle(cycle_data, cell_times, amp_genes, phi_genes, sigma_genes)
#'
#'  @export
sin_loglik_cellcycle <- function(cycle_data, cell_times, amp, phi, sigma, freq=1)
{
  G <- dim(cycle_data)[2]
  numcells <- dim(cycle_data)[1]
  sum <- 0

  for(s in 1:numcells)
  {
    sum <- sum + sum(mapply(dnorm, cycle_data[s,],
                            amp * sin(cell_times[s] + phi), sigma, log=TRUE));
  }

  return(sum)
}
