#' @title Fitting sinusoidal model
#'
#' @description Compute sinusoidal model estimates: amplitude, phase, and estimated cell times.
#'
#' @param cycle_data a N x G matrix, where N is number of cells, G number of genes
#' @param celltime_levels single numeric value. This determines the number of time intervals
#'            over (0, 2/pi). In the algorithm, the cells are initially assigned to an
#'            arbitrarily time interval. Then, based on the amplitude and phase estimates,
#'            the time intervals are updated accordingly.
#' @param cell_time_iter N x 1 vector of cell times that are set manually (need not be ordered).
#' @param freq Frequency of sinusoidal curve. Default=1, i.e., (0, 2\pi).
#' @param fix.phase TRUE if cell_times is an ordered vector, i.e., the
#'                  relative position of cells in time is known, and FALSE otherwise.
#' @param phase_in G x 1 vector of phase values.
#' @param n_core Number of cores used in the computations. Default is 1: not parallelized.
#'          n_cores > 1 turns to parallel computing.
#'
#' @author Kushal K Dey, Chiaowen Joyce Hsiao
#'
#'
sin_cell_ordering_iter_2 <- function(cycle_data,
                                   celltime_levels,
                                   cell_times_iter, freq=1,
                                   fix.phase=FALSE, phase_in=NULL,
                                   n_cores=1)
{

  G <- dim(cycle_data)[2]
  numcells <- dim(cycle_data)[1]
  sigma <- array(0,G)
  amp <- array(0,G)
  phi <- array(0,G)

  # Fit linear models for each gene $g$ given the cell times [ linear model depends on fix.phase]

  if(!fix.phase){

    fit_lm_notfixed_phase <- function(g) {
      fit <- lm(cycle_data[,g]  ~ sin(freq*cell_times_iter) + cos(freq*cell_times_iter) -1)
      sigma <- sd(fit$residuals)
      beta1 <- fit$coefficients[1]
      beta2 <- fit$coefficients[2]
      if(beta1==0 & beta2==0){
        stop(paste0("You have a gene with all 0 counts at gene",g))
      }
      # out_amp <- sqrt(beta1^2 + beta2^2);
      # out_phi <- atan3(as.numeric(beta2), as.numeric(beta1));
      # ll <- list("out_amp"=out_amp, "out_phi"=out_phi, "out_sigma"=out_sigma)
      ll <- list(beta1=beta1, beta2=beta2, sigma=sigma)
      return(ll)
    }

    fit_lm_all <- parallel::mclapply(1:G, function(g) fit_lm_notfixed_phase(g),
                                     mc.cores=n_cores)
    beta1 <- as.numeric(unlist(lapply(1:G, function(g) return(fit_lm_all[[g]]$beta1))))
    beta2 <- as.numeric(unlist(lapply(1:G, function(g) return(fit_lm_all[[g]]$beta2))))
    sigma <- as.numeric(unlist(lapply(1:G, function(g) return(fit_lm_all[[g]]$sigma))))

    # amp <- as.numeric(unlist(lapply(1:G, function(g) return(fit_lm_all[[g]]$out_amp))))
    # phi <- as.numeric(unlist(lapply(1:G, function(g) return(fit_lm_all[[g]]$out_phi))))
    # sigma <- as.numeric(unlist(lapply(1:G, function(g) return(fit_lm_all[[g]]$out_sigma))))
  }

  if(fix.phase){
    phi <- phase_in

    fit_lm_fixed_phase <- function(g) {
      fit <- lm(cycle_data[,g]  ~ sin(freq*cell_times_iter+phi[g]) -1);
      out_sigma <- sd(fit$residuals);
      out_amp <- abs(fit$coefficients[1]);
      out_phi <- phi;
      ll <- list("out_amp"=out_amp, "out_phi"=out_phi, "out_sigma"=out_sigma)
      return(ll)
    }

    fit_lm_all <- parallel::mclapply(1:G, function(g) fit_lm_fixed_phase(g),
                                     mc.cores=n_cores)

    amp <- as.numeric(unlist(lapply(1:G, function(g) return(fit_lm_all[[g]]$out_amp))))
    phi <- as.numeric(unlist(lapply(1:G, function(g) return(fit_lm_all[[g]]$out_phi))))
    sigma <- as.numeric(unlist(lapply(1:G, function(g) return(fit_lm_all[[g]]$out_sigma))))
  }

  # compute predictive values
  cell_times_class <- seq(0, 2*pi, 2*pi/(celltime_levels-1))
  num_celltime_class <- length(cell_times_class)

  sin_class_times <- sin(freq*cell_times_class)
  cos_class_times <- cos(freq*cell_times_class)
  # sin_phi_genes <- sin(phi)
  # cos_phi_genes <- cos(phi)
  # sinu_signal <- cbind(sin_class_times, cos_class_times) %*% rbind(amp*cos_phi_genes, amp*sin_phi_genes)
  sinu_signal <- cbind(sin_class_times, cos_class_times) %*% rbind(beta1, beta2)

  # compute log liklihood from the model for each cell at each time interval
  options(digits=12)
  signal_intensity_per_class <- matrix(0, numcells, num_celltime_class)

  signal_intensity_per_class <- do.call(rbind, lapply(1:numcells, function(cell)
  {
    res_error <- sweep(sinu_signal,2,cycle_data[cell,]);
    res_error_adjusted <- -(res_error^2);
    res_error_adjusted <- sweep(res_error_adjusted, 2, 2*sigma^2, '/');
    out <- rowSums(sweep(res_error_adjusted,2,log(sigma)) - 0.5*log(2*pi));
    return(out)
  }))

  # estimate expected time interval based on log likelihood

  # for every cell, normalize the log likelihood with respect to
  # the maximum log likelihood across time intervals
  signal_intensity_class_exp <- do.call(rbind,
       lapply(1:dim(signal_intensity_per_class)[1], function(x)
  {
    out <- exp(signal_intensity_per_class[x,]- max(signal_intensity_per_class[x,]))
    return(out)
  }))

  cell_times <- cell_times_class[unlist(lapply(1:dim(signal_intensity_class_exp)[1], function(x)
  {
    temp <- signal_intensity_class_exp[x,]
    if(length(unique(signal_intensity_class_exp[x,]))==1)
      out <- sample(1:dim(signal_intensity_class_exp)[2],1)
    else
      out <- which(rmultinom(1,1,signal_intensity_class_exp[x,])==1)
    return(out)
  }))]

  loglik <- sin_loglik_cellcycle_2(cycle_data=cycle_data, cell_times=cell_times,
                                 beta1=beta1, beta2=beta2, sigma=sigma, freq=1)

  out <- list(cell_times_iter=cell_times,
              beta1_iter=beta1, beta2_iter=beta2, sigma_iter=sigma,
              signal_intensity_iter=signal_intensity_per_class,
              loglik_iter=loglik)
  return(out)
}
