#! /usr/bin/env Rscript

# get environment variables
STEPSIZE <- as.numeric(Sys.getenv('STEPSIZE'))
TASKID <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

# set defaults if nothing comes from environment variables
STEPSIZE[is.na(STEPSIZE)] <- 1
TASKID[is.na(TASKID)] <- 0

# get command lines arguments
args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 1){
  stop("Not enough arguments. Please use args 'listsize', 'prepare', 'run <itemsize>' or 'merge'")
}

##############################################
# !!!! set a directory for output to save !!!! 
##############################################
save_dir <- "~/"

# parameters
ns <- c(100, 200, 500, 1000)
seed <- 1:1000
trt_effect <- c(0, 0.5, 1)
dgp <- c(1, 3)
parm <- expand.grid(n = ns, seed = seed, 
                    trt_effect = trt_effect,
                    dgp = dgp)

# load packages
library(drord)

# function to make data for hospitalized COVID-19 patients
make_data_1 <- function(n, tx_effect = 0){
  unif_var <- runif(n)
  age_grp <- as.numeric(cut(unif_var, c(0, 0.003848326, 0.192681847, 0.354730472, 0.520380178, 0.745410702, 0.888355056, 1), include.lowest = TRUE))

  # treatment
  treat <- rbinom(n, 1, 1/2)

  # probability of death
  prob_death_by_age <- c(0.000000000, 0.008547009, 0.026262626, 0.079051383, 
                         0.105409154, 0.165919283, 0.371062992)
  prob_death <- prob_death_by_age[age_grp]

  # probability of icu
  prob_icu_given_age <- function(age_vec, tx, tx_effect = 0){
    prob_vent_by_age <- c(0.0000000, 0.1766382, 0.3191919, 0.3142292, 0.3730929, 0.4652466, 0.3474409)
    intercept_by_age <- qlogis(prob_vent_by_age)
    n <- length(age_vec)
    prob_vent <- plogis(intercept_by_age[age_vec] + tx_effect * tx)
    return(prob_vent)
  }
  prob_vent <- prob_icu_given_age(age_vec = age_grp, tx = treat, 
                                   tx_effect = tx_effect)

  prob_nothing <- 1 - prob_death - prob_vent

  prob_matrix <- cbind(prob_death, prob_vent, prob_nothing)
  outcome <- apply(prob_matrix, 1, function(probs){
    sample(1:3, size = 1, replace = TRUE, prob = probs)
  })
  
  return(list(out = outcome, covar = data.frame(age_grp = age_grp), treat = treat))
}

get_truth_1 <- function(n = 1e6, tx_effect){
  unif_var <- runif(n)
  age_grp <- as.numeric(cut(unif_var, c(0, 0.003848326, 0.192681847, 0.354730472, 0.520380178, 0.745410702, 0.888355056, 1), include.lowest = TRUE))

  # treatment
  treat <- rbinom(n, 1, 1/2)

  # probability of death
  prob_death_by_age <- c(0.000000000, 0.008547009, 0.026262626, 0.079051383, 
                         0.105409154, 0.165919283, 0.371062992)
  prob_death <- prob_death_by_age[age_grp]

  # probability of ventilator
  # age_vec a vector of age group 1:7
  prob_icu_given_age <- function(age_vec, tx, tx_effect = 0){
    prob_vent_by_age <- c(0.0000000, 0.1766382, 0.3191919, 0.3142292, 0.3730929, 0.4652466, 0.3474409)
    intercept_by_age <- qlogis(prob_vent_by_age)
    n <- length(age_vec)
    prob_vent <- plogis(intercept_by_age[age_vec] + tx_effect * tx)
    return(prob_vent)
  }
  prob_vent1 <- prob_icu_given_age(age_vec = age_grp, tx = rep(1, n), 
                                   tx_effect = tx_effect)
  prob_vent0 <- prob_icu_given_age(age_vec = age_grp, tx = rep(0, n), 
                                   tx_effect = tx_effect)

  prob_nothing1 <- 1 - prob_death - prob_vent1
  prob_nothing0 <- 1 - prob_death - prob_vent0

  # mann whitney
  F_0_kminus1 <- c(0, F0[-3])
  mannwhitney_truth <- sum((F_0_kminus1 + 1/2 * f_0) * f_1)

  # log-odds ratio
  logodds_truth <- mean(qlogis(F1[2:3]) - qlogis(F0[2:3]))

  # weighted mean
  mean_y1 <- sum((1:3)*f_1)
  mean_y0 <- sum((1:3)*f_0)
  wmean_truth <- mean_y1 - mean_y0
  
  return(list(mann_whitney = mannwhitney_truth,
              log_odds = logodds_truth,
              weighted_mean = wmean_truth))
}


# function to make data to match hospitalizedControlArmDistributionsByAgeGroup
make_data_3 <- function(n, tx_effect = 0){
  # cut points of std normal to get cov dist.
  norm_cuts <- c(-Inf, -1.64267998, -0.41766186 , 0.03327091, 0.49144538, 1.06010333, 1.56493129, Inf)
  norm_var <- rnorm(n)
  age_grp <- as.numeric(cut(norm_var, norm_cuts, include.lowest = TRUE))

  # treatment
  treat <- rbinom(n, 1, 1/2)

  # probability of death
  prob_death_by_age <- c(0.0000 ,0.0015 ,0.0065 ,0.0200 ,0.0380 ,0.0740 ,0.1885)
  prob_death <- prob_death_by_age[age_grp]

  # probability of ventilator
  # age_vec a vector of age group 1:7
  prob_icu_given_age <- function(age_vec, tx, tx_effect = 0){
    prob_vent_by_age <- c( 0.0205, 0.1755, 0.2475, 0.2530, 0.3605, 0.4460, 0.5080)
    intercept_by_age <- qlogis(prob_vent_by_age)
    n <- length(age_vec)
    prob_vent <- plogis(intercept_by_age[age_vec] + tx_effect * tx)
    return(prob_vent)
  }
  prob_vent <- prob_icu_given_age(age_vec = age_grp, tx = treat, 
                                   tx_effect = tx_effect)

  prob_nothing <- 1 - prob_death - prob_vent

  prob_matrix <- cbind(prob_death, prob_vent, prob_nothing)
  outcome <- apply(prob_matrix, 1, function(probs){
    sample(1:3, size = 1, replace = TRUE, prob = probs)
  })
  
  return(list(out = outcome, covar = data.frame(age_grp = age_grp), treat = treat))
}

get_truth_3 <- function(n = 1e6, tx_effect = 0){
    # cut points of std normal to get cov dist.
  norm_cuts <- c(-Inf, -1.64267998, -0.41766186 , 0.03327091, 0.49144538, 1.06010333, 1.56493129, Inf)
  norm_var <- rnorm(n)
  age_grp <- as.numeric(cut(norm_var, norm_cuts, include.lowest = TRUE))

  # treatment
  treat <- rbinom(n, 1, 1/2)

  # probability of death
  prob_death_by_age <- c(0.0000 ,0.0015 ,0.0065 ,0.0200 ,0.0380 ,0.0740 ,0.1885)
  prob_death <- prob_death_by_age[age_grp]

  # probability of ventilator
  # age_vec a vector of age group 1:7
  prob_icu_given_age <- function(age_vec, tx, tx_effect = 0){
    prob_vent_by_age <- c( 0.0205, 0.1755, 0.2475, 0.2530, 0.3605, 0.4460, 0.5080)
    intercept_by_age <- qlogis(prob_vent_by_age)
    n <- length(age_vec)
    prob_vent <- plogis(intercept_by_age[age_vec] + tx_effect * tx)
    return(prob_vent)
  }
  prob_vent1 <- prob_icu_given_age(age_vec = age_grp, tx = 1, 
                                   tx_effect = tx_effect)
  prob_vent0 <- prob_icu_given_age(age_vec = age_grp, tx = 0, 
                                   tx_effect = tx_effect)

  prob_nothing1 <- 1 - prob_death - prob_vent1
  prob_nothing0 <- 1 - prob_death - prob_vent0
  
  F0 <- c(mean(prob_death), mean(prob_death + prob_vent0), 1)
  F1 <- c(mean(prob_death), mean(prob_death + prob_vent1), 1)
  f_0 <- c(mean(prob_death), mean(prob_vent0), mean(prob_nothing0))
  f_1 <- c(mean(prob_death), mean(prob_vent1), mean(prob_nothing1))

  # mann whitney
  F_0_kminus1 <- c(0, F0[-3])

  # estimate
  mannwhitney_truth <- sum((F_0_kminus1 + 1/2 * f_0) * f_1)

  # log-odds ratio
  logodds_truth <- mean(qlogis(F1[2:3]) - qlogis(F0[2:3]))

  # weighted mean
  mean_y1 <- sum((1:3)*f_1)
  mean_y0 <- sum((1:3)*f_0)
  wmean_truth <- mean_y1 - mean_y0
  
  return(list(mann_whitney = mannwhitney_truth,
              log_odds = logodds_truth,
              weighted_mean = wmean_truth))
}


# get the list size #########
if (args[1] == 'listsize') {
  cat(nrow(parm))
}

# execute prepare job ##################
if (args[1] == 'prepare') {
  print('Nothing to prepare')
}

# execute parallel job #################################################
if (args[1] == 'run') {
  if (length(args) < 2) {
    stop("Not enough arguments. 'run' needs a second argument 'id'")
  }
  id <- as.numeric(args[2])
  print(paste(Sys.time(), "arrid:" , id, "TASKID:",
              TASKID, "STEPSIZE:", STEPSIZE))
  for (i in (id+TASKID):(id+TASKID+STEPSIZE-1)) {
    print(paste(Sys.time(), "i:" , i))
    print(parm[i,])
    
    # set seed
    set.seed(parm$seed[i])

    if(parm$dgp[i] == 1){
      get_tx_eff <- function(n, trt_effect = 1){
        if(trt_effect == 1){
          if(n == 100){
            return(-1.73)
          }else if(n == 200){            
            return(-1.84)
          }else if(n == 500){            
            return(-1.0)
          }else{            
            return(-0.68)
          }
        }else if(trt_effect == 0.5){
          if(n == 100){            
            return(-1.15)
          }else if(n == 200){            
            return(-1.20)
          }else if(n == 500){            
            return(-0.68)
          }else{            
            return(-0.46)
          }
        }
      }

      if(parm$trt_effect[i] == 0){
        tx_eff <- 0
      }else{
        tx_eff <- get_tx_eff(parm$n[i])
      }
      dat <- make_data_1(n = parm$n[i], tx_effect = tx_eff)

      # get truth
      set.seed(1234)
      truth <- get_truth_1(n = 1e6, tx_effect = tx_eff)
    }else{
      get_tx_eff <- function(n, trt_effect){
        if(trt_effect == 1){
          if(n == 100){
            return(-1.58)
          }else if(n == 200){
            return(-1.58)
          }else if(n == 500){
            return(-0.83)
          }else{
            return(-0.57)
          }
        }else if(trt_effect == 0.5){
          if(n == 100){
            return(-1.17)
          }else if(n == 200){
            return(-0.82)
          }else if(n == 500){
            return(-0.47)
          }else{
            return(-0.36)
          }
        }
      }

      if(parm$trt_effect[i] == 0){
        tx_eff <- 0
      }else{
        tx_eff <- get_tx_eff(parm$n[i])
      }
      dat <- make_data_3(n = parm$n[i], tx_effect = tx_eff)

      # get truth
      set.seed(1234)
      truth <- get_truth_3(n = 1e6, tx_effect = tx_eff)
    }

    # get output for adjusted estimators
    rslt_adj <- drord(
      out = dat$out,
      treat = dat$treat,
      covar = dat$covar,
      out_form = "age_grp",
      out_model = "clm",
      param = c("weighted_mean", "mann_whitney", "log_odds"),
      ci = c("bca", "wald"),
      nboot = 1e3,
      est_dist = FALSE
    )

    # get output for unadjusted estimators
    rslt_unadj <- drord(
      out = dat$out,
      treat = dat$treat,
      covar = dat$covar,
      out_form = "1",
      treat_form = "1",
      param = c("weighted_mean", "mann_whitney", "log_odds"),
      ci = c("bca", "wald"),
      nboot = 1e3, 
      stratify = TRUE,
      est_dist = FALSE
    )

    in_ci <- function(ci, truth){
      truth < max(ci) & truth > min(ci)
    }

    format_out <- function(rslt){
      # format output
      logodds_out <- c(rslt$log_odds$est[3],
                       rslt$log_odds$ci$wald[3,],
                       in_ci(ci = rslt$log_odds$ci$wald[3,], truth = truth$log_odds),
                       rslt$log_odds$ci$bca[3,],
                       in_ci(ci = rslt$log_odds$ci$bca[3,], truth = truth$log_odds))

      mannwhitney_out <- c(rslt$mann_whitney$est[1],
                       rslt$mann_whitney$ci$wald,
                      in_ci(ci = rslt$mann_whitney$ci$wald, truth = truth$mann_whitney),
                       rslt$mann_whitney$ci$bca[1,],
                       in_ci(ci = rslt$mann_whitney$ci$bca[1,], truth = truth$mann_whitney))
      
      wmean_out <- c(rslt$weighted_mean$est$est[1] - rslt$weighted_mean$est$est[2],
                     rslt$weighted_mean$ci$wald[3,],
                     in_ci(ci = rslt$weighted_mean$ci$wald[3,], truth = truth$weighted_mean),
                     rslt$weighted_mean$ci$bca[3,],
                     in_ci(ci = rslt$weighted_mean$ci$bca[3,], truth = truth$weighted_mean))
      return(c(logodds_out, mannwhitney_out, wmean_out))
    }
   

    out <- c(as.numeric(parm[i,]), unlist(truth), 
             format_out(rslt = rslt_adj),
             format_out(rslt = rslt_unadj))

    # save output
    save(out, file = paste0(save_dir, "fit_", i, ".RData"))
  }
}

# merge job ###########################
if (args[1] == 'merge'){
  out_names <- c(
    "n", "seed", "trt_effect", "dgp",
    "mannwhitney_truth", "logodds_truth","wmean_truth",
    c(t(outer(c("logodds_adj","mannwhitney_adj","wmean_adj"), c("_est", "_wald_cil", "_wald_ciu",
                                                 "_wald_cover","_bca_cil","_bca_ciu",
                                                 "_bca_cover"), paste0))),
    c(t(outer(c("logodds_unadj","mannwhitney_unadj","wmean_unadj"), c("_est", "_wald_cil", "_wald_ciu",
                                                 "_wald_cover","_bca_cil","_bca_ciu",
                                                 "_bca_cover"), paste0)))
  )

  n_out <- length(out_names)

  rslt <- matrix(NA, nrow = nrow(parm), ncol = n_out)
  for(i in 1:nrow(parm)){
    tmp <- tryCatch({
      load(paste0(save_dir, "fit_", i, ".RData"))
      out
    }, error = function(e){
      rep(NA, n_out)
    })
    rslt[i, ] <- tmp
  }

  out <- data.frame(rslt)
  colnames(out) <- out_names
 
  make_output_table <- function(out, estimand, ci, scale = TRUE){
    null_val <- if(estimand == "mannwhitney"){
      0.5
    }else if(estimand == "logodds"){
      0
    }else if(estimand == "wmean"){
      0
    }
    out$logodds_unadj_est[out$logodds_unadj_est == Inf] <- NA
    out$logodds_unadj_est[out$logodds_unadj_est == NaN] <- NA
    out$logodds_unadj_est[out$logodds_unadj_est == -Inf] <- NA
    out$logodds_adj_est[out$logodds_adj_est == -Inf] <- NA
    out$logodds_adj_est[out$logodds_adj_est == Inf] <- NA
    out$logodds_adj_est[out$logodds_adj_est == NaN] <- NA

    sample_size <- sort(rep(c(100, 200, 500, 1000), 6))
    est_type <- rep(c("Unadjusted", "Adjusted"), 12)

    effect_col <- unlist(by(out, out$n, function(x){
      c(rep(x[x$trt_effect == 0, paste0(estimand, "_truth")][1],2),
        rep(x[x$trt_effect == 0.5, paste0(estimand, "_truth")][1],2),
        rep(x[x$trt_effect == 1, paste0(estimand, "_truth")][1],2))
    }), use.names = FALSE)

    power_col <- unlist(by(out, out$n, function(x){
      c(mean(na.rm = TRUE, as.numeric(x[x$trt_effect == 0, paste0(estimand, "_unadj_", ci, "_cil")] > null_val |
                                      x[x$trt_effect == 0, paste0(estimand, "_unadj_", ci, "_ciu")] < null_val)),
        mean(na.rm = TRUE, as.numeric(x[x$trt_effect == 0, paste0(estimand, "_adj_", ci, "_cil")] > null_val |
                                      x[x$trt_effect == 0, paste0(estimand, "_adj_", ci, "_ciu")] < null_val)),
        mean(na.rm = TRUE, as.numeric(x[x$trt_effect == 0.5, paste0(estimand, "_unadj_", ci, "_cil")] > null_val |
                                      x[x$trt_effect == 0.5, paste0(estimand, "_unadj_", ci, "_ciu")] < null_val)),
        mean(na.rm = TRUE, as.numeric(x[x$trt_effect == 0.5, paste0(estimand, "_adj_", ci, "_cil")] > null_val |
                                      x[x$trt_effect == 0.5, paste0(estimand, "_adj_", ci, "_ciu")] < null_val)),
        mean(na.rm = TRUE, as.numeric(x[x$trt_effect == 1, paste0(estimand, "_unadj_", ci, "_cil")] > null_val |
                                      x[x$trt_effect == 1, paste0(estimand, "_unadj_", ci, "_ciu")] < null_val)),
        mean(na.rm = TRUE, as.numeric(x[x$trt_effect == 1, paste0(estimand, "_adj_", ci, "_cil")] > null_val |
                                      x[x$trt_effect == 1, paste0(estimand, "_adj_", ci, "_ciu")] < null_val)))    
      }), use.names = FALSE)

    mse_col <- unlist(by(out, out$n, function(x){
      if(scale){
        scale_factor <- x$n[1]
      }else{
        scale_factor <- 1
      }
      c(scale_factor*mean(na.rm = TRUE, (x[x$trt_effect == 0, paste0(estimand, "_unadj_est")] - x[x$trt_effect == 0, paste0(estimand, "_truth")])^2),
        scale_factor*mean(na.rm = TRUE, (x[x$trt_effect == 0, paste0(estimand, "_adj_est")] - x[x$trt_effect == 0, paste0(estimand, "_truth")])^2),
        scale_factor*mean(na.rm = TRUE, (x[x$trt_effect == 0.5, paste0(estimand, "_unadj_est")] - x[x$trt_effect == 0.5, paste0(estimand, "_truth")])^2),
        scale_factor*mean(na.rm = TRUE, (x[x$trt_effect == 0.5, paste0(estimand, "_adj_est")] - x[x$trt_effect == 0.5, paste0(estimand, "_truth")])^2),
        scale_factor*mean(na.rm = TRUE, (x[x$trt_effect == 1, paste0(estimand, "_unadj_est")] - x[x$trt_effect == 1, paste0(estimand, "_truth")])^2),
        scale_factor*mean(na.rm = TRUE, (x[x$trt_effect == 1, paste0(estimand, "_adj_est")] - x[x$trt_effect == 1, paste0(estimand, "_truth")])^2))
    }), use.names = FALSE)

    bias_col <- unlist(by(out, out$n, function(x){
      if(scale){
        scale_factor <- sqrt(x$n[1])
      }else{
        scale_factor <- 1
      }
      c(scale_factor*mean(na.rm = TRUE, (x[x$trt_effect == 0, paste0(estimand, "_unadj_est")] - x[x$trt_effect == 0, paste0(estimand, "_truth")])),
        scale_factor*mean(na.rm = TRUE, (x[x$trt_effect == 0, paste0(estimand, "_adj_est")] - x[x$trt_effect == 0, paste0(estimand, "_truth")])),
        scale_factor*mean(na.rm = TRUE, (x[x$trt_effect == 0.5, paste0(estimand, "_unadj_est")] - x[x$trt_effect == 0.5, paste0(estimand, "_truth")])),
        scale_factor*mean(na.rm = TRUE, (x[x$trt_effect == 0.5, paste0(estimand, "_adj_est")] - x[x$trt_effect == 0.5, paste0(estimand, "_truth")])),
        scale_factor*mean(na.rm = TRUE, (x[x$trt_effect == 1, paste0(estimand, "_unadj_est")] - x[x$trt_effect == 1, paste0(estimand, "_truth")])),
        scale_factor*mean(na.rm = TRUE, (x[x$trt_effect == 1, paste0(estimand, "_adj_est")] - x[x$trt_effect == 1, paste0(estimand, "_truth")])))
    }), use.names = FALSE)
    
    var_col <- unlist(by(out, out$n, function(x){
      if(scale){
        scale_factor <- x$n[1]
      }else{
        scale_factor <- 1
      }
      c(scale_factor * var(na.rm = TRUE, x[x$trt_effect == 0, paste0(estimand, "_unadj_est")]),
        scale_factor * var(na.rm = TRUE, x[x$trt_effect == 0, paste0(estimand, "_adj_est")]),
        scale_factor * var(na.rm = TRUE, x[x$trt_effect == 0.5, paste0(estimand, "_unadj_est")]),
        scale_factor * var(na.rm = TRUE, x[x$trt_effect == 0.5, paste0(estimand, "_adj_est")]),
        scale_factor * var(na.rm = TRUE, x[x$trt_effect == 1, paste0(estimand, "_unadj_est")]),
        scale_factor * var(na.rm = TRUE, x[x$trt_effect == 1, paste0(estimand, "_adj_est")]))
    }), use.names = FALSE)
    rel_eff <- rep(NA, length(mse_col))
    rel_eff[seq(1, length(mse_col), by = 2)] <- 1
    rel_eff[seq(2, length(mse_col), by = 2)] <- mse_col[seq(2, length(mse_col), 2)] / mse_col[seq(1, length(mse_col), 2)]

    tab <- data.frame(sample_size, est_type, effect_col, power_col, mse_col, bias_col, var_col, rel_eff)

    return(tab)
  }

  library(xtable)
  # hospitalized
  tab_parm <- expand.grid(estimand = c("mannwhitney","logodds","wmean"),
                          ci = c("wald", "bca"))
  tab_list_1 <- vector(mode = "list", length = nrow(tab_parm))
  for(i in seq_len(nrow(tab_parm))){
    tab_list_1[[i]] <- make_output_table(out[out$dgp == 1,], estimand = tab_parm$estimand[i], ci = tab_parm$ci[i])
  }

  # DIM BCA
  print(xtable(tab_list_1[[6]], digits = c(0,0,rep(3,7))), 
        include.rownames = FALSE,
        hline.after = c(4, 8, 12))
  # MW BCA
  print(xtable(tab_list_1[[4]], digits = c(0,0,rep(3,7))), 
        include.rownames = FALSE,
        hline.after = c(4, 8, 12))
  # LOR BCA
  print(xtable(tab_list_1[[5]], digits = c(0,0,rep(3,7))), 
        include.rownames = FALSE,
        hline.after = c(4, 8, 12))

  # supplement Wald-based inference
  print(xtable(tab_list_1[[3]], digits = c(0,0,rep(3,7))), 
        include.rownames = FALSE,
        hline.after = c(4, 8, 12))
  print(xtable(tab_list_1[[1]], digits = c(0,0,rep(3,7))), 
        include.rownames = FALSE,
        hline.after = c(4, 8, 12))
  print(xtable(tab_list_1[[2]], digits = c(0,0,rep(3,7))), 
        include.rownames = FALSE,
        hline.after = c(4, 8, 12))
  
  # non-hospitalized results
  tab_list_3 <- vector(mode = "list", length = nrow(tab_parm))
  for(i in seq_len(nrow(tab_parm))){
    tab_list_3[[i]] <- make_output_table(out[out$dgp == 3,], estimand = tab_parm$estimand[i], ci = tab_parm$ci[i])
  }


  # DIM
  print(xtable(tab_list_3[[6]], digits = c(0,0,rep(3,7))), 
        include.rownames = FALSE,
        hline.after = c(4, 8, 12))
  print(xtable(tab_list_3[[3]], digits = c(0,0,rep(3,7))), 
        include.rownames = FALSE,
        hline.after = c(4, 8, 12))

  # MW
  print(xtable(tab_list_3[[5]], digits = c(0,0,rep(3,7))), 
        include.rownames = FALSE,
        hline.after = c(4, 8, 12))
  print(xtable(tab_list_3[[2]], digits = c(0,0,rep(3,7))), 
        include.rownames = FALSE,
        hline.after = c(4, 8, 12))
  # LOR
  print(xtable(tab_list_3[[4]], digits = c(0,0,rep(3,7))), 
        include.rownames = FALSE,
        hline.after = c(4, 8, 12))
  print(xtable(tab_list_3[[1]], digits = c(0,0,rep(3,7))), 
        include.rownames = FALSE,
        hline.after = c(4, 8, 12))
}