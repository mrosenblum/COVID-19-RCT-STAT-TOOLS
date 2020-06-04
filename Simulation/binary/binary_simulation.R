#! /usr/bin/env Rscript

# get environment variables
MYSCRATCH <- Sys.getenv('MYSCRATCH')
RESULTDIR <- Sys.getenv('RESULTDIR')
STEPSIZE <- as.numeric(Sys.getenv('STEPSIZE'))
TASKID <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

# set defaults if nothing comes from environment variables
MYSCRATCH[is.na(MYSCRATCH)] <- '.'
RESULTDIR[is.na(RESULTDIR)] <- '.'
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
trt_effect <- c(0, 1)
dgp <- c(1, 3)
parm <- expand.grid(n = ns, seed = seed, 
                    trt_effect = trt_effect,
                    dgp = dgp)

# load packages
library(drord)

# function to make data to match hospitalizedControlArmDistributionsByAgeGroup
make_data_1 <- function(n, tx_effect = 0){
  # cut points of std normal to get cov dist.
  norm_cuts <- c(-Inf, -2.31, -1.30, -0.77, -0.39, 0.06, 0.66, Inf)
  norm_var <- rnorm(n)
  age_grp <- as.numeric(cut(norm_var, norm_cuts, include.lowest = TRUE))

  # treatment
  treat <- rbinom(n, 1, 1/2)

  # probability of death
  prob_death_by_age <- c(0.000000000, 0.008547009, 0.026262626, 0.079051383, 
                         0.105409154, 0.165919283, 0.371062992)
  prob_death <- prob_death_by_age[age_grp]

  # probability of ventilator
  # age_vec a vector of age group 1:7
  prob_vent_given_age <- function(age_vec, tx, tx_effect = 0){
    prob_vent_by_age <- c(0.0000000, 0.1766382, 0.3191919, 0.3142292, 0.3730929, 0.4652466, 0.3474409)
    intercept_by_age <- qlogis(prob_vent_by_age)
    n <- length(age_vec)
    prob_vent <- plogis(intercept_by_age[age_vec] + tx_effect * tx)
    return(prob_vent)
  }
  prob_vent <- prob_vent_given_age(age_vec = age_grp, tx = treat, 
                                   tx_effect = tx_effect)

  prob_nothing <- 1 - prob_death - prob_vent

  prob_matrix <- cbind(prob_death, prob_vent, prob_nothing)
  outcome <- apply(prob_matrix, 1, function(probs){
    sample(1:3, size = 1, replace = TRUE, prob = probs)
  })
  
  return(list(out = outcome, covar = data.frame(age_grp = age_grp), treat = treat))
}

get_truth_1 <- function(n = 1e6, tx_effect){
  # cut points of std normal to get cov dist.
  norm_cuts <- c(-Inf, -2.31, -1.30, -0.77, -0.39, 0.06, 0.66, Inf)
  norm_var <- rnorm(n)
  age_grp <- as.numeric(cut(norm_var, norm_cuts, include.lowest = TRUE))

  # treatment
  treat <- rbinom(n, 1, 1/2)

  # probability of death
  prob_death_by_age <- c(0.000000000, 0.008547009, 0.026262626, 0.079051383, 
                         0.105409154, 0.165919283, 0.371062992)
  prob_death <- prob_death_by_age[age_grp]

  # probability of ventilator
  # age_vec a vector of age group 1:7
  prob_vent_given_age <- function(age_vec, tx, tx_effect = 0){
    prob_vent_by_age <- c(0.0000000, 0.1766382, 0.3191919, 0.3142292, 0.3730929, 0.4652466, 0.3474409)
    intercept_by_age <- qlogis(prob_vent_by_age)
    n <- length(age_vec)
    prob_vent <- plogis(intercept_by_age[age_vec] + tx_effect * tx)
    return(prob_vent)
  }
  prob_vent1 <- prob_vent_given_age(age_vec = age_grp, tx = rep(1, n), 
                                   tx_effect = tx_effect)
  prob_vent0 <- prob_vent_given_age(age_vec = age_grp, tx = rep(0, n), 
                                   tx_effect = tx_effect)

  prob_nothing1 <- 1 - prob_death - prob_vent1
  prob_nothing0 <- 1 - prob_death - prob_vent0

  # weighted mean -- modified for binarization
  mean_y1 <- sum(c(1,1,0)*f_1)
  mean_y0 <- sum(c(1,1,0)*f_0)
  wmean_truth <- mean_y1 - mean_y0
  
  return(list(mann_whitney = mannwhitney_truth))
}


# get_power_ttest1 <- function(n, tx_effect, B){
#   do_one_ttest <- function(n, tx_effect){
#     dat <- make_data_1(n = n, tx_effect = tx_effect)
#     ttest <- t.test(x = as.numeric(dat$out[dat$treat == 1] %in% c(1,2)), 
#                     y = as.numeric(dat$out[dat$treat == 0] %in% c(1,2)))
#     return(ttest$p.value <= 0.05)
#   }
#   many_ttests <- replicate(B, do_one_ttest(n = n, tx_effect = tx_effect))
#   pow <- mean(many_ttests)
#   return(pow)
# }

# # calibrate
# effect_seq <- seq(-1.5, -2, length = 25)
# ttest_rslt <- sapply(c(100), function(n){
#   sapply(effect_seq, function(tx_eff, n){
#     get_power_ttest1(n = n, tx_effect = tx_eff, B = 500)
#   }, n = n)
# })



# function to make data to match hospitalizedControlArmDistributionsByAgeGroup
make_data_3 <- function(n, tx_effect = 0){
  # cut points of std normal to get cov dist.
  norm_cuts <- c(-Inf, -1.64267998, -0.41766186 , 0.03327091, 0.49144538, 1.06010333, 1.56493129, Inf)
  norm_var <- rnorm(n)
  age_grp <- as.numeric(cut(norm_var, norm_cuts, include.lowest = TRUE))

  # treatment
  treat <- rbinom(n, 1, 1/2)

  # probability of death
  # hospitalizedControlArmDistributionsByAgeGroup$ProbDeath
  prob_death_by_age <- c(0.0000 ,0.0015 ,0.0065 ,0.0200 ,0.0380 ,0.0740 ,0.1885)
  prob_death <- prob_death_by_age[age_grp]

  # probability of ventilator
  # prob_vent_by_age <- hospitalizedControlArmDistributionsByAgeGroup$ProbVentilator
  # age_vec a vector of age group 1:7
  prob_vent_given_age <- function(age_vec, tx, tx_effect = 0){
    prob_vent_by_age <- c( 0.0205, 0.1755, 0.2475, 0.2530, 0.3605, 0.4460, 0.5080)
    intercept_by_age <- qlogis(prob_vent_by_age)
    n <- length(age_vec)
    prob_vent <- plogis(intercept_by_age[age_vec] + tx_effect * tx)
    return(prob_vent)
  }
  prob_vent <- prob_vent_given_age(age_vec = age_grp, tx = treat, 
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
  prob_vent_given_age <- function(age_vec, tx, tx_effect = 0){
    prob_vent_by_age <- c( 0.0205, 0.1755, 0.2475, 0.2530, 0.3605, 0.4460, 0.5080)
    intercept_by_age <- qlogis(prob_vent_by_age)
    n <- length(age_vec)
    prob_vent <- plogis(intercept_by_age[age_vec] + tx_effect * tx)
    return(prob_vent)
  }
  prob_vent1 <- prob_vent_given_age(age_vec = age_grp, tx = 1, 
                                   tx_effect = tx_effect)
  prob_vent0 <- prob_vent_given_age(age_vec = age_grp, tx = 0, 
                                   tx_effect = tx_effect)

  prob_nothing1 <- 1 - prob_death - prob_vent1
  prob_nothing0 <- 1 - prob_death - prob_vent0
  
  F0 <- c(mean(prob_death), mean(prob_death + prob_vent0), 1)
  F1 <- c(mean(prob_death), mean(prob_death + prob_vent1), 1)
  f_0 <- c(mean(prob_death), mean(prob_vent0), mean(prob_nothing0))
  f_1 <- c(mean(prob_death), mean(prob_vent1), mean(prob_nothing1))

  # weighted mean
  mean_y1 <- sum(c(1,1,0)*f_1)
  mean_y0 <- sum(c(1,1,0)*f_0)
  wmean_truth <- mean_y1 - mean_y0
  
  return(list(mann_whitney = mannwhitney_truth))
}

# get the list size #########
if (args[1] == 'listsize') {
  cat(nrow(parm))
}

# execute prepare job ##################
if (args[1] == 'prepare') {
  print(paste0('initial datasets saved to: ~/drinf/scratch/dataList ... .RData'))
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
      # make data set
      get_tx_eff <- function(n){
        if(n == 100){
          return(-1.83)
        }else if(n == 200){
          return(-1.13)
        }else if(n == 500){
          return(-0.625)
        }else{
          return(-0.438)
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
      # make data set
      get_tx_eff <- function(n){
        if(n == 100){
          return(-1.85)
        }else if(n == 200){
          return(-1.18)
        }else if(n == 500){
          return(-0.725)
        }else{
          return(-0.443)
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

    # binarize outcome
    out_bin <- dat$out
    out_bin[dat$out == 3] <- 0
    out_bin[dat$out == 2] <- 1

    rslt_adj <- drord(
      out = out_bin,
      treat = dat$treat,
      covar = dat$covar,
      out_form = "age_grp",
      param = c("weighted_mean"),
      ci = c("bca", "wald"),
      nboot = 1e3,
      est_dist = FALSE
    )

    # get output for unadjusted mw and weighted-mean
    unadj_0 <- mean(out_bin[dat$treat == 0])
    unadj_1 <- mean(out_bin[dat$treat == 1])
    pval <- chisq.test(table(out_bin, dat$treat))$p.value

    in_ci <- function(ci, truth){
      truth < max(ci) & truth > min(ci)
    }

    format_out <- function(rslt){
      # format output
      wmean_out <- c(rslt$weighted_mean$est$est[1] - rslt$weighted_mean$est$est[2],
                     rslt$weighted_mean$ci$wald[3,],
                     in_ci(ci = rslt$weighted_mean$ci$wald[3,], truth = truth$weighted_mean),
                     rslt$weighted_mean$ci$bca[3,],
                     in_ci(ci = rslt$weighted_mean$ci$bca[3,], truth = truth$weighted_mean))
      return(c(wmean_out))
    }

    out <- c(as.numeric(parm[i,]), truth$weighted_mean, 
             format_out(rslt = rslt_adj),
             c(unadj_1 - unadj_0, pval))

    # save output
    save(out, file = paste0(save_dir, "fit_", row.names(parm)[i], ".RData"))
  }
}

# merge job ###########################
if (args[1] == 'merge'){
  out_names <- c(
    "n", "seed", "trt_effect", "dgp",
    "wmean_truth", 
    c(t(outer(c("wmean_adj"), c("_est", "_wald_cil", "_wald_ciu",
                                                 "_wald_cover","_bca_cil","_bca_ciu",
                                                 "_bca_cover"), paste0))),
    c("wmean_unadj_est", "wmean_unadj_pval")
  )

  n_out <- length(out_names)

  rslt <- matrix(NA, nrow = nrow(parm), ncol = n_out)
  for(i in 1:nrow(parm)){
    tmp <- tryCatch({
      load(paste0(save_dir, "final_fit_binary_", i, ".RData"))
      out
    }, error = function(e){
      rep(NA, n_out)
    })
    rslt[i, ] <- tmp
  }

  out <- data.frame(rslt)
  colnames(out) <- out_names

  make_output_table <- function(out, estimand = "wmean", ci){
    null_val <- 0
    
    sample_size <- sort(rep(c(100, 200, 500, 1000), 4))
    est_type <- rep(c("Unadjusted", "Adjusted"), 8)

    effect_col <- unlist(by(out, out$n, function(x){
      c(rep(x[x$trt_effect == 0, paste0(estimand, "_truth")][1],2),
        rep(x[x$trt_effect == 1, paste0(estimand, "_truth")][1],2))
    }), use.names = FALSE)

    power_col <- unlist(by(out, out$n, function(x){
      c(mean(na.rm = TRUE, as.numeric(x$wmean_unadj_pval[x$trt_effect == 0] < 0.05)),
        mean(na.rm = TRUE, as.numeric(x[x$trt_effect == 0, paste0(estimand, "_adj_", ci, "_cil")] > null_val |
                                      x[x$trt_effect == 0, paste0(estimand, "_adj_", ci, "_ciu")] < null_val)),
        mean(na.rm = TRUE, as.numeric(x$wmean_unadj_pval[x$trt_effect == 1] < 0.05)),
        mean(na.rm = TRUE, as.numeric(x[x$trt_effect == 1, paste0(estimand, "_adj_", ci, "_cil")] > null_val |
                                      x[x$trt_effect == 1, paste0(estimand, "_adj_", ci, "_ciu")] < null_val)))    
      }), use.names = FALSE)

    mse_col <- unlist(by(out, out$n, function(x){
      c(mean(na.rm = TRUE, (x[x$trt_effect == 0, paste0(estimand, "_unadj_est")] - x[x$trt_effect == 0, paste0(estimand, "_truth")])^2),
        mean(na.rm = TRUE, (x[x$trt_effect == 0, paste0(estimand, "_adj_est")] - x[x$trt_effect == 0, paste0(estimand, "_truth")])^2),
        mean(na.rm = TRUE, (x[x$trt_effect == 1, paste0(estimand, "_unadj_est")] - x[x$trt_effect == 1, paste0(estimand, "_truth")])^2),
        mean(na.rm = TRUE, (x[x$trt_effect == 1, paste0(estimand, "_adj_est")] - x[x$trt_effect == 1, paste0(estimand, "_truth")])^2))
    }), use.names = FALSE)

    bias_col <- unlist(by(out, out$n, function(x){
      c(mean(na.rm = TRUE, (x[x$trt_effect == 0, paste0(estimand, "_unadj_est")] - x[x$trt_effect == 0, paste0(estimand, "_truth")])),
        mean(na.rm = TRUE, (x[x$trt_effect == 0, paste0(estimand, "_adj_est")] - x[x$trt_effect == 0, paste0(estimand, "_truth")])),
        mean(na.rm = TRUE, (x[x$trt_effect == 1, paste0(estimand, "_unadj_est")] - x[x$trt_effect == 1, paste0(estimand, "_truth")])),
        mean(na.rm = TRUE, (x[x$trt_effect == 1, paste0(estimand, "_adj_est")] - x[x$trt_effect == 1, paste0(estimand, "_truth")])))
    }), use.names = FALSE)
    
    var_col <- unlist(by(out, out$n, function(x){
      c(var(na.rm = TRUE, x[x$trt_effect == 0, paste0(estimand, "_unadj_est")]),
        var(na.rm = TRUE, x[x$trt_effect == 0, paste0(estimand, "_adj_est")]),
        var(na.rm = TRUE, x[x$trt_effect == 1, paste0(estimand, "_unadj_est")]),
        var(na.rm = TRUE, x[x$trt_effect == 1, paste0(estimand, "_adj_est")]))
    }), use.names = FALSE)
    rel_eff <- rep(NA, length(mse_col))
    rel_eff[seq(1, length(mse_col), by = 2)] <- 1
    rel_eff[seq(2, length(mse_col), by = 2)] <- mse_col[seq(2, length(mse_col), 2)] / mse_col[seq(1, length(mse_col), 2)]

    tab <- data.frame(sample_size, est_type, effect_col, power_col, mse_col, bias_col, var_col, rel_eff)

    return(tab)
  }

  library(xtable)
  # dgp 1
  tab_parm <- expand.grid(estimand = c("wmean"),
                          ci = c("wald", "bca"))
  tab_list_1 <- vector(mode = "list", length = nrow(tab_parm))
  for(i in seq_len(nrow(tab_parm))){
    tab_list_1[[i]] <- make_output_table(out[out$dgp == 1,], estimand = tab_parm$estimand[i], ci = tab_parm$ci[i])
  }

  # BCA
  print(xtable(tab_list_1[[2]], digits = c(0,0,rep(3,7))), 
        include.rownames = FALSE,
        hline.after = c(4, 8, 12))
  # WALD
  print(xtable(tab_list_1[[1]], digits = c(0,0,rep(3,7))), 
        include.rownames = FALSE,
        hline.after = c(4, 8, 12))

  # dgp 3
  tab_list_3 <- vector(mode = "list", length = nrow(tab_parm))
  for(i in seq_len(nrow(tab_parm))){
    tab_list_3[[i]] <- make_output_table(out[out$dgp == 3,], estimand = tab_parm$estimand[i], ci = tab_parm$ci[i])
  }

  # BCA
  print(xtable(tab_list_1[[2]], digits = c(0,0,rep(3,7))), 
      include.rownames = FALSE,
      hline.after = c(4, 8, 12))
  # WALD
  print(xtable(tab_list_1[[1]], digits = c(0,0,rep(3,7))), 
        include.rownames = FALSE,
        hline.after = c(4, 8, 12))
}












