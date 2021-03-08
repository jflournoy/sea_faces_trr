#Run the old ICC on parcel averages using brms
library(brms)
library(data.table)
later:::ensureInitialized()

task_id <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
setDTthreads(4)

if(is.na(task_id)){
  message('Not running as slurm array, setting ROI to some value')
  this_ROI <- 1 #TO BE SET BY SLURM ENV VAR
} else {
  this_ROI <- task_id
}

all_d <- readRDS('sch400_cope_fGTc.RDS')
this_d <- all_d[ROI == this_ROI]

chains <- 4
iterations <- 2000

get_icc <- function(d, sess_tag){
  icc_form <- bf(y ~ 1 + sess + (1 | id))
  compiled_icc_mod <- brm(icc_form, data = d,
                          family = 'gaussian',
                          cores = 1, chains = 1, iter = 1,
                          control = list(adapt_delta = 0.99, max_treedepth = 15), 
                          file = 'fits/compiled_icc')
  fm <- update(compiled_icc_mod, formula. = icc_form, newdata = d, 
               family = 'gaussian',
               chains = chains, cores = chains, iter = iterations, 
               control = list(adapt_delta = 0.99, max_treedepth = 15), 
               file = paste0('fits/fit-ROI', this_ROI, '-', sess_tag, '-icc'))         

  vc <- VarCorr(fm, summary = FALSE)
  
  sigma_lambda_2 <- vc$id$sd[, 'Intercept']^2
  sigma_eps_2 <- vc$residual__$sd^2
  ICC_posterior <- sigma_lambda_2 / (sigma_lambda_2 + sigma_eps_2)
  dimnames(ICC_posterior)$parameters <- 'ICC'
  # bayesplot::mcmc_areas(as.mcmc(ICC_posterior), prob = .95, prob_outer = .99)
  r <- as.data.table(posterior_summary(ICC_posterior))
  r[, sess := sess_tag]
  r[, est_logit := mean(arm::logit(ICC_posterior))]
  r[, sd_logit := sd(arm::logit(ICC_posterior))]
  return(r)
}
unique_sessions <- sort(unique(this_d$sess))

retvals <- lapply(1:(length(unique_sessions) - 1), function(sess_i){
  d <- this_d[sess %in% unique_sessions[sess_i:(sess_i+1)]]
  sess_tag <- gsub('month', '', paste(unique_sessions[sess_i:(sess_i+1)], collapse = '_'))
  return(get_icc(d, sess_tag = sess_tag))
})

meta_icc_dt <- rbindlist(retvals)
meta_icc_form <- bf(est_logit | se(sd_logit, sigma = TRUE) ~ 1)
compiled_meta_fit <- brm(meta_icc_form, data = meta_icc_dt,
                family = 'gaussian',
                chains = 1, cores = 1, iter = 1, 
                control = list(adapt_delta = 0.99, max_treedepth = 15),
                file = 'fits/copiled_icc_meta')
meta_fit <- update(compiled_meta_fit, formula. = meta_icc_form, newdata = meta_icc_dt,
                family = 'gaussian',
                chains = chains, cores = chains, iter = iterations, 
                control = list(adapt_delta = 0.99, max_treedepth = 15),
                file = paste0('fits/fit-ROI', this_ROI, '-meta_icc'))

meta_sum <- data.table(posterior_summary(arm::invlogit(posterior_samples(meta_fit, pars = 'b_Intercept'))))
meta_sum[, sess := 'meta']
all_sess_icc <- get_icc(this_d, sess_tag = 'all')

sess_iccs <- rbindlist(retvals)[, c('est_logit', 'sd_logit') := NULL]
all_sess_icc[, c('est_logit', 'sd_logit') := NULL]
retvals_dt <- rbindlist(list(sess_iccs, meta_sum, all_sess_icc))
retvals_dt[, ROI := this_ROI]
saveRDS(retvals_dt, paste0('fits/oldicc-', this_ROI, '.RDS'))

##Using the new data

fd <- readRDS('sch400_RML_fit_processed.RDS')
fd <- d[cond %in% c('Fear', 'Calm')]

faccols <- c('sess', 'cond', 'id')
fd[, (faccols) := lapply(.SD, as.factor), .SDcols = faccols]
fd[, cond_code := fifelse(cond == levels(cond)[1], -0.5, 0.5)]

this_fd <- fd[ROI %in% this_ROI]

get_icc_fd <- function(d, sess_tag){
  icc_form <- bf(est|se(se, sigma = TRUE) ~ 1 + cond_code + sess + cond_code:sess + (1 + cond_code | id))
  compiled_icc_mod <- brm(icc_form, data = d,
                          family = 'gaussian',
                          cores = 1, chains = 1, iter = 1,
                          control = list(adapt_delta = 0.99, max_treedepth = 15), 
                          file = 'fits/compiled_icc_fulldat')
  fm <- update(compiled_icc_mod, formula. = icc_form, newdata = d, 
               family = 'gaussian',
               chains = chains, cores = chains, iter = iterations, 
               control = list(adapt_delta = 0.99, max_treedepth = 15), 
               file = paste0('fits/fit-ROI', this_ROI, '-', sess_tag, '-icc-fulldat'))         
  
  vc <- VarCorr(fm, summary = FALSE)
  
  sigma_lambda_2 <- vc$id$sd[, 'cond_code']^2
  sigma_lambda_2_int <- vc$id$sd[, 'Intercept']^2
  sigma_eps_2 <- vc$residual__$sd^2
  ICC_posterior <- sigma_lambda_2 / (sigma_lambda_2 + sigma_lambda_2_int + sigma_eps_2)
  dimnames(ICC_posterior)$parameters <- 'ICC'
  # bayesplot::mcmc_areas(as.mcmc(ICC_posterior), prob = .95, prob_outer = .99)
  r <- as.data.table(posterior_summary(ICC_posterior))
  r[, sess := sess_tag]
  r[, est_logit := mean(arm::logit(ICC_posterior))]
  r[, sd_logit := sd(arm::logit(ICC_posterior))]
  return(r)
}
unique_sessions_fd <- sort(unique(this_fd$sess))

retvals_fd <- lapply(1:(length(unique_sessions_fd) - 1), function(sess_i){
  d <- this_fd[sess %in% unique_sessions_fd[sess_i:(sess_i+1)]]
  sess_tag <- paste(unique_sessions_fd[sess_i:(sess_i+1)], collapse = '_')
  return(get_icc_fd(d, sess_tag = sess_tag))
})

meta_icc_fd_dt <- rbindlist(retvals_fd)
meta_fit_fd <- update(compiled_meta_fit, formula. = meta_icc_form, newdata = meta_icc_fd_dt,
                   family = 'gaussian',
                   chains = chains, cores = chains, iter = iterations, 
                   control = list(adapt_delta = 0.99, max_treedepth = 15),
                   file = paste0('fits/fit-ROI', this_ROI, '-meta_icc-fulldat'))

meta_fd_sum <- data.table(posterior_summary(arm::invlogit(posterior_samples(meta_fit_fd, pars = 'b_Intercept'))))
meta_fd_sum[, sess := 'meta']
all_sess_icc_fd <- get_icc_fd(this_fd, sess_tag = 'all')

sess_iccs_fd <- rbindlist(retvals_fd)[, c('est_logit', 'sd_logit') := NULL]
# all_sess_icc_fd[, c('est_logit', 'sd_logit') := NULL]
retvals_fd_dt <- rbindlist(list(sess_iccs_fd, meta_fd_sum))
retvals_fd_dt[, ROI := this_ROI]
saveRDS(retvals_dt, paste0('fits/oldicc-', this_ROI, '-fulldat.RDS'))
