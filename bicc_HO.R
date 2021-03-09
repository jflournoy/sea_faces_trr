#Much taken from: https://github.com/afni/afni/blob/2e428852abd096340cf8133f014e68093b65d773/src/R_scripts/TRR.R
library(brms)
library(coda)
library(data.table)
library(ggplot2) 

task_id <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
setDTthreads(4)

if(is.na(task_id)){
  message('Not running as slurm array, setting ROI to some value')
  this_ROI <- 243 #TO BE SET BY SLURM ENV VAR
} else {
  this_ROI <- task_id
}

d <- readRDS('sch400_RML_fit_processed.RDS')
d[, ROI := as.numeric(ROI) + 1] #account for 0 indexing

if(! identical(range(as.numeric(d$ROI)), c(1, 400)) ){
  stop("Error in ROI Indexing: \n```\nd <- readRDS('sch400_RML_fit_processed.RDS')\nd[, ROI := as.numeric(ROI) + 1]\n```")
}

d <- d[cond %in% c('Fear', 'Calm')]

faccols <- c('sess', 'cond', 'id')
d[, (faccols) := lapply(.SD, as.factor), .SDcols = faccols]
d[, cond_code := fifelse(cond == levels(cond)[1], -0.5, 0.5)]

uq_sess <- unique(d$sess)
sess_pairs <- lapply(1:(length(uq_sess)-1), function(i){
  pair <- c(as.character(uq_sess[[i]]), as.character(uq_sess[[i + 1]]))
  name <- paste(gsub(pattern = 'month', replacement = '', x = pair), collapse = '_')
  return(list(pair = pair, name = name))
})

iterations <- 3000
chains <- 4
grouping_term <- 'id'
modelForm <- as.formula(est|se(se, sigma = TRUE) ~ 0 + sess + sess:cond_code + (0 + sess | id) + (0 + sess:cond_code | id))
bigCovarModelForm <- bf(est|se(se, sigma = TRUE) ~ 0 + sess + sess:cond_code + (0 + sess | id) + (0 + sess:cond_code | id))

### Big Full Covar model:
d_bigmodel <- d[ROI == this_ROI]
compiled_bigcovmodel<- brm(bigCovarModelForm, data = d_bigmodel, 
                           family = 'gaussian',
                           chains = 1, cores = 1,
                           iter = 1, control = list(adapt_delta = 0.99, max_treedepth = 15), 
                           file = 'fits/compiled_bigcov')
# summary(compiled_bigcovmodel)
fm_bigcov <- update(compiled_bigcovmodel, formula. = bigCovarModelForm, newdata = d_bigmodel, 
                    family = 'gaussian',
                    chains = 4, cores = 4,
                    iter = 2000, control = list(adapt_delta = 0.99, max_treedepth = 15), 
                    file = paste0('fits/fit-ROI', this_ROI, '-bigcov'))

summary(fm_bigcov)
parnames(fm_bigcov)
fm_bigcov_arr <- as.array(fm_bigcov, pars = c(sprintf('cor_id__sess%02d:cond_code__sess%02d:cond_code', 1:9, 2:10),
                                              sprintf('cor_id__sess%02d__sess%02d', 1:9, 2:10)))
dimnames(fm_bigcov_arr)$parameters <- c(sprintf('cor(sess%02d,sess%02d)_cont', 1:9, 2:10),
                                        sprintf('cor(sess%02d,sess%02d)_avg', 1:9, 2:10))
p1 <- bayesplot::mcmc_areas(fm_bigcov_arr, pars = sprintf('cor(sess%02d,sess%02d)_cont', 1:9, 2:10), 
                             prob = .95, prob_outer = .99) 
p2 <- bayesplot::mcmc_areas(fm_bigcov_arr, pars = sprintf('cor(sess%02d,sess%02d)_avg', 1:9, 2:10), 
                             prob = .95, prob_outer = .99)
ggsave(paste0('plots/plot-cov-ROI_', this_ROI, '-con.pdf'), plot = p1)
ggsave(paste0('plots/plot-cov-ROI_', this_ROI, '-ave.pdf'), plot = p2)

#Extract correlations from this model to save
##metanalysis of fisher-z transformed.

fm_bigcov_post_samp <- posterior_samples(fm_bigcov, pars = c(sprintf('cor_id__sess%02d:cond_code__sess%02d:cond_code', 1:9, 2:10),
                                      sprintf('cor_id__sess%02d__sess%02d', 1:9, 2:10)))
fm_bigcov_post_sum <- t(sapply(fm_bigcov_post_samp, function(x) {
  c(trr = mean(x), 
    trr_lower = quantile(x, probs = .025)[[1]], trr_upper = quantile(x, probs = .975)[[1]],
    z = mean(atanh(x)), z_se = sd(atanh(x)))
}))

fm_bigcov_post_sum_rn <- rownames(fm_bigcov_post_sum)
fm_bigcov_sum_dt <- data.table(fm_bigcov_post_sum)
fm_bigcov_sum_dt[, effect := fifelse(grepl('cond_code', fm_bigcov_post_sum_rn),
                                 'contrast', 'average')]
fm_bigcov_sum_dt[, rep_name := gsub('cor_id__sess(\\d{2})(?::cond_code)*__sess(\\d{2})(?::cond_code)*', '\\1_\\2', fm_bigcov_post_sum_rn)]

meta_form <- bf(z|se(z_se, sigma = TRUE) ~ 1)
meta_compiled <- brms::brm(meta_form, data = fm_bigcov_sum_dt,
                  cores = 1, chains = 1, iter = 1,
                  control = list(adapt_delta = 0.99, max_treedepth = 15),
                  file = 'fits/compiled_meta')
meta_fit_cont <- update(meta_compiled, formula. = meta_form, newdata = fm_bigcov_sum_dt[effect == 'contrast'],
                        family = 'gaussian',
                        chains = chains, cores = chains, iter = iterations, 
                        control = list(adapt_delta = 0.99, max_treedepth = 15),
                        file = paste0('fits/fit-ROI', this_ROI, '-meta_cont'))
meta_fit_avg <- update(meta_compiled, formula. = meta_form, newdata = fm_bigcov_sum_dt[effect == 'average'],
                        family = 'gaussian',
                        chains = chains, cores = chains, iter = iterations, 
                        control = list(adapt_delta = 0.99, max_treedepth = 15),
                        file = paste0('fits/fit-ROI', this_ROI, '-meta_avg'))
p1_comb <- bayesplot::mcmc_areas(meta_fit_cont, pars = 'b_Intercept', 
                            prob = .95, prob_outer = .99) 
p2_comb <- bayesplot::mcmc_areas(meta_fit_avg, pars = 'b_Intercept', 
                            prob = .95, prob_outer = .99)
ggsave(paste0('plots/plot-covmeta-ROI_', this_ROI, '-con.pdf'), plot = p1_comb)
ggsave(paste0('plots/plot-covmeta-ROI_', this_ROI, '-avg.pdf'), plot = p2_comb)

#Maybe do AR model from above model


#each session pair separately.
d_compile <- d[sess %in% c("01", "02") &
                 ROI == this_ROI]

compiled_model <- brm(modelForm, data = d_compile, 
                      family = 'gaussian',
                      chains = 1, cores = 1,
                      iter = 1, control = list(adapt_delta = 0.99, max_treedepth = 15), 
                      file = 'fits/compiled')

r <- lapply(sess_pairs, function(sess_pair){
  apair <- sess_pair[['pair']]
  apair_name <- sess_pair[['name']]
  
  d_chunk <- d[sess %in% apair &
                 ROI == this_ROI]
  fm <- update(compiled_model, formula. = modelForm, newdata = d_chunk, 
               family = 'gaussian',
               chains = chains, cores = chains,
               iter = iterations, control = list(adapt_delta = 0.99, max_treedepth = 15), 
               file = paste0('fits/fit-ROI', this_ROI, '-',  apair_name))
  rs <- summary(fm)
  return(list(rs = rs, fm = fm, name = apair_name))
})

sum_r <- lapply(r, function(r_){
  fm <- r_[['fm']]
  r_name <- r_[['name']]
  
  pE <- fixef(fm, summary = FALSE)
  
  pars <- parnames(fm)
  corsessname <- pars[grepl('cor_\\w+__sess\\d{2}__sess\\d{2}$', pars)]
  corsesscondname <- pars[grepl('cor_\\w+__sess\\d{2}:cond_code__sess\\d{2}:cond_code$', pars)]
  ll <- posterior_samples(fm, pars = corsessname)[[1]]
  mm <- posterior_samples(fm, pars = corsesscondname)[[1]]
  
  m <- median
  
  lp <- m(ll); mp <- m(mm)
  
  lq <- HPDinterval(as.mcmc(ll)); mq <- HPDinterval(as.mcmc(mm))
  
  d_res <- data.frame(
    rep_name = r_name,
    effect = c('average', 'contrast'),
    mean = c(mean((pE[,1]+pE[,2])/2), mean((pE[,3]+pE[,4])/2)),
    sd = c(sd((pE[,1]+pE[,2])/2), sd((pE[,3]+pE[,4])/2)),
    trr = c(lp, mp),
    trr_lower = c(lq[1], mq[1]),
    trr_upper = c(lq[2], mq[2]))
  
  p1 <- ggplot(data.frame(x=ll), aes(x=x)) + theme_bw() + geom_density(size=2) +
    geom_vline(aes(xintercept= lp), color="blue", linetype="dashed", size=2) +
    geom_area(data = subset(data.frame(x = density(ll)$x, y = density(ll)$y),
                            x >= lq[1] & x <= lq[2]), aes(x=x,y=y, color='gray'), alpha=0.15) +
    xlab("TRR") +
    theme(legend.position="none", axis.text=element_text(size=18), axis.title=element_text(size=18))
  # ggsave(paste0('plots/plot-ROI_', this_ROI, '-',  r_name, '-ave.pdf'), plot = p1)
  # 
  p2 <- ggplot(data.frame(x=mm), aes(x=x)) + theme_bw() + geom_density(size=2) +
    geom_vline(aes(xintercept= mp), color="blue", linetype="dashed", size=2) +
    geom_area(data = subset(data.frame(x = density(mm)$x, y = density(mm)$y),
                            x >= mq[1] & x <= mq[2]), aes(x=x,y=y, color='gray'), alpha=0.15) +
    xlab("TRR") +
    theme(legend.position="none", axis.text=element_text(size=18), axis.title=element_text(size=18))
  # ggsave(paste0('plots/plot-ROI_', this_ROI, '-',  r_name, '-con.pdf'), plot = p2)
  
  return(list(p_ave = p1, p_con = p2, summary = d_res))
})

big_df <- data.table::rbindlist(lapply(sum_r, `[[`, 'summary'))
# readr::write_csv(big_df, paste0('csv/TRR-ROI_', this_ROI, '.csv'))
comp_plot <- ggplot(merge(fm_bigcov_sum_dt, big_df, by = c('effect', 'rep_name'), suffix = c('_cov', '_sep')),
       aes(x = trr_cov, y = trr_sep)) + 
  geom_abline(intercept = 0, slope = 1, color = 'red') + 
  geom_errorbar(aes(xmin = trr_lower_cov, xmax = trr_upper_cov), alpha = .2) + 
  geom_errorbar(aes(ymin = trr_lower_sep, ymax = trr_upper_sep), alpha = .2) + 
  geom_point() + 
  coord_cartesian(x = c(-1, 1), y = c(-1, 1)) + 
  facet_grid(effect ~ .)
ggsave(paste0('plots/plot-comp-ROI_', this_ROI, '.pdf'), plot = comp_plot)

####--------------
### Harvard 0xford
##----------------

if(is.na(task_id)){
  message('Not running as slurm array, setting ROI to some value')
  this_ho_ROI <- 10 #TO BE SET BY SLURM ENV VAR
} else {
  if(task_id < 1 | task_id > 21) stop(sprintf('HO ROI task-id out of range: %d', task_id))
  this_ho_ROI <- task_id
}

dh <- readRDS('HO_RML_fit_processed.RDS')
dh[, ROI := as.numeric(ROI) + 1] #account for 0 indexing

if(! identical(range(as.numeric(dh$ROI)), c(1, 21)) ){
  stop("Error in HO ROI Indexing")
}

dh <- dh[cond %in% c('Fear', 'Calm')]

faccols <- c('sess', 'cond', 'id')
dh[, (faccols) := lapply(.SD, as.factor), .SDcols = faccols]
dh[, cond_code := fifelse(cond == levels(cond)[1], -0.5, 0.5)]

### Big Full Covar model:
d_ho_bigmodel <- dh[ROI == this_ho_ROI]
fm_ho_bigcov <- update(compiled_bigcovmodel, formula. = bigCovarModelForm, newdata = d_ho_bigmodel, 
                       family = 'gaussian',
                       chains = chains, cores = chains,
                       iter = iterations, control = list(adapt_delta = 0.99, max_treedepth = 15), 
                       file = paste0('fits/fit-HO-ROI', this_ho_ROI, '-bigcov'))

# summary(fm_ho_bigcov)
# parnames(fm_ho_bigcov)
fm_ho_bigcov_arr <- as.array(fm_ho_bigcov, pars = c(sprintf('cor_id__sess%02d:cond_code__sess%02d:cond_code', 1:9, 2:10),
                                              sprintf('cor_id__sess%02d__sess%02d', 1:9, 2:10)))
dimnames(fm_ho_bigcov_arr)$parameters <- c(sprintf('cor(sess%02d,sess%02d)_cont', 1:9, 2:10),
                                        sprintf('cor(sess%02d,sess%02d)_avg', 1:9, 2:10))
p1 <- bayesplot::mcmc_areas(fm_ho_bigcov_arr, pars = sprintf('cor(sess%02d,sess%02d)_cont', 1:9, 2:10), 
                            prob = .95, prob_outer = .99) 
p2 <- bayesplot::mcmc_areas(fm_ho_bigcov_arr, pars = sprintf('cor(sess%02d,sess%02d)_avg', 1:9, 2:10), 
                            prob = .95, prob_outer = .99)
ggsave(paste0('plots/plot-cov-HO-ROI_', this_ho_ROI, '-con.pdf'), plot = p1)
ggsave(paste0('plots/plot-cov-HO-ROI_', this_ho_ROI, '-ave.pdf'), plot = p2)

#Extract correlations from this model to save
##metanalysis of fisher-z transformed.

fm_ho_bigcov_post_samp <- posterior_samples(fm_ho_bigcov, pars = c(sprintf('cor_id__sess%02d:cond_code__sess%02d:cond_code', 1:9, 2:10),
                                                             sprintf('cor_id__sess%02d__sess%02d', 1:9, 2:10)))
fm_ho_bigcov_post_sum <- t(sapply(fm_ho_bigcov_post_samp, function(x) {
  c(trr = mean(x), 
    trr_lower = quantile(x, probs = .025)[[1]], trr_upper = quantile(x, probs = .975)[[1]],
    z = mean(atanh(x)), z_se = sd(atanh(x)))
}))

fm_ho_bigcov_post_sum_rn <- rownames(fm_ho_bigcov_post_sum)
fm_ho_bigcov_sum_dt <- data.table(fm_ho_bigcov_post_sum)
fm_ho_bigcov_sum_dt[, effect := fifelse(grepl('cond_code', fm_ho_bigcov_post_sum_rn),
                                     'contrast', 'average')]
fm_ho_bigcov_sum_dt[, rep_name := gsub('cor_id__sess(\\d{2})(?::cond_code)*__sess(\\d{2})(?::cond_code)*', '\\1_\\2', fm_ho_bigcov_post_sum_rn)]

meta_fit_ho_cont <- update(meta_compiled, formula. = meta_form, newdata = fm_ho_bigcov_sum_dt[effect == 'contrast'],
                        family = 'gaussian',
                        chains = chains, cores = chains, iter = iterations, 
                        control = list(adapt_delta = 0.99, max_treedepth = 15),
                        file = paste0('fits/fit-HO-ROI', this_ROI, '-meta_cont'))
meta_fit_ho_avg <- update(meta_compiled, formula. = meta_form, newdata = fm_ho_bigcov_sum_dt[effect == 'average'],
                       family = 'gaussian',
                       chains = chains, cores = chains, iter = iterations, 
                       control = list(adapt_delta = 0.99, max_treedepth = 15),
                       file = paste0('fits/fit-HO-ROI', this_ROI, '-meta_avg'))
p1_comb <- bayesplot::mcmc_areas(meta_fit_ho_cont, pars = 'b_Intercept', 
                                 prob = .95, prob_outer = .99) 
p2_comb <- bayesplot::mcmc_areas(meta_fit_ho_avg, pars = 'b_Intercept', 
                                 prob = .95, prob_outer = .99)
ggsave(paste0('plots/plot-covmeta-HO-ROI_', this_ROI, '-con.pdf'), plot = p1_comb)
ggsave(paste0('plots/plot-covmeta-HO-ROI_', this_ROI, '-avg.pdf'), plot = p2_comb)

if(FALSE){
  ### Big model + 2nd stage AR approach:
  bigModelForm <- bf(est|se(se, sigma = TRUE) ~ 0 + sess + sess:cond_code + (0 + sess || id) + (0 + sess:cond_code || id))
  
  compiled_bigmodel<- brm(bigModelForm, data = d_bigmodel, 
                          family = 'gaussian',
                          chains = 1, cores = 1,
                          iter = 1, control = list(adapt_delta = 0.99, max_treedepth = 15), 
                          file = 'fits/compiled_big')
  
  fm_big <- update(compiled_bigmodel, formula. = bigModelForm, newdata = d_bigmodel, 
                   family = 'gaussian',
                   chains = 4, cores = 4,
                   iter = 2000, control = list(adapt_delta = 0.99, max_treedepth = 15), 
                   file = paste0('fits/fit-ROI', this_ROI, '-big'))
  
  summary(fm_big)
  big_RE <- ranef(fm_big, summary =  FALSE) #check that these are right
  dim(big_RE$id)
  sess_contrast_ids <- grep('cond', dimnames(big_RE$id)[[3]])
  mean_dt <- melt(as.data.table(apply(big_RE$id[,,sess_contrast_ids], c(2,3), mean), keep.rownames = 'id'),
                  value.name = 'mean', id.vars = 'id')
  sd_dt <- melt(as.data.table(apply(big_RE$id[,,sess_contrast_ids], c(2,3), sd), keep.rownames = 'id'),
                value.name = 'sd', id.vars = 'id')
  
  re_dt <- mean_dt[sd_dt, on = c('id', 'variable')] 
  re_dt[, time := as.numeric(gsub('sess(\\d{2}):cond_code', '\\1', variable))]
  ar_model <- bf(mean | se(sd, sigma = TRUE) ~ ar(time = time, gr = id, cov = TRUE))
  brms::get_prior(ar_model, data = re_dt)
  
  compiled_ARmodel<- brm(ar_model, data = re_dt, 
                         family = 'gaussian',
                         chains = 1, cores = 1,
                         iter = 1, control = list(adapt_delta = 0.99, max_treedepth = 15), 
                         file = 'fits/compiled_AR')
  fm_bigAR <- update(compiled_ARmodel, formula. = ar_model, newdata = re_dt, 
                     family = 'gaussian',
                     chains = 4, cores = 4,
                     iter = 2000, control = list(adapt_delta = 0.99, max_treedepth = 15), 
                     file = paste0('fits/fit-ROI', this_ROI, '-bigAR'))
  summary(fm_bigAR)
}