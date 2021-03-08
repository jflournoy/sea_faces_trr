#Largely taken from: https://github.com/afni/afni/blob/2e428852abd096340cf8133f014e68093b65d773/src/R_scripts/TRR.R
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
bigModelForm <- bf(est|se(se, sigma = TRUE) ~ 0 + sess + sess:cond_code + (0 + sess || id) + (0 + sess:cond_code || id))
bigCovarModelForm <- bf(est|se(se, sigma = TRUE) ~ 0 + sess + sess:cond_code + (0 + sess | id) + (0 + sess:cond_code | id))

d_compile <- d[sess %in% c("01", "02") &
                 ROI == this_ROI]

compiled_model <- brm(modelForm, data = d_compile, 
          family = 'gaussian',
          chains = 1, cores = 1,
          iter = 1, control = list(adapt_delta = 0.99, max_treedepth = 15), 
          file = 'fits/compiled')

outfileRDS <- paste0('raw/ROI_', this_ROI, '.RDS')
if(!file.exists(outfileRDS)){
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
  
  saveRDS(r, outfileRDS)  
} else {
  r <- readRDS(outfileRDS)
}

sum_r <- lapply(r, function(r_){
  fm <- r_[['fm']]
  r_name <- r_[['name']]
  print(r_name)
  
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
  # ggsave(paste0('plots/plot-', r_name, '-ave.pdf'), plot = p1)
  
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
readr::write_csv(big_df, paste0('csv/TRR-ROI_', this_ROI, '.csv'))

##ADD metanalysis of fisher-z transformed.

TRR_z_cont_est_se <- data.table::rbindlist(lapply(r, function(r_){
  fm <- r_[['fm']]

  pars <- parnames(fm)
  corsesscondname <- pars[grepl('cor_\\w+__sess\\d{2}:cond_code__sess\\d{2}:cond_code$', pars)]
  mm <- posterior_samples(fm, pars = corsesscondname)[[1]]

  z <- atanh(mm)
  zp <- mean(z)
  zse <- sd(z)
  return(data.frame(z_hat = zp, zse = zse, sid = r_[['name']]))
}))

priors <- c(prior(student_t(3, 0, .5), class = Intercept))
get_prior(z_hat|se(zse, sigma = FALSE) ~ 1 + (1 | sid), data = TRR_z_cont_est_se)
meta <- brms::brm(z_hat|se(zse, sigma = FALSE) ~ 1 + (1 | sid), data = TRR_z_cont_est_se,
                  cores = chains, chains = chains, iter = 4000,
                  control = list(adapt_delta = 0.99, max_treedepth = 15),
                  prior = priors,
                  file = paste0('fits/meta-', this_ROI))
# pp_check(meta)
# summary(meta)

# parnames(meta)
# bayesplot::mcmc_rank_overlay(meta, pars = c('b_Intercept', 'sd_sid__Intercept'))

# TRR_z_posteriors <- data.table::rbindlist(lapply(r, function(r_){
#   fm <- r_[['fm']]
# 
#   pars <- parnames(fm)
#   corsesscondname <- pars[grepl('cor_\\w+__sess\\d{2}:cond_code__sess\\d{2}:cond_code$', pars)]
#   mm <- posterior_samples(fm, pars = corsesscondname)[[1]]
# 
#   z <- atanh(mm)
# 
#   return(data.frame(z = z, pair = r_[['name']]))
# }))

# mean(TRR_z_posteriors$z)
# sd(TRR_z_posteriors$z)
# 
# get_prior(z ~ 1 + (1 | pair), data = TRR_z_posteriors)
# meta_2 <- brms::brm(z ~ 1 + (1 | pair), data = TRR_z_posteriors,
#                     cores = chains, chains = chains, iter = 2000,
#                     file = paste0('fits/meta_2-', this_ROI))
# pp_check(meta_2)
# summary(meta_2)

### Big model + 2nd stage AR approach:
d_bigmodel <- d[ROI == this_ROI]

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


### Big Full Covar model:
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
bayesplot::mcmc_areas_ridges(fm_bigcov, pars = sprintf('cor_id__sess%02d:cond_code__sess%02d:cond_code', 1:9, 2:10), 
                             prob = .95, prob_outer = .99)

