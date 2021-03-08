library(brms)
library(data.table)
library(parallel)
### Note: ROIs in processed data are off by -1 relative to Sch400 mask voxel
#values. In other words, all fits are labeled with ROI 0-399 and all Sch400 mask
#voxel values are 1-400.

### Note: Sch400 ROIs for dACC: 

meta_fn <- dir('fits', pattern = 'meta-.*rds', full.names = TRUE)
AR_fn <- dir('fits', pattern = 'fit.*bigAR\\.rds', full.names = TRUE)
old_icc_fn <- dir('fits', pattern = 'oldicc-.*RDS', full.names = TRUE)

# anAR <-  readRDS(AR_fn[grepl('311', AR_fn)])
# summary(ameta)
# mcmc_areas(ameta, pars = 'b_Intercept', prob = .95, prob_outer = .99, 
#            transformations = 'tanh') + 
#   ggplot2::coord_cartesian(x = c(-1, 1)) + 
# mcmc_areas(anAR, pars = 'ar[1]', prob = .95, prob_outer = .99) + 
#   ggplot2::coord_cartesian(x = c(-1, 1))

cl <- parallel::makePSOCKcluster(4)
parallel::clusterExport(cl = cl, varlist = c('meta_fn'))
nada <- parallel::clusterEvalQ(cl, {library(data.table); library(brms)})

ma_ests <- parallel::parLapply(cl = cl, X = (1:length(meta_fn)), fun = function(i){
  ameta <- readRDS(meta_fn[[i]])
  ROI <- as.numeric(gsub('fits/meta-(\\d{1,3}).rds', '\\1', meta_fn[[i]]))
  r <- data.table(t(fixef(ameta)['Intercept', c('Estimate', 'Q2.5', 'Q97.5')]))
  r$ROI <- ROI
  return(r)
})

parallel::stopCluster(cl)

ma_ests_dt <- rbindlist(ma_ests)

cl <- parallel::makePSOCKcluster(4)
parallel::clusterExport(cl = cl, varlist = c('AR_fn'))
nada <- parallel::clusterEvalQ(cl, {library(data.table); library(brms)})

AR_ests <- parallel::parLapply(cl = cl, X = (1:length(AR_fn)), fun = function(i){
  anAR <- readRDS(AR_fn[[i]])
  ROI <- as.numeric(gsub('fits/fit-ROI(\\d{1,3})-bigAR.rds', '\\1', AR_fn[[i]]))
  parnames(anAR)
  r <- data.table(t(brms::posterior_summary(anAR, pars = 'ar\\[1\\]')['ar[1]', c('Estimate', 'Q2.5', 'Q97.5')]))
  r$ROI <- ROI
  return(r)
})

parallel::stopCluster(cl)

AR_ests_dt <- rbindlist(AR_ests)

cl <- parallel::makePSOCKcluster(4)
parallel::clusterExport(cl = cl, varlist = c('old_icc_fn'))
nada <- parallel::clusterEvalQ(cl, {library(data.table)})

oldicc_ests <- parallel::parLapply(cl = cl, X = (1:length(old_icc_fn)), fun = function(i){
  anoldicc <- readRDS(old_icc_fn[[i]])
  return(anoldicc)
})

parallel::stopCluster(cl)

oldicc_ests_dt <- rbindlist(oldicc_ests)

#Put the meta-analytic results back on r scale rather than z scale
cols <- c('Estimate', 'Q2.5', 'Q97.5')
ma_ests_dt[, (cols) := lapply(.SD, tanh), .SDcols = cols]
#Ensure ROI ids match up with actual sch400 1-400 IDs
ma_ests_dt[, ROI := ROI + 1]
AR_ests_dt[, ROI := ROI + 1]

setorder(ma_ests_dt, Estimate)

library(ggplot2)
ma_ests_dt[, ROI_fac := factor(ROI, levels = ROI[order(Estimate)])]
ggplot(ma_ests_dt, aes(x = Estimate, y = ROI_fac)) + 
  geom_errorbarh(aes(xmin = Q2.5, xmax = Q97.5), alpha = .2) + 
  geom_point(alpha = .5) + 
  scale_y_discrete(breaks = c(366, 330, 332, 182))

oldicc_ests_meta_dt <- oldicc_ests_dt[sess == 'meta']
oldicc_ests_meta_dt[, ROI_fac := factor(ROI, levels = ROI[order(Estimate)])]
ggplot(oldicc_ests_meta_dt, aes(x = Estimate, y = ROI_fac)) + 
  geom_errorbarh(aes(xmin = Q2.5, xmax = Q97.5), alpha = .2) + 
  geom_point(alpha = .5) + 
  scale_y_discrete(breaks = c(366, 330, 332, 182))

sch400_dump <- as.data.table(readr::read_delim('sch400_maskdump.tsv', delim = ' ',  col_names = c('i', 'j', 'k', 'ROI')))
fear.icc.np_dump <- as.data.table(readr::read_delim('fear.icc_MEAN_maskdump.tsv', delim = ' ',  col_names = c('i', 'j', 'k', 'fear.icc')))

ROI_out_dt <- ma_ests_dt[sch400_dump, on = 'ROI']
ROI_out_dt[is.na(Estimate), c('Estimate', 'Q2.5', 'Q97.5') := 0]

readr::write_delim(ROI_out_dt[, c('i', 'j', 'k', 'Estimate')], 'TRR_meta_est.txt', delim = ' ', col_names = FALSE)
readr::write_delim(ROI_out_dt[, c('i', 'j', 'k', 'Q2.5')], 'TRR_meta_ll.txt', delim = ' ', col_names = FALSE)
readr::write_delim(ROI_out_dt[, c('i', 'j', 'k', 'Q97.5')], 'TRR_meta_uu.txt', delim = ' ', col_names = FALSE)

ROI_oldicc_out_dt <- oldicc_ests_meta_dt[sch400_dump, on = 'ROI']
ROI_oldicc_out_dt[is.na(Estimate), c('Estimate', 'Q2.5', 'Q97.5') := 0]

readr::write_delim(ROI_oldicc_out_dt[, c('i', 'j', 'k', 'Estimate')], 'oldICC_meta_est.txt', delim = ' ', col_names = FALSE)
readr::write_delim(ROI_oldicc_out_dt[, c('i', 'j', 'k', 'Q2.5')], 'oldICC_meta_ll.txt', delim = ' ', col_names = FALSE)
readr::write_delim(ROI_oldicc_out_dt[, c('i', 'j', 'k', 'Q97.5')], 'oldICC_meta_uu.txt', delim = ' ', col_names = FALSE)

comparison_DT <- ROI_out_dt[fear.icc.np_dump, on = c('i', 'j', 'k')]
comparison_DT_sum <- comparison_DT[, list(mean_icc = mean(fear.icc), 
                                          meta_trr = unique(Estimate)), by = ROI]
  
ggplot(comparison_DT, aes(x = Estimate, y = fear.icc)) + 
  geom_hex(aes(color = ..count.., fill = ..count..), binwidth = c(.001, .01)) + 
  geom_smooth(method = 'lm') + 
  geom_point(data = comparison_DT_sum, aes(x = meta_trr, y = mean_icc), color = 'yellow')

psych::corr.test(comparison_DT_sum[, c('mean_icc', 'meta_trr')])

oldICC_comparison_DT <- ROI_oldicc_out_dt[fear.icc.np_dump, on = c('i', 'j', 'k')]
oldICC_comparison_DT_sum <- oldICC_comparison_DT[, list(mean_icc = mean(fear.icc), 
                                                        meta_icc = unique(Estimate)), by = ROI]

ggplot(oldICC_comparison_DT, aes(x = Estimate, y = fear.icc)) + 
  geom_hex(aes(color = ..count.., fill = ..count..), binwidth = c(.001, .01)) + 
  geom_smooth(method = 'lm') + 
  geom_point(data = oldICC_comparison_DT_sum, aes(x = meta_icc, y = mean_icc), color = 'yellow')

psych::corr.test(comparison_DT_sum[, c('mean_icc', 'meta_trr')])

TRR_old_icc_comp_DT <- merge(ma_ests_dt, oldicc_ests_meta_dt, by = 'ROI', suffixes = c('.trr', '.icc'))
ggplot(TRR_old_icc_comp_DT, aes(x = Estimate.trr, y = Estimate.icc)) + 
  geom_errorbar(aes(ymin = Q2.5.icc, ymax = Q97.5.icc), alpha = .1) +
  geom_errorbarh(aes(xmin = Q2.5.trr, xmax = Q97.5.trr), alpha = .1) +
  geom_point(alpha = .8, size = .5) + 
  geom_smooth(method = 'lm')
