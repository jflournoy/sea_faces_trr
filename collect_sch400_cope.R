library(data.table)
setDTthreads(4)
fn <- dir('sch400_cope/', full.names = TRUE)
d <- rbindlist(lapply(fn, fread))
d[, c('V2', 'V3') := NULL]
regex <- '/mnt/stressdevlab/stress_pipeline/(\\d{4})/(month\\d{2})/faceReactivity/FaceReactivity1.feat/reg_standard/stats/cope4.*'
d[, c('id', 'sess', 'name') := 
    list(gsub(regex, '\\1', name),
         gsub(regex, '\\2', name),
         NULL)]

d_l <- melt(d, id.vars = c('id', 'sess'))
d_l[, c('ROI', 'variable') := 
      list(gsub('Mean_(\\d{1,3})', '\\1', variable),
           NULL)]
setnames(d_l, 'value', 'y')
saveRDS(d_l, 'sch400_cope_fGTc.RDS')
