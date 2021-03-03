library(brms)
library(data.table)

f <- dir('sch400/', pattern = '.*tsv', full.names = TRUE)
d <- rbindlist(lapply(f, fread))

# #For testing
# d_ <- copy(d)
# #d <- copy(d_)

contrast_map <- c("cope10" = "calm",
                  "cope11" = "fear",
                  "cope12" = "happy",
                  "varcope10" = "calm",
                  "varcope11" = "fear",
                  "varcope12" = "happy")

d[, c('id', 'session', 'est') := 
    list(gsub('.*/(\\d{4})/.*', '\\1', name),
         gsub('.*/(month\\d{2})/.*', '\\1', name),
         gsub('.*/stats/(.*).nii.gz.*', '\\1', name))]
d[, contrast := contrast_map[est]]
d[, est := fifelse(grepl('var', est), 'se', 'y')]
d[, c('name', 'V2', 'V3') := NULL]

d_l <- data.table::melt(d, id.vars = c('id', 'session', 'est', 'contrast'))
d_l[, c('ROI', 'variable') := 
      list(gsub('Mean_', '', variable), NULL)]

d_w <- data.table::dcast(d_l, ... ~ est, value.var = 'value')
d_w[, se := se^.5] #varcopes are actually the variance error not the standard error so we need to sqrt
saveRDS(d_w, 'sch400_data.RDS')
