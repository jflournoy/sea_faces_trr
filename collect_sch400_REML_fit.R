library(data.table)
setDTthreads(5)

col.names_ <- expand.grid(stat = c('est', 't-score'), basis = c('t', 'dt'), block = 1:3, cond = c('Calm', 'Fear', 'Happy'))
col.names_ <- col.names_[, dim(col.names_)[[2]]:1]
col.names <- c('F', apply(col.names_, 1, paste, collapse = '_'))

f <- dir('afni_out/', pattern = '.*REML.*.1D', full.names = TRUE)
fn <- 'sch400_REML_fit.RDS'
if(!file.exists(fn)){
  d <- rbindlist(lapply(f, fread, skip = 9, header = FALSE, nrows = 1, col.names = col.names))
  saveRDS(d, fn)
} else {
  d <- readRDS(fn)
  d$file <- f
}

fn_regex <- 'afni_out//(\\d{4})_(\\d{2})_(\\d+)_REMLfit.1D'

d <- d[, F := NULL]
d[, c('id', 'sess', 'ROI') := 
    list(gsub(fn_regex, '\\1', file),
         gsub(fn_regex, '\\2', file),
         gsub(fn_regex, '\\3', file))]
d[, file := NULL]

d_l <- melt(d, id.vars = c('id', 'sess', 'ROI'))

var_regex <- '(\\w+)_([123])_(\\w{1,2})_(est|t-score)'
d_l[, c('cond', 'block', 'basis', 'stat') :=
      list(gsub(var_regex, '\\1', variable),
           gsub(var_regex, '\\2', variable),
           gsub(var_regex, '\\3', variable),
           gsub(var_regex, '\\4', variable))]
d_l[, variable := NULL]
d_w <- data.table::dcast(d_l, ... ~ stat, value.var = 'value')

# t = est/se
# t*se = est
# se = est/t

d_w[, se := est / `t-score`]
saveRDS(d_w, 'sch400_RML_fit_processed.RDS')
