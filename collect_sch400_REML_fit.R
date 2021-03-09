library(data.table)
setDTthreads(5)

collect_data <- function(filelist, outname, fn_regex){
  if(!file.exists(outname)){
    col.names_ <- expand.grid(stat = c('est', 't-score'), basis = c('t', 'dt'), block = 1:3, cond = c('Calm', 'Fear', 'Happy'))
    col.names_ <- col.names_[, dim(col.names_)[[2]]:1]
    col.names <- c('F', apply(col.names_, 1, paste, collapse = '_'))
    
    d <- rbindlist(lapply(filelist, fread, skip = 9, header = FALSE, nrows = 1, col.names = col.names))
    d$file <- filelist
    readr::write_rds(d, outname)
  } else {
    d <- readRDS(outname)
  }
  
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
  return(d_w)
}

f <- dir('afni_out/', pattern = '^\\d{4}.*REML.*.1D', full.names = TRUE)
fn <- 'sch400_REML_fit.RDS'
f_ho <- dir('afni_out/', pattern = '^HO_\\d{4}.*REML.*.1D', full.names = TRUE)
fn_ho <- 'HO_REML_fit.RDS'
fn_regex <- 'afni_out//(\\d{4})_(\\d{2})_(\\d+)_REMLfit.1D'
fn_regex_ho <- 'afni_out//HO_(\\d{4})_(\\d{2})_(\\d+)_REMLfit.1D'

d_w_sch <- collect_data(f, fn, fn_regex)
d_w_ho <- collect_data(f_ho, fn_ho, fn_regex_ho)

readr::write_rds(d_w, 'sch400_RML_fit_processed.RDS')
readr::write_rds(d_w_ho, 'HO_RML_fit_processed.RDS')
