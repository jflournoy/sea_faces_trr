library(data.table)
ev_files <- unlist(lapply(
  unlist(lapply(dir('/ncf/mclaughlin/stressdevlab/stress_pipeline/', pattern = '1\\d{3}', full.names = T), 
                dir, pattern = 'month\\d{2}', full.names = T)),
  function(basedir){
    dir(file.path(basedir, 'faceReactivity/evfiles/'), pattern = '*txt', full.names = T)
  }))

files_r1 <- ev_files[grepl('1_[CHF].txt$', ev_files)]
names(files_r1) <- files_r1
d <- data.table::rbindlist(lapply(files_r1, data.table::fread, 
                                  col.names = c('onset', 'duration', 'V3'), 
                                  header = FALSE), 
                           idcol = 'file')

file_reg <- '.*(\\d{4}).*(month\\d{2}).*FaceReactivity1_([CHF]).txt'
d[, c('id', 'sess', 'cond') :=
     list(gsub(file_reg, '\\1', file),
          gsub(file_reg, '\\2', file),
          gsub(file_reg, '\\3', file))]


d[, newfname := file.path('')]