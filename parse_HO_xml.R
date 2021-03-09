library(XML)
x <- xmlParse('HarvardOxford-Subcortical.xml')
x_l <- xmlToList(x)
d <- data.frame(label = unlist(lapply(x_l$data, `[[`, 'text')),
           index = unlist(lapply(x_l$data, function(l){
             as.numeric(l$.attrs[['index']]) + 1
           })))
readr::write_csv(d, 'HO_sub.csv')

