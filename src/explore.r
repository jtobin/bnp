#!/usr/bin/Rscript

data_dir       = paste(HOME, 'projects/bnp/data', sep = '/')
input_data_dir = paste(data_dir, 'input', sep = '/')

gdp_data = paste(input_data_dir, 'gdp.csv', sep = '/')

d = read.csv(gdp_data, header = T, colClasses = c('factor', 'numeric'))

require(ggplot2)

g = ggplot(d[,2], aes(rate))
