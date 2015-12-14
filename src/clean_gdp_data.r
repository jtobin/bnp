#!/usr/bin/Rscript

HOME = Sys.getenv('HOME')

wb_gdp_data_file = 'ny.gdp.mktp.kd.zg_Indicator_en_csv_v2.csv'

data_dir         = paste(HOME, 'projects/bnp/data', sep = '/')
working_data_dir = paste(data_dir, 'working', sep = '/')
input_data_dir   = paste(data_dir, 'input', sep = '/')
wb_gdp_data      = paste(working_data_dir, wb_gdp_data_file, sep = '/')

prune = function(data) {
    pruned = data[,c('Country.Name', 'X2014')]
    names(pruned)  = c('country', 'rate')
    completes_only = pruned[complete.cases(pruned),]
    return(completes_only)
  }

d = read.csv(wb_gdp_data, header = T, skip = 4)

write.csv(
    prune(d)
  , paste(input_data_dir, 'gdp.csv', sep = '/')
  , row.names = F
  )

