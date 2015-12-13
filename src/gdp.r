project_dir = '/Users/jtobin/projects/bnp/'
data_dir = paste(project_dir, 'data/input', sep = '/')

wb_gdp_data = paste(
    data_dir
  , 'ny.gdp.mktp.kd.zg_Indicator_en_csv_v2.csv'
  , sep = '/'
  )

d = read.csv(
    wb_gdp_data
  , header = T
  , skip = 4
  )


