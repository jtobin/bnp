PROJECT_DIR    = $(HOME)/projects/bnp
DATA_DIR       = $(PROJECT_DIR)/data/
INPUT_DATA_DIR = $(PROJECT_DIR)/data/input
RAW_DATA_DIR   = $(DATA_DIR)/raw

# http://data.worldbank.org/indicator/NY.GDP.MKTP.KD.ZG
RAW_GDP_DATA_URL = 'http://api.worldbank.org/v2/en/indicator/ny.gdp.mktp.kd.zg?downloadformat=csv'

# get raw data
$(RAW_DATA_DIR)/ny.gdp.mktp.kd.zg_Indicator_en_csv_v2.zip:
	curl $(RAW_GDP_DATA_URL) > $@

# clean raw data (FIXME do more here)
ny.gdp.mktp.kd.zg_Indicator_en_csv_v2.csv: \
		$(RAW_GDP_DATA)
	unzip $< $@ && mv $@ $(INPUT_DATA_DIR)

