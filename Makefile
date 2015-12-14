PROJECT_DIR    = $(HOME)/projects/bnp

DATA_DIR         = $(PROJECT_DIR)/data
RAW_DATA_DIR     = $(DATA_DIR)/raw
WORKING_DATA_DIR = $(PROJECT_DIR)/data/working
INPUT_DATA_DIR   = $(PROJECT_DIR)/data/input

SRC_DIR = $(PROJECT_DIR)/src
CLEAN_DATA_SCRIPT = $(SRC_DIR)/clean_gdp_data.r

# http://data.worldbank.org/indicator/NY.GDP.MKTP.KD.ZG
RAW_GDP_DATA_URL = 'http://api.worldbank.org/v2/en/indicator/ny.gdp.mktp.kd.zg?downloadformat=csv'
RAW_GDP_DATA_BASE = 'ny.gdp.mktp.kd.zg_Indicator_en_csv_v2'

$(RAW_DATA_DIR)/$(RAW_GDP_DATA_BASE).zip:
	curl $(RAW_GDP_DATA_URL) > $@

$(WORKING_DATA_DIR)/%.csv: \
		$(RAW_DATA_DIR)/%.zip
	unzip $< $(notdir $@) -d $(WORKING_DATA_DIR)

$(INPUT_DATA_DIR)/%.csv: \
		$(WORKING_DATA_DIR)/%.csv
	$(CLEAN_DATA_SCRIPT)

all: $(INPUT_DATA_DIR)/$(RAW_GDP_DATA_BASE).csv

